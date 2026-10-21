#----------------------------------------#
#   literally a comparison plot betwn    #
#   particle spray, nbody thin stream,   #
#   nbody thin stream + cocoon  `        #
#----------------------------------------#
# %%
import sys, os
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np
import h5py

from scipy.stats import binned_statistic, norm, multivariate_normal
from scipy.optimize import curve_fit, minimize, Bounds
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline
from scipy.special import expit, logit, logsumexp # inverse-logit / logit, for the f_1 reparameterization


# import astropy.coordinates as coord
from astropy.coordinates import Galactocentric, ICRS, CartesianRepresentation,CartesianDifferential
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.constants as const
from astropy.table import Table

import gala.coordinates as gc
import gala.dynamics as gd
import gala.potential as gp
from gala.dynamics import mockstream as ms
from gala.units import galactic
from gala.coordinates import reflex_correct

import matplotlib.pyplot as plt
# %matplotlib inline
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse

plt.style.use(script_path+'/vedant.mplstyle')
# %config InlineBackend.figure_format='retina'
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
import cmasher as cm

from tqdm import tqdm

sys.path.append(script_path)
sys.path.append(repo_path)
from results import case_label, gmm_selection, observed, observed_errors, p_cocoon_from_table #<--- TODO make results.py so that not all the stuff runs when i import functions. 
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
import inspect_new_sims as simspect 
import pickle
import gmm
import noise as n
# %%
datapath = '/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
table_path = "/n/home02/amphillips/p27_nbody/data/gmm_tables/"
grid_info = paf.extended_grid_info(scratch=False) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors[:-1] + lm_colors[::-1]
cc = reordered_colors[:-1]
ccc = cc[1:]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')

# ordering decided here. 
orbits = np.array(['gd1','aau','pa5','jet','c19'])
# %%
#--------------------------------------------------------#
#       FIGURE 1: every rvir0 = 0.75 pc stream at        #
#       present day, galactocentric z vs x.              #
#       left: thin streams only (cocoons removed);       #
#       right: everything unbound from the progenitor.   #
#--------------------------------------------------------#

#--------------------------------------------------------#
#   Chen particle spray streams, one per orbit, to        #
#   match the trajectory of each rvir0 = 0.75 pc, hm       #
#   N-body stream: same t=0 phase-space position and       #
#   potential (BovyMWPotential2014) as petar, with the     #
#   DF's progenitor mass following the N-body bound mass   #
#   track from data.tidal. saved as hdf5 (every frame_dt   #
#   snapshot, not just the present day) to data/particle_  #
#   sprays/ so they don't have to be regenerated. adapted   #
#   from animations/particle_spray_vs_nbody.py.             #
#--------------------------------------------------------#
spray_dir = repo_path+"/data/particle_sprays/"
os.makedirs(spray_dir, exist_ok=True)

fig1_stellar_pop = 'hm'
fig1_rvir_index = 0 #<-- 0.75 pc, matches fig1_rvir below
fig1_N_spray = 5000 #<-- total particles at the present day, both tails
fig1_spray_dt = 1 #Myr
fig1_frame_dt = 10 #Myr, matches petar's snapshot cadence

def find_finished_copy(orbit, stellar_pop, rvir_index, copy_options):
    """
    retrieve_sim_info()'s `copy` just picks a subdirectory; not every copy ran to the present
    day. try them in copy_options order (the same order final_datasets.py uses) and take the
    first whose data.tidal reaches the present day -- cheap, unlike prepare_nbody_data_anycopy.
    """
    for copy in copy_options:
        path, apo, age, init_displacement = grid_info.retrieve_sim_info(
            orbit=orbit, stellar_pop=stellar_pop, rvir_index=rvir_index, copy=copy)
        tidal = paf.load_tidal(path)
        if tidal.time.max() >= age:
            return path, apo, age, init_displacement, copy
    raise FileNotFoundError("no finished copy for %s rvir_index=%i"%(orbit, rvir_index))

def nbody_bound_mass_track(path, t):
    """bound mass from petar's data.tidal (10 Myr cadence), linearly interpolated onto t [Myr]."""
    tidal = paf.load_tidal(path)
    return np.interp(t, tidal.time, tidal.mass)*u.Msun

def make_spray(path, init_displacement, age, rvir0, spray_file,
               dt=fig1_spray_dt, output_every=int(fig1_frame_dt/fig1_spray_dt), N=fig1_N_spray):
    """
    Chen-type particle spray from the same t=0 progenitor phase-space position as petar, in
    the same potential (BovyMWPotential2014). writes every frame_dt snapshot to spray_file.
    """
    n_steps = int(age/dt)
    t = np.arange(n_steps+1)*dt
    n_per_tail = int(N/2)
    release_steps = np.round(np.linspace(0, n_steps, n_per_tail)).astype(int)
    assert len(np.unique(release_steps))==n_per_tail, "more releases than steps; lower dt"
    n_particles = np.zeros(n_steps+1, dtype=int)
    n_particles[release_steps] = 1

    prog_mass = nbody_bound_mass_track(path, t)
    plummer_b = rvir0*3*np.pi/16 #<-- virial radius of a plummer sphere is (16/3pi) b

    prog_w0 = gd.PhaseSpacePosition(init_displacement[:3]*u.kpc, init_displacement[3:]*(u.km/u.s))
    prog_pot = gp.PlummerPotential(m=prog_mass[0], b=plummer_b, units=galactic)
    df = ms.ChenStreamDF(lead=True, trail=True)
    H = gp.Hamiltonian(gp.BovyMWPotential2014(units=galactic))
    gen = ms.MockStreamGenerator(df, H, progenitor_potential=prog_pot)

    stream, prog = gen.run(prog_w0, prog_mass, dt=dt*u.Myr, n_steps=n_steps,
                           n_particles=n_particles,
                           output_every=output_every, output_filename=spray_file, overwrite=True,
                           progress=True)
    print("wrote %s: %i particles at the present day"%(spray_file, len(stream.x)))

def load_spray_present_day(spray_file):
    """galactocentric (x,y,z) [kpc] of the released particles at the last saved snapshot."""
    with h5py.File(spray_file, 'r') as f:
        stream_pos = f['stream']['pos'][:, -1, :] #<-- kpc
    released = np.all(np.isfinite(stream_pos), axis=0)
    return stream_pos[:, released]

pericenters_kpc = []

orbit_titles = ['GD-1', 'AAU','Pal 5','Jet','C-19']
for ii, orbit in enumerate(tqdm(orbits)): #<--- this i can do later i think. 

    pr = prog_tab[prog_tab['name']==orbit] #<-- prog row
    stream_age = pr['age_gyr'][0]

    stream_age = np.round(stream_age+0.1, 2)


    disp = [pr['x'][0], pr['y'][0], pr['z'][0], pr['vx'][0], pr['vy'][0], pr['vz'][0]]
    orbit_obj = paf.integrate_prog_orbit(disp, steps=10000, dt=1*u.Myr)

    peri = orbit_obj.pericenter().to(u.kpc).value
    pericenters_kpc.append(peri)

orbit_titles = np.array(orbit_titles)

pericenters_kpc = np.array(pericenters_kpc)

# %%
fit_keys = ['phi2','v_phi1','v_phi2','v_gsr'] #<-- gmm.py's keys; also the order of the M_*/S_* columns

fig1_rvir = 0.75
fig1_noise, fig1_include_binaries = None, False #<-- the gmm.py fit that decides who's in the cocoon
# right panel: True = every unbound star, including the ones trim_new dropped before the fit
#   (phi1 tails, distance outliers, early escapers), so left vs right is cocoon + trim.
#   False = only the stars the gmm saw, so left vs right is exactly the cocoon.
fig1_include_trimmed = True
tt_fig1 = Table.read(table_path+case_label(fig1_noise, fig1_include_binaries)+".fits", format="fits")

orbit_title_map = {o:t for o, t in zip(orbits, orbit_titles)}
reordered = np.argsort(pericenters_kpc) #<-- same as the main result plot

#### load everything first, then plot.
#   loop in pericenter order so ccc[ii] is the same color per orbit as in the summary plots.
fig1_x, fig1_y, fig1_z = [], [], [] # galactocentric, kpc
fig1_n_sim = [] # which stream each star belongs to (index into orbits[reordered])
fig1_thin = [] # in the gmm fit AND p_cocoon < 0.5
fig1_fit = [] # in the gmm fit at all

fig1_spray_x, fig1_spray_y, fig1_spray_z = [], [], [] # galactocentric, kpc
fig1_spray_n_sim = [] # same indexing as fig1_n_sim

for ii, orbit in enumerate(tqdm(orbits[reordered])):
    # if orbit=="pa5":
    #     continue
    filename = datapath+"%s_%.2f.pickle"%(orbit, fig1_rvir)
    with open(filename, 'rb') as handle:
        data_dict = pickle.load(handle)

    #### Chen particle spray on the same trajectory, for the third panel
    copy_options = [4,3,2,1,0] if orbit=='aau' else [0,1,2,3,4] #<-- same order final_datasets.py uses
    spray_path, spray_apo, spray_age, spray_init_displacement, spray_copy = find_finished_copy(
        orbit, fig1_stellar_pop, fig1_rvir_index, copy_options)

    spray_file = spray_dir+"spray_%s_%s_%.2f_%i.h5"%(orbit, fig1_stellar_pop, fig1_rvir, spray_copy)
    if not os.path.exists(spray_file):
        print("DUDE THIS IS NOT ONE WON'T WORK, NO SPRAY FilE")
        continue
        # print("particle spraying %s"%orbit)
        # make_spray(spray_path, spray_init_displacement, spray_age, fig1_rvir*u.pc, spray_file)

    xs, ys, zs = load_spray_present_day(spray_file)
    fig1_spray_x.append(xs)
    fig1_spray_y.append(ys)
    fig1_spray_z.append(zs)
    fig1_spray_n_sim.append(np.full(len(xs), fill_value=ii))

    use, sc_fig1, noise_fig1 = gmm_selection(data_dict, noise=fig1_noise, include_binaries=fig1_include_binaries)
    row = tt_fig1[(tt_fig1['orbit']==orbit) & (tt_fig1['Rvir0']==fig1_rvir)]
    x_data = np.column_stack([observed(k, sc_fig1, noise_fig1, use) for k in fit_keys])

    p_cocoon = np.full(len(use), fill_value=np.nan) #<-- nan for stars that weren't in the fit (bound, or outside the trim)
    p_cocoon[use] = p_cocoon_from_table(x_data, row, observed_errors(data_dict, noise_fig1, fig1_noise, use))
    thin = use & (p_cocoon<0.5) #<-- nan < 0.5 is False

    #### the progenitor (2 rtid) is dropped from both panels; everything else unbound goes in the right one
    unbound = data_dict['unbound']
    x, y, z = data_dict['CoM']['pos'].to(u.kpc).value[unbound].T

    fig1_x.append(x)
    fig1_y.append(y)
    fig1_z.append(z)
    fig1_n_sim.append(np.full(len(x), fill_value=ii))
    fig1_thin.append(thin[unbound])
    fig1_fit.append(use[unbound])

    print("%s: %i unbound, %i in the fit, %i thin stream"%(orbit, unbound.sum(), use.sum(), thin.sum()))

fig1_x, fig1_y, fig1_z = np.concatenate(fig1_x), np.concatenate(fig1_y), np.concatenate(fig1_z)
fig1_n_sim = np.concatenate(fig1_n_sim)
fig1_thin = np.concatenate(fig1_thin)
fig1_fit = np.concatenate(fig1_fit)

fig1_spray_x = np.concatenate(fig1_spray_x)
fig1_spray_y = np.concatenate(fig1_spray_y)
fig1_spray_z = np.concatenate(fig1_spray_z)
fig1_spray_n_sim = np.concatenate(fig1_spray_n_sim)


n_orbits = len(orbits)
cmap_fig1 = mcolors.ListedColormap(ccc[:n_orbits])
# discrete norm: one color band per orbit, boundaries on the half-integers
norm_fig1 = mcolors.BoundaryNorm(np.arange(n_orbits+1)-0.5, cmap_fig1.N)

fig, axs = plt.subplots(1, 3, figsize=[24, 8], sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.05)


for ax in axs:
    for width, height in [(2, 2), (3, 2), (0.5, 0.5), (30, 0.5)]:
        ax.add_patch(Ellipse(xy=[0, 0], width=width, height=height, angle=0,
                            facecolor='none', edgecolor='black', linewidth=1))

#### leftmost panel: Chen particle spray on the matching trajectory (no progenitor to remove)
order = np.argsort(fig1_spray_x)
axs[0].scatter(fig1_spray_y[order], fig1_spray_z[order],
               c=fig1_spray_n_sim[order], cmap=cmap_fig1, norm=norm_fig1,
               s=1, rasterized=True)
axs[0].text(0.05, 0.95, 'Particle spray', ha='left', va='top', transform=axs[0].transAxes)

panel_selections = [fig1_thin, np.ones(len(fig1_thin), dtype=bool) if fig1_include_trimmed else fig1_fit]
panel_titles = ['Thin stream only', 'All stars']
for ax, sel, title in zip(axs[1:], panel_selections, panel_titles):
    # order = np.argsort(fig1_y[sel]) #<-- most negative y drawn first (bottom), most positive y last (top)
    order = np.argsort(fig1_x[sel])
    ax.scatter(fig1_y[sel][order], fig1_z[sel][order],
               c=fig1_n_sim[sel][order], cmap=cmap_fig1, norm=norm_fig1,
               s=1, rasterized=True)
    ax.text(0.05, 0.95, title, ha='left', va='top', transform=ax.transAxes)

for ax in axs:
    ax.set_xlabel(r'$y~[\rm kpc]$')
    ax.set_aspect('equal')

axs[0].set_xlim(-35, 35)
axs[0].set_ylim(-35, 35)
axs[0].set_ylabel(r'$z~[\rm kpc]$')



handles = [Line2D([], [], ls='', marker='o', markersize=8, color=ccc[ii], label=orbit_title_map[orbit])
           for ii, orbit in enumerate(orbits[reordered])]
axs[2].legend(handles=handles, loc='lower right', fontsize=15)

plt.savefig("plots/fig1_motiv_notrim_with_spray.pdf", dpi=300, bbox_inches='tight')