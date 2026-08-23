#--------------------------------------------------------------#
#   in inspect_new_sims.py i defined a bunch of functions.     #
#   that will help me define stream frames and make good cuts  # 
#   on the data; here i will use those functions to measure    # 
#   cocoon fractions for the different orbits and cluster      #
#   densities.                                                 #
#--------------------------------------------------------------#
# %%
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline


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
plt.style.use(script_path+'/vedant.mplstyle')
# %config InlineBackend.figure_format='retina'
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap

from tqdm import tqdm

sys.path.append(script_path)
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
import inspect_new_sims as simspect # lol idk. but i need basically all the funcitons. 

# first import... very beefy...
from sklearn.mixture import GaussianMixture
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty
# %%
# I'm going to define get_cocoon_selection here, so that I can easily edit it within this notebook 
def get_cocoon_selection(coords, cuts):
    """
    assume that coords and cuts are in the same order,
    i.e., the convention of 
    phi2, pmphi1, pmphi2, vgsr
    """
    coords_to_cut = {"phi2":coords['phi2'],"pm_phi1":coords['pm_phi1'],"pm_phi2":coords['pm_phi2'], 'v_gsr':coords['v_gsr']}

    selections = []
    for ii, key in enumerate(coords_to_cut.keys()):
        selection = np.abs(coords_to_cut[key])>cuts[ii]
        selections.append(selection)
    cocoon_selection = np.logical_or.reduce(selections) #<-- this should be true if any of the selections are true. 

    return cocoon_selection

# %%
# define stuff about the grid:
grid_info = paf.extended_grid_info(scratch=True) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')


orbits = ['gd1','aau','pa5','jet','m3','c19']
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
copy_options = [0,1,2,3,4]

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']
# do cuts as phi2, pmphi1, pmphi2, vgsr
gd1_cuts = [0.75, 0.15,0.04, 1.5]
aau_cuts = [0.5, 0.05, 0.0125, 2]
pa5_cuts = [0.3, 0.02, 0.04, 4.5]
jet_cuts = [0.4, 0.013, 0.005, 2.]
m3_cuts = [1.8, 0.4, 0.15, 4.75]
c19_cuts = [0.32, 0.023, 0.018, 2.75]


orbit_cuts = [
    gd1_cuts, aau_cuts, pa5_cuts, jet_cuts, m3_cuts, c19_cuts
]




dicts = []
sf_coords_obs = [] # before straightening
straight_sf_coords_obs = [] 

f_cocoons = []
vgsr_dispersions = []
phi2_dispersions = []


pericenters_kpc = []
init_displacements = [
    grid_info.gd1_init_displacement, 
    grid_info.aau_init_displacement,
    grid_info.pa5_init_displacement,
    grid_info.jet_init_displacement,
    grid_info.m3_init_displacement,
    grid_info.c19_init_displacement
]

for ii, orbit in enumerate(tqdm(orbits)):

    ## do the orbit-wise check -- integrate prog orbit and find the pericenter. 
    init_displacement = init_displacements[ii]
    orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
    peri = orbit_obj.pericenter()
    pericenters_kpc.append(peri.to(u.kpc).value)

    ### eventually eventually another inner loop will go here for masses.
    f_cocoons_this_orbit = []
    vgsr_dispersions_this_orbit = []
    phi2_dispersions_this_orbit = []

    mass_index = 1
    for rvir_index in range(4):
        (core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
            simspect.prepare_nbody_data_anycopy(
                orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
                include_photometry=False
            )

        # dicts.append(data_dict)

        coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab)
        # sf_coords_obs.append(sf_coords_obs)

        # straightened coords
        sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab)



        unbound = ~CMdict['in_rtid']
        unbound = unbound[inMW][trim]

        trimmed_sc = simspect.clip_coords(sc, [inMW, trim]) #<-- this applies inMW, trim to the coordinate dictionary
        sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction
        ol_clip = simspect.outlier_clip( #<-- avoid biasing the cocoon dispersion with a few crazy outliers. 
            sc_straighter['v_gsr'], sc_straighter['pm_phi1'], sc_straighter['pm_phi2'] 
        )
        cocoon_clips = orbit_cuts[ii]
        cocoon_selection = get_cocoon_selection(sc_straighter, cocoon_clips)     #<-- so now i'll want to index ol_clip & unbound & cocoon_selection

        use = ol_clip & unbound
        cocoon_selection = use & cocoon_selection
        thinStream_selection = use & ~cocoon_selection

        phi2 = sc_straighter['phi2']
        vgsr = sc_straighter['v_gsr']

        f_cocoon = len(phi2[cocoon_selection]) / len(phi2[use])

        sigma_phi2 = paf.dispersion_5_95(phi2[cocoon_selection])/2
        sigma_vgsr = paf.dispersion_5_95(vgsr[cocoon_selection])/2

        # sigma_phi2 = np.std(phi2[cocoon_selection])
        # sigma_vgsr = np.std(vgsr[cocoon_selection])

        f_cocoons_this_orbit.append(f_cocoon)
        vgsr_dispersions_this_orbit.append(sigma_vgsr)
        phi2_dispersions_this_orbit.append(sigma_phi2)

    f_cocoons.append(f_cocoons_this_orbit)
    vgsr_dispersions.append(vgsr_dispersions_this_orbit)
    phi2_dispersions.append(phi2_dispersions_this_orbit)
# %%
pericenters_kpc = np.array(pericenters_kpc)

reordered = np.argsort(pericenters_kpc)
pericenters_kpc = pericenters_kpc[reordered]

f_cocoons = np.array(f_cocoons)[reordered]
vgsr_dispersions = np.array(vgsr_dispersions)[reordered]
phi2_dispersions = np.array(phi2_dispersions)[reordered]
orbits = np.array(orbits)[reordered]

# ccc = plt.cm.magma(np.linspace(0.9, 0, len(pericenters_kpc)))
ccc = cc[1:]

fig, axs = plt.subplots(1,3,figsize=[21,7], sharex=True)

### iterate through orbits
for ii, orbit in enumerate(tqdm(orbits)):
    # if orbit=='m3':
    #     continue

    if orbit not in ['gd1','jet','c19','pa5','aau']:
        continue

    f_cocoons_this_orbit = f_cocoons[ii]
    vgsr_dispersions_this_orbit = vgsr_dispersions[ii]
    phi2_dispersions_this_orbit = phi2_dispersions[ii]

    x = rvirs
    axs[0].plot(x, f_cocoons_this_orbit, label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[ii],
                marker='o', color=ccc[ii], markersize=10)
    axs[1].plot(x, phi2_dispersions_this_orbit,
                marker='o', color=ccc[ii], markersize=10)
    axs[2].plot(x, vgsr_dispersions_this_orbit,
                marker='o', color=ccc[ii], markersize=10)

axs[0].legend(loc='upper right')
for ax in axs:
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$R_{\rm vir, 0}~[\rm kpc]$')
axs[0].set_ylabel(r'$f_{\rm cocoon}$')
axs[1].set_ylabel(r'$\sigma_{\phi_2, \rm cocoon}~[\degree]$')
axs[2].set_ylabel(r'$\sigma_{v_{\rm GSR, cocoon}}~[\rm km~s^{-1}]$')
axs[0].set_ylim(0, 0.13)
axs[1].set_ylim(0, 2.1)
axs[1].set_yticks([0, 0.5, 1, 1.5, 2])
axs[2].set_ylim(0, 23/2)


plt.savefig(repo_path+'/plots/probeCombination_workshop/prog_summary_all_%s.pdf'%masses[mass_index])
# plt.savefig(repo_path+'/plots/probeCombination_workshop/prog_summary_gd1_jet_c19_%s.pdf'%masses[mass_index])
# plt.savefig(repo_path+'/plots/probeCombination_workshop/prog_summary_gd1_%s.pdf'%masses[mass_index])



# plt.savefig('/n/home02/amphillips/p27_nbody/plots/cocoon_separation/%s/prog_summary.pdf'%masses[mass_index],
#              dpi=300, bbox_inches='tight')



