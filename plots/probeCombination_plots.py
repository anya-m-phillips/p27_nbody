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
# %%
orbit='gd1'
mass = 'hm'
rvir_index=3
init_displacement = grid_info.circ_init_displacement
(core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
    simspect.prepare_nbody_data_anycopy(
        orbit, stellar_pop=mass, rvir_index=rvir_index, copies=[0,1,2,3,4],
        include_photometry=False
    )

params1, _ = paf.get_orbital_params(paths=[path], times=[int(age/10)], n=0, rng=np.random.default_rng(42))
params2, _ = paf.get_orbital_params(paths=[path], times=[int(age/10)], n=0, rng=np.random.default_rng(43)) #<-- sneakily add 2x binaries lol
params = np.vstack([params1, params2])
params
obstimes = np.zeros((len(params), 1,1))
binary_rvs = paf.get_rvs(params, obstimes)[:,0][:,0] * u.km/u.s
# %%
coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab)
sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab)

trimmed_sc = simspect.clip_coords(sc, [inMW, trim]) #<-- this applies inMW, trim to the coordinate dictionary
sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction
# ol_clip = simspect.outlier_clip( #<-- avoid biasing the cocoon dispersion with a few crazy outliers. 
#     sc_straighter['v_gsr'], sc_straighter['pm_phi1'], sc_straighter['pm_phi2'] 
# )

# %%
def vkick(phi1, b, r0,
           wpar, wz, wr=0*u.km/u.s):
    Msh=1e6*u.Msun
    rs = 0.1 * u.kpc
    phi1_rad = phi1*u.degree.to(u.radian)
    wperp = np.sqrt(wr**2 + wz**2)
    wrel = np.sqrt(wpar**2 + wperp**2)

    num = 2*const.G*Msh*phi1_rad
    denom_1 = wrel*r0
    denom_2 = phi1_rad**2 + ((b**2+rs**2)*wrel**2 / (r0**2*wperp**2))
    denom = denom_1 * denom_2
    dv = -num/denom
    return dv.to(u.km/u.s)
    

dPhi1 = -45 #<-- where we will put the perturbation
r0= apo*u.kpc

x = sc_straighter['phi1']-dPhi1
y = sc_straighter['v_gsr']*u.km/u.s

y_new = y+vkick(x, 
                wpar=0*u.km/u.s,
                wz=100*u.km/u.s, wr=0*u.km/u.s,
                b=0.1*u.kpc, r0=r0)

### lying with the indexing a little bit...
# nsingles = data_dict['nsingles']
# nbinaries = data_dict['nbinaries']
nbinaries = int(2 * data_dict['nbinaries'])
step = 5


fig, ax = plt.subplots(figsize=[15,7])
ax.set_ylim(-10,10)
ax.set_xlim(-20,20)
ax.set_aspect('equal')
ax.set_xlabel(r'$\phi_1~[\degree]$')
ax.set_ylabel(r'$\Delta v_{\rm GSR}~[\rm km~s^{-1}]$')
# step 1:
if step==1:
    ax.scatter(x, y, c='k', s=.3)

# step 2:
if step==2:
    ax.scatter(x, y_new, c='k', s=.3)

# step 3:
if step==3:
    ax.scatter(x[:-nbinaries], y_new[:-nbinaries], c='k', s=.3)
    ax.scatter(x[-nbinaries:], y_new[-nbinaries:], c=cc[-6], s=.3)

# step 4
if step==4:
    ax.scatter(x[:-nbinaries], y_new[:-nbinaries], c='k', s=.3)
    ax.scatter(x[-nbinaries:], y_new[-nbinaries:]+binary_rvs, c='k', s=0.3)

if step<5:
    plt.savefig(repo_path+'/plots/probeCombination_workshop/binary_demo%i.pdf'%step, dpi=300, bbox_inches='tight')


if step==5:
    plt.close() # <-- get rid of other panel. 
    obstimes_multi = np.array(
        [[np.arange(0, 50.25, 0.25)]]*len(params)
    )
    binary_rvs_multi = paf.get_rvs(params, obstimes_multi)[:,0] * u.km/u.s

    for ii in tqdm(range(len(binary_rvs_multi[0]))):
        brvs = binary_rvs_multi[:,ii]
        fig, axx = plt.subplots(figsize=[15,7])

        axx.set_ylim(-10,10)
        axx.set_xlim(-20,20)
        axx.set_aspect('equal')
        axx.set_xlabel(r'$\phi_1~[\degree]$')
        axx.set_ylabel(r'$\Delta v_{\rm GSR}~[\rm km~s^{-1}]$')
        axx.scatter(x[:-nbinaries], y_new[:-nbinaries], c='k', s=.3)
        # axx.scatter(x[-nbinaries:], y_new[-nbinaries:]+brvs, c=cc[-6], s=5)
        axx.scatter(x[-nbinaries:], y_new[-nbinaries:]+binary_rvs, c='k', s=0.3)

        axx.text(0.38, 0.38, 'day %i'%(obstimes_multi[0][0][ii]), transform=axx.transAxes)

        plt.savefig(repo_path+f'/plots/probeCombination_workshop/binary_demo_mov/frame_{ii:05d}.png', dpi=300, bbox_inches='tight')
        plt.close()
# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%




# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
unbound = ~CMdict['in_rtid']
unbound = unbound[inMW][trim]
trimmed_sc = simspect.clip_coords(sc, [inMW, trim]) #<-- this applies inMW, trim to the coordinate dictionary
sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction
ol_clip = simspect.outlier_clip( #<-- avoid biasing the cocoon dispersion with a few crazy outliers. 
    sc_straighter['v_gsr'], sc_straighter['pm_phi1'], sc_straighter['pm_phi2'] 
)

# %%

# %%

# %%

# %%

# %%

# %%

# %%

orbits = ['gd1','aau','pa5','jet','m3','c19']
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
copy_options = [0,1,2,3,4]

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']
# cuts in phi2, pmphi1, pmphi2, v_gsr off-track to define cocoon. 
gd1_cuts = [0.75, 0.12,0.04, 1.5]
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