#--------------------------------------------------------------#
#   in inspect_new_sims.py i defined a bunch of functions.     #
#   that will help me define stream frames and make good cuts  # 
#   on the data; here i will use those functions and 
#   develop cocoon separation using a GMM instead
#   of hard cuts. inspired by jarvis+2026, we will start
#   with a maximized likelihood of the GMM, then use 
#   emcee or similar to sample posteriors. hopefully this
#   will be more robust/useful for down-sampled "observed" 
#   data.                                                      #
#--------------------------------------------------------------#
# %%
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

from scipy.stats import binned_statistic, norm, multivariate_normal
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
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty
# %%
# define stuff about the grid:
grid_info = paf.extended_grid_info(scratch=True) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')


orbits = ['gd1','aau','pa5','jet','m3','c19']
init_displacements = [
    grid_info.gd1_init_displacement, 
    grid_info.aau_init_displacement,
    grid_info.pa5_init_displacement,
    grid_info.jet_init_displacement,
    grid_info.m3_init_displacement,
    grid_info.c19_init_displacement
]
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
# copy_options = [0,1,2,3,4]
copy_options = [4,3,2,1,0]

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']
# %%
def gmm_likelihood_multivariate(x_data, 
                   f_1, Mu_1, Sigma_1, 
                   Mu_2, Sigma_2
                   ):
    """
    Mu_1, Mu_2 should be 4-vectors (phi2, transverse+radial velocities)
    Sigma_1, Sigma_2 should I guess be covariance matrices. yikes. 
    """
    Q_1 = f_1
    Q_2 = 1 - f_1 
    component_1 = Q_1 * multivariate_normal.pdf(x_data, mean=Mu_1, cov=Sigma_1)
    component_2 = Q_2 * multivariate_normal.pdf(x_data, mean=Mu_2, cov=Sigma_2)
    li = component_1 + component_2

    ln_li = np.log(li)
    ln_L = np.sum(ln_li)

    return ln_L # idk how to test. 
# except that idk what a covariance matrix is so let's take a different approach:
def gmm_likelihood_simple(x_data,
                          f_1, mu_1, sigma_1, mu_2, sigma_2):
    """
    now mu's are still 4-vectors, but so are sigmas. 
    """
    Q_1 = f_1
    Q_2 = 1 - f_1 

    L_1, L_2 = 1,1 #<-- we will add to these. 
    for k in range(len(mu_1)): #<-- this is a loop over phase space dimensions
        p1_k = norm.logpdf(x_data[:,k], loc=mu_1[k], scale=sigma_1[k]) #<-- might have to change how data is indexed
        L_1 += p1_k # <-- numerically dubious

        p2_k = norm.logpdf(x_data[:,k], loc=mu_2[k], scale=sigma_2[k])
        L_2 += p2_k # <-- numerically dubious

    # sum log likelihoods for each component
    L1 = np.sum(L_1)
    L2 = np.sum(L_2)



    L = Q_1*np.exp(L1) + Q_2*np.exp(L2) #<-- ;-; is there a better way to keep this in log space.. 

    return np.log(L) #<-- eventually will want to minimize the negative log likelihood ig. 
    
# %%
ii = 0 # <--- gd1 i think. 
orbit = orbits[ii]
## do the orbit-wise check -- integrate prog orbit and find the pericenter. 
init_displacement = init_displacements[ii]
orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
peri = orbit_obj.pericenter()

mass_index = 1
rvir_index=0
(core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
    simspect.prepare_nbody_data_anycopy(
        orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
        include_photometry=False
    )

coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab) # observed frame
sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) # straight coords

# select unbound. 
unbound = ~CMdict['in_rtid']
unbound = unbound[inMW][trim]

# clip the straight coords before putting into polynomial straightening step
trimmed_sc = simspect.clip_coords(sc, [inMW, trim]) #<-- this applies inMW, trim to the coordinate dictionary
sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction
# %%
### for right now i am blindly putting these into the gmm, eventually will need 
# to convert proper motions -> transverse velocities and _then_ straighten before 
# doing the mixture modeling step. 
