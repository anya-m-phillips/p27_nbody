#----------------------------------------------------------#
#   scratch aside on formatting of bonaca/price-whelan26   #
#   want to check how many gaia members these things have  #
#   and to call it a lower limit for what via/desi should  #
#   see. also does the jarvis catalog have more/fewer      #
#   members than gaia for gd1 ? maybe3 more because of     #
#   spectroscopic confirmation? anyway let's do it         #
#----------------------------------------------------------#
# %%
#--------------------------------------------------------------#
#  save dictionaries that are [final data dict] plus observed  #
#   coordinates, both with and without noise/binary motions    #
#--------------------------------------------------------------#
# %%
import sys
soft_path = '/n/home02/amphillips/software/'
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

from scipy.stats import binned_statistic, norm, multivariate_normal
from scipy.optimize import curve_fit, minimize, Bounds
# from scipy.ndimage import gaussian_filter1d
# from scipy.interpolate import CubicSpline
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
plt.style.use(script_path+'/vedant.mplstyle')
# %config InlineBackend.figure_format='retina'
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap

from tqdm import tqdm
import pickle

sys.path.append(script_path)
sys.path.append(repo_path)
sys.path.append(soft_path+'viamock')
from streamframe import StreamFrame
import artpop
import viamock
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty

import PETAR_ANALYSIS_FUNCTIONS as paf
import inspect_new_sims as simspect 
import noise as noise #<-- for like isochrone stuff. 

import argparse
# %%
dp = repo_path+'/data/bpw25_catalogs/'
# sns = ['ATLAS-Aliqa Uma.fits', 'C-19.fits','GD-1.fits','Jet.fits','Palomar 5.fits']
sns = ['GD-1.fits','Palomar 5.fits','C-19.fits', 'ATLAS-Aliqa Uma.fits', 'Jet.fits']

fig, ax = plt.subplots()
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
ccc = cc[1:]

for ii, cat in enumerate(tqdm(sns)):
    t = Table.read(dp+cat, format='fits')
    n = len(t)
    ax.scatter(t['bp_rp'], t['phot_g_mean_mag'], c=ccc[ii],
               label=cat[:-5]+"; n=%i"%n)

ax.invert_yaxis()
ax.set_xlabel(r'$G_{B_p}-G_{R_p}$')
ax.set_ylabel(r'$G$')
ax.legend(loc='upper left', bbox_to_anchor=[1,1])

# %%
t_jarvis = Table.read(repo_path+"/data/jarvis26_Table7.fits", format='fits')
print("N* from DESI catalog:", len(t_jarvis))