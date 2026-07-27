# %%
################ import packages
import sys
script_path = "/n/home02/amphillips/p27_nbody/scripts" # for cannon

import petar
import numpy as np

from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

import astropy.coordinates as coord
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
# from get_escapers import rms, make_key, dedupe_true_copies

# first import... very beefy...
from sklearn.mixture import GaussianMixture
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty

sys.path.append("/n/home02/amphillips/p26/scripts/")
from gen_king_clust import rotation_matrix_from_vectors, random_unit_vector

# %%

### DEFINE STUFF ABOUT THE GRID: 
grid_info = paf.extended_grid_info(scratch=True) 

# main helper function: 
def prepare_nbody_data(path = grid_info.gd1_lm_paths[0]+"0/", 
                       include_photometry=False, 
                       i=grid_info.gd1_age,
                       apo = grid_info.gd1_apo,
                       init_displacement = grid_info.gd1_init_displacement,
                       ):
    """
    this is expecting N-body data that outputs every 10 Myr. 
    """
    core = paf.load_core(path)
    file_index = int(i/10)

    data_dict = paf.intrinsic_stream_data_v3(
        path, i, core, apo, init_displacement,
        use_core=False,
        binary_treatments=['CoM','companions','luminous']
    )

    CMdict = data_dict['CoM']
    lumdict = data_dict['luminous']
    inMW, trim = CMdict['inMW'], CMdict['trim']


    if include_photometry==True:
        # do the calculations -- ORDERED AS ALL_PARTICLES
        primary_Ls = lumdict["L"][inMW][trim] #* u.Lsun
        primary_Rs = lumdict['R'][inMW][trim] #* u.Rsun
        primary_Teffs = (primary_Ls / (4*np.pi*primary_Rs**2 * const.sigma_sb))**(1/4)

        R_cgs = primary_Rs.cgs.value
        Teff_cgs = primary_Teffs.cgs.value
        G, BP, RP, z = [],[],[],[]
        for Tval, Rval in tqdm(zip(Teff_cgs, R_cgs)):
            Gval, BPval, RPval = paf.get_gaia_photometry(Tval, Rval, 10) # 10 pc. 
            G.append(Gval)
            BP.append(BPval)
            RP.append(RPval)

            zval = paf.get_z_photometry(Tval, Rval, 10) # 10pc
            z.append(zval)
        G = np.array(G)
        BP = np.array(BP)
        RP = np.array(RP)
        z = np.array(z)

        return core, data_dict, CMdict, lumdict, inMW, trim, G, BP, RP, z
    
    else:
        return core, data_dict, CMdict, lumdict, inMW, trim

# %%
orbits = ['circ','gd1','aau','pa5','jet','m3','c19']

dicts = []
for ii, orbit in tqdm(enumerate(orbits)):
    path, apo, age, init_displacement = grid_info.retrieve_sim_info(
        orbit=orbit, stellar_pop='lm', rvir_index=0, copy=0
    )

    core, data_dict, CMdict, lumdict, inMW, trim = prepare_nbody_data(
        path = path,
        include_photometry=False,
        i=age if orbit != "circ" else 30000,
        init_displacement = init_displacement
    )
    dicts.append(data_dict)
# %%
def rotation_matrix(a,b,c):
    """
    rotate angles a,b,c about x,y,z axes
    https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
    """
    Rx = np.array([
        [1,0,0],
        [0, np.cos(a), -np.sin(a)],
        [0, np.sin(a), np.cos(a)]
    ])
    Ry = np.array([
        [np.cos(b), 0, np.sin(b)],
        [0,1,0],
        [-np.sin(b), 0, np.cos(b)]
    ])
    Rz = np.array([
        [np.cos(c), -np.sin(c), 0],
        [np.sin(c), np.cos(c), 0],
        [0,0,1]
    ])

    # R = Rx @ Ry @ Rz
    R = Rz @ Ry @ Rx # <-- do x rotation first, _then_ z rotation.. mostly for movie purposes. 
    return R
# %%
lm_colors, hm_colors, simcolors = paf.define_simcolors()
# cc = simcolors[3:]
reordered_colors = lm_colors+hm_colors
cc = reordered_colors[1:]
### generate a rotation matrix here, so that I can view the streams from another angle!

# theta_list = np.arange(0, 361, 1)*u.degree.to(u.radian)
# it = 0
# for theta in tqdm(theta_list):
theta = 10*u.degree.to(u.radian)
phi =  10*u.degree.to(u.radian)
R = rotation_matrix(a=phi, b=0, c=theta)

fig, ax = plt.subplots()
for ii, data_dict in enumerate(dicts):

    x, y, z = data_dict['CoM']['pos'].T.to(u.kpc)
    pos = data_dict['CoM']['pos'].to(u.kpc)

    rotated_pos = (R@pos.T).T
    ax.scatter(rotated_pos[:,0], rotated_pos[:,2], 
            c=cc[ii], 
            rasterized=True, 
            s=1)
    
pos_earth = np.array([8,0,0])
rotated_pos_earth = np.dot(R, pos_earth)
ax.scatter(rotated_pos_earth[0], rotated_pos_earth[2], c='k', marker='*', s=100)
ax.set_xlim(-40,40)
ax.set_ylim(-40,40)
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_xticks([])
ax.set_xticklabels([])

    # dir='/n/netscratch/conroy_lab/Lab/amphillips/movies/grid_rotation/'
    # filename = f"frame_{it:05d}.png"
    # plt.savefig(dir+filename)
    # plt.close()
    # it+=1
# %%
# and do we have streamframes for everyone???
# - [x] GD1 (koposov; in gala)
# - [x] pa5 (price-whelan; in gala)
# - [ ] c19 (ibata+24; not in gala; see Nasser's paper)
# - [ ] M3 (??) lmao, does exist: Such a rotation can be easily done by some tools, e.g., Gala 7 (Price-Whelan 2017). -- Yang+23
# - [ ] AAU: https://iopscience.iop.org/article/10.3847/1538-4357/abeb18#apjabeb18app1 (Li+21; not in gala)
# - [ ] jet: Do+26 sec 3
# %%