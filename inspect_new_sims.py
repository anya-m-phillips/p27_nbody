# %%
################ import packages
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

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
# from get_escapers import rms, make_key, dedupe_true_copies

# first import... very beefy...
from sklearn.mixture import GaussianMixture
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty
# %%

### DEFINE STUFF ABOUT THE GRID: 
grid_info = paf.extended_grid_info(scratch=True) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]

# %%
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
# check data_dict outputs...
_, data_dict, CMdict, lumdict, inMW, trim = prepare_nbody_data(
    path = grid_info.m3_lm_paths[0]+"0/",
    i=grid_info.m3_age,
    apo=grid_info.m3_apo,
    init_displacement=grid_info.m3_init_displacement
) 
# %%
# example astropy coorddinate xforms
pos, vel = CMdict['pos'], CMdict['vel']
coords_ICRS = paf.galcen_to_ICRS(pos, vel)
coords_ICRS = reflex_correct(coords_ICRS)

# fig, ax = plt.subplots()
# ax.scatter(coords_ICRS.ra, coords_ICRS.dec)

# ### quick check to make sure i set up my orbits right. 
# m3_catalog = Table.read('/n/home02/amphillips/stream_catalogs/M3.fits')
# ax.scatter(m3_catalog['ra'], m3_catalog['dec'])
# m3_catalog

endpoints = SkyCoord(
    ra=[186.45,197.20]*u.degree,
    dec=[19.06,27.76]*u.degree
)

# choose the origin based on the ICRS coordinate of the progenitor 
# at present day. might be smart to read these in from data/FINAL_ics_nolmc.csv rather than hard coding here.
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')
row = prog_tab[prog_tab['name']=='M3'] 
prog_pos = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.kpc
prog_vel = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.km/u.s
prog_coord = paf.galcen_to_ICRS(prog_pos, prog_vel) #<-- debug and print htis out to see if its between the endpoints...
origin = SkyCoord(
    ra=prog_coord.ra, dec=prog_coord.dec #<-- choose an origin... in ra
)
M3Yang23 = gc.GreatCircleICRSFrame.from_endpoints(
    endpoints[0], endpoints[1], origin=origin
)
M3Yang23
pole = M3Yang23.pole
ra0 = M3Yang23.origin.ra
M3Yang23_tryagain = gc.GreatCircleICRSFrame.from_pole_ra0(
    pole=pole, ra0=ra0, origin_disambiguate=origin
)

# %%
origin
# %%
sc = coords_ICRS.tranform_to(M3Yang23_tryagain) #<--- why does the .from_endpoints variant not work but the from_pole_ra0 does?

# sc = coords_ICRS.transform_to(C19Ibata24)
# %%

# plt.scatter(coords_ICRS.ra[inMW][trim], coords_ICRS.pm_dec[inMW][trim])
fig, ax = plt.subplots(figsize=[8,3])
ax.scatter(sc.phi1[inMW][trim], sc.phi2[inMW][trim], c='k', s=.1)
# # ax.set_ylim(-3, 6)
# ax.set_xlim(-15,15)


 
# %%
def streamframe_coords_observed(orbit, data_dict):
    """
    orbit options are:
    gd1, pa5, jet, aau, m3, c19, circ
    
    (input as strings), will put into a great-circle coordinate frame 
    based on ICRS ra/dec. this will allow accounting for projection 
    effects. also will correct radial velocity + proper motion for solar
    reflex motion, using a distance that for now is assumed to 
    be known exactly. For real streams this would be based on 
    a dynamically-fit distance track and hopefully would also be 
    pretty well-constrained....

    data_dict should be one of the three choicdes of subdictionary (CoM, luminous, companions)
    so that it has position and velocity info. those should have astropy units. 
    """
    pos, vel = data_dict['pos'], data_dict['vel']
    coords_ICRS = paf.galcen_to_ICRS(pos, vel) 
    coords_ICRS = reflex_correct(coords_ICRS) #<-- correct for solar reflex motion

    matrix_orbits = ['aau','jet'] #<-- these are given as rotation matrices in the literature. will not use gala. 
    
    if orbit not in matrix_orbits: # if orbit=='gd1' or orbit=='pa5' or orbit=='c19':
        coords_stream = {}
        if orbit=='gd1':
            sc=coords_ICRS.transform_to(gc.GD1Koposov10)
        if orbit=='pa5':
            sc = coords_ICRS.transform_to(gc.Pal5PriceWhelan18)

        if orbit=='c19':
            # as defined in ibata+24 and explained in mohammed 26
            alpha_0 = 354.356*u.degree #<-- sets the phi1 zero point

            pole = SkyCoord( #<-- defines the pole
                ra= 81.45*u.degree,
                dec = -6.346*u.degree
            )
            prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')
            row = prog_tab[prog_tab['name']=='c19'] 
            prog_pos = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.kpc
            prog_vel = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.km/u.s
            prog_coord = paf.galcen_to_ICRS(prog_pos, prog_vel) #<-- debug and print htis out to see if its between the endpoints...
            origin = SkyCoord(
                ra=prog_coord.ra, dec=prog_coord.dec #<-- choose an origin... in ra
            )

            C19Ibata24 = gc.GreatCircleICRSFrame.from_pole_ra0(
                pole=pole, ra0=alpha_0, origin_disambiguate=origin #[coordinate that the origin should be closest to?]
            )
            sc = coords_ICRS.transform_to(C19Ibata24)

        if orbit=='m3':
            # from Yang+2023, section 4.4. Endpoints are given but origin is not. 
            endpoints = SkyCoord(
                ra=[186.45,197.20]*u.degree,
                dec=[19.06,27.76]*u.degree
            )

            # choose the origin based on the ICRS coordinate of the progenitor 
            # at present day. might be smart to read these in from data/FINAL_ics_nolmc.csv rather than hard coding here.
            prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')
            row = prog_tab[prog_tab['name']=='M3'] 
            prog_pos = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.kpc
            prog_vel = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.km/u.s
            prog_coord = paf.galcen_to_ICRS(prog_pos, prog_vel) #<-- debug and print htis out to see if its between the endpoints...
            origin = SkyCoord(
                ra=prog_coord.ra, dec=prog_coord.dec #<-- choose an origin... in ra
            )
            M3Yang23 = gc.GreatCircleICRSFrame.from_endpoints(
                endpoints[0], endpoints[1], origin=origin
            )
            sc = coords_ICRS.tranform_to(M3Yang23)

        coords_stream['phi1'] = sc.phi1.to(u.degree).value
        coords_stream['phi2'] = sc.phi2.to(u.degree).value
        coords_stream['pm_phi1'] = sc.pm_phi1_cosphi2.to(u.mas/u.yr).value
        coords_stream['pm_phi2'] = sc.pm_phi2.to(u.mas/u.yr).value
        coords_stream['distance'] = sc.distance.to(u.kpc).value
        coords_stream['v_gsr'] = sc.radial_velocity.to(u.km/u.s).value

    if orbit in matrix_orbits:
        if orbit=='aau':
            # rotation matrix from shipp+2019, Li+2021
            M = np.array([
                [0.83697865, 0.29481904, -0.4610298], 
                [0.51616778, -0.70514011, 0.4861566],
                [0.18176238, 0.64487142, 0.74236331]
            ])
        if orbit=='jet':
            # Do+26
            M = np.array([
                [-0.69798645, 0.61127501, -0.37303856], 
                [-0.62615889, -0.26819784, 0.73211677],
                [0.34747655, 0.74458900, 0.56995374]
            ])
        ### rotate ICRS coords by M. see paf.xform_to_koposov_coords(ra, dec) for how to do this. 
        # create the coords_stream dictionary

        pass


    return coords_stream #<-- a dictionary, like Jake's 

def desi_RVerr(zmag, feh=-2.0):
    """
    get RV error for desi data model, which 
    depends on z magnitude and metallicity
    """
    log_err = -0.47 + 0.27*(zmag-16) - 0.23*feh
    return 10**log_err

### define Via RVerr(mag, metallicity)

def add_noise(icrs_coords, survey='DESI'):
    """
    will want to have "DESI" and "Via" options 
    for the survey velocity errors. 
    For Desi these will come from a Koposov paper, 
    for Via they will come from viamock. 
    """
    return
# %%
# if __name__=='__main__': #### stuff i won't want to run when i import functions to other scripts below. 

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


# %%

# %%
# and do we have streamframes for everyone???
# - [x] GD1 (koposov; in gala)
# - [x] pa5 (price-whelan; in gala)
# - [x ] c19 (ibata+24; not in gala; see Nasser's paper)
# - [ ] M3 (??) lmao, does exist: Such a rotation can be easily done by some tools, e.g., Gala 7 (Price-Whelan 2017). -- Yang+23
# - [x ] AAU: https://iopscience.iop.org/article/10.3847/1538-4357/abeb18#apjabeb18app1 (Li+21; not in gala)
# - [x ] jet: Do+26 sec 3
# %%