#--------------------------------------------------------------#
#  add noise to data and fit the gaussian mixture model again  #
#                                                              #
#--------------------------------------------------------------#
# %%
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

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
plt.style.use(script_path+'/vedant.mplstyle')
# %config InlineBackend.figure_format='retina'
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap

from tqdm import tqdm

sys.path.append(script_path)
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
import inspect_new_sims as simspect 
import pickle
import gmm as gmm #<-- functions from gmm
import read_mist_models
sys.path.append(repo_path+"/old/")
# %%
# np.log10(12e9) #<-- print the log(age) isochrone i want. 
# 10**10.07918 / 1e9

# np.log10(2700e6)
# 10**( 9.43136) / 1e9 #<-- roughly 2.7 Gyr. 
# %%
### 12 Gyr isochrone
isocmd = read_mist_models.ISOCMD('/n/home02/amphillips/data/MIST_12Gyr_fehn2_ubvraplus/MIST_iso_6ab1420039622.iso.UBVRIplus')
age_ind = isocmd.age_index(10.07918)


### 2700 Myr (dynamical age) isochrone
# isocmd = read_mist_models.ISOCMD('/n/home02/amphillips/data/MIST_2700Myr_fehn2_ubvraplus/MIST_iso_6ab14f0410dce.iso.UBVRIplus')
# age_ind = isocmd.age_index(9.43136)

G_iso = isocmd.isocmds[age_ind]['Gaia_G_EDR3']
BP_iso = isocmd.isocmds[age_ind]['Gaia_BP_EDR3']
RP_iso = isocmd.isocmds[age_ind]['Gaia_RP_EDR3']
mass_iso = isocmd.isocmds[age_ind]['star_mass']
Teff_iso = 10**isocmd.isocmds[age_ind]['log_Teff']
L_iso = 10**isocmd.isocmds[age_ind]['log_L']


plt.hist(isocmd.isocmds[age_ind]['initial_mass'], bins=10)
plt.xlabel(r'$M_{ini}$')
# print info
    # print(isocmd.photo_sys)
    # print(isocmd.ages)
    # print(isocmd.hdr_list)
# fig, ax = plt.subplots()
# ax.scatter(BP_iso - RP_iso, G_iso, c=np.log10(mass_iso))
# ax.invert_yaxis()

### load Jarvis+26 table 7
tt = Table.read('/n/home02/amphillips/data/jarvis26_Table7.fits', format='fits')
# tt.colnames
# %%
### photometry functions to convert from legacy survey z and gaia G
def gaia_g_to_lsst_z(G, bp_rp):
    """
    Convert Gaia G magnitude to LSST ComCam z magnitude.

    Uses the polynomial transformation from RTN-099 (Section 1.3.4):
        z - G = +0.034*(BP-RP)^2 - 0.747*(BP-RP) + 0.416
    RMS residual: 0.014 mag.

    Parameters
    ----------
    G : float or array-like
        Gaia G magnitude(s).
    bp_rp : float or array-like
        Gaia BP-RP colour(s).

    Returns
    -------
    float or ndarray
        LSST z magnitude(s).
    """
    G     = np.asarray(G,     dtype=float)
    bp_rp = np.asarray(bp_rp, dtype=float)

    # out_of_range = (bp_rp < _BPRP_MIN) | (bp_rp > _BPRP_MAX)
    # if np.any(out_of_range):
    #     warnings.warn(
    #         f"{np.sum(out_of_range)} BP-RP value(s) outside the valid range "
    #         f"[{_BPRP_MIN}, {_BPRP_MAX}]. Transformation may be unreliable."
    #     )
 
    return G + 0.034 * bp_rp**2 - 0.747 * bp_rp + 0.416

def lsst_z_to_gaia_g(z, bp_rp):
    """
    Convert LSST ComCam z magnitude to Gaia G magnitude.

    Exact inverse of :func:`gaia_g_to_lsst_z`, using the same polynomial
    from RTN-099 (Section 1.3.4):
        G - z = -0.034*(BP-RP)^2 + 0.747*(BP-RP) - 0.416
    RMS residual: 0.014 mag.

    Parameters
    ----------
    z : float or array-like
        LSST z magnitude(s).
    bp_rp : float or array-like
        Gaia BP-RP colour(s).

    Returns
    -------
    float or ndarray
        Gaia G magnitude(s).
    """
    z     = np.asarray(z,     dtype=float)
    bp_rp = np.asarray(bp_rp, dtype=float)

    # out_of_range = (bp_rp < _BPRP_MIN) | (bp_rp > _BPRP_MAX)
    # if np.any(out_of_range):
    #     warnings.warn(
    #         f"{np.sum(out_of_range)} BP-RP value(s) outside the valid range "
    #         f"[{_BPRP_MIN}, {_BPRP_MAX}]. Transformation may be unreliable."
    #     )

    return z - 0.034 * bp_rp**2 + 0.747 * bp_rp - 0.416

def get_Teff(R, L):
    return (L / (4*np.pi*R**2 * const.sigma_sb))**(1/4)

def desi_RVerr(zmag, feh=-2.0):
    """
    get RV error for desi data model, which 
    depends on z magnitude and metallicity
    """
    log_err = -0.47 + 0.27*(zmag-16) - 0.23*feh
    return 10**log_err

nus = paf.define_photometric_bands()
nu_G_min, nu_G_max, nu_BP_min, nu_BP_max, nu_RP_min, nu_RP_max, nu_z_min, nu_z_max = nus
def get_gaia_photometry(Teff, Radius, distance):
    # frequencies = [nu_G.to(u.Hz).value, nu_BP.to(u.Hz).value, nu_RP.to(u.Hz).value]
    f_min = [nu_G_min.to(u.Hz).value, nu_BP_min.to(u.Hz).value, nu_RP_min.to(u.Hz).value]
    f_max = [nu_G_max.to(u.Hz).value, nu_BP_max.to(u.Hz).value, nu_RP_max.to(u.Hz).value]

    mags = []
    for nu_min, nu_max in zip(f_min, f_max):
        mag = paf.integrated_mag(nu_min, nu_max, Teff, Radius, distance)
        mags.append(mag)
    return mags

### add the BB curve integration step: 
def g_phot(T, R, dpc):
    ### T, R just need astropy units. 
    R_cgs = R.cgs.value
    Teff_cgs = T.cgs.value
    G, BP, RP = [],[],[]
    for Tval, Rval in tqdm(zip(Teff_cgs, R_cgs)):
        Gval, BPval, RPval = get_gaia_photometry(Tval, Rval, dpc) # 10 pc. 
        G.append(Gval)
        BP.append(BPval)
        RP.append(RPval)   
    return np.array(G), np.array(BP), np.array(RP)    

### test whether z CMD looks reasonable. 
z_iso = gaia_g_to_lsst_z(G_iso, BP_iso-RP_iso)
# fig, ax = plt.subplots()
# ax.scatter(BP_iso - RP_iso, z_iso, c=np.log10(mass_iso))
# ax.invert_yaxis()
# %%
### okay... now i guess go about painting mags onto my stars...
#   a ~12 Gyr isochrone is not quite fair, since my simulated stars are not 
#   that old. hmmm..... idk 
grid_info = paf.extended_grid_info(scratch=False) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')

# ordering decided here. 
orbits = ['gd1','aau','pa5','jet','m3','c19']
init_displacements = [
    grid_info.gd1_init_displacement, 
    grid_info.aau_init_displacement,
    grid_info.pa5_init_displacement,
    grid_info.jet_init_displacement,
    grid_info.m3_init_displacement,
    grid_info.c19_init_displacement]
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
copy_options = [0,1,2,3,4]
# copy_options = [4,3,2,1,0]

# keys = ['phi2','pm_phi1','pm_phi2','v_gsr'] #<--- EDITING TO INCLUDE DISTANCE ?
keys = ['d_phi2','v_phi1','v_phi2','v_gsr'] #<-- i guess like why not do this

ii=0
orbit = orbits[ii]

rvir_index=3
rvir = rvirs[rvir_index]

mass_index=1


(core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
    simspect.prepare_nbody_data_anycopy(
        orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
        include_photometry=False
    )

sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) #<-- sc is returned as a DICTIONARY! 

unbound = ~CMdict['in_rtid']
# unbound = unbound[inMW][trim] # don't care about this. 

stellar_masses = lumdict['mass'].to(u.Msun).value
stellar_types = lumdict['type']
nonrem = stellar_types < 10 # <-- i think 10 starts to be WD. 

L, R = lumdict['L'].to(u.Lsun), lumdict['R'].to(u.Rsun)
Teff = get_Teff(R, L).to(u.K)
G, BP, RP = g_phot(Teff, R, dpc=10)

### NEW scheme for trimming the stream just dropped, no 'inMW' necessary now. 
inMW_na = np.ones(len(sc['phi1']), dtype=bool) #<-- i don't actually want to do a "inMW" trim here. 
trim_new = gmm.trim_obstream_percentile(sc) # & ((sc['phi1']>5) & (sc['phi1']<15))

usePhot = nonrem & unbound & trim_new
# %%
BP_RP = BP - RP
z = gaia_g_to_lsst_z(G, BP_RP)
fig, ax = plt.subplots()
## CMD
ax.scatter(BP_RP[usePhot], z[usePhot], c=stellar_masses[usePhot],
           edgecolor='k', lw=.5, s=30,vmin=0.1, vmax=2)

iso_cutoff = -700
ax.scatter(BP_iso[:iso_cutoff]-RP_iso[:iso_cutoff], z_iso[:iso_cutoff],
           c=mass_iso[:iso_cutoff], vmin=0.1, vmax=2, zorder=0)
ax.set_xlabel(r'$G_{\rm BP} - G_{\rm RP}$')
ax.set_ylabel(r'$z$')
ax.invert_yaxis()



########## Kiel diagram
# ax.scatter(Teff.to(u.K)[usePhot], L[usePhot],
#            edgecolor='k', lw=.5, s=30, 
#            c= stellar_masses[usePhot], vmin=0.01, vmax=2)
#         #    c='cornflowerblue')
# ax.scatter(Teff_iso[:iso_cutoff], L_iso[:iso_cutoff], zorder=0, 
#            c = mass_iso[:iso_cutoff], vmin=0.01, vmax=2)
#         #    c='k')
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xticks([1e4, 1e3])
# ax.invert_xaxis()
# # sel = Teff_iso<1e4
# ax.set_xlabel(r'$T_{\rm eff}~[\rm K]$')
# ax.set_ylabel(r'$L~[L_{\odot}]$')



# bins = np.linspace(0, 2, 50)
# ax.hist(mass_iso, bins=bins, color='k', alpha=0.2)
# ax.hist(stellar_masses[nonrem], bins=bins, color='k', histtype='step', lw=3)
# %%
