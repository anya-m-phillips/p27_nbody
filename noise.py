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
m0_iso = isocmd.isocmds[age_ind]['initial_mass']
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

### ------------------------------------------------------------------------ ###
###  isochrone interpolation: ZAMS mass --> synthetic Gaia photometry        ###
###  (replaces the blackbody + top-hat scheme; swap the isochrone above to   ###
###   alter the stellar population artificially)                             ###
### ------------------------------------------------------------------------ ###
### PRECISION NOTES, since the evolved sequence is nearly degenerate in mass:
###   - everything stays float64 and is interpolated against the RAW
###     `initial_mass` column. no rounding, no binning, no float32 anywhere.
###     for the 12 Gyr isochrone the whole post-turnoff sequence (RGB + CHeB +
###     AGB) spans only 0.0168 Msun, and the tightest node spacing is
###     1.6e-11 Msun -- that is ~1e5 x float64 eps at 0.8 Msun, so linear
###     interpolation still resolves it, but there is no headroom to waste.
###   - LINEAR interpolation, deliberately, not a spline: on the near-vertical
###     RGB/AGB segments dM/dmag is ~1e-9, so a cubic overshoots enormously
###     between nodes. linear is monotonic and cannot invent points off the
###     isochrone.
###   - using `initial_mass` (NOT `star_mass`) is what makes it invertible:
###     initial_mass is strictly increasing along the isochrone, while
###     star_mass turns over once winds kick in. this also matches
###     `lumdict['m0_zams']` -- the genuine birth mass -- and not
###     `m0_effective` (BSE's fitting coordinate) or `mass` (current mass).
###   - interpolating in mass means an evolutionary phase is sampled in
###     proportion to the initial-mass interval it occupies, which is exactly
###     its lifetime x IMF weight. so the giants come out rare automatically;
###     nothing extra is needed to get the relative numbers right.
ISO_BANDS = ('Gaia_G_EDR3', 'Gaia_BP_EDR3', 'Gaia_RP_EDR3')

def build_isochrone_table(iso, bands=ISO_BANDS, max_phase=5, mass_col='initial_mass'):
    """
    Build a strictly-increasing initial-mass grid and the matching magnitudes,
    ready for `np.interp`.

    Parameters
    ----------
    iso : structured array
        one age's rows, i.e. `isocmd.isocmds[age_ind]`.
    bands : sequence of str
        isochrone magnitude columns to carry along.
    max_phase : int or None
        drop rows with MIST `phase` above this. default 5 keeps MS(0), RGB(2),
        CHeB(3), EAGB(4) and TPAGB(5) and drops post-AGB(6), which is the
        proto-WD tail -- those stars are remnants in the simulation and are
        already cut by `nonrem`. set None to keep everything.
    mass_col : str
        must be a monotonically increasing column; see the precision notes.

    Returns
    -------
    m0_grid : (M,) float64, strictly increasing
    table : dict band -> (M,) float64
    """
    m0 = np.asarray(iso[mass_col], dtype=np.float64)

    keep = np.ones(len(m0), dtype=bool)
    if max_phase is not None:
        keep &= np.asarray(iso['phase'], dtype=np.float64) <= max_phase
    idx = np.flatnonzero(keep)

    # `np.interp` REQUIRES increasing xp and does not check -- it silently
    # returns garbage otherwise (same trap as straighten_stream_orbit_interp).
    # so sort, then drop any non-increasing node rather than trusting the file.
    idx = idx[np.argsort(m0[idx], kind='stable')]
    strictly_increasing = np.ones(len(idx), dtype=bool)
    strictly_increasing[1:] = np.diff(m0[idx]) > 0.0
    n_dropped = np.count_nonzero(~strictly_increasing)
    if n_dropped:
        print(f'build_isochrone_table: dropped {n_dropped} non-increasing '
              f'{mass_col} node(s)')
    idx = idx[strictly_increasing]

    m0_grid = m0[idx]
    table = {b: np.asarray(iso[b], dtype=np.float64)[idx] for b in bands}
    return m0_grid, table


def isochrone_photometry(m0_query, m0_grid, table):
    """
    Interpolate isochrone magnitudes at the ZAMS masses `m0_query`.

    Magnitudes are ABSOLUTE (MIST isochrones are), i.e. the same convention as
    `g_phot(..., dpc=10)`, so downstream code needs no change.

    Stars outside the isochrone's initial-mass range get NaN and are flagged
    `False` in `on_iso` -- the isochrone is *old*, so its upper limit is the
    turnoff-ish 0.8 Msun and every more massive N-body star has to be thrown
    out before interpolating. NaN rather than a clamped edge value, so a
    forgotten mask shows up as a hole in the CMD instead of a fake pile-up at
    the tip of the AGB.

    Returns
    -------
    phot : dict band -> (N,) float64, NaN off the isochrone
    on_iso : (N,) bool, True where the star was interpolated
    """
    m0_query = np.asarray(m0_query, dtype=np.float64)
    on_iso = (m0_query >= m0_grid[0]) & (m0_query <= m0_grid[-1])

    phot = {}
    for band, y in table.items():
        vals = np.full(m0_query.shape, np.nan, dtype=np.float64)
        vals[on_iso] = np.interp(m0_query[on_iso], m0_grid, y)
        phot[band] = vals
    return phot, on_iso


def gaia_from_isochrone(m0_query, iso, bands=ISO_BANDS, max_phase=5):
    """
    Convenience wrapper: ZAMS masses -> (G, BP, RP, on_iso), absolute mags.
    """
    m0_grid, table = build_isochrone_table(iso, bands=bands, max_phase=max_phase)
    phot, on_iso = isochrone_photometry(m0_query, m0_grid, table)
    return phot[bands[0]], phot[bands[1]], phot[bands[2]], on_iso


### test whether z CMD looks reasonable.
z_iso = gaia_g_to_lsst_z(G_iso, BP_iso-RP_iso)
# fig, ax = plt.subplots()
# ax.scatter(BP_iso - RP_iso, z_iso, c=np.log10(mass_iso))
# ax.invert_yaxis()
# %%


###### MAIN PROGRAM BELOW 
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

keys = ['d_phi2','v_phi1','v_phi2','v_gsr'] 

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

# get distances also:
coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab)
distances = coords_obs.distance
sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) #<-- sc is returned as a DICTIONARY! 

unbound = ~CMdict['in_rtid']
# unbound = unbound[inMW][trim] # don't care about this. 

stellar_masses = lumdict['mass'].to(u.Msun).value
stellar_types = lumdict['type']
nonrem = stellar_types < 10 # <-- i think 10 starts to be WD. 

### the genuine birth mass -- NOT m0_effective (BSE's fitting coordinate) and
### NOT mass (current mass). see the star.mass0 section of the readme.
m0s = np.asarray(lumdict['m0_zams'], dtype=np.float64)

L, R = lumdict['L'].to(u.Lsun), lumdict['R'].to(u.Rsun)
Teff = get_Teff(R, L).to(u.K)

USE_ISOCHRONE = True   # False falls back to the blackbody + top-hat scheme
ISO_MAX_PHASE = 5      # drop post-AGB(6); see build_isochrone_table

if USE_ISOCHRONE:
    ### paint photometry by interpolating the isochrone in INITIAL mass.
    m0_grid, iso_table = build_isochrone_table(isocmd.isocmds[age_ind],
                                               max_phase=ISO_MAX_PHASE)
    min_m0, max_m0 = m0_grid[0], m0_grid[-1]

    iso_phot, on_iso = isochrone_photometry(m0s, m0_grid, iso_table)
    G  = iso_phot['Gaia_G_EDR3']
    BP = iso_phot['Gaia_BP_EDR3']
    RP = iso_phot['Gaia_RP_EDR3']

    ### the isochrone is old, so it stops at the turnoff: everything more
    ### massive than max_m0 has to be thrown out. (a handful of sim stars can
    ### also sit below the isochrone's low-mass end, so this is two-sided.)
    alive = on_iso
    n_hi = int(np.count_nonzero(m0s > max_m0))
    n_lo = int(np.count_nonzero(m0s < min_m0))
    print(f'isochrone m0 range [{min_m0:.10f}, {max_m0:.10f}] Msun, '
          f'{len(m0_grid)} nodes (min spacing {np.diff(m0_grid).min():.3e})')
    print(f'painted {alive.sum()}/{len(m0s)} stars; dropped {n_hi} above and '
          f'{n_lo} below the isochrone')
else:
    max_m0 = m0_iso.max()
    min_m0 = m0_iso.min()
    alive = (m0s <= max_m0) & (m0s >= min_m0)
    G, BP, RP = g_phot(Teff, R, dpc=10)

BP_RP = BP - RP

### NEW scheme for trimming the stream just dropped, no 'inMW' necessary now. 
inMW_na = np.ones(len(sc['phi1']), dtype=bool) #<-- i don't actually want to do a "inMW" trim here. 
trim_new = gmm.trim_obstream_percentile(sc) # & ((sc['phi1']>5) & (sc['phi1']<15))

#### apply the trim to all of the coordinates:
trimmed_sc = simspect.clip_coords(sc, [inMW_na, trim_new]) #<-- this applies inMW, trim to the coordinate dictionary
sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction

usePos = nonrem & unbound & alive #<-- sc_straighter already has trim_new applied. 
usePhot = nonrem & unbound & trim_new & alive




z = gaia_g_to_lsst_z(G, BP_RP)

mz = paf.m_from_M(z, dist=distances)# 10*u.kpc)
mG = paf.m_from_M(G, dist=distances)# 10*u.kpc)

# mz_iso = paf.m_from_M(z_iso, dist=10*u.kpc)

### CHECK CMD
fig, ax = plt.subplots()
ax.scatter(BP_RP[usePhot], mG[usePhot], c=np.log10(m0s[usePhot]),
           edgecolor='k', lw=.5, s=30,vmin=-1, vmax=np.log10(max_m0))

iso_cutoff = -700
# ax.scatter(BP_iso[:iso_cutoff]-RP_iso[:iso_cutoff], mz_iso[:iso_cutoff],
#            c=np.log10(m0_iso[:iso_cutoff]), vmin=-1, vmax=np.log10(max_m0), zorder=0)
ax.set_xlabel(r'$G_{\rm BP} - G_{\rm RP}$')
ax.set_ylabel(r'$G$')
ax.invert_yaxis()
ax.set_ylim(bottom=21o)
ax.set_xlim(right=1.5)


# %%