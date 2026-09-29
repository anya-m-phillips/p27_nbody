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
sys.path.append("/n/home02/amphillips/software/viamock")
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
import inspect_new_sims as simspect 
import pickle
# import gmm as gmm #<-- functions from gmm
import read_mist_models

### stuff to add errors. 
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty
import viamock
import artpop

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
    TODO: add Teff dependence ??? 
    """
    log_err = -0.47 + 0.27*(zmag-16) - 0.23*feh
    return 10**log_err

############# decided i don't need a wrapper function for via errors. 
def via_RVerr(G, log_Teff, feh=-2.0, exptime_s = 3600, nexp=1, **kwargs):
    """
    maybe a little silly to wrap tbh. 
    kwargs might be seeing, moon %, airmass, idk. 
    """
    erv, _, _ = viamock.get_viaspec_errors(
        G = G, 
        feh = feh,
        logteff = log_Teff, 
        exptime_s = exptime_s,
        nexp=nexp,
        **kwargs
    )
    return erv


########## also deciding against blackbody integration... 
#    but keeping the functions perhaps to illustrate the difference in stellar populations. 
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
###########################


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
ISO_BANDS = ('Gaia_G_EDR3', 'Gaia_BP_EDR3', 'Gaia_RP_EDR3', 'log_Teff') #<-- lol i ... think this will work...

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


def load_gaia_catalog(orbit):
    catalog = Table.read(repo_path+'/data/bpw25_catalogs/%s.fits'%orbit)
    return catalog


### ------------------------------------------------------------------------ ###
###  CMD-space lookup: (G, BP-RP) --> nearest point on the isochrone         ###
###  (for real catalog photometry, which never sits exactly on the track)    ###
### ------------------------------------------------------------------------ ###
### NOTES:
###   - the track is ordered by EEP, NOT initial_mass. consecutive EEPs are
###     consecutive points along the sequence, and on the TP-AGB many nodes
###     share an identical initial_mass (build_isochrone_table drops those as
###     non-increasing, which is right for mass interpolation but would punch
###     holes in the CMD track).
###   - the artpop isochrone labels the post-AGB / WD-cooling tail as phase 5,
###     not 6, so `max_phase` can't remove it. cut on EEP instead: MIST primary
###     EEP 1409 is the start of post-AGB (Dotter 2016, table 2). without the
###     cut, faint blue catalog stars snap onto a WD track at logTeff ~ 5.
###   - "nearest" is measured in units of (sigma_color, sigma_G), and it is the
###     nearest point on the piecewise-linear track (projection onto each
###     segment), not the nearest node -- node spacing is very uneven.
POST_AGB_EEP = 1409

def isochrone_cmd_track(iso, max_eep=POST_AGB_EEP, G_col='Gaia_G_EDR3',
                        BP_col='Gaia_BP_EDR3', RP_col='Gaia_RP_EDR3',
                        carry=('log_Teff',)):
    """
    The isochrone as an ordered polyline in (BP-RP, G), for nearest_isochrone_point.

    Parameters
    ----------
    iso : structured array / table with an 'EEP' column (artpop or MIST).
    max_eep : int or None
        keep EEP < max_eep. default drops post-AGB + WD cooling. None keeps all.
    carry : sequence of str
        extra isochrone columns to interpolate at the matched point.

    Returns
    -------
    track : dict with 'EEP', 'G', 'BP_RP' and each of `carry`, (M,) float64,
        ordered along the sequence. magnitudes are absolute.
    """
    eep = np.asarray(iso['EEP'], dtype=np.float64)
    keep = np.ones(len(eep), dtype=bool) if max_eep is None else eep < max_eep
    idx = np.flatnonzero(keep)
    idx = idx[np.argsort(eep[idx], kind='stable')]

    track = {
        'EEP':   eep[idx],
        'G':     np.asarray(iso[G_col], dtype=np.float64)[idx],
        'BP_RP': (np.asarray(iso[BP_col], dtype=np.float64)
                  - np.asarray(iso[RP_col], dtype=np.float64))[idx],
    }
    for col in carry:
        track[col] = np.asarray(iso[col], dtype=np.float64)[idx]
    return track


def nearest_isochrone_point(G, BP_RP, track, sigma_G=0.2, sigma_color=0.05,
                            chunk=1024):
    """
    For each (G, BP-RP), find the closest point on the isochrone track and
    interpolate every track column there.

    distance is sqrt((dcolor/sigma_color)^2 + (dG/sigma_G)^2), minimized over
    each straight segment of the track. the defaults weight colour 4x more than
    G: the absolute G of a painted catalog star carries the N-body distance
    (~0.2 mag of spread along the stream), while BP-RP is what actually
    constrains Teff. the scales only matter relative to each other.

    Parameters
    ----------
    G, BP_RP : (N,) array-like, ABSOLUTE G and colour. NaNs pass through.
    track : dict from isochrone_cmd_track.
    sigma_G, sigma_color : float, scale of each CMD axis in the metric.
    chunk : int, query stars per vectorized block (memory ~ chunk x M).

    Returns
    -------
    nearest : dict with every track key ('EEP', 'G', 'BP_RP', 'log_Teff', ...)
        evaluated at the matched point, plus 'dist' (in the scaled units above;
        e.g. dist > 3 ~ "not on this isochrone": blue stragglers, BHB stars
        bluer than the model HB, contaminants). all (N,), NaN where the input
        was NaN. an interpolated 'EEP' tells you which phase it landed on.
    """
    G = np.asarray(G, dtype=np.float64)
    c = np.asarray(BP_RP, dtype=np.float64)

    # scaled coordinates, so one unit is "one sigma" along either axis
    x_iso = track['BP_RP'] / sigma_color
    y_iso = track['G'] / sigma_G
    x0, y0 = x_iso[:-1], y_iso[:-1]
    dx, dy = np.diff(x_iso), np.diff(y_iso)
    L2 = dx**2 + dy**2
    degenerate = L2 == 0.0
    L2_safe = np.where(degenerate, 1.0, L2)

    nearest = {k: np.full(G.shape, np.nan) for k in list(track) + ['dist']}
    good = np.flatnonzero(np.isfinite(G) & np.isfinite(c))

    for start in range(0, len(good), chunk):
        q = good[start:start+chunk]
        px = (c[q] / sigma_color)[:, None]
        py = (G[q] / sigma_G)[:, None]

        # projection parameter along each segment, clipped to the segment
        t = ((px - x0)*dx + (py - y0)*dy) / L2_safe
        t = np.clip(np.where(degenerate, 0.0, t), 0.0, 1.0)
        d2 = (x0 + t*dx - px)**2 + (y0 + t*dy - py)**2

        j = np.argmin(d2, axis=1)
        rows = np.arange(len(q))
        tj = t[rows, j]
        for k, v in track.items():
            nearest[k][q] = v[j] + tj*(v[j+1] - v[j])
        nearest['dist'][q] = np.sqrt(d2[rows, j])

    return nearest


def Teff_from_gaia_isochrone(track, G, BP_RP, **kwargs):
    """
    log Teff of stars given their (absolute) G, BP-RP, read off the nearest
    point on the isochrone. kwargs go to nearest_isochrone_point.
    """
    return nearest_isochrone_point(G, BP_RP, track, **kwargs)['log_Teff']


def trim_obstream_percentile(sc, p=[1,99], 
                            trim_keys=['phi1','d_phi2','v_phi1','v_phi2','v_gsr','distance']):
    criteria = []
    for key in trim_keys:
        key_low, key_high = np.percentile(sc[key], q=p)
        key_crit = (sc[key]<=key_high) & (sc[key]>=key_low)
        criteria.append(key_crit)

    trim_criteria = np.logical_and.reduce(criteria)
    return trim_criteria
# %%
# SCRATCH::::

isocmd = artpop.fetch_mist_iso_cmd(
    log_age=np.log10(12e9),
    feh=-2.0,
    phot_system='UBVRIplus',
    #v_over_vcrit=0.0 #<-- idk
)
m0_grid, iso_table = build_isochrone_table(isocmd) #<-- will do gaia bands + teff automatically, have max_phase=5 (remove post agb evolution)
min_m0, max_m0 = m0_grid[0], m0_grid[-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')
# %%
orbit='gd1'
rvir_index=0
masses=['lm','hm']
mass_index=1
tidal_boundary=2.0
copy_options=[0,1,2,3,4]
(core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
    simspect.prepare_nbody_data_anycopy(
        orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
        include_photometry=False, N_rtid_boundary = tidal_boundary, #<--- not sure what i'm going to use for tthis: 
        verbose=True
    )
# %%
coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab)
sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) #<-- sc is returned as a DICTIONARY! 
trim_new = trim_obstream_percentile(sc=sc)

distances = coords_obs.distance[trim_new]

m0s = np.asarray(lumdict['m0_zams'], dtype=np.float64)
iso_phot, on_iso = isochrone_photometry(m0s, m0_grid, iso_table)
G  = iso_phot['Gaia_G_EDR3'][trim_new]
BP = iso_phot['Gaia_BP_EDR3'][trim_new] #<-- i think i never need these since they are getting replaced completely... 
RP = iso_phot['Gaia_RP_EDR3'][trim_new] #<-- i think i never need these since they are getting replaced completely... 

mG = paf.m_from_M(G, dist=distances)

iso_sort = np.argsort(mG)

mG_sorted = mG[iso_sort]
BP_sorted, RP_sorted, G_sorted = BP[iso_sort], RP[iso_sort], G[iso_sort]

distances_sorted = distances[iso_sort]



t = Table.read(repo_path+"/data/bpw25_catalogs/%s.fits"%orbit, format='fits')
mG_cat = t['phot_g_mean_mag']
cat_sort = np.argsort(mG_cat)
N_jarvis = 679
N = N_jarvis


t_use = t[cat_sort][:N]
G_cat = paf.M_from_m(t_use['phot_g_mean_mag'], distances_sorted[:N])
BP_RP_cat = t_use['bp_rp']

# get nearest analog in iso_table??
iso_track = isochrone_cmd_track(isocmd) #<-- EEP-ordered, post-AGB/WD tail dropped
nearest = nearest_isochrone_point(G_cat, BP_RP_cat, iso_track) #<-- default sigma_G=0.2, sigma_color=0.05
logTeff_cat = nearest['log_Teff']
print('match dist (sigma units) 50/90/99:', np.nanpercentile(nearest['dist'], [50, 90, 99]))

fig, ax = plt.subplots()
ax.plot(iso_track['BP_RP'], iso_track['G'], 'k-', lw=0.8)
ax.plot(np.vstack([BP_RP_cat, nearest['BP_RP']]), np.vstack([G_cat, nearest['G']]), 'r-', lw=0.3)
im = ax.scatter(BP_RP_cat, G_cat, c=logTeff_cat, s=5)
plt.colorbar(im, ax=ax, label=r'$\log T_{\rm eff}$')
ax.set_xlim(0, 2); ax.set_ylim(8, -3)

# check that the distance translation for the gaia data doesn't make
# things look absolutely crazy... I'm pretty happy with it i would say. 
fig, axs = plt.subplots(1,2)
axs[0].scatter(t_use["bp_rp"], t_use["phot_g_mean_mag"], s=5)
axs[1].scatter(BP_RP_cat, G_cat, s=5)
axs[1].scatter(BP_sorted[:N]-RP_sorted[:N], G_sorted[:N])
axs[0].scatter(BP_sorted[:N]-RP_sorted[:N], mG_sorted[:N])
for ax in axs:
    ax.invert_yaxis()

# %%


# %%


def assign_photometry_from_catalog(gaia_iso, distances,
                                   catalog, N=None):
    """
    gaia_iso should be [G, BP, RP] from the isochrone. 

    returns photometry row matched brightest to faintest. 
    below the catalog limit (or after N bright stars) rest of photometry 
    is nans. 
    will return mG, BP_RP matched to catalog data. 

    catalog should be an astropy table object. can load with the function above. 
    """
    N_jarvis = 679 #<-- # stars in the jarvis catalog. 
    # STEP 1: 
    #   compute apparent magnitudes given the gaia isochrones and distances
    #
    # STEP 2: 
    #   order the iso tables and the 
    #   catalog by aparent G-band magnitude. for the
    #   catalog, the column will be named phot_g_mean_mag
    #
    # STEP 3:
    #   assign the brightest N stars from the isochrone to the
    #   brightest N stars in the catalog. if N is None, assign the 
    #   brightest len(catalog) stars from the isochrone to the mags 
    #   from the catalog. I care about keeping the catalog apparent G mags
    #   and BP-RP colors. 
    #   make all other rows of the isochrone table should be nan. 
    #
    # STEP 4:
    #   use the distances to convert the catalog apparent magnitudes
    #   back to aboslute magnitudes. we now have G, BP-RP (absolute) 
    #   from catalogs.
    # 
    # STEP 5: 
    #   look up what the Teffs should be given G, BP, RP using the provided
    #   function Teff_from_gaia_isochrone. Because this is in real data world
    #   the points won't lie exactly on the isochrone, so find the nearest point
    #   on like idk the interpolated isochrone or something and assign the Teff
    #   that way?
    #
    # STEP 6:
    #   put the row-matched table back in the original order, return a flag
    #   for stuff that got successfully matched to photometry. 
    return


# %%
# ###### MAIN PROGRAM BELOW 
# if __name__=="__main__":
#     ############# LOAD ISOCHRONES
#     #   old stuff that was determining iso ages to download from like mist.com. using artpop instead. 
#     # np.log10(12e9) #<-- print the log(age) isochrone i want. 
#     # 10**10.07918 / 1e9
#     # np.log10(2700e6)
#     # 10**( 9.43136) / 1e9 #<-- roughly 2.7 Gyr. 

#     ### 2700 Myr (dynamical age) isochrone
#     # isocmd = read_mist_models.ISOCMD('/n/home02/amphillips/data/MIST_2700Myr_fehn2_ubvraplus/MIST_iso_6ab14f0410dce.iso.UBVRIplus')
#     # age_ind = isocmd.age_index(9.43136)

#     ### 12 Gyr isochrone
#     # isocmd = read_mist_models.ISOCMD('/n/home02/amphillips/data/MIST_12Gyr_fehn2_ubvraplus/MIST_iso_6ab1420039622.iso.UBVRIplus')
#     # age_ind = isocmd.age_index(10.07918)

#     # G_iso = isocmd.isocmds[age_ind]['Gaia_G_EDR3']
#     # BP_iso = isocmd.isocmds[age_ind]['Gaia_BP_EDR3']
#     # RP_iso = isocmd.isocmds[age_ind]['Gaia_RP_EDR3']
#     # z_iso = gaia_g_to_lsst_z(G_iso, BP_iso-RP_iso)

#     # mass_iso = isocmd.isocmds[age_ind]['star_mass']
#     # m0_iso = isocmd.isocmds[age_ind]['initial_mass']
#     # Teff_iso = 10**isocmd.isocmds[age_ind]['log_Teff']
#     # L_iso = 10**isocmd.isocmds[age_ind]['log_L']

#     isocmd = artpop.fetch_mist_iso_cmd(
#         log_age=np.log10(12e9),
#         feh=-2.0,
#         phot_system='UBVRIplus',
#         #v_over_vcrit=0.0
#     )
#     G_iso = isocmd['Gaia_G_EDR3']
#     BP_iso = isocmd['Gaia_BP_EDR3']
#     RP_iso = isocmd['Gaia_RP_EDR3']
#     z_iso = gaia_g_to_lsst_z(G_iso, BP_iso-RP_iso)

#     mass_iso = isocmd['star_mass']
#     m0_iso = isocmd['initial_mass']
#     Teff_iso = 10**isocmd['log_Teff']
#     L_iso = 10**isocmd['log_L']


#     ##### plots showing what masses we have... 
#     # plt.hist(isocmd.isocmds[age_ind]['initial_mass'], bins=10)
#     # plt.xlabel(r'$M_{ini}$')
#     # print info
#         # print(isocmd.photo_sys)
#         # print(isocmd.ages)
#         # print(isocmd.hdr_list)
#     # fig, ax = plt.subplots()
#     # ax.scatter(BP_iso - RP_iso, G_iso, c=np.log10(mass_iso))
#     # ax.invert_yaxis()

#     ### load Jarvis+26 table 7
#     # tt = Table.read('/n/home02/amphillips/data/jarvis26_Table7.fits', format='fits')
#     # tt.colnames

#     ### test whether z CMD looks reasonable.
#     # fig, ax = plt.subplots()
#     # ax.scatter(BP_iso - RP_iso, z_iso, c=np.log10(mass_iso))
#     # ax.invert_yaxis()


#     ### okay... now i guess go about painting mags onto my stars...
#     #   a ~12 Gyr isochrone is not quite fair, since my simulated stars are not 
#     #   that old. hmmm..... idk 
#     grid_info = paf.extended_grid_info(scratch=False) 
#     lm_colors, hm_colors, simcolors = paf.define_simcolors()
#     reordered_colors = hm_colors + lm_colors[::-1]
#     cc = reordered_colors[:-1]
#     prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')

#     # ordering decided here. 
#     orbits = ['gd1','aau','pa5','jet','m3','c19']
#     init_displacements = [
#         grid_info.gd1_init_displacement, 
#         grid_info.aau_init_displacement,
#         grid_info.pa5_init_displacement,
#         grid_info.jet_init_displacement,
#         grid_info.m3_init_displacement,
#         grid_info.c19_init_displacement]
#     masses = ['lm','hm']
#     rvirs = [0.75, 1.5, 3, 6]
#     copy_options = [0,1,2,3,4]
#     # copy_options = [4,3,2,1,0]

#     keys = ['d_phi2','v_phi1','v_phi2','v_gsr'] 

#     ii=0
#     orbit = orbits[ii]

#     rvir_index=0
#     rvir = rvirs[rvir_index]

#     mass_index=1


#     (core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
#         simspect.prepare_nbody_data_anycopy(
#             orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
#             include_photometry=False, N_rtid_boundary=2.0 #<-- increase tidal boundary... will this work?
#         )

#     # get distances also:
#     coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab)
#     distances = coords_obs.distance
#     sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) #<-- sc is returned as a DICTIONARY! 

#     unbound = ~CMdict['in_rtid']
#     # unbound = unbound[inMW][trim] # don't care about this. 

#     stellar_masses = lumdict['mass'].to(u.Msun).value
#     stellar_types = lumdict['type']
#     nonrem = stellar_types < 10 # <-- i think 10 starts to be WD. 

#     ### the genuine birth mass -- NOT m0_effective (BSE's fitting coordinate) and
#     ### NOT mass (current mass). see the star.mass0 section of the readme.
#     m0s = np.asarray(lumdict['m0_zams'], dtype=np.float64)

#     L, R = lumdict['L'].to(u.Lsun), lumdict['R'].to(u.Rsun)
#     Teff = get_Teff(R, L).to(u.K)

#     USE_ISOCHRONE = True   # False falls back to the blackbody + top-hat scheme
#     ISO_MAX_PHASE = 5      # drop post-AGB(6); see build_isochrone_table

#     if USE_ISOCHRONE:
#         ### paint photometry by interpolating the isochrone in INITIAL mass.
#         m0_grid, iso_table = build_isochrone_table(isocmd,
#                                                 max_phase=ISO_MAX_PHASE)
#         min_m0, max_m0 = m0_grid[0], m0_grid[-1]

#         iso_phot, on_iso = isochrone_photometry(m0s, m0_grid, iso_table)
#         G  = iso_phot['Gaia_G_EDR3']
#         BP = iso_phot['Gaia_BP_EDR3']
#         RP = iso_phot['Gaia_RP_EDR3']
#         log_Teff = iso_phot['log_Teff']


#         ### the isochrone is old, so it stops at the turnoff: everything more
#         ### massive than max_m0 has to be thrown out. (a handful of sim stars can
#         ### also sit below the isochrone's low-mass end, so this is two-sided.)
#         alive = on_iso
#         n_hi = int(np.count_nonzero(m0s > max_m0))
#         n_lo = int(np.count_nonzero(m0s < min_m0))
#         print(f'isochrone m0 range [{min_m0:.10f}, {max_m0:.10f}] Msun, '
#             f'{len(m0_grid)} nodes (min spacing {np.diff(m0_grid).min():.3e})')
#         print(f'painted {alive.sum()}/{len(m0s)} stars; dropped {n_hi} above and '
#             f'{n_lo} below the isochrone')
#     else:
#         max_m0 = m0_iso.max()
#         min_m0 = m0_iso.min()
#         alive = (m0s <= max_m0) & (m0s >= min_m0)
#         G, BP, RP = g_phot(Teff, R, dpc=10)



#     ### NEW scheme for trimming the stream just dropped, no 'inMW' necessary now. 
#     inMW_na = np.ones(len(sc['phi1']), dtype=bool) #<-- i don't actually want to do a "inMW" trim here. 
#     trim_new = gmm.trim_obstream_percentile(sc) # & ((sc['phi1']>5) & (sc['phi1']<15))

#     #### apply the trim to all of the coordinates:
#     trimmed_sc = simspect.clip_coords(sc, [inMW_na, trim_new]) #<-- this applies inMW, trim to the coordinate dictionary
#     sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction



#     ### distance things, get RV uncertainties
#     BP_RP = BP - RP
#     z = gaia_g_to_lsst_z(G, BP_RP)

#     mz = paf.m_from_M(z, dist=distances)# 10*u.kpc)
#     mG = paf.m_from_M(G, dist=distances)# 10*u.kpc)

#     rverr_desi = desi_RVerr(zmag=mz, feh=-2.0)
#     rverr_via = viamock.get_viaspec_errors(
#         G = mG, 
#         log_Teff=log_Teff, 
#         feh=-2.0, 
#         exptime_s=3600,
#         nexp=1 #<--- ** WHETHER TO BREAK UP THE EXPOSURES IS A CHOICE TO BE MADE; IT DOES MAKE A SLIGHT DiFFERENCE IN THE NOISE 
#         # could add seeing, moon %, ...
#     )


#     pm_err = total_proper_motion_uncertainty(mG, 'dr3') / np.sqrt(2) #<-- we'll add some in two dimensions
#     pos_err = total_position_uncertainty(mG, 'dr3') / np.sqrt(2)

#     rng = np.random.default_rng(seed=42)
#     pmphi1_noise = (rng.normal(0, pm_err) * u.microarcsecond / u.yr).to(u.mas/u.yr)
#     pmphi2_noise = (rng.normal(0, pm_err) * u.microarcsecond / u.yr).to(u.mas/u.yr)
#     phi1_noise = (rng.normal(0, pos_err) * u.microarcsecond).to(u.degree)
#     phi2_noise = (rng.normal(0, pos_err) * u.microarcsecond).to(u.degree)




#     usePos = nonrem[trim_new] & unbound[trim_new] & alive[trim_new] #<-- sc_straighter already has trim_new applied. 
#     usePhot = nonrem & unbound & alive & trim_new


#     ### DECIDE ON AN RV THRESHOLD AND WHICH SURVEY TO EMULATE
#     RV_thresh = 1. # km/s
#     errs_used = rverr_via # <rverr_desi or via
#     vgsr_noise = rng.normal(0, errs_used)

#     noise_dict = {
#         'phi1': phi1_noise.to(u.degree).value,
#         'phi2': phi2_noise.to(u.degree).value,
#         'pm_phi1': pmphi1_noise.to(u.mas/u.yr).value,
#         'pm_phi2': pmphi2_noise.to(u.mas/u.yr).value,
#         'v_gsr': vgsr_noise #<-- already in km/s ig. 
#     }


#     noisey_selection_pos = errs_used[trim_new] < RV_thresh
#     noisey_selection_phot = errs_used < RV_thresh


#     ### idK guys.... .....
#     fig, ax = plt.subplots()
#     # ax.hist(rverr_via, bins=30)
#     ax.hist(mG[usePhot & noisey_selection_phot])



#     keys = ['phi2','pm_phi1','pm_phi2','v_gsr']
#     fig, axs = plt.subplots(len(keys), figsize=[8, 10])

#     plt.subplots_adjust(hspace=0.03, wspace=0.03)

#     key_labels = [
#         r'$\phi_2~[\degree]$',
#         r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$',
#         r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$',
#         r'$v_{\rm GSR}~[\rm km~s^{-1}]$'
#     ]


#     for jj, key in enumerate(keys):
#         # cocoon_selection = p_thin<0.5

#         ydata = sc_straighter[key][usePos] 
#         std = np.std(ydata)
#         cut = 3*std #cocoon_sigmas[jj]

#         ax = axs[jj]#,0]

#         ax.scatter(sc_straighter['phi1'][usePos & noisey_selection_pos] + noise_dict['phi1'][usePhot & noisey_selection_phot],  # plot cocoon on top. 
#                 sc_straighter[key][usePos & noisey_selection_pos] + noise_dict[key][usePhot & noisey_selection_phot], # plot cocoon on top. 
#                 # x_data[:,ii],
#                     c='k', s=5, 
#                     rasterized=True) 
#         # ax.set_ylim(-cut,cut)


#         # ax.set_ylim(-3*cut, 3*cut)
#         ax.set_ylabel(key_labels[jj], fontsize=15)
#         ax.set_xlim(-80, 0) #<-- gd1
#         # ax.set_xlim(-20,10) #<-- jet



#         # too annoying to get the limits to work out. being unrigorous for now...
#         # ax.set_xticks([])
#         # ax.set_xticklabels([])
#         # ax.set_xlim(0, 0.4)
#         # ax.set_yticks([])
#         # ax.set_yticklabels([])

#         if jj<3:
#             # print("REMOVING TICK LABLES>>>>>")
#             axs[jj].set_xticklabels([])


#     axs[-1].set_xlabel(r'$\phi_1~[\degree]$')


#     ### CHECK CMD
#     # fig, ax = plt.subplots()
#     # time_cmap = paf.define_time_cmap()
#     # ax.scatter(BP_RP[usePhot], mG[usePhot], 
#     #            c=np.log10(m0s[usePhot]),
#     #         #    c=log_Teff[usePhot],
#     #         cmap=time_cmap.reversed(),
#     #            edgecolor='k', lw=.5, s=30)
#     # # iso_cutoff = -700
#     # # ax.scatter(BP_iso[:iso_cutoff]-RP_iso[:iso_cutoff], mz_iso[:iso_cutoff],
#     # #            c=np.log10(m0_iso[:iso_cutoff]), vmin=-1, vmax=np.log10(max_m0), zorder=0)
#     # ax.set_xlabel(r'$G_{\rm BP} - G_{\rm RP}$')
#     # ax.set_ylabel(r'$G$')
#     # ax.invert_yaxis()
#     # %%