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

# %%
### import functions from the desi comparison script:
sys.path.append('/n/home02/amphillips/p27_nbody/old')
# from DESI_comparison import get_gaia_photometry, get_z_photometry, unpack_data_dict, prepare_nbody_data
import DESI_comparison as dc
# %%
print("doing setup for petar data...")
paths = paf.define_paths()
tdis_vals, t_peri_vals, t_apo_vals = paf.get_tdis_tplot(paths) ### later. 

apocenters = paf.define_apocenters()
init_displacements = paf.define_init_displacements() 

rvir0_values = np.array([0.75,0.75, 1.5,1.5,3.,3.,6.,6.,]*3)
# %%
print("defining some helpers")
def get_streamcoords_nbody(data_dict, lumdict, subdict, inMW, trim, core, 
                           stream_frame, i=10000):
    """
    a more general version of get_GD1_coords_nbody() from old/DESI_comparison.py
    streamframe is:
        gc.GD1Koposiv10
        gc.Pal5PriceWhelan18
        or
        'circular'
    """
    init_displacement = data_dict['init_displacement']
    pos, vel = subdict['pos'], subdict['vel']
    coords = subdict['coords']

    primary_types = lumdict["type"][inMW][trim] #<-- that said, always store the type (and eventually stellar parameters) from the luminous companion. since I order the center of mass dictionary and the "luminous" dictionary as [singles, binaries] this ordering should be the same. 
    nonrem = primary_types<10

    ### step 1: transform galactocentric position+velocity to ICRS
    icrs_coords = paf.galcen_to_ICRS(pos, vel)
    ### step 2: correct for solar reflex motion
    icrs_coords = reflex_correct(icrs_coords)

    ### step 3: transform to great circle frame. 
    if stream_frame=='circular':
        return "NOT READY FOR THIS CASE YET"
    else:
        stream_coords = icrs_coords.transform_to(stream_frame)

    ### step 4 catalog these. note v_gsr should come from the reflex-corrected icrs coords. 
    phi1 = stream_coords.phi1[inMW][trim].to(u.degree)
    phi2 = stream_coords.phi2[inMW][trim].to(u.degree)
    v_gsr = icrs_coords.radial_velocity[inMW][trim].to(u.km/u.s)
    pmphi1 = stream_coords.pm_phi1_cosphi2[inMW][trim].to(u.mas/u.yr)
    pmphi2 = stream_coords.pm_phi2[inMW][trim].to(u.mas/u.yr)
    dist = icrs_coords.distance[inMW][trim].to(u.kpc)

    # these are already inMW + clip clipped, but not straightened or "outlier clipped."
    nbody_coords = {
        "phi1":phi1,
        "phi2":phi2,
        "dist":dist,
        "v_gsr":v_gsr,
        "pm_phi1":pmphi1,
        "pm_phi2":pmphi2,
    }

    ### also want to integrate progenitor orbit and xform to icrs...
    sc, interp_orbit, orbit_coords, orbit_w, itr = paf.straighten_stream_orbit_interp(
        coords, 'phi2', core, i, 
        trim_criteria=[inMW, trim], Dt_start=240,
        return_orbit_chunk=True, #<-- mostly performing this function to get orbit_w.
        use_core=False,
        init_displacement=init_displacement
    )

    orbit_pos, orbit_vel = orbit_w[:,:3] * u.kpc, orbit_w[:,3:]*u.kpc/u.Myr
    orbit_icrs_coords = paf.galcen_to_ICRS(orbit_pos, orbit_vel)
    orbit_icrs_coords = reflex_correct(orbit_icrs_coords)
    orbit_gd1_coords = orbit_icrs_coords.transform_to(gc.GD1Koposov10)
    orbit_phi1 = orbit_gd1_coords.phi1
    orbit_phi2 = orbit_gd1_coords.phi2
    orbit_v_gsr = orbit_icrs_coords.radial_velocity
    orbit_pmphi1 = orbit_gd1_coords.pm_phi1_cosphi2
    orbit_pmphi2 = orbit_gd1_coords.pm_phi2
    orbit_dist = orbit_icrs_coords.distance.to(u.kpc)

    trim_orbit = (orbit_phi1<=np.max(phi1)) & (orbit_phi1>=np.min(phi1))
    orbit_sort = np.argsort(orbit_phi1[trim_orbit])

    orbit_coords = {
        "phi1":orbit_phi1[trim_orbit][orbit_sort],
        "phi2":orbit_phi2[trim_orbit][orbit_sort],
        "dist":orbit_dist[trim_orbit][orbit_sort],
        "v_gsr":orbit_v_gsr[trim_orbit][orbit_sort],
        "pm_phi1":orbit_pmphi1[trim_orbit][orbit_sort],
        "pm_phi2":orbit_pmphi2[trim_orbit][orbit_sort]
    }

    def subtract_orbit_interp(yval):
        x = nbody_coords['phi1']
        y = nbody_coords[yval]
        xp, yp = orbit_coords['phi1'], orbit_coords[yval]
        return y - np.interp(x, xp, yp)


    nbody_coords_straight = {key:subtract_orbit_interp(key) for key in list(orbit_coords.keys())[1:]}
    nbody_coords_straight["phi1"] = nbody_coords["phi1"]

    return nbody_coords_straight, nbody_coords, nonrem

def stream_plot(phi1_, phi2, vr, vphi1, vphi2, memprobs=None,
                phi1_lim = 100, phi2_lim = 2, vr_lim = 6, vphi1_lim = 20, vphi2_lim = 30,
                # colors = ["#87CBAC", "#156064",], 
                colors = ["#8DCB87", "#156064"],
                all_panels=True):
    """
    !!! apply any outlier clips before this ! ! ! 


    """
    cmap = LinearSegmentedColormap.from_list('cmap', colors)
    if memprobs is not None:
        probs=memprobs  
        ts_flag = probs>0.5 
        ordering = np.argsort(probs)
    else:
        probs = np.zeros(len(phi1_))
        ordering = np.arange(0, len(phi1_), 1)
        ts_flag = np.ones(len(phi1_)).astype(bool)


    phi2_bins = np.linspace(-phi2_lim, phi2_lim, 50)
    vr_bins = np.linspace(-vr_lim, vr_lim, 50)
    vphi1_bins = np.linspace(-vphi1_lim, vphi1_lim, 50)
    vphi2_bins = np.linspace(-vphi2_lim, vphi2_lim, 50)

    if all_panels==True:
        fig, axs = plt.subplots(4,2,width_ratios=[2,1], figsize=[12,14])
        axs[0,1].set_ylabel(r'$\phi_2~[\degree]$')
        axs[1,1].set_ylabel(r'$v_r~\rm[km~s^{-1}]$')
        axs[2,1].set_ylabel(r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$')
        axs[3,1].set_ylabel(r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$')


    if all_panels==False:
        fig, axs = plt.subplots(2,2, width_ratios=[2,1], figsize=[10,6])
        plt.subplots_adjust(wspace=0.03, hspace=0.03)
        axs[1,0].set_xlabel(r'$\phi_1~[\rm \degree]$')
        axs[1,1].set_xlabel("PDF")


    axs[0,0].scatter(phi1_[ts_flag], phi2[ts_flag], c=colors[1], s=4, rasterized=True)
    axs[0,0].scatter(phi1_[~ts_flag], phi2[~ts_flag], c=colors[0],
        edgecolor='k', lw=0.5, s=20, rasterized=True)

    
    axs[0,0].set_ylim(-phi2_lim, phi2_lim)
    axs[0,0].set_xlim(-phi1_lim, phi1_lim)
    axs[0,1].hist(phi2[ts_flag], bins=phi2_bins, density=True, histtype='step', color=colors[1], orientation='horizontal', lw=2);
    axs[0,1].hist(phi2[~ts_flag], bins=phi2_bins, density=True, histtype='step', color=colors[0], orientation='horizontal', lw=2);
    axs[0,0].set_ylabel(r'$\phi_2~[\degree]$')


    axs[1,0].scatter(phi1_[ts_flag], vr[ts_flag], c=colors[1], s=4, rasterized=True)
    axs[1,0].scatter(phi1_[~ts_flag], vr[~ts_flag], c=colors[0],
        edgecolor='k', lw=0.5, s=20, rasterized=True)    
    axs[1,0].set_ylim(-vr_lim, vr_lim)
    axs[1,0].set_xlim(-phi1_lim, phi1_lim)
    axs[1,0].set_ylabel(r'$v_r~\rm[km~s^{-1}]$')

    axs[1,1].hist(vr[ts_flag], bins=vr_bins, density=True, histtype='step', color=colors[1], orientation='horizontal', lw=3);
    axs[1,1].hist(vr[~ts_flag], bins=vr_bins, density=True, histtype='step', color=colors[0], orientation='horizontal', lw=3);

    if all_panels==True:
        axs[2,0].scatter(phi1_[ts_flag], vphi1[ts_flag], c=colors[1], s=4, rasterized=True)
        axs[2,0].scatter(phi1_[~ts_flag], vphi1[~ts_flag], c=colors[0],
        edgecolor='k', lw=0.5, s=20, rasterized=True)            
        axs[2,0].set_ylim(-vphi1_lim, vphi1_lim)
        axs[2,0].set_xlim(-phi1_lim, phi1_lim)
        axs[2,1].hist(vphi1[ts_flag], bins=vphi1_bins, density=True, histtype='step', color=colors[1], orientation='horizontal');
        axs[2,1].hist(vphi1[~ts_flag], bins=vphi1_bins, density=True, histtype='step', color=colors[0], orientation='horizontal');
        axs[2,0].set_ylabel(r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$')

        axs[3,0].scatter(phi1_[ts_flag], vphi2[ts_flag], c=colors[1], s=4 if all_panels==False else 1, rasterized=True)
        axs[3,0].scatter(phi1_[~ts_flag], vphi2[~ts_flag], c=colors[0],
        edgecolor='k', lw=0.5, s=20, rasterized=True)            
        axs[3,0].set_ylim(-vphi2_lim, vphi2_lim)
        axs[3,0].set_xlim(-phi1_lim, phi1_lim)
        axs[3,1].hist(vphi2[ts_flag], bins=vphi2_bins, density=True, histtype='step', color=colors[1], orientation='horizontal');
        axs[3,1].hist(vphi2[~ts_flag], bins=vphi2_bins, density=True, histtype='step', color=colors[0], orientation='horizontal');
        axs[3,0].set_ylabel(r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$')

        axs[3,0].set_xlabel(r'$\phi_1~[\rm \degree]$')
        axs[3,1].set_xlabel("PDF")
    

    return fig, axs

def get_outlier_clip(vr, pmphi1, pmphi2):
    ### this avoids biasing the dispersion of the cocoon
    # with just a very few outliers. calibrated to GD1
    outlier_clip = (np.abs(pmphi1)<1.5) & (np.abs(pmphi2)<1.5) & (np.abs(vr)<100)
    return outlier_clip

def poly_straightening(nbody_coords, nbody_for_straightening, deg=5,
                       dostraight=True):
    """
    returns a unitless dictionary of straightened stream coordinates
    first coords are to be straightened, 
    second coords are for fitting the straightening polynomial. 
    """
    # fit and subtract polynomials from the residual curves. 
    nbody_coords_straighter = {"phi1":nbody_coords['phi1']}
    units = [u.degree, u.kpc, u.km/u.s, u.mas/u.yr, u.mas/u.yr] # phi2, dist, vr, pms
    if dostraight==True:
        for k, key in enumerate(list(nbody_coords.keys())[:-1]):
            y = nbody_coords[key].to(units[k]).value
            _, _, poly, fit = paf.straighten_stream_polynomial(nbody_for_straightening['phi1'].to(u.degree).value, 
                                                            nbody_for_straightening[key].to(units[k]).value,
                                                            degree=deg,
                                                            trim_criteria = [np.ones(len(nbody_for_straightening[key])).astype(bool),np.ones(len(nbody_for_straightening[key])).astype(bool)],
                                                            return_poly_fn=True)
            y -= poly(nbody_coords_straighter['phi1'].to(u.degree).value, *fit)
            nbody_coords_straighter[key] = y
    else:
        for k, key in enumerate(list(nbody_coords.keys())[:-1]):
            y = nbody_coords[key].to(units[k]).value
            nbody_coords_straighter[key] = y
    return nbody_coords_straighter
# %%
print("potentially looping...")

phi2_line_GD1, vr_line_GD1, vphi1_line_GD1, vphi2_line_GD1 = [1.5,2.6,0.15,0.1] # GD1


phi2_line_pal5, vr_line_pal5, vphi1_line_pal5, vphi2_line_pal5 = [2.5,25.,0.15,0.4] # pa5

phi2_line_circ, vr_line_circ, vphi1_line_circ, vphi2_line_circ = [0.18,2.2,0.03,0.01] # circ
n_list = np.arange(8, 16, 1)
# n_list = [8,16, 10,18, 12,20, 14,22]
gd1_n = [8,9,10,11,12,13,14,15]
pal5_n = [16,17,18,19,20,21,22,23]
circ_n = [0,1,2,3,4,5,6,7]
### lists to store stream widths and cocoon fraction/dispersions
cocoon_fractions_cm = []
sigvrs_ts_cm, sigphi2s_ts_cm = [],[]
sigvrs_c_cm, sigphi2s_c_cm = [],[]

cocoon_fractions_lum = []
sigvrs_ts_lum, sigphi2s_ts_lum = [],[]
sigvrs_c_lum, sigphi2s_c_lum = [],[]


# n_list = [8,10,12,14,16,18,20,22]

n_list = [8]
for j, n in tqdm(enumerate(n_list)):
    if n in gd1_n:
        d=5
        phi2_line, vr_line, vphi1_line, vphi2_line = phi2_line_GD1, vr_line_GD1, vphi1_line_GD1, vphi2_line_GD1
    if n in pal5_n:
        d=8
        phi2_line, vr_line, vphi1_line, vphi2_line = phi2_line_pal5, vr_line_pal5, vphi1_line_pal5, vphi2_line_pal5
    if n in circ_n:
        d=0
        phi2_line, vr_line, vphi1_line, vphi2_line = phi2_line_circ, vr_line_circ, vphi1_line_circ, vphi2_line_circ


    core, data_dict, CMdict, lumdict, inMW, trim = dc.prepare_nbody_data(n, include_photometry=False)

    if n in gd1_n:
        d=8
        nbody_coords_cm, _, _ = get_streamcoords_nbody(data_dict, lumdict, subdict=CMdict, inMW=inMW, trim=trim, core=core, 
                                                        stream_frame=gc.GD1Koposov10, i=10000)
        nbody_coords_lum, _, _ = get_streamcoords_nbody(data_dict, lumdict, subdict=lumdict, inMW=inMW, trim=trim, core=core, 
                                                        stream_frame=gc.GD1Koposov10, i=10000)

    if n in pal5_n:
        nbody_coords_cm, _, _ = get_streamcoords_nbody(data_dict, lumdict, subdict=CMdict, inMW=inMW, trim=trim, core=core, 
                                                        stream_frame=gc.Pal5PriceWhelan18, i=10000)
        nbody_coords_lum, _, _ = get_streamcoords_nbody(data_dict, lumdict, subdict=lumdict, inMW=inMW, trim=trim, core=core, 
                                                        stream_frame=gc.Pal5PriceWhelan18, i=10000)
        
    if n in circ_n:
        ## in this case just use the orbit-subtracted Jake coordinates. 
        nbody_coords_cm = {
            'phi2':CMdict['phi2_straight']*u.degree,
            'r':CMdict['r_straight']*u.kpc,
            'vr':CMdict['vr_straight']*u.km/u.s,
            'pm_phi1':CMdict['pm_phi1_straight']*u.mas/u.yr,
            'pm_phi2':CMdict['pm_phi2_straight']*u.mas/u.yr,
            'phi1':CMdict['coords']['phi1']*u.degree
        }
        nbody_coords_lum = {
            'phi2':lumdict['phi2_straight']*u.degree,
            'r':lumdict['r_straight']*u.kpc,
            'vr':lumdict['vr_straight']*u.km/u.s,
            'pm_phi1':lumdict['pm_phi1_straight']*u.mas/u.yr,
            'pm_phi2':lumdict['pm_phi2_straight']*u.mas/u.yr,
            'phi1':lumdict['coords']['phi1']*u.degree
        }

    ### subtract any residual wiggles. 
    # if n in circ_n:
    #     nbody_coords_straight_cm = nbody_coords_cm
    #     nbody_coords_straight_lum = nbody_coords_lum
    
    # else:
    nbody_coords_straight_cm = poly_straightening(nbody_coords_cm, nbody_coords_cm, deg=d, dostraight=False if n in circ_n else True)
    nbody_coords_straight_lum = poly_straightening(nbody_coords_lum, nbody_coords_cm, deg=d, dostraight=False if n in circ_n else True)

    phi1_cm, phi2_cm, dist_cm, v_gsr_cm, pm_phi1_cm, pm_phi2_cm = [nbody_coords_straight_cm[key] for key in nbody_coords_straight_cm.keys()]
    phi1_lum, phi2_lum, dist_lum, v_gsr_lum, pm_phi1_lum, pm_phi2_lum = [nbody_coords_straight_lum[key] for key in nbody_coords_straight_lum.keys()]




    ### woaoh woah woah sorry ... also need to get outlier clips lmao...
    oc_cm = get_outlier_clip(v_gsr_cm, pm_phi1_cm, pm_phi2_cm)
    oc_lum = get_outlier_clip(v_gsr_lum, pm_phi1_lum, pm_phi2_lum)
    phi2_cm = phi2_cm[oc_cm]
    phi2_lum = phi2_lum[oc_lum]

    v_gsr_cm = v_gsr_cm[oc_cm]
    v_gsr_lum = v_gsr_lum[oc_lum]

    pm_phi1_cm = pm_phi1_cm[oc_cm]
    pm_phi1_lum = pm_phi1_lum[oc_lum]

    pm_phi2_cm = pm_phi2_cm[oc_cm]
    pm_phi2_lum = pm_phi2_lum[oc_lum]

    ##### 1. CM data:
    cocoon_selection = (np.abs(phi2_cm)>phi2_line) | (np.abs(v_gsr_cm)>vr_line) | (np.abs(pm_phi1_cm)>vphi1_line) | (np.abs(pm_phi2_cm)>vphi2_line)
    ts_selection = ~cocoon_selection
    cocoon, stream = np.zeros(len(phi2_cm)), np.ones(len(phi2_cm))
    memberships = np.where(cocoon_selection, cocoon, stream)
    # memberships = np.where(ts_selection, cocoon, stream)

    ########## COMMENT BELOW IF I ACTUALLY WANT TO RUN THE FULL GRID SUMMARY ! 
    if n==8:
        ### plot options. 
        fig, axs = stream_plot(phi1_cm[oc_cm], phi2_cm, v_gsr_cm, pm_phi1_cm, pm_phi2_cm,
                            memprobs=memberships,
                            phi2_lim = 2*phi2_line,
                            vr_lim= 3*vr_line,
                            vphi1_lim = 2*vphi1_line,
                            vphi2_lim = 2*vphi2_line,
                            # phi1_lim=150,
                            all_panels=False) 
        
        if n in gd1_n:
            for ax in axs[:,0]:
                ax.set_xlim(-150, 75)
        if n in pal5_n:
            for ax in axs[:,0]:
                ax.set_xlim(-50,50)
       

        axs[0,0].set_xticklabels([])
        axs[0,1].set_yticklabels([])
        axs[0,1].set_xticklabels([])
        axs[1,1].set_yticklabels([])
        axs[0,0].axhline(-phi2_line,c='k', lw=1, ls='--')
        axs[0,0].axhline(phi2_line,c='k', lw=1, ls='--')
        axs[1,0].axhline(-vr_line,c='k', lw=1, ls='--')
        axs[1,0].axhline(vr_line,c='k', lw=1, ls='--')

        # axs[2,0].axhline(-vphi1_line,c='k', lw=1, ls='--')
        # axs[2,0].axhline(vphi1_line,c='k', lw=1, ls='--')
        # axs[3,0].axhline(-vphi2_line,c='k', lw=1, ls='--')
        # axs[3,0].axhline(vphi2_line,c='k', lw=1, ls='--')


        # plt.subplots_adjust(wspace=0.3)
        # plt.savefig("fig/example_cocoon_separation.pdf", dpi=300, bbox_inches='tight')
        plt.savefig("plots/example_cocoon_separation.pdf", dpi=300, bbox_inches='tight')
        break


    cocoon_fraction = len(phi2_cm[cocoon_selection])/len(phi2_cm)
    sigvr_ts, sigphi2_ts = np.std(v_gsr_cm[ts_selection]), np.std(phi2_cm[ts_selection])
    sigvr_cocoon, sigphi2_cocoon = np.std(v_gsr_cm[cocoon_selection]), np.std(phi2_cm[cocoon_selection])

    sigvrs_ts_cm.append(sigvr_ts)
    sigphi2s_ts_cm.append(sigphi2_ts)
    sigvrs_c_cm.append(sigvr_cocoon)
    sigphi2s_c_cm.append(sigphi2_cocoon)
    cocoon_fractions_cm.append(cocoon_fraction)


    ##### 2. data with binaries:
    cocoon_selection = (np.abs(phi2_lum)>phi2_line) | (np.abs(v_gsr_lum)>vr_line) | (np.abs(pm_phi1_lum)>vphi1_line) | (np.abs(pm_phi2_lum)>vphi2_line)
    ts_selection = ~cocoon_selection
    cocoon, stream = np.zeros(len(phi2_lum)), np.ones(len(phi2_lum))
    memberships = np.where(cocoon_selection, cocoon, stream)

    cocoon_fraction = len(phi2_lum[cocoon_selection])/len(phi2_lum)
    sigvr_ts, sigphi2_ts = np.std(v_gsr_lum[ts_selection]), np.std(phi2_lum[ts_selection])
    sigvr_cocoon, sigphi2_cocoon = np.std(v_gsr_lum[cocoon_selection]), np.std(phi2_lum[cocoon_selection])

    sigvrs_ts_lum.append(sigvr_ts)
    sigphi2s_ts_lum.append(sigphi2_ts)
    sigvrs_c_lum.append(sigvr_cocoon)
    sigphi2s_c_lum.append(sigphi2_cocoon)
    cocoon_fractions_lum.append(cocoon_fraction)


sigvrs_ts_cm = np.array(sigvrs_ts_cm)
sigphi2s_ts_cm = np.array(sigphi2s_ts_cm)
sigvrs_c_cm = np.array(sigvrs_c_cm)
sigphi2s_c_cm = np.array(sigphi2s_c_cm)
cocoon_fractions_cm = np.array(cocoon_fractions_cm)

sigvrs_ts_lum = np.array(sigvrs_ts_lum)
sigphi2s_ts_lum = np.array(sigphi2s_ts_lum)
sigvrs_c_lum = np.array(sigvrs_c_lum)
sigphi2s_c_lum = np.array(sigphi2s_c_lum)
cocoon_fractions_lum = np.array(cocoon_fractions_lum)
# %%
