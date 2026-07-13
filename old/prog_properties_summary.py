# Demonstrating why we might suspect no
#  external heating to create a cocoon 
# structure in the GD-1 stream.
# *** NOTE: this will be replaced with an updated simulation grid (more streams, better dynamical ages)
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
# %%
print("doing setup for petar data...")
paths = paf.define_paths()
tdis_vals, t_peri_vals, t_apo_vals = paf.get_tdis_tplot(paths) ### later. 

apocenters = paf.define_apocenters()
init_displacements = paf.define_init_displacements() 

rvir0_values = np.array([0.75,0.75, 1.5,1.5,3.,3.,6.,6.,]*3)
# %%
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


    axs[0,0].scatter(phi1_[ordering][::-1], phi2[ordering][::-1], c=probs[ordering][::-1], cmap=cmap, s=4 if all_panels==False else 1, rasterized=True)
    axs[0,0].set_ylim(-phi2_lim, phi2_lim)
    axs[0,0].set_xlim(-phi1_lim, phi1_lim)
    axs[0,1].hist(phi2[ts_flag], bins=phi2_bins, density=True, histtype='step', color=colors[1], orientation='horizontal', lw=2);
    axs[0,1].hist(phi2[~ts_flag], bins=phi2_bins, density=True, histtype='step', color=colors[0], orientation='horizontal', lw=2);
    axs[0,0].set_ylabel(r'$\phi_2~[\degree]$')


    axs[1,0].scatter(phi1_[ordering][::-1], vr[ordering][::-1], c=probs[ordering][::-1], cmap=cmap, s=4 if all_panels==False else 1, rasterized=True)
    axs[1,0].set_ylim(-vr_lim, vr_lim)
    axs[1,0].set_xlim(-phi1_lim, phi1_lim)
    axs[1,0].set_ylabel(r'$v_r~\rm[km~s^{-1}]$')

    axs[1,1].hist(vr[ts_flag], bins=vr_bins, density=True, histtype='step', color=colors[1], orientation='horizontal', lw=2);
    axs[1,1].hist(vr[~ts_flag], bins=vr_bins, density=True, histtype='step', color=colors[0], orientation='horizontal', lw=2);

    if all_panels==True:
        axs[2,0].scatter(phi1_[ordering][::-1], vphi1[ordering][::-1], c=probs[ordering][::-1], cmap=cmap, s=1, rasterized=True)
        axs[2,0].set_ylim(-vphi1_lim, vphi1_lim)
        axs[2,0].set_xlim(-phi1_lim, phi1_lim)
        axs[2,1].hist(vphi1[ts_flag], bins=vphi1_bins, density=True, histtype='step', color=colors[1], orientation='horizontal');
        axs[2,1].hist(vphi1[~ts_flag], bins=vphi1_bins, density=True, histtype='step', color=colors[0], orientation='horizontal');
        axs[2,0].set_ylabel(r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$')

        axs[3,0].scatter(phi1_[ordering][::-1], vphi2[ordering][::-1], c=probs[ordering][::-1], cmap=cmap, s=1, rasterized=True)
        axs[3,0].set_ylim(-vphi2_lim, vphi2_lim)
        axs[3,0].set_xlim(-phi1_lim, phi1_lim)
        axs[3,1].hist(vphi2[ts_flag], bins=vphi2_bins, density=True, histtype='step', color=colors[1], orientation='horizontal');
        axs[3,1].hist(vphi2[~ts_flag], bins=vphi2_bins, density=True, histtype='step', color=colors[0], orientation='horizontal');
        axs[3,0].set_ylabel(r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$')

        axs[3,0].set_xlabel(r'$\phi_1~[\rm \degree]$')
        axs[3,1].set_xlabel("PDF")
    

    return fig, axs

def get_outlier_clip(vr, pmphi1, pmphi2):
    ### this maybe avoids biasing the dispersion of the cocoon ? ?
    outlier_clip = (np.abs(pmphi1)<1.5) & (np.abs(pmphi2)<1.5) & (np.abs(vr)<100)
    return outlier_clip

def poly_straightening(nbody_coords, nbody_for_straightening):
    """
    returns a unitless dictionary of straightened stream coordinates
    first coords are to be straightened, 
    second coords are for fitting the straightening polynomial. 
    """
    # fit and subtract polynomials from the residual curves. 
    nbody_coords_straighter = {"phi1":nbody_coords['phi1']}
    units = [u.degree, u.kpc, u.km/u.s, u.mas/u.yr, u.mas/u.yr] # phi2, dist, vr, pms
    for k, key in enumerate(list(nbody_coords.keys())[:-1]):
        y = nbody_coords[key].to(units[k]).value
        _, _, poly, fit = paf.straighten_stream_polynomial(nbody_for_straightening['phi1'].to(u.degree).value, nbody_for_straightening[key].to(units[k]).value,
                                                        trim_criteria = [np.ones(len(nbody_for_straightening[key])).astype(bool),np.ones(len(nbody_for_straightening[key])).astype(bool)],
                                                        return_poly_fn=True)
        y -= poly(nbody_coords_straighter['phi1'].to(u.degree).value, *fit)
        nbody_coords_straighter[key] = y
    return nbody_coords_straighter
# %%
#--------------------------------------------------------------------#
#                 THE LOOP IS HERE                                   #
#--------------------------------------------------------------------#

phi2_line_GD1, vr_line_GD1, vphi1_line_GD1, vphi2_line_GD1 = [1.5,2.5,0.15,0.1] # GD1


phi2_line_pal5, vr_line_pal5, vphi1_line_pal5, vphi2_line_pal5 = [2.5,6,0.15,0.1] # pa5

n_list = np.arange(8, 16, 1)
# n_list = [8,16, 10,18, 12,20, 14,22]
gd1_n = [8,9,10,11,12,13,14,15]
pal5_n = [16,17,18,19,20,21,22,23]

### lists to store stream widths and cocoon fraction/dispersions
cocoon_fractions_cm = []
sigvrs_ts_cm, sigphi2s_ts_cm = [],[]
sigvrs_c_cm, sigphi2s_c_cm = [],[]

cocoon_fractions_lum = []
sigvrs_ts_lum, sigphi2s_ts_lum = [],[]
sigvrs_c_lum, sigphi2s_c_lum = [],[]


# n_list=[16]
for j, n in tqdm(enumerate(n_list)):
    if n in gd1_n:
        phi2_line, vr_line, vphi1_line, vphi2_line = phi2_line_GD1, vr_line_GD1, vphi1_line_GD1, vphi2_line_GD1
    if n in pal5_n:
        phi2_line, vr_line, vphi1_line, vphi2_line = phi2_line_pal5, vr_line_pal5, vphi1_line_pal5, vphi2_line_pal5


    core, data_dict, CMdict, lumdict, inMW, trim = dc.prepare_nbody_data(n, include_photometry=False)
    nbody_coords_cm, _, _ = dc.get_GD1_coords_nbody(data_dict, lumdict, subdict=CMdict, inMW=inMW, trim=trim, core=core, i=10000)
    nbody_coords_lum, _, _ = dc.get_GD1_coords_nbody(data_dict, lumdict, subdict=lumdict, inMW=inMW, trim=trim, core=core, i=10000)


    nbody_coords_straight_cm = poly_straightening(nbody_coords_cm, nbody_coords_cm)
    nbody_coords_straight_lum = poly_straightening(nbody_coords_lum, nbody_coords_cm)

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

    ########### COMMENT BELOW IF I ACTUALLY WANT TO RUN THE FULL GRID SUMMARY ! 
    # if n==8:
    #     ### plot options. 
    #     fig, axs = stream_plot(phi1_cm[oc_cm], phi2_cm, v_gsr_cm, pm_phi1_cm, pm_phi2_cm,
    #                         memprobs=memberships,
    #                         phi2_lim = 2*phi2_line,
    #                         vr_lim=3*vr_line,
    #                         vphi1_lim = 2*vphi1_line,
    #                         vphi2_lim = 2*vphi2_line,
    #                         # phi1_lim=150,
    #                         all_panels=True) 
        
    #     if n in gd1_n:
    #         for ax in axs[:,0]:
    #             ax.set_xlim(-150, 75)
    #     if n in pal5_n:
    #         for ax in axs[:,0]:
    #             ax.set_xlim(0,)
       

    #     axs[0,0].set_xticklabels([])
    #     axs[0,1].set_yticklabels([])
    #     axs[0,1].set_xticklabels([])
    #     axs[1,1].set_yticklabels([])
    #     axs[0,0].axhline(-phi2_line,c='k', lw=1, ls='--')
    #     axs[0,0].axhline(phi2_line,c='k', lw=1, ls='--')
    #     axs[1,0].axhline(-vr_line,c='k', lw=1, ls='--')
    #     axs[1,0].axhline(vr_line,c='k', lw=1, ls='--')

    #     # axs[2,0].axhline(-vphi1_line,c='k', lw=1, ls='--')
    #     # axs[2,0].axhline(vphi1_line,c='k', lw=1, ls='--')
    #     # axs[3,0].axhline(-vphi2_line,c='k', lw=1, ls='--')
    #     # axs[3,0].axhline(vphi2_line,c='k', lw=1, ls='--')


    #     # plt.subplots_adjust(wspace=0.3)
    #     # plt.savefig("fig/example_cocoon_separation.pdf", dpi=300, bbox_inches='tight')
    #     # plt.savefig("plots/example_cocoon_separation.pdf", dpi=300, bbox_inches='tight')
    #     break


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
def make_summary_plot(cocoon_fractions,
                     sigphi2s_ts, sigphi2s_c,
                      sigvrs_ts, sigvrs_c,
                      add_legend=True,
                      add_hiOB=True,
                      add_ts=True):
    #### plotting the summarized result here. 
    rvir0_list = np.array([0.75,0.75, 1.5,1.5, 3,3, 6,6])

    lm_index = np.array([True, False]*4)
    markers = np.array(['o', '^']*4)

    # ii=lm_index

    # ii = np.ones(len(rvir0_list)).astype(bool)
    cc = "#729B79"
    ss=70

    fig, axs = plt.subplots(1, 3, figsize=[16, 4.5])
    plt.subplots_adjust(wspace=0.3)

    ax=axs[0]
    ii=lm_index
    ax.scatter(rvir0_list[ii], cocoon_fractions[ii], c='k', marker = markers[ii][0], s=ss, edgecolor='k')
    
    if add_hiOB==True:
        ii=~lm_index
        ax.scatter(rvir0_list[ii], cocoon_fractions[ii], c='k', marker = markers[ii][0], s=ss, edgecolor='k')
    
    ax.set_ylabel(r'$f_{\rm cocoon}$')

    ax = axs[1]
    ii=lm_index
    if add_ts==True:
        ax.scatter(rvir0_list[ii], sigphi2s_ts[ii], label='thin stream', c=cc, marker=markers[ii][0], s=ss, edgecolor='k')
    ax.scatter(rvir0_list[ii], sigphi2s_c[ii], label='cocoon', c='k', marker = markers[ii][0], s=ss, edgecolor='k')
    if add_hiOB==True:
        ii=~lm_index
        if add_ts==True:
            ax.scatter(rvir0_list[ii], sigphi2s_ts[ii], label='thin stream', c=cc, marker=markers[ii][0], s=ss, edgecolor='k')
        ax.scatter(rvir0_list[ii], sigphi2s_c[ii], label='cocoon', c='k', marker = markers[ii][0], s=ss, edgecolor='k')
    ax.set_ylabel(r'$\sigma_{\phi_2}~[\degree]$')


    ax = axs[2]
    ii=lm_index
    if add_ts==True:
        ax.scatter(rvir0_list[ii], sigvrs_ts[ii], label='thin stream, lo-OB', c=cc, marker=markers[ii][0], s=ss, edgecolor='k')
    ax.scatter(rvir0_list[ii], sigvrs_c[ii], label='cocoon, lo-OB', c='k', marker = markers[ii][0], s=ss, edgecolor='k')

    if add_hiOB==True:
        ii=~lm_index
        if add_ts==True:
            ax.scatter(rvir0_list[ii], sigvrs_ts[ii], label='thin stream, hi-OB', c=cc, marker=markers[ii][0], s=ss, edgecolor='k')
        
        ax.scatter(rvir0_list[ii], sigvrs_c[ii], label='cocoon, hi-OB', c='k', marker = markers[ii][0], s=ss, edgecolor='k')

    ax.set_ylabel(r'$\sigma_{v_r}~[\rm km~s^{-1}]$')


    if add_legend==True:
        ax.legend(loc='upper left', bbox_to_anchor=[1,1])

    for ax in axs:
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
        ax.minorticks_off()
        ax.set_xticks([0.75, 1.5, 3., 6.])
        # ax.set_xticklabels(["0.75", "1.5", "3", "6"], fontsize=15)
        ax.tick_params(labelsize=15)

    return fig, axs
# %%
#### no binaries:
fig, axs = make_summary_plot(cocoon_fractions_cm,
                             sigphi2s_ts_cm, 
                             sigphi2s_c_cm,
                             sigvrs_ts_cm, 
                             sigvrs_c_cm,
                             add_hiOB=False,add_legend=False, add_ts=False)

axs[0].set_ylim(0, 0.07)
axs[1].set_ylim(0, 2.9)
axs[2].set_ylim(0, 9)

# plt.savefig('plots/prog_properties_cocoon_fraction.pdf', dpi=300, bbox_inches='tight')

# fig.suptitle("``system'' data")
# plt.savefig("fig/summary_noiseless_CM.pdf", dpi=300, bbox_inches='tight')

# fig, axs = make_summary_plot(cocoon_fractions_lum,
#                              sigphi2s_ts_lum, sigphi2s_c_lum,
#                              sigvrs_ts_lum, sigvrs_c_lum)
# fig.suptitle("single epoch with binaries")
# plt.savefig("fig/summary_noiseless_lum.pdf", dpi=300, bbox_inches='tight')
# %%
