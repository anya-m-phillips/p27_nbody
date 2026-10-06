#---------------------------------------------------#
#       bold decisions re the naming of             #
#       this notebook                               #
#---------------------------------------------------#
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
from matplotlib.patches import Ellipse

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
import gmm
import noise as n
# %%
datapath = '/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
table_path = "/n/home02/amphillips/p27_nbody/data/gmm_tables/"

# %%
#-----------------------------------------------#
#   opening all of the data and checking out    #
#   the observability given distances to all    #
#   the streams, which noise i use, etc.        #
#-----------------------------------------------#


make_plots=False
constrain_widths=False
noise = 'via' #<-- None or 'via' or 'desi'
include_binaries=False

if noise is None:
    cd0 = 'noiseless'
if noise is not None:
    cd0 = noise+"_noise" #< so via_noise or desi_noise

if include_binaries==True:
    cd1 = 'binaries'
if include_binaries==False:
    cd1 = 'CoM'

if constrain_widths==True:
    cd2 = '_constrained'
if constrain_widths==False:
    cd2 = ''

case_name = cd0+"_"+cd1+cd2
sfx = '_primaries' if include_binaries else '' #<-- picks the coords, trim, and catalog products, as in gmm.py

print("running case:", case_name)


grid_info = paf.extended_grid_info(scratch=False) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors[:-1] + lm_colors[::-1]
cc = reordered_colors[:-1]
ccc = cc[1:]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')

# ordering decided here. 
orbits = ['gd1','aau','pa5','jet','c19']
phi1_lims = [
    (-100, 10),
    (-15,15),
    (-15,15),
    (-20,20),
    (-20,10)
]
orbit_lim_map = {o:l for o, l in zip(orbits, phi1_lims)}

init_displacements = [
    grid_info.gd1_init_displacement, 
    grid_info.aau_init_displacement,
    grid_info.pa5_init_displacement,
    grid_info.jet_init_displacement,
    grid_info.c19_init_displacement]
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
copy_options = [0,1,2,3,4]


pericenters_kpc = []
apocenters_kpc = []
present_rs = []
med_distances = []

# f_cocoons, cocoon_sigvgsrs, cocoon_sigphi2s, thin_sigvgsrs, thin_sigphi2s = [], [], [], [], []

# fig, ax = plt.subplots()
figg, axx = plt.subplots()

bins = np.linspace(0, 40, 100)


orbit_titles = ['GD-1', 'AAU','Pal 5','Jet','C-19']
for ii, orbit in enumerate(tqdm(orbits)): #<--- this i can do later i think. 

    pr = prog_tab[prog_tab['name']==orbit] #<-- prog row
    disp = [pr['x'][0], pr['y'][0], pr['z'][0], pr['vx'][0], pr['vy'][0], pr['vz'][0]]
    orbit_obj = paf.integrate_prog_orbit(disp, steps=10000, dt=1*u.Myr)

    peri = orbit_obj.pericenter().to(u.kpc).value
    pericenters_kpc.append(peri)
    apo = orbit_obj.apocenter().to(u.kpc).value
    apocenters_kpc.append(apo)

    x,y,z = disp[:3]
    r = np.sqrt(x**2 + y**2 + z**2)
    present_rs.append(r)

    ##### let's open all of the dictionaries also to get like a median distance. use the most diffuse guy.
    filename = datapath+"%s_%.2f.pickle"%(orbit, 6.00)
    with open(filename, 'rb') as handle:
        data_dict = pickle.load(handle)

    coords_obs = data_dict['coords_obs']
    sc = data_dict['sc_straighter'+sfx]

    distances = coords_obs.distance.to(u.kpc).value
    med_distances.append(np.median(distances))
    trim_new = data_dict['trim_new'+sfx]
    unbound = data_dict['unbound']

    nonrem = data_dict['nonrem']
    phot = data_dict['catalog_photometry'+sfx]

    phot_iso = data_dict['phot']

    cut = data_dict['cut_for_catalog_photometry'+sfx]
    matched_flag = data_dict['matched_to_catalog_photometry'+sfx]
    matched_full = np.full(len(cut), fill_value = False)
    matched_full[cut] = matched_flag


    noise_dict = data_dict['noise_catalog_photometry'+sfx]
    noise_dict['v_gsr'] = noise_dict['v_gsr_'+noise] #<-- ie tack on 'via' or 'desi to get the right key here
    rverr = noise_dict['rverr_'+noise] #<-- this is the RV uncertainty. the above is the noise sampled from a gaussian of width rverr_[survey]
    alive = data_dict['alive']
    acceptable_G = data_dict['acceptable_G']
    cf = (u.microarcsecond/u.yr).to(u.mas/u.yr)
    good_pm = noise_dict['pm_err_gaia']*cf < 0.5 #<-- mas/yr. avoid crazy cocoon inflation due to bad gaia pms.


    if noise=='desi':
        N_jarvis = 679 #<-- length of jarvis catalog. 
        particle_IDs = np.arange(0, len(phot['mG']), 1).astype(int)

        
        bright_ordering = np.argsort(phot['mG']) #<-- nans have moved to the end
        ordered_IDs = particle_IDs[bright_ordering]
        used_IDs = ordered_IDs[:N_jarvis]

        top_N_jarvis = np.isin(particle_IDs, used_IDs)

        good_RV = (top_N_jarvis) & (rverr<10.) #km/s       

    else:
        good_RV = rverr<5.0 #km/s #<--- pretty happy with how this mag distribution comes out...


    use = cut & matched_full & good_pm & good_RV #& nonrem

    # lw=3
    # ax.hist(distances[use], bins=40,# bins=bins, 
    #         color=ccc[ii], histtype='step', lw=lw,
    #         label=orbit)
    # if orbit=='gd1'
    x = coords_obs.phi1
    if orbit=='gd1':
        x+=40*u.degree

    axx.scatter(x[unbound & trim_new], distances[unbound & trim_new], c='0.9', s=5, rasterized=True)
    axx.scatter(x[use], distances[use], c=ccc[ii], s=20, rasterized=True,
                edgecolor='k', lw=0.2)
    axx.set_xlabel(r'$\phi_1~[\degree]$')
    axx.set_ylabel(r'Distance [kpc]')
    # axx.set_title(noise)
    # axx.set_ylim(0, 20)

    # fig_cmd, cmd = plt.subplots()
    # cmd.scatter(phot_iso['BP_RP'][unbound & trim_new], phot_iso['mG'][unbound & trim_new], c='0.9', zorder=0,
    #             label='all isochrone photometry', rasterized=True)
    # cmd.scatter(phot_iso['BP_RP'][use], phot_iso['mG'][use], c='0.7', label='matched stars', rasterized=True)
    # cmd.scatter(phot['BP_RP'][use], phot['mG'][use], c=ccc[ii], label='catalog photometry', rasterized=True)
    # cmd.set_xlabel(r'$G_{B_P}-G_{R_P}$')
    # cmd.set_ylabel(r'$G$')
    # cmd.invert_yaxis()
    # cmd.set_title(orbit_titles[ii])
    # # cmd.set_lim(bottom=)
    # if orbit=='gd1':
    #     cmd.legend(loc='lower left')
    # fig_cmd.savefig(repo_path+"/plots/photometry_explainer/%s.pdf"%orbit, dpi=300, bbox_inches='tight')

# figg.savefig(repo_path+"/plots/photometry_explainer/distance_phi1.pdf", dpi=300, bbox_inches='tight')


med_distances = np.array(med_distances)
present_rs = np.array(present_rs)
pericenters_kpc = np.array(pericenters_kpc)
apocenters_kpc = np.array(apocenters_kpc)
### roughly ~amount of the way through orbit
orbital_phases = (present_rs - pericenters_kpc) / (apocenters_kpc - pericenters_kpc)
orbits = np.array(orbits)

eccentricities = (apocenters_kpc - pericenters_kpc) / (apocenters_kpc + pericenters_kpc)


# tt = Table.read(table_path+case_name+".fits", format="fits")
# okay_mean = np.abs(tt['M_c']) < 0.5*np.asarray(tt['S_c'])
# rvir_cut = tt['Rvir0']<6.

# reordered = np.argsort(pericenters_kpc )

# ccc = cc[1:]
# fig, axs = plt.subplots(2,3,figsize=[21,14])
# plt.subplots_adjust(wspace=0.2, hspace=0.2)
# for ii, orbit in enumerate(tqdm(orbits[reordered])):
#     orbsel = tt['orbit'] == orbit

#     selection=orbsel & np.logical_and.reduce(okay_mean.T) #& rvir_cut

#     f_cocoons_this_orbit = tt['f_cocoon'][selection]

    
#     x = tt['Rvir0'][selection]
#     axs[0,0].plot(x, f_cocoons_this_orbit, 
#                 label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[reordered][ii],
#                 # label = orbit+r'; $\varphi_{\rm orb} =%.2f$'%orbital_phases[reordered][ii],
#                 # label = orbit+r'; $e=%.2f$'%eccentricities[reordered][ii],
#                 marker='o', color=ccc[ii], markersize=10)


#     ### cocoon ! ! !
#     phi2_dispersions_this_orbit = tt['S_c'][:,0][selection] #<-- TODO: translate back to angle from distance. 
#     # phi2_dispersions_this_orbit = (dphi2_dispersions_this_orbit / med_distances[ii]) * u.radian.to(u.degree)
    
#     vgsr_dispersions_this_orbit = tt['S_c'][:,-1][selection]

#     axs[0,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')
#     axs[0,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')

#     ### thin ! ! !
#     phi2_dispersions_this_orbit = tt['S_ts'][:,0][selection] #<-- TODO: translate back to angle from distance. 
#     vgsr_dispersions_this_orbit = tt['S_ts'][:,-1][selection]
#     axs[1,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--')
#     axs[1,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--')

# for ax in np.concatenate([axs[0], axs[1]]):
#     ax.set_ylim(bottom=0)
#     ax.set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
#     ax.minorticks_off()
#     ax.set_xticks([.75, 1.5, 3., 6.])



# axs[0,0].set_ylabel(r'$f_{\rm cocoon}$')
# axs[0,0].set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
# # axs[0,0].set_ylim(0.0, 0.2)

# axs[0,1].set_ylabel(r'$\sigma_{\phi_2, \rm cocoon}~[\degree]$')
# axs[1,1].set_ylabel(r'$\sigma_{\phi_2, \rm thin}~[\degree]$')

# axs[0,2].set_ylabel(r'$\sigma_{v_{\rm GSR, cocoon}}~[\rm km~s^{-1}]$')
# axs[1,2].set_ylabel(r'$\sigma_{v_{\rm GSR, thin}}~[\rm km~s^{-1}]$')

# axs[0,0].legend(loc='upper center', bbox_to_anchor=[0.5,-0.25], fontsize=25)

# axs[1,0].remove()

# plot_filename = "summary_"+case_name
# plt.savefig(repo_path+"/plots/"+plot_filename, dpi=300, bbox_inches='tight')
# %%
#----------------------------#
# !!!!!! THE MAIN RESULT PLOT W MANY PANELS 
#----------------------------#

### getting jarvis values to paste on:
jarvis_svgsr_ts = 2.49
jarvis_svgsr_ts_err = 0.28
jarvis_sphi2_ts = 0.23
jarvis_sphi2_ts_err = 0.01

jarvis_svgsr_c = 6.13
jarvis_svgsr_c_err = 0.75
jarvis_sphi2_c = 2.18
jarvis_sphi2_c_err = 0.17

jarvis_fcocoon = 0.3 #<-- really it's more like 0.33 once you subtract the background but whatever. 
jarvis_fcocoon_err = 0.02 

caseI = 'noiseless_CoM'
caseII = 'noiseless_binaries'


tt = Table.read(table_path+caseI+".fits", format="fits")
okay_mean = np.abs(tt['M_c']) < 0.5*np.asarray(tt['S_c'])
# rvir_cut = tt['Rvir0']<6.



tt2 = Table.read(table_path+caseII+".fits", format='fits')
okay_mean2 = np.abs(tt2['M_c']) < 0.5*np.asarray(tt['S_c'])


reordered = np.argsort(pericenters_kpc)

ccc = cc[1:]
fig, axs = plt.subplots(2,3,figsize=[21,14])


axs[0,0].axhspan(jarvis_fcocoon-jarvis_fcocoon_err, jarvis_fcocoon+jarvis_fcocoon_err, 
                 color='k', alpha=0.1, label = "Values from Jarvis+26")
axs[0,1].axhspan(jarvis_sphi2_c-jarvis_sphi2_c_err, jarvis_sphi2_c+jarvis_sphi2_c_err, 
                 color='k', alpha=0.1)
axs[1,1].axhspan(jarvis_sphi2_ts-jarvis_sphi2_ts_err, jarvis_sphi2_ts+jarvis_sphi2_ts_err, 
                 color='k', alpha=0.1)
axs[0,2].axhspan(jarvis_svgsr_c-jarvis_svgsr_c_err, jarvis_svgsr_c+jarvis_svgsr_c_err, 
                 color='k', alpha=0.1)
axs[1,2].axhspan(jarvis_svgsr_ts-jarvis_svgsr_ts_err, jarvis_svgsr_ts+jarvis_svgsr_ts_err, 
                 color='k', alpha=0.1)

plt.subplots_adjust(wspace=0.2, hspace=0.2)
for ii, orbit in enumerate(tqdm(orbits[reordered])):
    orbsel = tt['orbit'] == orbit
    orbsel2 = tt2['orbit'] == orbit

    selection=orbsel #& np.logical_and.reduce(okay_mean.T) #& rvir_cut
    selection2 = orbsel2 #& np.logical_and.reduce(okay_mean2.T)


    f_cocoons_this_orbit = tt['f_cocoon'][selection]
    f_cocoons_this_orbit2 = tt2['f_cocoon'][selection2]
    
    x = tt['Rvir0'][selection]
    x2 = tt2['Rvir0'][selection2]
    axs[0,0].plot(x, f_cocoons_this_orbit, 
                label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[reordered][ii],
                # label = orbit+r'; $\varphi_{\rm orb} =%.2f$'%orbital_phases[reordered][ii],
                # label = orbit+r'; $e=%.2f$'%eccentricities[reordered][ii],
                marker='o', color=ccc[ii], markersize=10, zorder=0)

    axs[0,0].plot(x2, f_cocoons_this_orbit2, 
                # label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[reordered][ii],
                # label = orbit+r'; $\varphi_{\rm orb} =%.2f$'%orbital_phases[reordered][ii],
                # label = orbit+r'; $e=%.2f$'%eccentricities[reordered][ii],
                marker='o', color=ccc[ii], markersize=10,
                ls=':', zorder=0)


    ### cocoon ! ! !
    phi2_dispersions_this_orbit = tt['S_c'][:,0][selection] #<-- TODO: translate back to angle from distance. 
    phi2_dispersions_this_orbit2 = tt2['S_c'][:,0][selection2]

    vgsr_dispersions_this_orbit = tt['S_c'][:,-1][selection]
    vgsr_dispersions_this_orbit2 = tt2['S_c'][:,-1][selection]

    axs[0,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-', zorder=0)
    axs[0,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-', zorder=0)

    axs[0,1].plot(x2, phi2_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':', zorder=0)
    axs[0,2].plot(x2, vgsr_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':', zorder=0)
    # axs[0,2].set_yscale('log')

    ### thin ! ! !
    phi2_dispersions_this_orbit = tt['S_ts'][:,0][selection] #<-- TODO: translate back to angle from distance. 
    vgsr_dispersions_this_orbit = tt['S_ts'][:,-1][selection]

    phi2_dispersions_this_orbit2 = tt2['S_ts'][:,0][selection2] #<-- TODO: translate back to angle from distance. 
    vgsr_dispersions_this_orbit2 = tt2['S_ts'][:,-1][selection2]

    axs[1,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--', zorder=0)
    axs[1,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--', zorder=0)


    # axs[1,1].plot(x2, phi2_dispersions_this_orbit2, marker='o', color='k', markersize=10, ls='-')
    # axs[1,2].plot(x2, vgsr_dispersions_this_orbit2, marker='o', color='k', markersize=10, ls='-')
    axs[1,1].plot(x2, phi2_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':', zorder=0)
    axs[1,2].plot(x2, vgsr_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':', zorder=0)



axs[0,0].plot([None],[None], c='k', marker='o', markersize=10, label='system CMs')
axs[0,0].plot([None],[None], c='k', marker='o', markersize=10, ls=':', label='with binary orbits')

axs[0,0].set_ylabel(r'$f_{\rm cocoon}$')
axs[0,0].set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
# axs[0,0].set_ylim(0.0, 0.2)

axs[0,1].set_ylabel(r'$\sigma_{\phi_2, \rm cocoon}~[\degree]$')
axs[1,1].set_ylabel(r'$\sigma_{\phi_2, \rm thin}~[\degree]$')

axs[0,2].set_ylabel(r'$\sigma_{v_{\rm GSR, cocoon}}~[\rm km~s^{-1}]$')
axs[1,2].set_ylabel(r'$\sigma_{v_{\rm GSR, thin}}~[\rm km~s^{-1}]$')

for ax in np.concatenate([axs[0], axs[1]]):
    ax.minorticks_off()
    ax.set_xticks([.75, 1.5, 3., 6.])


axs[1,0].remove()


#   ADD A DATA POINT FOR THE FIT WITH DESI-LIKE NOISE AND DOWN-SAMPLING
files = [
    'desi_noise_binaries.fits',
    'desi_noise_CoM.fits',
    'via_noise_binaries.fits',
    'via_noise_CoM.fits'
]
markers = ['v', 'v', 'o', 'o']
labels=[
    'DESI noise with binaries',
    'Desi noise CM',
    'Via Y1 with binaries',
    'Via CM'
]

for ii, f in enumerate(files):
    t_obs = Table.read(table_path+f, format='fits')
    row = t_obs[(t_obs['Rvir0']==0.75) & (t_obs['orbit']=='gd1')] # only making one comparison. 

    axs[0,0].scatter(row['Rvir0'], row['f_cocoon'], 
                label=labels[ii],
                marker=markers[ii], color='none' if 'binaries' in labels[ii] else 'k', 
                edgecolor='k', lw=1.5, s=500)
    axs[0,1].scatter(row['Rvir0'], row['S_c'][:,0], marker=markers[ii], color='none' if 'binaries' in labels[ii] else 'k', 
                            edgecolor='k', lw=1., s=500)
    axs[0,2].scatter(row['Rvir0'], row['S_c'][:,-1], marker=markers[ii], color='none' if 'binaries' in labels[ii] else 'k', 
                            edgecolor='k', lw=1., s=500)
    axs[1,1].scatter(row['Rvir0'], row['S_ts'][:,0], marker=markers[ii], color='none' if 'binaries' in labels[ii] else 'k',
                             edgecolor='k', lw=1., s=500)
    axs[1,2].scatter(row['Rvir0'], row['S_ts'][:,-1], marker=markers[ii], color='none' if 'binaries' in labels[ii] else 'k', 
                            edgecolor='k', lw=1., s=500)


axs[0,0].legend(loc='upper center', bbox_to_anchor=[0.5,-0.25], fontsize=20)
for ax in np.concatenate([axs[0], axs[1]]):
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')

plt.savefig("plots/summary_cocoon_ts_properties.pdf", dpi=300, bbox_inches='tight')
# %%
#### not clearly resolving the thin stream. are there 
#   significant differences in the rv err distribution between
#   my data and the jarvis catalog? 
orbit='gd1'
noise='desi'
rvir_index=0

##### let's open all of the dictionaries also to get like a median distance. use the most diffuse guy.
filename = datapath+"%s_%.2f.pickle"%(orbit, rvirs[rvir_index])
with open(filename, 'rb') as handle:
    data_dict = pickle.load(handle)

coords_obs = data_dict['coords_obs']
sc = data_dict['sc_straighter']

distances = coords_obs.distance.to(u.kpc).value
trim_new = data_dict['trim_new']
unbound = data_dict['unbound']

nonrem = data_dict['nonrem']
phot = data_dict['catalog_photometry']

phot_iso = data_dict['phot']

cut = data_dict['cut_for_catalog_photometry']
matched_flag = data_dict['matched_to_catalog_photometry']
matched_full = np.full(len(cut), fill_value = False)
matched_full[cut] = matched_flag


noise_dict = data_dict['noise_catalog_photometry']
noise_dict['v_gsr'] = noise_dict['v_gsr_'+noise] #<-- ie tack on 'via' or 'desi to get the right key here
rverr = noise_dict['rverr_'+noise] #<-- this is the RV uncertainty. the above is the noise sampled from a gaussian of width rverr_[survey]
alive = data_dict['alive']
acceptable_G = data_dict['acceptable_G']
cf = (u.microarcsecond/u.yr).to(u.mas/u.yr)
good_pm = noise_dict['pm_err_gaia']*cf < 0.5 #<-- mas/yr. avoid crazy cocoon inflation due to bad gaia pms. 

N_jarvis = 679 #<-- length of jarvis catalog. 
particle_IDs = np.arange(0, len(phot['mG']), 1).astype(int)
bright_ordering = np.argsort(phot['mG']) #<-- nans have moved to the end
ordered_IDs = particle_IDs[bright_ordering]
used_IDs = ordered_IDs[:N_jarvis]
top_N_jarvis = np.isin(particle_IDs, used_IDs)
good_RV = (top_N_jarvis) & (rverr<10.) #km/s       

use = cut & matched_full & good_pm & good_RV #& nonrem


### jarvis catalog...
jt = Table.read(repo_path+'/data/jarvis26_Table7.fits',format='fits')
plt.hist(jt['V_ERR'], bins=np.arange(0, 10, 0.25), histtype='step', lw=3)
plt.hist(rverr[use], bins=np.arange(0, 10, 0.25), histtype='step', lw=3)
### a little bit not good idk. 
# %%
#-----------------------------------------------------#
#   binary fractions in thin stream vs cocoon         # 
#-----------------------------------------------------#
case = "noiseless_CoM"
tt = Table.read(table_path+case+".fits", format="fits")
# c_labels = ["#CCC9E7", "#2F2F2F"]
c_labels = ['orange','white','midnightblue']
cocoon_cmap = LinearSegmentedColormap.from_list('cocoon_cmap', c_labels)

rvir0s = [0.75, 1.5, 3.0, 6.0]

fig, ax = plt.subplots()

# for ii, orbit in enumerate(tqdm(orbits[reordered])): #<--- this i can do later i think. 

fbins_ts = []
fbls_ts = []
fbus_ts = []

fbins_c = []
fbls_c = []
fbus_c = []

for rvir in rvir0s:
    orbit = 'gd1' #<-- idk. 
    #### extract info about the mixture modele from the saved table:
    print(orbit, rvir)
    row = tt[(tt['orbit']==orbit) & (tt['Rvir0']==rvir)]

    means_fit = np.array([
        row['M_ts'][0], row['M_c'][0]
    ])
    sigmas_fit = np.array([
        row['S_ts'][0], row['S_c'][0]
    ])
    fracs_fit = np.array([
        1-row['f_cocoon'][0], row['f_cocoon'][0]
    ])

    print(sigmas_fit[-1])

    ##### let's open all of the dictionaries also to get like a median distance. use the most diffuse guy.
    filename = datapath+"%s_%.2f.pickle"%(orbit, rvir) #<-- go for the 
    with open(filename, 'rb') as handle:
        data_dict = pickle.load(handle)

    coords_obs = data_dict['coords_obs']
    kk = 'sc_straighter'
    tn = 'trim_new'
    if 'binaries' in case:
        kk+='_primaries'
        tn+='_primaries' #<-- binaries cases are trimmed with binary orbital motion included, as in gmm.py
    sc_straighter = data_dict[kk]

    distances = coords_obs.distance.to(u.kpc).value

    nsingles = data_dict['nsingles']

    # all_IDs = data_dict['IDs']

    unbound, trim_new = data_dict['unbound'], data_dict[tn]

    #### put things into thin stream vs cocoon: 
    keys = ['phi2','v_phi1','v_phi2','v_gsr']
    use = unbound & trim_new
    bound = ~unbound

    x_data = np.column_stack([
        sc_straighter[k][use] for k in keys
    ])

    pcs = np.full(len(bound), fill_value = np.nan)

    ncomponents = 2
    p1, p2 = [gmm.component_membership_probability(x_data, fracs_fit[:-1], means_fit, sigmas_fit, component=cc_i) for cc_i in range(ncomponents)]
    p_thin = p1
    ts = p1>0.5
    p_cocoon = 1-p_thin

    pcs[use] = p_cocoon

    # single/binary info is contained in the ordering. [:nsingles] are single and [nsingles:] are the binaries. 
    n_bound = len(pcs[bound])
    binaries_bound = np.sum(bound.astype(int)[nsingles:])
    fbin_bound = binaries_bound/n_bound

    p_cocoon_singles = pcs[:nsingles]
    p_cocoon_binaries = pcs[nsingles:]

    n_ts = len(pcs[(pcs<0.5) & (~np.isnan(pcs))])
    n_c = len(pcs[(pcs>0.5) & (~np.isnan(pcs))])



    nbin_c = len(p_cocoon_binaries[(p_cocoon_binaries>0.5) & (~np.isnan(p_cocoon_binaries))])
    nbin_ts = len(p_cocoon_binaries[(p_cocoon_binaries<0.5) & (~np.isnan(p_cocoon_binaries))])



    print(orbit)
    fbin_ts = nbin_ts / n_ts
    fbin_c = nbin_c / n_c
    print("\tthin stream binary fraction: ", nbin_ts / n_ts)
    print("\tcocoon binary fraction: ", nbin_c / n_c)
    print("\tbound fbin", fbin_bound)


    p_ts, P_ts = paf.Prob_of_frac(nbin_ts, n_ts)
    lo_ts,m_ts,up_ts = paf.percentile(p_ts, P_ts)
    xerr_l, xerr_u = m_ts-lo_ts, up_ts-m_ts
    fbins_ts.append(m_ts)
    fbls_ts.append(xerr_l)
    fbus_ts.append(xerr_u)

    p_c, P_c = paf.Prob_of_frac(nbin_c, n_c)
    lo_c,m_c,up_c = paf.percentile(p_c, P_c)
    yerr_l, yerr_u = m_c-lo_c, up_c-m_c
    fbins_c.append(m_c)
    fbls_c.append(yerr_l)
    fbus_c.append(yerr_u)


    # ax.errorbar([m_ts], [m_c], c=ccc[ii],
    #             xerr=[[xerr_l], [xerr_u]], yerr=[[yerr_l], [yerr_u]], 
    #             markersize=5, marker='o', markeredgecolor='k', ls='', capsize=3,
    #            label=orbit)#+r"; $f_{\rm bin, bound}=%.2f$"%fbin_bound)

ax.errorbar(rvir0s, fbins_c, yerr=[fbls_c, fbus_c], label='cocoon', capsize=3, marker='o', 
            markeredgecolor='k')
ax.errorbar(rvir0s, fbins_ts, yerr=[fbls_ts, fbus_ts], label='thin stream', capsize=3, marker='o', 
            markeredgecolor='k')
ax.set_xlabel(r'$R_{\rm vir,0}~[\rm pc]$')
ax.set_ylabel(r'$f_{\rm bin}$')
ax.set_title(orbit)
# ax.legend(loc='lower right')
l = [0,0.4]
# ax.set_xlim(l)
# ax.set_ylim(l)
# ax.plot(l,l, c='k', lw=1)
# ax.set_xlabel(r'$f_{\rm bin, ts}$')
# ax.set_ylabel(r'$f_{\rm bin, c}$')
# ax.set_title(r'$R_{\rm vir, 0}=%.2f$'%rvir)

# plt.savefig("plots/binary_fractions.pdf", dpi=300, bbox_inches='tight')
# %%
#--------------------------------------------------------#
#       plots to demonstrate the need for matching       #
#   to catalog photometry and the associated binary      #
#   selection effects.                                   #
#--------------------------------------------------------#
orbit = 'gd1'
rvir_index=0
filename = datapath+"%s_%.2f.pickle"%(orbit, rvirs[rvir_index])
with open(filename, 'rb') as handle:
    data_dict = pickle.load(handle)

phot = data_dict['phot']
catphot = data_dict['catalog_photometry']
alive = data_dict['alive']


ttgd1 = n.load_gaia_catalog(orbit)

bins = np.arange(12, 23, 1)

fig, ax = plt.subplots()#1,2,figsize=[14,7])
N_jarvis = 679

mGphot_cat_argsort = np.argsort(ttgd1['phot_g_mean_mag'])
mGphot_cat = np.sort(ttgd1['phot_g_mean_mag'])

mgPhot_argsort = np.argsort(phot['mG'][alive])
mGphot_brightsort = np.sort(phot['mG'][alive])

# ax = axs[0]
ax.scatter(ttgd1['bp_rp'][mGphot_cat_argsort][:N_jarvis], 
           ttgd1['phot_g_mean_mag'][mGphot_cat_argsort][:N_jarvis], c=ccc[-2],
           label='brightest stars\nin catalog', zorder=1, rasterized=True)
ax.scatter(phot['BP_RP'][alive][mgPhot_argsort][:N_jarvis],
           phot['mG'][alive][mgPhot_argsort][:N_jarvis], c='0.2', 
           label='brightest stars\nfrom isochrone\nmapping', zorder=2, rasterized=True)
ax.scatter(phot['BP_RP'], phot['mG'], c='0.8', 
           label='all stars\nfrom isochrone\nmapping', zorder=0, rasterized=True)


ax.invert_yaxis()
ax.set_xlabel(r'$G_{B_P} - G_{R_P}$')
ax.set_ylabel(r'$G$')
ax.legend(loc='lower left', fontsize=15)
plt.savefig("plots/matching_catalog_mags.pdf", dpi=300, bbox_inches='tight')

# ax = axs[1]


# ax.hist(mGphot_cat[:N_jarvis], bins=bins,
#         # histtype='step',
#         color='k', alpha=0.2, 
#         lw=3, label='brightest stars \nin catalog');



# ax.hist(mGphot_brightsort[:N_jarvis], bins=bins,
#         histtype='step', lw=3, 
#         color='k', 
#         label="brightest stars\nfrom isochrone \nmapping");

# ax.legend(loc='upper left')

# ax.set_xlabel(r'$G$')
# ax.set_ylabel(r'Count')


# %%
#


#### helpers: rebuild exactly what gmm.py fit, for any case, and evaluate a
#   saved mixture model on it. used by the demo panels and figure 1 below.
fit_keys = ['phi2','v_phi1','v_phi2','v_gsr'] #<-- gmm.py's keys; also the order of the M_*/S_* columns

def case_label(noise=None, include_binaries=False, constrain_widths=False):
    '''the gmm.py case_name for a set of flags, e.g. desi_noise_CoM'''
    cd0 = 'noiseless' if noise is None else noise+"_noise"
    cd1 = 'binaries' if include_binaries else 'CoM'
    cd2 = '_constrained' if constrain_widths else ''
    return cd0+"_"+cd1+cd2


def gmm_selection(data_dict, noise=None, include_binaries=False):
    '''
    the stars that went into the gmm.py fit for one case.

    returns use, sc_straighter, noise_dict:
      use           -- full-length bool mask (CoM/luminous ordering)
      sc_straighter -- the straightened coord dict that case fits
      noise_dict    -- the saved catalog-photometry noise draw, with v_gsr,
                       v_phi1, v_phi2 filled in for this survey (None if noiseless)
    '''
    sfx = '_primaries' if include_binaries else '' #<-- coords, trim, and catalog products all follow it
    sc_straighter = data_dict['sc_straighter'+sfx]
    unbound, nonrem = data_dict['unbound'], data_dict['nonrem']

    if noise is None:
        use = unbound & data_dict['trim_new'+sfx]
        return use, sc_straighter, None

    #### catalog-based photometry + noise, same selection as gmm.py
    phot = data_dict['catalog_photometry'+sfx]
    cut = data_dict['cut_for_catalog_photometry'+sfx]
    matched_full = np.full(len(cut), fill_value=False)
    matched_full[cut] = data_dict['matched_to_catalog_photometry'+sfx]

    noise_dict = dict(data_dict['noise_catalog_photometry'+sfx]) #<-- copy, so the survey keys don't leak back into data_dict
    noise_dict['v_gsr'] = noise_dict['v_gsr_'+noise] #<-- ie tack on 'via' or 'desi to get the right key here
    rverr = noise_dict['rverr_'+noise]

    distances = data_dict['coords_obs'].distance
    noise_dict['v_phi1'] = distances.to(u.km).value * (noise_dict['pm_phi1']*u.mas/u.yr).to(u.radian/u.s).value # km/s
    noise_dict['v_phi2'] = distances.to(u.km).value * (noise_dict['pm_phi2']*u.mas/u.yr).to(u.radian/u.s).value # km/s

    cf = (u.microarcsecond/u.yr).to(u.mas/u.yr)
    good_pm = noise_dict['pm_err_gaia']*cf < 0.5 #<-- mas/yr

    if noise in ['via', 'via10hr']:
        good_RV = rverr<5.0 #km/s
    if noise=='desi':
        N_jarvis = 679 #<-- length of jarvis catalog.
        top_N_jarvis = np.zeros(len(phot['mG']), dtype=bool)
        top_N_jarvis[np.argsort(phot['mG'])[:N_jarvis]] = True #<-- nans sort to the end
        good_RV = top_N_jarvis & (rverr<10.) #km/s

    use = cut & matched_full & good_pm & good_RV & nonrem
    return use, sc_straighter, noise_dict


def observed(key, sc_straighter, noise_dict, use):
    '''one coordinate for the `use` stars, with the saved noise draw added if there is one'''
    y = sc_straighter[key][use]
    if noise_dict is not None:
        y = y + noise_dict[key][use]
    return y


def p_cocoon_from_table(x_data, row):
    '''p_cocoon per star under the mixture model saved in one row of a gmm table'''
    means_fit = np.array([row['M_ts'][0], row['M_c'][0]])
    sigmas_fit = np.array([row['S_ts'][0], row['S_c'][0]])
    fracs_fit = np.array([1-row['f_cocoon'][0], row['f_cocoon'][0]])
    p_thin = gmm.component_membership_probability(x_data, fracs_fit[:-1], means_fit, sigmas_fit, component=0)
    return 1-p_thin


def demo_data(orbit, rvir, tt, noise=None, include_binaries=False):
    '''
    rebuild the data vector for one simulation exactly as gmm.py does, and
    evaluate the saved mixture model on it.

    returns phi1, ydict, p_cocoon, row -- phi1 and every ydict entry are
    already cut down to the stars that went into the fit (`use`), and
    ydict holds the plotted keys with the saved noise draw added.
    '''
    row = tt[(tt['orbit']==orbit) & (tt['Rvir0']==rvir)]

    filename = datapath+"%s_%.2f.pickle"%(orbit, rvir)
    with open(filename, 'rb') as handle:
        data_dict = pickle.load(handle)

    use, sc_straighter, noise_dict = gmm_selection(data_dict, noise=noise, include_binaries=include_binaries)

    #### the fit is in transverse velocities; the plot shows proper motions
    x_data = np.column_stack([observed(k, sc_straighter, noise_dict, use) for k in fit_keys])
    p_cocoon = p_cocoon_from_table(x_data, row)

    ydict = {k: observed(k, sc_straighter, noise_dict, use) for k in ['phi2','pm_phi1','pm_phi2','v_gsr']}
    return sc_straighter['phi1'][use], ydict, p_cocoon, row

# %%
#
#--------------------------------------------------------#
#       plots to demonstrate the GMM results             #
#       for individuaal simulations (success and fails)  #
#--------------------------------------------------------#

# c_labels = ["#FBBA72","#F5AE66","#EFA15A","#E9944E","#E38741","#DD7A35","#D76D29","#D1601D","#CA5310"]
noise = 'desi' #<-- None, 'via', 'via10hr', 'desi'
include_binaries = False
case_name = case_label(noise, include_binaries) #<-- e.g. desi_noise_CoM; can't disagree with the flags any more
tt = Table.read(table_path+case_name+".fits", format="fits")
c_labels = ["#CCC9E7", "#2F2F2F"]
# c_labels = ['orange','white','midnightblue']
cocoon_cmap = LinearSegmentedColormap.from_list('cocoon_cmap', c_labels)

orbit = 'gd1'
rvir_indices = [0] #<-- any subset of range(len(rvirs)), one column each

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']
key_labels = [
    r'$\phi_2~[\degree]$',
    r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$',
    r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$',
    r'$v_{\rm GSR}~[\rm km~s^{-1}]$'
]

# the noisy samples are ~1e2-1e3 stars, the noiseless ones ~1e4: bigger, outlined points for "observed" data
point_kw = dict(s=5) if noise is None else dict(s=50, edgecolor='k', lw=.5)

ncols = len(rvir_indices)
fig = plt.figure(figsize=[5*ncols+1, 8])
# outer grid: the panels, plus a thin column on the far right for the colorbar
outer = fig.add_gridspec(1, 2, width_ratios=[ncols, 0.05], wspace=0.03)
axs = outer[0].subgridspec(len(keys), ncols, hspace=0.03, wspace=0.03).subplots(
    sharex=True, sharey='row', squeeze=False)
cax = fig.add_subplot(outer[1])

norm_pc = Normalize(vmin=0, vmax=1)
ylims = np.zeros(len(keys)) #<-- per row, the widest of the columns

for ii, rvir_index in enumerate(rvir_indices):
    rvir = rvirs[rvir_index]
    phi1, ydict, p_cocoon, row = demo_data(orbit, rvir, tt, noise=noise, include_binaries=include_binaries)
    order = np.argsort(p_cocoon) # plot cocoon on top.
    cocoon_selection = p_cocoon>0.5

    for jj, key in enumerate(keys):
        y = ydict[key]
        ax = axs[jj, ii]
        ax.scatter(phi1[order], y[order],
                   c=p_cocoon[order], cmap=cocoon_cmap, norm=norm_pc,
                   rasterized=True, **point_kw)

        # 3 sigma of the cocoon stars (or of everything, if there's no cocoon)
        ysel = y[cocoon_selection] if cocoon_selection.sum()>1 else y
        ylims[jj] = max(ylims[jj], 3*np.std(ysel))

    if noise is None:
        axs[0, ii].set_title(r'$R_{\rm vir,0}=%.2f$ pc, $f_{\rm c}=%.2f$'%(rvir, row['f_cocoon'][0]), fontsize=13)

    if noise is not None:
        axs[0, ii].set_title(
            r'$f_{\rm c}=%.2f, \sigma_{\phi_2, c} = %.2f~\degree, \sigma_{v_{\rm GSR}, c} = %.2f~\mathrm{km~s^{-1}} $'%(row['f_cocoon'][0], row['S_c'][:,0][0], row['S_c'][:,-1][0]))


    axs[-1, ii].set_xlabel(r'$\phi_1~[\degree]$')

for jj in range(len(keys)):
    axs[jj, 0].set_ylim(-ylims[jj], ylims[jj]) #<-- shared across the row
    axs[jj, 0].set_ylabel(key_labels[jj], fontsize=15)
axs[0, 0].set_xlim(orbit_lim_map[orbit]) #<-- sharex, so this sets every panel

cb = fig.colorbar(ScalarMappable(norm=norm_pc, cmap=cocoon_cmap), cax=cax)
cb.set_label(r'$p_{\rm cocoon}$', fontsize=15)

# keep the old filename for the noiseless CoM panels; every other case gets its own file
plot_name = "demo_cocoon_separation_rvir0" if case_name=='noiseless_CoM' else "demo_cocoon_separation_"+case_name
if orbit!='gd1':
    plot_name += "_"+orbit
plt.savefig("plots/"+plot_name+".pdf", dpi=300, bbox_inches='tight')
# %%
#--------------------------------------------------------#
#       FIGURE 1: every rvir0 = 0.75 pc stream at        #
#       present day, galactocentric z vs x.              #
#       left: thin streams only (cocoons removed);       #
#       right: everything unbound from the progenitor.   #
#--------------------------------------------------------#
fig1_rvir = 0.75
fig1_noise, fig1_include_binaries = None, False #<-- the gmm.py fit that decides who's in the cocoon
# right panel: True = every unbound star, including the ones trim_new dropped before the fit
#   (phi1 tails, distance outliers, early escapers), so left vs right is cocoon + trim.
#   False = only the stars the gmm saw, so left vs right is exactly the cocoon.
fig1_include_trimmed = False
tt_fig1 = Table.read(table_path+case_label(fig1_noise, fig1_include_binaries)+".fits", format="fits")

orbit_title_map = {o:t for o, t in zip(orbits, orbit_titles)}
reordered = np.argsort(pericenters_kpc) #<-- same as the main result plot

#### load everything first, then plot.
#   loop in pericenter order so ccc[ii] is the same color per orbit as in the summary plots.
fig1_x, fig1_y, fig1_z = [], [], [] # galactocentric, kpc
fig1_n_sim = [] # which stream each star belongs to (index into orbits[reordered])
fig1_thin = [] # in the gmm fit AND p_cocoon < 0.5
fig1_fit = [] # in the gmm fit at all

for ii, orbit in enumerate(tqdm(orbits[reordered])):
    filename = datapath+"%s_%.2f.pickle"%(orbit, fig1_rvir)
    with open(filename, 'rb') as handle:
        data_dict = pickle.load(handle)

    use, sc_fig1, noise_fig1 = gmm_selection(data_dict, noise=fig1_noise, include_binaries=fig1_include_binaries)
    row = tt_fig1[(tt_fig1['orbit']==orbit) & (tt_fig1['Rvir0']==fig1_rvir)]
    x_data = np.column_stack([observed(k, sc_fig1, noise_fig1, use) for k in fit_keys])

    p_cocoon = np.full(len(use), fill_value=np.nan) #<-- nan for stars that weren't in the fit (bound, or outside the trim)
    p_cocoon[use] = p_cocoon_from_table(x_data, row)
    thin = use & (p_cocoon<0.5) #<-- nan < 0.5 is False

    #### the progenitor (2 rtid) is dropped from both panels; everything else unbound goes in the right one
    unbound = data_dict['unbound']
    x, y, z = data_dict['CoM']['pos'].to(u.kpc).value[unbound].T

    fig1_x.append(x)
    fig1_y.append(y)
    fig1_z.append(z)
    fig1_n_sim.append(np.full(len(x), fill_value=ii))
    fig1_thin.append(thin[unbound])
    fig1_fit.append(use[unbound])

    print("%s: %i unbound, %i in the fit, %i thin stream"%(orbit, unbound.sum(), use.sum(), thin.sum()))

fig1_x, fig1_y, fig1_z = np.concatenate(fig1_x), np.concatenate(fig1_y), np.concatenate(fig1_z)
fig1_n_sim = np.concatenate(fig1_n_sim)
fig1_thin = np.concatenate(fig1_thin)
fig1_fit = np.concatenate(fig1_fit)


n_orbits = len(orbits)
cmap_fig1 = mcolors.ListedColormap(ccc[:n_orbits])
# discrete norm: one color band per orbit, boundaries on the half-integers
norm_fig1 = mcolors.BoundaryNorm(np.arange(n_orbits+1)-0.5, cmap_fig1.N)

fig, axs = plt.subplots(1, 2, figsize=[16, 8], sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.05)


for ax in axs:
    for width, height in [(2, 2), (3, 2), (0.5, 0.5), (30, 0.5)]:
        ax.add_patch(Ellipse(xy=[0, 0], width=width, height=height, angle=0,
                            facecolor='none', edgecolor='black', linewidth=1))

panel_selections = [fig1_thin, np.ones(len(fig1_thin), dtype=bool) if fig1_include_trimmed else fig1_fit]
panel_titles = ['Thin stream only', 'All stars']
for ax, sel, title in zip(axs, panel_selections, panel_titles):
    # order = np.argsort(fig1_y[sel]) #<-- most negative y drawn first (bottom), most positive y last (top)
    order = np.argsort(fig1_x[sel])
    ax.scatter(fig1_y[sel][order], fig1_z[sel][order],
               c=fig1_n_sim[sel][order], cmap=cmap_fig1, norm=norm_fig1,
               s=1, rasterized=True)
    # ax.set_title(title)
    ax.set_xlabel(r'$y~[\rm kpc]$')
    ax.set_aspect('equal')

axs[0].set_xlim(-35, 35)
axs[0].set_ylim(-35, 35)
axs[0].set_ylabel(r'$z~[\rm kpc]$')



handles = [Line2D([], [], ls='', marker='o', markersize=8, color=ccc[ii], label=orbit_title_map[orbit])
           for ii, orbit in enumerate(orbits[reordered])]
axs[1].legend(handles=handles, loc='lower right', fontsize=15)

plt.savefig("plots/fig1_motiv_trim.pdf", dpi=300, bbox_inches='tight')
# %%
