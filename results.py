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

# %%

##### finish later. 
# def get_p_cocoon(data_dict, use, row, noise ):
#     keys = ['phi2','v_phi1','v_phi2','v_gsr']
#     sc_straighter = data_dict['sc_straighter']
#     row = tt[(tt['orbit']==orbit) & (tt['Rvir0']==rvir)]

#     means_fit = np.array([
#         row['M_ts'][0], row['M_c'][0]
#     ])
#     sigmas_fit = np.array([
#         row['S_ts'][0], row['S_c'][0]
#     ])
#     fracs_fit = np.array([
#         1-row['f_cocoon'][0], row['f_cocoon'][0]
#     ])
#     keys = ['phi2','v_phi1','v_phi2','v_gsr']
#     x_data = np.column_stack([
#         sc_straighter[k][use] for k in keys
#     ])

#     ncomponents = 2
#     p1, p2 = [gmm.component_membership_probability(x_data, fracs_fit[:-1], means_fit, sigmas_fit, component=cc_i) for cc_i in range(ncomponents)]
#     p_thin = p1
#     ts = p1>0.5
#     p_cocoon = 1-p_thin
#     return p_cocoon

# def get_cocoon_flag(data_dict, use, means, sigmas):

#     return p_cocoon>0.5

# %%
datapath = '/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
table_path = "/n/home02/amphillips/p27_nbody/data/gmm_tables/"


#-----------------------------------------------#
#   opening all of the data and checking out    #
#   the observability given distances to all    #
#   the streams, which noise i use, etc.        #
#-----------------------------------------------#


make_plots=True
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
    (-80, 10),
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

for ii, orbit in enumerate(tqdm(orbits)): #<--- this i can do later i think. 
    init_displacement = init_displacements[ii]
    orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
    peri = orbit_obj.pericenter()
    apo = orbit_obj.apocenter()
    pericenters_kpc.append(peri.to(u.kpc).value)
    apocenters_kpc.append(apo.to(u.kpc).value)

    x,y,z = init_displacement[:3]
    r = np.sqrt(x**2 + y**2 + z**2)
    present_rs.append(r) # kpc

    ##### let's open all of the dictionaries also to get like a median distance. use the most diffuse guy.
    filename = datapath+"%s_%.2f.pickle"%(orbit, 6.00)
    with open(filename, 'rb') as handle:
        data_dict = pickle.load(handle)

    coords_obs = data_dict['coords_obs']
    sc = data_dict['sc_straighter']

    distances = coords_obs.distance.to(u.kpc).value
    med_distances.append(np.median(distances))
    trim_new = data_dict['trim_new']
    unbound = data_dict['unbound']

    nonrem = data_dict['nonrem']
    phot = data_dict['phot']
    noise_dict = data_dict['noise']
    noise_dict['v_gsr'] = noise_dict['v_gsr_'+noise] #<-- ie tack on 'via' or 'desi to get the right key here
    rverr = noise_dict['rverr_'+noise] #<-- this is the RV uncertainty. the above is the noise sampled from a gaussian of width rverr_[survey]
    alive = data_dict['alive']
    acceptable_G = data_dict['acceptable_G']
    cf = (u.microarcsecond/u.yr).to(u.mas/u.yr)
    good_pm = noise_dict['pm_err_gaia']*cf < 0.5 #<-- mas/yr. avoid crazy cocoon inflation due to bad gaia pms. 

    if noise=='via':
        good_RV = rverr<1.0 #km/s #<--- pretty happy with how this mag distribution comes out...
    else:
        good_RV = rverr<10. #km/s

    good_RV5 = rverr<5.0



    use = trim_new & unbound & nonrem & alive & acceptable_G & nonrem & good_RV
    use5 = trim_new & unbound & nonrem & alive & acceptable_G & nonrem & good_RV5

    # lw=3
    # ax.hist(distances[use], bins=40,# bins=bins, 
    #         color=ccc[ii], histtype='step', lw=lw,
    #         label=orbit)
    # if orbit=='gd1'
    x = coords_obs.phi1
    if orbit=='gd1':
        x+=40*u.degree

    axx.scatter(x[unbound & trim_new], distances[unbound & trim_new], c='0.9', s=5, rasterized=True)
    axx.scatter(x[use5], distances[use5], c='0.7', s=5, rasterized=True)
    axx.scatter(x[use], distances[use], c=ccc[ii], s=20, rasterized=True,
                edgecolor='k', lw=0.2)
    axx.set_xlabel(r'$\phi_1~[\degree]$')
    axx.set_ylabel(r'Distance [kpc]')
    axx.set_title(noise)
    # axx.set_ylim(0, 20)

    fig, cmd = plt.subplots()
    cmd.scatter(phot['BP_RP'][use5], phot['mG'][use5], c='0.7' )# ccc[ii])
    cmd.scatter(phot['BP_RP'][use], phot['mG'][use], c=ccc[ii])

    cmd.scatter(phot['BP_RP'][unbound & trim_new], phot['mG'][unbound & trim_new], c='0.9', zorder=0)
    cmd.invert_yaxis()
    cmd.set_title(orbit)

# ax.legend()


med_distances = np.array(med_distances)
present_rs = np.array(present_rs)
pericenters_kpc = np.array(pericenters_kpc)
apocenters_kpc = np.array(apocenters_kpc)
### roughly ~amount of the way through orbit
orbital_phases = (present_rs - pericenters_kpc) / (apocenters_kpc - pericenters_kpc)
orbits = np.array(orbits)

eccentricities = (apocenters_kpc - pericenters_kpc) / (apocenters_kpc + pericenters_kpc)






# %%
tt = Table.read(table_path+case_name+".fits", format="fits")
okay_mean = np.abs(tt['M_c']) < 0.5*np.asarray(tt['S_c'])
rvir_cut = tt['Rvir0']<6.

reordered = np.argsort(pericenters_kpc )

ccc = cc[1:]
fig, axs = plt.subplots(2,3,figsize=[21,14])
plt.subplots_adjust(wspace=0.2, hspace=0.2)
for ii, orbit in enumerate(tqdm(orbits[reordered])):
    orbsel = tt['orbit'] == orbit

    selection=orbsel & np.logical_and.reduce(okay_mean.T) #& rvir_cut

    f_cocoons_this_orbit = tt['f_cocoon'][selection]

    
    x = tt['Rvir0'][selection]
    axs[0,0].plot(x, f_cocoons_this_orbit, 
                label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[reordered][ii],
                # label = orbit+r'; $\varphi_{\rm orb} =%.2f$'%orbital_phases[reordered][ii],
                # label = orbit+r'; $e=%.2f$'%eccentricities[reordered][ii],
                marker='o', color=ccc[ii], markersize=10)


    ### cocoon ! ! !
    phi2_dispersions_this_orbit = tt['S_c'][:,0][selection] #<-- TODO: translate back to angle from distance. 
    # phi2_dispersions_this_orbit = (dphi2_dispersions_this_orbit / med_distances[ii]) * u.radian.to(u.degree)
    
    vgsr_dispersions_this_orbit = tt['S_c'][:,-1][selection]

    axs[0,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')
    axs[0,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')

    ### thin ! ! !
    phi2_dispersions_this_orbit = tt['S_ts'][:,0][selection] #<-- TODO: translate back to angle from distance. 
    vgsr_dispersions_this_orbit = tt['S_ts'][:,-1][selection]
    axs[1,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--')
    axs[1,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--')

for ax in np.concatenate([axs[0], axs[1]]):
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
    ax.minorticks_off()
    ax.set_xticks([.75, 1.5, 3., 6.])



axs[0,0].set_ylabel(r'$f_{\rm cocoon}$')
axs[0,0].set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
# axs[0,0].set_ylim(0.0, 0.2)

axs[0,1].set_ylabel(r'$\sigma_{\phi_2, \rm cocoon}~[\degree]$')
axs[1,1].set_ylabel(r'$\sigma_{\phi_2, \rm thin}~[\degree]$')

axs[0,2].set_ylabel(r'$\sigma_{v_{\rm GSR, cocoon}}~[\rm km~s^{-1}]$')
axs[1,2].set_ylabel(r'$\sigma_{v_{\rm GSR, thin}}~[\rm km~s^{-1}]$')

axs[0,0].legend(loc='upper center', bbox_to_anchor=[0.5,-0.25], fontsize=25)

axs[1,0].remove()

plot_filename = "summary_"+case_name
# plt.savefig(repo_path+"/plots/"+plot_filename, dpi=300, bbox_inches='tight')
# %%
#----------------------------#
# a version of this table    #
# that includes the no noise #
# + binary motions gmm fits. #
#----------------------------#

caseI = 'noiseless_CoM'
caseII = 'noiseless_binaries'


tt = Table.read(table_path+caseI+".fits", format="fits")
okay_mean = np.abs(tt['M_c']) < 0.5*np.asarray(tt['S_c'])
# rvir_cut = tt['Rvir0']<6.
tt


tt2 = Table.read(table_path+caseII+".fits", format='fits')
okay_mean2 = np.abs(tt2['M_c']) < 0.5*np.asarray(tt['S_c'])


reordered = np.argsort(pericenters_kpc)

ccc = cc[1:]
fig, axs = plt.subplots(2,3,figsize=[21,14])
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
                marker='o', color=ccc[ii], markersize=10)

    axs[0,0].plot(x2, f_cocoons_this_orbit2, 
                # label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[reordered][ii],
                # label = orbit+r'; $\varphi_{\rm orb} =%.2f$'%orbital_phases[reordered][ii],
                # label = orbit+r'; $e=%.2f$'%eccentricities[reordered][ii],
                marker='o', color=ccc[ii], markersize=10,
                ls=':')


    ### cocoon ! ! !
    phi2_dispersions_this_orbit = tt['S_c'][:,0][selection] #<-- TODO: translate back to angle from distance. 
    phi2_dispersions_this_orbit2 = tt2['S_c'][:,0][selection2]

    vgsr_dispersions_this_orbit = tt['S_c'][:,-1][selection]
    vgsr_dispersions_this_orbit2 = tt2['S_c'][:,-1][selection]

    axs[0,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')
    axs[0,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')

    axs[0,1].plot(x2, phi2_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':')
    axs[0,2].plot(x2, vgsr_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':')


    ### thin ! ! !
    phi2_dispersions_this_orbit = tt['S_ts'][:,0][selection] #<-- TODO: translate back to angle from distance. 
    vgsr_dispersions_this_orbit = tt['S_ts'][:,-1][selection]

    phi2_dispersions_this_orbit2 = tt2['S_ts'][:,0][selection2] #<-- TODO: translate back to angle from distance. 
    vgsr_dispersions_this_orbit2 = tt2['S_ts'][:,-1][selection2]

    axs[1,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--')
    axs[1,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='--')


    # axs[1,1].plot(x2, phi2_dispersions_this_orbit2, marker='o', color='k', markersize=10, ls='-')
    # axs[1,2].plot(x2, vgsr_dispersions_this_orbit2, marker='o', color='k', markersize=10, ls='-')
    axs[1,1].plot(x2, phi2_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':')
    axs[1,2].plot(x2, vgsr_dispersions_this_orbit2, marker='o', color=ccc[ii], markersize=10, ls=':')

for ax in np.concatenate([axs[0], axs[1]]):
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
    ax.minorticks_off()
    ax.set_xticks([.75, 1.5, 3., 6.])

axs[0,0].plot([None],[None], c='k', marker='o', markersize=10, label='system CMs')
axs[0,0].plot([None],[None], c='k', marker='o', markersize=10, ls=':', label='with binary orbits')

axs[0,0].set_ylabel(r'$f_{\rm cocoon}$')
axs[0,0].set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
# axs[0,0].set_ylim(0.0, 0.2)

axs[0,1].set_ylabel(r'$\sigma_{\phi_2, \rm cocoon}~[\degree]$')
axs[1,1].set_ylabel(r'$\sigma_{\phi_2, \rm thin}~[\degree]$')

axs[0,2].set_ylabel(r'$\sigma_{v_{\rm GSR, cocoon}}~[\rm km~s^{-1}]$')
axs[1,2].set_ylabel(r'$\sigma_{v_{\rm GSR, thin}}~[\rm km~s^{-1}]$')

axs[0,0].legend(loc='upper center', bbox_to_anchor=[0.5,-0.25], fontsize=25)

axs[1,0].remove()

plot_filename = "summary_"+case_name

# plt.savefig("plots/summary_CMs_binaries_combined_noiseless.pdf")

# %%
#-----------------------------------------------------#
#   binary fractions in thin stream vs cocoon         # 
#-----------------------------------------------------#
case = "noiseless_binaries"
tt = Table.read(table_path+case+".fits", format="fits")
# c_labels = ["#CCC9E7", "#2F2F2F"]
c_labels = ['orange','white','midnightblue']
cocoon_cmap = LinearSegmentedColormap.from_list('cocoon_cmap', c_labels)

rvir=0.75

fig, ax = plt.subplots()

for ii, orbit in enumerate(tqdm(orbits[reordered])): #<--- this i can do later i think. 
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

    ##### let's open all of the dictionaries also to get like a median distance. use the most diffuse guy.
    filename = datapath+"%s_%.2f.pickle"%(orbit, rvir) #<-- go for the 
    with open(filename, 'rb') as handle:
        data_dict = pickle.load(handle)

    coords_obs = data_dict['coords_obs']
    kk = 'sc_straighter'
    tt = 'trim_new'
    if 'binaries' in case:
        kk+='_primaries'
        tt+='_primaries'
    sc_straighter = data_dict[kk]

    distances = coords_obs.distance.to(u.kpc).value

    nsingles = data_dict['nsingles']

    # all_IDs = data_dict['IDs']

    unbound, trim_new = data_dict['unbound'], data_dict[tt]

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

    p_c, P_c = paf.Prob_of_frac(nbin_c, n_c)
    lo_c,m_c,up_c = paf.percentile(p_c, P_c)
    yerr_l, yerr_u = m_c-lo_c, up_c-m_c



    ax.errorbar([m_ts], [m_c], c=ccc[ii],
                xerr=[[xerr_l], [xerr_u]], yerr=[[yerr_l], [yerr_u]], 
                markersize=5, marker='o', markeredgecolor='k', ls='', capsize=3,
               label=orbit+r"; $f_{\rm bin, bound}=%.2f$"%fbin_bound)

ax.legend(loc='lower right')
l = [0,0.4]
ax.set_xlim(l)
ax.set_ylim(l)
ax.plot(l,l, c='k', lw=1)
ax.set_xlabel(r'$f_{\rm bin, ts}$')
ax.set_ylabel(r'$f_{\rm bin, c}$')

# plt.savefig("plots/binary_fractions.pdf", dpi=300, bbox_inches='tight')
# %%
#
#
#--------------------------------------------------------#
#       plots to demonstrate the GMM results             #
#       for individuaal simulations (success and fails)  #
#--------------------------------------------------------#
# %%
# c_labels = ["#FBBA72","#F5AE66","#EFA15A","#E9944E","#E38741","#DD7A35","#D76D29","#D1601D","#CA5310"]
include_binaries=True
case_name = 'desi_noise_CoM'
tt = Table.read(table_path+case_name+".fits", format="fits")
# c_labels = ["#CCC9E7", "#2F2F2F"]
c_labels = ['orange','white','midnightblue']
cocoon_cmap = LinearSegmentedColormap.from_list('cocoon_cmap', c_labels)

orbit = 'gd1'
rvir_index=0
rvir = rvirs[rvir_index]


#### extract info about the mixture modele from the saved table:
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


#### open the N-body data
filename = datapath+"%s_%.2f.pickle"%(orbit, rvir)
with open(filename, 'rb') as handle:
    data_dict = pickle.load(handle)


if include_binaries==False:
    sc_straighter = data_dict['sc_straighter'] #<-- dictionary
if include_binaries==True:
    sc_straighter = data_dict['sc_straighter_primaries'] #<-- dictionary

unbound = data_dict['unbound']
trim_new = data_dict['trim_new']

keys = ['phi2','v_phi1','v_phi2','v_gsr']


if noise is not None:
    nonrem = data_dict['nonrem']


    phot = data_dict['phot']
    noise_dict = data_dict['noise']
    noise_dict['v_gsr'] = noise_dict['v_gsr_'+noise] #<-- ie tack on 'via' or 'desi to get the right key here

    rverr = noise_dict['rverr_'+noise] #<-- this is the RV uncertainty. the above is the noise sampled from a gaussian of width rverr_[survey]
    alive = data_dict['alive']
    acceptable_G = data_dict['acceptable_G']

    cf = (u.microarcsecond/u.yr).to(u.mas/u.yr)
    good_pm = noise_dict['pm_err_gaia']*cf < 0.5 #<-- mas/yr. avoid crazy cocoon inflation due to bad gaia pms. 

    distances = data_dict['coords_obs'].distance

    noise_dict['v_phi1'] = distances.to(u.km).value * (noise_dict['pm_phi1']*u.mas/u.yr).to(u.radian/u.s).value # km/s
    noise_dict['v_phi2'] = distances.to(u.km).value * (noise_dict['pm_phi2']*u.mas/u.yr).to(u.radian/u.s).value # km/s
    noise_dict['d_phi2'] = distances.to(u.kpc).value * (noise_dict['phi2']*u.degree).to(u.radian).value


    if noise=='via':
        good_RV = rverr<1.0 #km/s #<--- pretty happy with how this mag distribution comes out...
        use = trim_new & unbound & alive & nonrem & acceptable_G & good_pm & good_RV
    else:
        good_RV = rverr<10. #km/s
        use = trim_new & unbound & alive & nonrem & good_pm & good_RV#<-- no acceptable G range for DESI errors. 
    
    # keys = ['phi2','v_phi1','v_phi2','v_gsr_'+noise]

    x_data = np.column_stack(
        [sc_straighter[k][use]+noise_dict[k][use] for k in keys]
    )


else:
    use = unbound & trim_new
    x_data = np.column_stack([
        sc_straighter[k][use] for k in keys
    ])



ncomponents = 2
p1, p2 = [gmm.component_membership_probability(x_data, fracs_fit[:-1], means_fit, sigmas_fit, component=cc_i) for cc_i in range(ncomponents)]
p_thin = p1
ts = p1>0.5
p_cocoon = 1-p_thin
order = np.argsort(p_cocoon)

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']
fig, axs = plt.subplots(len(keys), 2, figsize=[10, 10], width_ratios = [4,1])

plt.subplots_adjust(hspace=0.03, wspace=0.03)

key_labels = [
    r'$\phi_2~[\degree]$',
    r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$',
    r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$',
    r'$v_{\rm GSR}~[\rm km~s^{-1}]$'
]


for jj, key in enumerate(keys):
    cocoon_selection = p_thin<0.5

    ydata = sc_straighter[key][use][cocoon_selection]
    if noise is not None:
        ydata += noise_dict[key][use][cocoon_selection]

    std = np.std(ydata)
    cut = 3*std #cocoon_sigmas[jj]

    ax = axs[jj,0]

    y = sc_straighter[key][use][order]
    if noise is not None:
        y+=noise_dict[key][use][order]

    ax.scatter(sc_straighter['phi1'][use][order],  # plot cocoon on top. 
            y, # plot cocoon on top. 
            # x_data[:,ii],
                c=p_cocoon[order], s=50, edgecolor='k', lw=.5,
                cmap=cocoon_cmap,
                rasterized=True) 
    ax.set_ylim(-cut,cut)
    # ax.set_xlim(orbit_lim_map[orbit])

    # ax.set_ylim(-3*cut, 3*cut)
    ax.set_ylabel(key_labels[jj], fontsize=15)


    ax = axs[jj,1]
    bins = np.linspace(-cut, cut, 50)



    ax.hist(y[~cocoon_selection], 
            alpha=1., density=False, weights = np.zeros_like(sc_straighter[key][use][~cocoon_selection])+1/sc_straighter[key][use][~cocoon_selection].size, 
            color=c_labels[0],orientation='horizontal',
            bins=bins)
    ax.hist(y[cocoon_selection],
            histtype='step', density=False, weights = np.zeros_like(sc_straighter[key][use][cocoon_selection])+1/sc_straighter[key][use][cocoon_selection].size, 
            lw=2, 
            color=c_labels[-1],orientation='horizontal',
            bins=bins)
    ax.set_ylim(-cut,cut)



    # too annoying to get the limits to work out. being unrigorous for now...
    # ax.set_xticks([])
    # ax.set_xticklabels([])
    ax.set_xlim(0, 0.4)
    ax.set_yticks([])
    ax.set_yticklabels([])

    if jj<3:
        # print("REMOVING TICK LABLES>>>>>")
        axs[jj,0].set_xticklabels([])
        axs[jj,1].set_xticklabels([])


axs[-1,0].set_xlabel(r'$\phi_1~[\degree]$')
axs[-1,1].set_xlabel(r'fraction in bin', fontsize=15)

# plt.savefig("plots/demo_cocoon_separation_%s.pdf"%orbit, dpi=300, bbox_inches='tight')
# %%
