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
# %%
datapath = '/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
table_path = "/n/home02/amphillips/p27_nbody/data/gmm_tables/"

make_plots=True
constrain_widths=False
noise = None #<-- None or 'via' or 'desi'
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

f_cocoons, cocoon_sigvgsrs, cocoon_sigphi2s, thin_sigvgsrs, thin_sigphi2s = [], [], [], [], []

use_constrained = False
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
    distances = coords_obs.distance.to(u.kpc).value
    med_distances.append(np.median(distances))

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
plt.savefig(repo_path+"/plots/"+plot_filename, dpi=300, bbox_inches='tight')

# %%
#
#
#--------------------------------------------------------#
#       plots to demonstrate the GMM results             #
#       for individuaal simulations (success and fails)  #
#--------------------------------------------------------#
import gmm
# %%
# c_labels = ["#FBBA72","#F5AE66","#EFA15A","#E9944E","#E38741","#DD7A35","#D76D29","#D1601D","#CA5310"]
tt = Table.read(table_path+case_name+".fits", format="fits")


c_labels = ["#CCC9E7", "#2F2F2F"]
cocoon_cmap = LinearSegmentedColormap.from_list('cocoon_cmap', c_labels)

orbit = 'pa5'
rvir_index=2
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


# if include_binaries==False:
sc_straighter = data_dict['sc_straighter'] #<-- dictionary
# if include_binaries==True:
#     sc_straighter = data_dict['sc_straighter_primaries'] #<-- dictionary

unbound = data_dict['unbound']
trim_new = data_dict['trim_new']
use = unbound & trim_new

keys = ['phi2','v_phi1','v_phi2','v_gsr']
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
    std = np.std(ydata)
    cut = 3*std #cocoon_sigmas[jj]

    ax = axs[jj,0]

    ax.scatter(sc_straighter['phi1'][use][order],  # plot cocoon on top. 
            sc_straighter[key][use][order], # plot cocoon on top. 
            # x_data[:,ii],
                c=p_cocoon[order], s=5, cmap=cocoon_cmap,
                rasterized=True) 
    ax.set_ylim(-cut,cut)
    ax.set_xlim(orbit_lim_map[orbit])

    # ax.set_ylim(-3*cut, 3*cut)
    ax.set_ylabel(key_labels[jj], fontsize=15)


    ax = axs[jj,1]
    bins = np.linspace(-cut, cut, 50)



    ax.hist(sc_straighter[key][use][~cocoon_selection], 
            alpha=1., density=False, weights = np.zeros_like(sc_straighter[key][use][~cocoon_selection])+1/sc_straighter[key][use][~cocoon_selection].size, 
            color=c_labels[0],orientation='horizontal',
            bins=bins)
    ax.hist(sc_straighter[key][use][cocoon_selection],
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
