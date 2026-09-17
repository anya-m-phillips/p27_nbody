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

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']

pericenters_kpc = []
apocenters_kpc = []
present_rs = []

f_cocoons, cocoon_sigvgsrs, cocoon_sigphi2s, thin_sigvgsrs, thin_sigphi2s = [], [], [], [], []

use_constrained = True
for ii, orbit in enumerate(tqdm(orbits)): #<--- this i can do later i think. 
    # if orbit=='pa5' or orbit=='m3':
    #     continue
    # do the orbit-wise check -- integrate prog orbit and find the pericenter. 
    init_displacement = init_displacements[ii]
    orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
    peri = orbit_obj.pericenter()
    apo = orbit_obj.apocenter()
    pericenters_kpc.append(peri.to(u.kpc).value)
    apocenters_kpc.append(apo.to(u.kpc).value)

    x,y,z = init_displacement[:3]
    r = np.sqrt(x**2 + y**2 + z**2)
    present_rs.append(r) # kpc


    mass_index = 1 # <-- high mass stellar population... 

    f_cocoons_this_orbit = []
    cocoon_sig_vgsr_this_orbit = []
    cocoon_sig_phi2_this_orbit = []

    thin_sig_vgsr_this_orbit = []
    thin_sig_phi2_this_orbit = []

    for rvir_index in range(4):
        rvir = rvirs[rvir_index]

        ### load dictionary
        data_path = "/n/home02/amphillips/p27_nbody/data/data_dicts/"
        if use_constrained==True:
            data_path+='constrained/'
        else:
            data_path+='unconstrained/'
        filename = data_path+"%s_%.2f.pickle"%(orbit, rvir)
        with open(filename, 'rb') as handle:
            data_dict = pickle.load(handle)

        cocoon_dict = data_dict['cocoon_info']

        mu_thin, sigma_thin = cocoon_dict['mu_thin'], cocoon_dict['sigma_thin']
        mu_cocoon, sigma_cocoon = cocoon_dict['mu_cocoon'], cocoon_dict['sigma_cocoon']
        f_cocoon = cocoon_dict['f_cocoon']

        f_cocoons_this_orbit.append(f_cocoon)

        cocoon_sig_phi2_this_orbit.append(sigma_cocoon[0])
        cocoon_sig_vgsr_this_orbit.append(sigma_cocoon[-1])

        thin_sig_phi2_this_orbit.append(sigma_thin[0])
        thin_sig_vgsr_this_orbit.append(sigma_thin[-1])


        ### for now don't care about this. 
        # ol_clip = cocoon_dict['ol_clip']
        # unbound = cocoon_dict['unbound']
        # use = ol_clip & unbound

        # p_thin = cocoon_dict['p_thin']
        # sc_straighter = cocoon_dict['sc_straighter']


    f_cocoons.append(f_cocoons_this_orbit)
    cocoon_sigvgsrs.append(cocoon_sig_vgsr_this_orbit)
    cocoon_sigphi2s.append(cocoon_sig_phi2_this_orbit)
    thin_sigvgsrs.append(thin_sig_vgsr_this_orbit)
    thin_sigphi2s.append(thin_sig_phi2_this_orbit) 

pericenters_kpc = np.array(pericenters_kpc)
apocenters_kpc = np.array(apocenters_kpc)
present_rs = np.array(present_rs)

f_cocoons = np.array(f_cocoons)
cocoon_sigvgsrs = np.array(cocoon_sigvgsrs)
cocoon_sigphi2s = np.array(cocoon_sigphi2s)
thin_sigvgsrs = np.array(thin_sigvgsrs)
thin_sigphi2s = np.array(thin_sigphi2s)


### roughly ~amount of the way through orbit
orbital_phases = (present_rs - pericenters_kpc) / (apocenters_kpc - pericenters_kpc)
orbits = np.array(orbits)

eccentricities = (apocenters_kpc - pericenters_kpc) / (apocenters_kpc + pericenters_kpc)
# %%
reordered = np.argsort(orbital_phases)

ccc = cc[1:]
fig, axs = plt.subplots(2,3,figsize=[21,14])
for ii, orbit in enumerate(tqdm(orbits[reordered])):
    # if orbit in ['m3','pa5']:
    #     continue

    f_cocoons_this_orbit = f_cocoons[reordered][ii]


    x = rvirs
    axs[0,0].plot(x, f_cocoons_this_orbit, 
                # label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[reordered][ii],
                label = orbit+r'; $\varphi_{\rm orb} =%.2f$'%orbital_phases[reordered][ii],
                # label = orbit+r'; $e=%.2f$'%eccentricities[reordered][ii],
                marker='o', color=ccc[ii], markersize=10)


    ### cocoon ! ! !
    phi2_dispersions_this_orbit = cocoon_sigphi2s[reordered][ii]
    vgsr_dispersions_this_orbit = cocoon_sigvgsrs[reordered][ii]
    axs[0,1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')
    axs[0,2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10, ls='-')

    ### thin ! ! !
    phi2_dispersions_this_orbit = thin_sigphi2s[reordered][ii]
    vgsr_dispersions_this_orbit = thin_sigvgsrs[reordered][ii]
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

filename="plots/gmm_summary"
if use_constrained==True:
    filename+="_constrained.pdf"
else:
    filename+="_unconstrained.pdf"
plt.savefig(filename, dpi=300, bbox_inches='tight')

# %%

# c_labels = ["#FBBA72","#F5AE66","#EFA15A","#E9944E","#E38741","#DD7A35","#D76D29","#D1601D","#CA5310"]
c_labels = ["#CCC9E7", "#2F2F2F"]
cocoon_cmap = LinearSegmentedColormap.from_list('cocoon_cmap', c_labels)

orbit = 'gd1'
rvir_index=0
rvir = rvirs[rvir_index]
filename = data_path+"%s_%.2f.pickle"%(orbit, rvir)
with open(filename, 'rb') as handle:
    data_dict = pickle.load(handle)

cocoon_dict = data_dict['cocoon_info']

## for now don't care about this. 
# ol_clip = cocoon_dict['ol_clip']
trim_new = cocoon_dict['trim_new']
unbound = cocoon_dict['unbound']
# use = ol_clip & unbound
use = unbound

p_thin = cocoon_dict['p_thin']
p_cocoon = 1-p_thin
order = np.argsort(p_cocoon)
sc_straighter = cocoon_dict['sc_straighter']

cocoon_sigmas = cocoon_dict['sigma_cocoon']

fig, axs = plt.subplots(len(keys), 2, figsize=[10, 10], width_ratios = [4,1])

plt.subplots_adjust(hspace=0.03, wspace=0.03)

key_labels = [
    r'$\phi_2~[\degree]$',
    r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$',
    r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$',
    r'$v_{\rm GSR}~[\rm km~s^{-1}]$'
]
for jj, key in enumerate(keys):
    cut = 3*cocoon_sigmas[jj]

    ax = axs[jj,0]

    ax.scatter(sc_straighter['phi1'][use][order],  # plot cocoon on top. 
            sc_straighter[key][use][order], # plot cocoon on top. 
            # x_data[:,ii],
                c=p_cocoon[order], s=5, cmap=cocoon_cmap,
                rasterized=True) 
    ax.set_ylim(-cut,cut)


    # ax.set_ylim(-3*cut, 3*cut)
    ax.set_ylabel(key_labels[jj], fontsize=15)
    ax.set_xlim(-100, 15) #<-- gd1
    # ax.set_xlim(-20,10) #<-- jet


    ax = axs[jj,1]
    bins = np.linspace(-cut, cut, 50)



    cocoon_selection = p_thin<0.5
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

plt.savefig("plots/demo_cocoon_separation.pdf", dpi=300, bbox_inches='tight')
# %%
