#--------------------------------------------------------------#
#  save dictionaries that are [final data dict] plus observed  #
#   coordinates, both with and without noise/binary motions    #
#--------------------------------------------------------------#
# %%
import sys
soft_path = '/n/home02/amphillips/software/'
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

from scipy.stats import binned_statistic, norm, multivariate_normal
from scipy.optimize import curve_fit, minimize, Bounds
# from scipy.ndimage import gaussian_filter1d
# from scipy.interpolate import CubicSpline
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
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.colors as mcolors
# import matplotlib.cm as cm
# from matplotlib.gridspec import GridSpec
# from matplotlib.lines import Line2D
# plt.style.use(script_path+'/vedant.mplstyle')
# # %config InlineBackend.figure_format='retina'
# from matplotlib.cm import ScalarMappable
# from matplotlib.colors import Normalize
# from matplotlib.colors import LinearSegmentedColormap

from tqdm import tqdm
import pickle

sys.path.append(script_path)
sys.path.append(repo_path)
sys.path.append(soft_path+'viamock')
from streamframe import StreamFrame
import artpop
import viamock
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty

import PETAR_ANALYSIS_FUNCTIONS as paf
import inspect_new_sims as simspect 
import noise as noise #<-- for like isochrone stuff. 

import argparse
# %%
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

# if __name__=="__main__":
print("rng prep")
rng = np.random.default_rng(seed=42)

#### load 12 Gyr feh=-2.0 isocshrone with artpop:
print("preparing CMD")
isocmd = artpop.fetch_mist_iso_cmd(
    log_age=np.log10(12e9),
    feh=-2.0,
    phot_system='UBVRIplus',
    #v_over_vcrit=0.0 #<-- idk
)
m0_grid, iso_table = noise.build_isochrone_table(isocmd) #<-- will do gaia bands + teff automatically, have max_phase=5 (remove post agb evolution)
min_m0, max_m0 = m0_grid[0], m0_grid[-1]


print("outside of loop -- loading info about simulation grid. ")
grid_info = paf.extended_grid_info(scratch=False) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')

# ordering decided here. M3 is getting IGNORED :P
orbits = ['gd1','aau','pa5','jet','c19']
init_displacements = [
    grid_info.gd1_init_displacement, 
    grid_info.aau_init_displacement,
    grid_info.pa5_init_displacement,
    grid_info.jet_init_displacement,
    grid_info.c19_init_displacement]
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
copy_options = [0,1,2,3,4] #<-- order in which to try out copies. in practice there are <=5x copies per sim.


print("beginning loop...")
for ii, orbit in enumerate(tqdm(orbits)): #<--- this i can do later i think. 
    # if orbit !='aau':
    #     continue

    mass_index = 1 # <-- HIGH mass stellar population... should maximize cocoon contributions from stellar evolution-related kicks i think. 
    

    ### determine orbital phase -- informs the rtid boundary.

    init_displacement = init_displacements[ii]
    orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
    peri = orbit_obj.pericenter().to(u.kpc).value
    apo = orbit_obj.apocenter().to(u.kpc).value


    x,y,z = init_displacement[:3]
    r = np.sqrt(x**2 + y**2 + z**2)

    orbital_phase = (r - peri) / (apo - peri)

    if orbital_phase<0.5:
        tidal_boundary = 2.0 #<-- if we're closer to pericenter, make the tidal boundary 2 rtid since things will be recaptured. this would be gd1 and jet.
        print("orbital phase = ", orbital_phase, "; using boundary of %.1f rtid"%tidal_boundary)

    if orbital_phase>0.5:
        tidal_boundary = 1.0
    if orbital_phase>1:
        print("something's wrong dawg")

    print("orbital phase = ", orbital_phase, "; using boundary of %.1f rtid"%tidal_boundary)


    for rvir_index in range(4):
        # if rvir_index!=3:
        #     continue

        (core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
            simspect.prepare_nbody_data_anycopy(
                orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
                include_photometry=False, N_rtid_boundary = tidal_boundary, #<--- not sure what i'm going to use for tthis: 
                verbose=False
            )


        coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab) #<-- i think i straight up never actually need these. 
        distances = coords_obs.distance
        data_dict['coords_obs'] = coords_obs

        # ** this is always going to be centers of mass. gets returned as an immutible coordinate object, idk. 
        coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab)


        stellar_types = lumdict['type']
        nonrem = stellar_types<10
        m0s = np.asarray(lumdict['m0_zams'], dtype=np.float64) #<-- use masses from beginning of simulations
        # L, R = lumdict['L'].to(u.Lsun), lumdict['R'].to(u.Rsun)
        # Teff = noise.get_Teff(R, L).to(u.K)

        # straightened coords. first with an orbit, 
        sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) #<-- sc is returned as a DICTIONARY! 

        # straightened coords with primaries 
        PM_treatment = 'CoM'#<-- DECISION ABOUT PROPER MOTIONS BEING MADE HERE!!! 'CoM' or 'primary
        sc_primaries = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab,
                                                                    lumdict=lumdict, 
                                                                    PM_treatment=PM_treatment 
                                                                    )


        unbound = ~CMdict['in_rtid'] #<-- flag what's unbound from the cluster. 
        data_dict['unbound'] = unbound
        data_dict['nonrem'] = nonrem

        ### quick check on how restrictive unbound is... for debugging purposes. 
        # print(len(distances), len(distances[unbound]))


        ### NEW scheme for trimming the stream just dropped, no 'inMW' necessary now. 
        inMW_na = np.ones(len(sc['phi1']), dtype=bool) #<-- i don't actually want to do a "inMW" trim here. keep everything true but make the mask so functions downstream still work. 
        trim_new = trim_obstream_percentile(sc) 
        trim_new_primaries = trim_obstream_percentile(sc_primaries)
        data_dict['inMW_na'] = inMW_na
        data_dict['trim_new'] = trim_new
        data_dict['trim_new_primaries'] = trim_new_primaries


        sc_straighter = simspect.poly_straightening(sc, tc=[inMW_na, trim_new]) #<-- provide tc (trim criteria) so that the fitter doesn't lock to outliers but they're still included in the dataset. can exclude them later. 
        sc_straighter_primaries = simspect.poly_straightening(sc_primaries, tc=[inMW_na, trim_new_primaries])



        ##### GETTING A HANDLE ON THE CM VS PRIMARY TREATMENT FOR WHAT PROPER MOTIONS 
        #   THE BINARIES GET: 
        # fig, ax = plt.subplots()
        # nsin = data_dict['nsingles']

        # ## plot binaries
        # ts = trim_new[:nsin] #<-- note that the "trim" criteria with primary treatment will be different from the CM treatment and therefore also change depending on the PM treatment. 
        # tb = trim_new[nsin:] #<-- note that the "trim" criteria with primary treatment will be different from the CM treatment and therefore also change depending on the PM treatment. 

        # if PM_treatment=='primary':
        #     ax.set_title("binaries have instantaneous proper motion")
        # if PM_treatment=='CoM':
        #     ax.set_title("binaries have cm proper motion")
        # ax.scatter(sc_straighter_primaries['v_gsr'][nsin:][tb], sc_straighter_primaries['pm_phi2'][nsin:][tb],
        #            c='k', label='binaries')

        # ## plot singles
        # ax.scatter(sc_straighter_primaries['v_gsr'][:nsin][ts], sc_straighter_primaries['pm_phi2'][:nsin][ts],
        #            c='tomato', label='single stars')
        # ax.set_xlabel(r'$\Delta v_{\rm GSR}~[\rm km~s^{-1}]$')
        # ax.set_ylabel(r'$\Delta \mu_{\phi_2}~[\rm mas~yr^{-1}]$')
        # ax.set_ylim(-0.45, 0.45)
        # ax.set_xlim(-20,20)
        # ax.legend(loc='upper left', bbox_to_anchor=(1,1))


        data_dict['sc_straighter'] = sc_straighter # <--- same length as coords_obs. 
        data_dict['sc_straighter_primaries'] = sc_straighter_primaries
        ###### get noise: 
        
        ### isochrone photometry as loaded above:
        iso_phot, on_iso = noise.isochrone_photometry(m0s, m0_grid, iso_table)
        G  = iso_phot['Gaia_G_EDR3']
        BP = iso_phot['Gaia_BP_EDR3']
        RP = iso_phot['Gaia_RP_EDR3']
        BP_RP = BP-RP
        z = noise.gaia_g_to_lsst_z(G, BP_RP)
        mz = paf.m_from_M(z, dist=distances)
        mG = paf.m_from_M(G, dist=distances)
        log_Teff = iso_phot['log_Teff']

        alive = on_iso #<-- removes high mass things that are remnants; should already be removed by remnants cut tho? 
        acceptable_G_viamock = (mG>5) & (mG<30)
        data_dict['alive'] = alive
        data_dict['acceptable_G'] = acceptable_G_viamock

        photdict = {
            "G":G,
            "BP_RP":BP_RP,
            "z":z,
            'mG':mG,
            'mz':mz,
            'log_Teff':log_Teff
        }
        data_dict['phot'] = photdict #<-- add a bunch more numbers

        #### RV, POSITION, and PM errors 
        rverr_desi = noise.desi_RVerr(zmag=mz, feh=-2.0)

        rverr_via = np.full(len(mG), fill_value = np.nan)
        rverr_via[acceptable_G_viamock] = noise.via_RVerr(
            G = mG[acceptable_G_viamock],
            feh = -2.0,
            log_Teff=log_Teff[acceptable_G_viamock],
            exptime_s = 3600,
            nexp=1 #<-- a choice to make
            ### could also add seeing, airmass, moon, etc. 
        ) #<-- if G is out of the viamock table range (5-30) I will just have nan values. 



        vgsr_noise_via = rng.normal(0, rverr_via)
        vgsr_noise_desi = rng.normal(0, rverr_desi)



        pm_err = total_proper_motion_uncertainty(mG, 'dr3') / np.sqrt(2) #<-- we'll add some in two dimensions
        pos_err = total_position_uncertainty(mG, 'dr3') / np.sqrt(2)

        rng = np.random.default_rng(seed=42)
        pmphi1_noise = (rng.normal(0, pm_err) * u.microarcsecond / u.yr).to(u.mas/u.yr)
        pmphi2_noise = (rng.normal(0, pm_err) * u.microarcsecond / u.yr).to(u.mas/u.yr)
        phi1_noise = (rng.normal(0, pos_err) * u.microarcsecond).to(u.degree)
        phi2_noise = (rng.normal(0, pos_err) * u.microarcsecond).to(u.degree)

        noise_dict = {
            'phi1': phi1_noise.to(u.degree).value,
            'phi2': phi2_noise.to(u.degree).value,
            'pm_phi1': pmphi1_noise.to(u.mas/u.yr).value,
            'pm_phi2': pmphi2_noise.to(u.mas/u.yr).value,
            'v_gsr_via': vgsr_noise_via,#<-- the sampled RV noise from a gausian of std rverr_[survey]
            'v_gsr_desi':vgsr_noise_desi,#<-- the sampled RV noise from a gausian of std rverr_[survey]
            'rverr_via': rverr_via, #<-- the rv uncertainties
            'rverr_desi': rverr_desi, #<-- the rv uncertainties
            'pm_err_gaia': pm_err #<-- pm uncertainty ( total / sqrt2 )
            }

        data_dict['noise'] = noise_dict



        ### dump everything in scratch until I figure out how large the files will be all together...
        datapath='/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
        rvir = rvirs[rvir_index]
        print("dumping to pkl file...")
        with open(datapath+'%s_%.2f.pickle'%(orbit, rvir), 'wb') as handle:
            pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

# %%
