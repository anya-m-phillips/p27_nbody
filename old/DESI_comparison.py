# creating a faithful comparison between Jarvis+26 GD-1 dataset from 
# DESI DR2 and N-body mocks of GD-1. Determining the requirement of 
# external heating for a cocoon structure. 
# %%
################ import packages
import sys
script_path = "/n/home02/amphillips/stream_velocity_structures/scripts" # for cannon

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

from tqdm import tqdm

sys.path.append(script_path)
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
from get_escapers import rms, make_key, dedupe_true_copies

# first import... very beefy...
from sklearn.mixture import GaussianMixture
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty
# %%
print("loading GD-1 data")
tt = Table.read("/n/home02/amphillips/data/jarvis26_Table7.fits", format="fits")
jarvis_coords = {
    "phi1":tt['phi1'],
    "phi2":tt["DELTA_PHI2"],
    "v_gsr":tt["DELTA_VGSR"],
    "pm_phi1":tt["DELTA_PM_PHI1"],
    "pm_phi2":tt["DELTA_PM_PHI2"]
}

coord_errors = {
    "phi1":np.zeros(len(tt['phi1'])),
    "phi2":np.zeros(len(tt['DELTA_PHI2'])),
    "v_gsr":tt['V_ERR']*u.km/u.s,
    'pm_phi1':tt['PM_PHI1_ERR']*u.mas/u.yr,
    "pm_phi2":tt['PM_PHI2_ERR']*u.mas/u.yr
}
# %%
### setup stuff:
print("doing setup for petar data...")
paths = paf.define_paths()
tdis_vals, t_peri_vals, t_apo_vals = paf.get_tdis_tplot(paths) ### later. 

apocenters = paf.define_apocenters()
init_displacements = paf.define_init_displacements()
escaperdict_path = "/n/home02/amphillips/stream_velocity_structures/data/"

rvir0_values = np.array([0.75,0.75, 1.5,1.5,3.,3.,6.,6.,]*3)

# %%
print("defining relevant functions and parameters...")
#------------------------------------------------------------------------------#
#           PHOTOMETRY SETUP                                                   #
#------------------------------------------------------------------------------#
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

def get_z_photometry(Teff, Radius, distance):
    f_min = nu_z_min.to(u.Hz).value
    f_max = nu_z_max.to(u.Hz).value
    mag = paf.integrated_mag(f_min, f_max, Teff, Radius, distance)
    return mag

#------------------------------------------------------------------------------#
#                       DESI measurement errors                                #
#------------------------------------------------------------------------------#
def desi_RVerr(zmag, feh=-2.0):
    """
    get RV error for desi data model, which 
    depends on z magnitude and metallicity
    """
    log_err = -0.47 + 0.27*(zmag-16) - 0.23*feh
    return 10**log_err

#------------------------------------------------------------------------------#
#                     data unpacking stuff...                                  #
#------------------------------------------------------------------------------#
def unpack_data_dict(data_dict, binary_criterion, clip_outliers=False, i=10000):
    in_dp = data_dict["init_displacement"]

    dd = data_dict[binary_criterion]
    init_displacement = data_dict['init_displacement']
    coords = dd['coords']

    phi1 = coords["phi1"]
    phi2_straight = dd["phi2_straight"]
    vr_straight = dd["vr_straight"]
    inMW, trim = dd["inMW"], dd["trim"]

    nsingles = data_dict["nsingles"]
    single_mask = np.zeros(len(phi1), dtype=bool)
    single_mask[:nsingles]=True

    single_mask = single_mask[inMW][trim]

    #### yikes hold up -- I want to make sure I use the center of mass frame to get polynomials for subtracting out trends over the orbits.
    if binary_criterion!="CoM":
        ddd = data_dict["CoM"]

        #introduce "fs" (for straightening) coordin=ates
        phi1_fs = ddd['coords']['phi1']
        vr_straight_fs = ddd['vr_straight']
        phi2_straight_fs = ddd['phi2_straight']
        coords_fs = ddd['coords']
        inMW_fs, trim_fs = ddd['inMW'], ddd['trim']
        #### treating the proper motions differently -- want to convert to velocities _then_ straighten with an interpolated orbit + polynomial.
        pmphi1_fs, pmphi2_fs = coords_fs['pm_phi1']*u.degree.to(u.radian)/u.Myr, coords_fs['pm_phi2']*u.degree.to(u.radian)/u.Myr # these are now in radians per Myr
        v_phi1_fs = pmphi1_fs * (coords_fs['r']*u.kpc).to(u.km)
        v_phi2_fs = pmphi2_fs * (coords_fs['r']*u.kpc).to(u.km)
        v_phi1_fs = v_phi1_fs.to(u.km/u.s).value
        v_phi2_fs = v_phi2_fs.to(u.km/u.s).value

        _, _, poly_phi2, fit_phi2 = paf.straighten_stream_polynomial(phi1_fs, phi2_straight_fs, return_poly_fn=True, trim_criteria=[inMW_fs, trim_fs])
        _, _, poly_vr, fit_vr = paf.straighten_stream_polynomial(phi1_fs, vr_straight_fs, return_poly_fn=True, trim_criteria=[inMW_fs, trim_fs])

        phi2_straighter = phi2_straight - poly_phi2(phi1, *fit_phi2)
        vr_straighter = vr_straight - poly_vr(phi1, *fit_vr)

    else:
        # after subtracting off the progenitor orbit from phi2, fit a polynomial to remove the additional curves. 
        _, phi2_straighter = paf.straighten_stream_polynomial(phi1, phi2_straight, trim_criteria = [inMW, trim], 
                                                            degree=5)
        _, vr_straighter = paf.straighten_stream_polynomial(phi1, vr_straight, trim_criteria = [inMW, trim],
                                                            degree=5)



    #### treating the proper motions differently -- want to convert to velocities _then_ straighten with an interpolated orbit + polynomial.
    pmphi1, pmphi2 = coords['pm_phi1']*u.degree.to(u.radian)/u.Myr, coords['pm_phi2']*u.degree.to(u.radian)/u.Myr # these are now in radians per Myr
    v_phi1 = pmphi1 * (coords['r']*u.kpc).to(u.km)
    v_phi2 = pmphi2 * (coords['r']*u.kpc).to(u.km)
    v_phi1 = v_phi1.to(u.km/u.s).value
    v_phi2 = v_phi2.to(u.km/u.s).value


    ### orbit interpolating to remove trends in pmphi1 and pmphi2 as well? 
    y, interp_orbit, orbit_coords, orbit_w, ii  = paf.straighten_stream_orbit_interp(coords, yval=['pm_phi1', 'pm_phi2'], core=None, i=i,
                                                            use_core=False, init_displacement=init_displacement, return_orbit_chunk=True)
    pmphi1_straight, pmphi2_straight = y

    ## get and interpolate the v phi1, v phi2
    pmphi1_orbit, pmphi2_orbit, r_orbit = orbit_coords['pm_phi1']*u.degree.to(u.radian)/u.Myr, orbit_coords['pm_phi2']*u.degree.to(u.radian)/u.Myr, orbit_coords['r']*u.kpc
    vphi1_orbit, vphi2_orbit = r_orbit*pmphi1_orbit, r_orbit*pmphi2_orbit
    vphi1_orbit = vphi1_orbit.to(u.km/u.s).value
    vphi2_orbit = vphi2_orbit.to(u.km/u.s).value
    phi1_orbit = orbit_coords['phi1']

    vphi1_straight = v_phi1 - np.interp(phi1, phi1_orbit, vphi1_orbit)
    vphi2_straight = v_phi2 - np.interp(phi1, phi1_orbit, vphi2_orbit)

    if binary_criterion!="CoM":
        vphi1_straight_fs = v_phi1_fs - np.interp(phi1_fs, phi1_orbit, vphi1_orbit)
        vphi2_straight_fs = v_phi2_fs - np.interp(phi1_fs, phi1_orbit, vphi2_orbit)

        _, _, poly_vphi1, fit_vphi1 = paf.straighten_stream_polynomial(phi1_fs, vphi1_straight_fs, return_poly_fn=True, trim_criteria=[inMW_fs, trim_fs])
        _, _, poly_vphi2, fit_vphi2 = paf.straighten_stream_polynomial(phi1_fs, vphi2_straight_fs, return_poly_fn=True, trim_criteria=[inMW_fs, trim_fs])

        vphi1_straighter = vphi1_straight - poly_vphi1(phi1, *fit_vphi1)
        vphi2_straighter = vphi2_straight - poly_vphi2(phi1, *fit_vphi2)

    else:
        _, vphi1_straighter = paf.straighten_stream_polynomial(phi1, vphi1_straight, trim_criteria=[inMW, trim], degree=5)
        _, vphi2_straighter = paf.straighten_stream_polynomial(phi1, vphi2_straight, trim_criteria=[inMW, trim], degree=5)


    ### get everything all together:
    phi1 = phi1[inMW][trim]
    phi2 = phi2_straighter[inMW][trim]
    vr = vr_straighter[inMW][trim]
    vphi1 = vphi1_straighter[inMW][trim]
    vphi2 = vphi2_straighter[inMW][trim]

    # trim outliers (this was helpful for the GMM, not so necessary for a "bespoke" cocoon separation.)
    OUTLIER_CLIP = (np.abs(vphi1_straighter)<100) & (np.abs(vphi2_straighter)<100) & (np.abs(vr_straighter)<100)
    om = OUTLIER_CLIP[inMW][trim]


    if clip_outliers:
        return phi1[om], phi2[om], vr[om], vphi1[om], vphi2[om], single_mask[om]#, om
    else:
        return phi1, phi2, vr, vphi1, vphi2, single_mask #, om
    
# %%
print("functions to load N-body data and prepare it in a dictionary, calculate photometry...")
# load info about the simulation: 


def prepare_nbody_data(n, include_photometry=True, i=10000):
    path = paths[n]
    apo = apocenters[n]
    core = paf.load_core(path)
    init_displacement = init_displacements[n]
    

    if n==10:
        file_index = int(i/10)
    else:
        file_index=None

    data_dict = paf.intrinsic_stream_data_v2(n, i, core, path, apo, 
                                             use_core=False, 
                                             init_displacement=init_displacement,
                                            file_index=file_index)

    CMdict = data_dict['CoM']
    lumdict = data_dict['luminous']
    inMW, trim = CMdict['inMW'], CMdict['trim']

    if include_photometry==True:
        # do the calculations -- ORDERED AS ALL_PARTICLES
        primary_Ls = lumdict["L"][inMW][trim] #* u.Lsun
        primary_Rs = lumdict['R'][inMW][trim] #* u.Rsun
        primary_Teffs = (primary_Ls / (4*np.pi*primary_Rs**2 * const.sigma_sb))**(1/4)

        R_cgs = primary_Rs.cgs.value
        Teff_cgs = primary_Teffs.cgs.value
        G, BP, RP, z = [],[],[],[]
        for Tval, Rval in tqdm(zip(Teff_cgs, R_cgs)):
            Gval, BPval, RPval = get_gaia_photometry(Tval, Rval, 10) # 10 pc. 
            G.append(Gval)
            BP.append(BPval)
            RP.append(RPval)

            zval = get_z_photometry(Tval, Rval, 10) # 10pc
            z.append(zval)
        G = np.array(G)
        BP = np.array(BP)
        RP = np.array(RP)
        z = np.array(z)

        return core, data_dict, CMdict, lumdict, inMW, trim, G, BP, RP, z
    
    else:
        return core, data_dict, CMdict, lumdict, inMW, trim

def get_GD1_coords_nbody(data_dict, lumdict, subdict, inMW, trim, core, i=10000):
    init_displacement = data_dict['init_displacement']
    pos, vel = subdict['pos'], subdict['vel']
    coords = subdict['coords']

    primary_types = lumdict["type"][inMW][trim] #<-- that said, always store the type (and eventually stellar parameters) from the luminous companion. since I order the center of mass dictionary and the "luminous" dictionary as [singles, binaries] this ordering should be the same. 
    nonrem = primary_types<10

    ### step 1: transform galactocentric position+velocity to ICRS
    icrs_coords = paf.galcen_to_ICRS(pos, vel)
    ### step 2: correct for solar reflex motion
    icrs_coords = reflex_correct(icrs_coords)

    ### step 3: transform to Koposov GD1 great circle frame
    gd1_coords = icrs_coords.transform_to(gc.GD1Koposov10)

    ### step 4 catalog these. note v_gsr should come from the reflex-corrected icrs coords. 
    phi1 = gd1_coords.phi1[inMW][trim].to(u.degree)
    phi2 = gd1_coords.phi2[inMW][trim].to(u.degree)
    v_gsr = icrs_coords.radial_velocity[inMW][trim].to(u.km/u.s)
    pmphi1 = gd1_coords.pm_phi1_cosphi2[inMW][trim].to(u.mas/u.yr)
    pmphi2 = gd1_coords.pm_phi2[inMW][trim].to(u.mas/u.yr)
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
        return_orbit_chunk=True, use_core=False,
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

def add_noise(z_dist, G_dist):
    rng = np.random.default_rng(seed=42) #<-- seed re-initialized every time I re-run the error finding cell.
    vgsr_errs = desi_RVerr(z_dist, feh=-2.0)
    vgsr_noise = rng.normal(0, vgsr_errs) * u.km/u.s
    ### not 100% sure this is accurate since this is technically a proper motion uncertainty in icrs, not in transformed streamframe coords. .
    pm_errs = total_proper_motion_uncertainty(G_dist, 'dr3') 
    pm_phi1_noise = (rng.normal(0, pm_errs) * u.microarcsecond / u.yr).to(u.mas/u.yr)
    pm_phi2_noise = (rng.normal(0, pm_errs) * u.microarcsecond / u.yr).to(u.mas/u.yr)
    pos_errs = total_position_uncertainty(G_dist, 'dr3')
    phi2_noise = (rng.normal(0, pos_errs) * u.microarcsecond).to(u.degree)
    
    nbody_noise_dict = {
        "phi2":phi2_noise,
        "v_gsr":vgsr_noise,
        "pm_phi1":pm_phi1_noise,
        "pm_phi2":pm_phi2_noise
    }
    # return phi2_noise, vgsr_noise, pm_phi1_noise, pm_phi2_noise
    return nbody_noise_dict

def get_mag_mask(G_dist):
    DESI_mags = tt['GMAG0']

    bins = np.linspace(min(DESI_mags), max(DESI_mags), 21)

    # match the CDF ... 
    uu = np.random.rand(len(DESI_mags))
    target_mags = np.quantile(DESI_mags, uu)
    sim_inds = np.argsort(G_dist)
    sim_sorted = G_dist[sim_inds]
    matched_inds = np.searchsorted(sim_sorted, target_mags)
    matched_inds = np.clip(matched_inds, 0, len(sim_sorted)-1)
    selected_sim = sim_inds[matched_inds]

    selected_sim = np.unique(selected_sim)
    return selected_sim

# %%
print("plotting stuff")
def make_comparison_plots(nbody_coords_straight,
                          nbody_noise_dict,
                          nonrem,
                          mag_mask,
                          add_noise=True
                          ):
    fig, axs = plt.subplots(4,1,figsize=[8,12], sharex=True)
    plt.subplots_adjust(wspace=0.03, hspace=0.03)

    yvals = ['phi2','v_gsr','pm_phi1','pm_phi2']
    y_labels = [r'$\phi_2~[\degree]$', r'$v_{\rm GSR}~[\rm km~s^{-1}]$', r'$\mu_{\rm{\phi_1}}~\rm[mas~yr^{-1}]$',r'$\mu_{\rm{\phi_1}}~\rm[mas~yr^{-1}]$']
    lims = [5, 45, 3, 3]
    for k, yval in enumerate(yvals):
        ax = axs[k]

        nbody_y = nbody_coords_straight[yval][nonrem].value

        ### making sure the straightening polynomiail is fit to the full dataset ([inMW][trim]-ed ofc)
        _, _, poly, fit = paf.straighten_stream_polynomial(nbody_coords_straight['phi1'].value, nbody_coords_straight[yval].value,
                                                        trim_criteria = [np.ones(len(nbody_coords_straight[yval])).astype(bool),np.ones(len(nbody_coords_straight[yval])).astype(bool)],
                                                        return_poly_fn=True)

        nbody_y -= poly(nbody_coords_straight['phi1'][nonrem].value, *fit)

        if add_noise==True:
            y_noise = nbody_noise_dict[yval][nonrem].value
        else:
            y_noise = np.zeros(len(nbody_noise_dict[yval][nonrem].value))
        ################################ PLOTTING HERE ################################################

        ax.set_ylim(-lims[k], lims[k])

        ### an option to match the cdf of the desi Gmag distribution
        ax.scatter(nbody_coords_straight['phi1'][nonrem][mag_mask].value, nbody_y[mag_mask]+y_noise[mag_mask], 
                    c='k', s=3)
        ax.scatter(jarvis_coords['phi1'], jarvis_coords[yval], s=5, zorder=0,  c=tt['P_THIN'], cmap='cool', edgecolor='k', lw=0.1)

        ax.set_ylabel(y_labels[k])

    axs[-1].set_xlabel(r'$\phi_1~[\degree]$')

    return fig, axs

def make_cmd(G, BP, RP):
    fig, ax = plt.subplots()
    ob = ax.scatter(BP-RP, G_dist, c=nbody_coords['phi1'].to(u.degree).value, cmap='coolwarm', s=5, rasterized=True)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(position='right', size='5%', pad=0.03)
    fig.colorbar(ob, cax, ax, label=r'$\phi_1$')

    ax.invert_yaxis()
    # kind of lit u guys .... not so so bad u guys.....
    ax.set_xlim(-0.75, 2)
    ax.set_ylim(25, 12)

    ax.set_xlabel(r'$G_{\rm bp}-G_{\rm rp}$')
    ax.set_ylabel(r'$G$')

    ax.axhline(max(tt['GMAG0']), c='k', lw=1, ls='--')
    return fig, ax 
# %%
#----------------------------------------------------------------#
#                            main program below                  #
#----------------------------------------------------------------#
if __name__=='__main__': #### basically I do not want this part to run when I import this script elsewhere. 
    # decide which simulation to use 
    n=8
    # step 1: prepare data -> outputs core, data_dict, CMdict, lumdict, inMW, trim, G, BP, RP, z
    core, data_dict, CMdict, lumdict, inMW, trim, G, BP, RP, z = prepare_nbody_data(n)
    # %%


    # step 2: put the data in straightened GD-1 coords from the literature; extract nonremnants
    nbody_coords_straight, nbody_coords, nonrem = get_GD1_coords_nbody(data_dict, lumdict,
                                                                    subdict=CMdict,
                                                                    inMW=inMW, trim=trim, core=core)

    # step 3: compute new mags, output an intermediate cmd?
    G_dist = paf.m_from_M(G, dist=nbody_coords['dist']) # NOT STRAIGHTENED ! 
    z_dist = paf.m_from_M(z, dist=nbody_coords['dist'])
    fig, ax = make_cmd(G_dist, BP, RP)
    # step 4: noise up the data.
    nbody_noise_dict = add_noise(z_dist, G_dist)
    # %%
    # step 5: determine the mag selection (reject repeated entries in the CDF.)
    mag_mask = get_mag_mask(G_dist[nonrem]) # <-- since mag mask will be indices, need to trim remnants first
    
    # step 6: run the plotting function. 
    fig, axs = make_comparison_plots(nbody_coords_straight, 
                                    nbody_noise_dict,
                                    nonrem,
                                    mag_mask)
    # plt.savefig("fig/noised_comparison_%i.pdf"%n, dpi=300, bbox_inches='tight')
# %%
#### alternatively, dw about noising. 
    alltrue_mask = np.ones(len(G_dist[nonrem]), dtype=bool)
    fig, axs = make_comparison_plots(nbody_coords_straight,
                                    nbody_noise_dict,
                                    nonrem,
                                    alltrue_mask,
                                    add_noise=False)
    plt.savefig("fig/noiseless_comparison_%i.pdf"%n, dpi=300, bbox_inches='tight')
# %%

# %%
