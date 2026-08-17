# %%
################ import packages
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline


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
# from get_escapers import rms, make_key, dedupe_true_copies

# first import... very beefy...
from sklearn.mixture import GaussianMixture
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty
# %%

### DEFINE STUFF ABOUT THE GRID: 
grid_info = paf.extended_grid_info(scratch=True) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')

# %%
# main helper function: 
def prepare_nbody_data(path = grid_info.gd1_lm_paths[0]+"0/", 
                       include_photometry=False, 
                       i=grid_info.gd1_age,
                       apo = grid_info.gd1_apo,
                       init_displacement = grid_info.gd1_init_displacement,
                       ):
    """
    this is expecting N-body data that outputs every 10 Myr. 
    """
    core = paf.load_core(path)
    file_index = int(i/10)

    data_dict = paf.intrinsic_stream_data_v3(
        path, i, core, apo, init_displacement,
        use_core=False,
        binary_treatments=['CoM','companions','luminous']
    )

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
            Gval, BPval, RPval = paf.get_gaia_photometry(Tval, Rval, 10) # 10 pc. 
            G.append(Gval)
            BP.append(BPval)
            RP.append(RPval)

            zval = paf.get_z_photometry(Tval, Rval, 10) # 10pc
            z.append(zval)
        G = np.array(G)
        BP = np.array(BP)
        RP = np.array(RP)
        z = np.array(z)

        return core, data_dict, CMdict, lumdict, inMW, trim, G, BP, RP, z
    
    else:
        return core, data_dict, CMdict, lumdict, inMW, trim

def prepare_nbody_data_anycopy(orbit, stellar_pop='lm', rvir_index=0, copies=range(5),
                               verbose=True, **kwargs):
    """
    retrieve_sim_info()'s `copy` only picks which subdirectory to look in, and not every
    copy in the grid ran all the way to the present day. an unfinished one still has a
    data.core (just a short one), so it gets as far as intrinsic_stream_data_v3 and then
    raises FileNotFoundError on the missing data.<file_index> snapshot.

    every (orbit, stellar_pop, rvir_index) has at least one copy that finished, so try
    them in order and take the first that loads.

    extra kwargs go straight to prepare_nbody_data (e.g. include_photometry=True).
    returns (whatever prepare_nbody_data returns), path, apo, age, init_displacement, copy
    """
    tried = []
    for copy in copies:
        path, apo, age, init_displacement = grid_info.retrieve_sim_info(
            orbit=orbit, stellar_pop=stellar_pop, rvir_index=rvir_index, copy=copy
        )
        if path in tried: #<-- circ has no copy subdir, so retrieve_sim_info ignores `copy`
            continue
        tried.append(path)

        try:
            out = prepare_nbody_data(
                path = path,
                i = age if orbit != "circ" else 30000,
                init_displacement = init_displacement,
                apo = apo, #<-- always pass apo so that inMW, trim is correct and not always based on GD-1 orbit.
                **kwargs
            )
        except FileNotFoundError as err:
            if verbose:
                print("  %s rvir_index=%i copy %i unfinished, missing %s"%(
                    orbit, rvir_index, copy, err.filename))
            continue

        if verbose:
            print("  %s rvir_index=%i: using copy %i"%(orbit, rvir_index, copy))
        return out, path, apo, age, init_displacement, copy

    raise FileNotFoundError("no finished copy of %s/%s/rvir_index=%i; tried %s"%(
        orbit, stellar_pop, rvir_index, list(copies)))

def rotation_matrix(a,b,c):
    """
    rotate angles a,b,c about x,y,z axes
    https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
    """
    Rx = np.array([
        [1,0,0],
        [0, np.cos(a), -np.sin(a)],
        [0, np.sin(a), np.cos(a)]
    ])
    Ry = np.array([
        [np.cos(b), 0, np.sin(b)],
        [0,1,0],
        [-np.sin(b), 0, np.cos(b)]
    ])
    Rz = np.array([
        [np.cos(c), -np.sin(c), 0],
        [np.sin(c), np.cos(c), 0],
        [0,0,1]
    ])

    # R = Rx @ Ry @ Rz
    R = Rz @ Ry @ Rx # <-- do x rotation first, _then_ z rotation.. mostly for movie purposes. 
    return R


# %%
def streamframe_coords_observed(orbit, data_dict, prog_tab):
    """
    orbit options are:
    gd1, pa5, jet, aau, m3, c19, circ
    
    (input as strings), will put into a great-circle coordinate frame 
    based on ICRS ra/dec. this will allow accounting for projection 
    effects. also will correct radial velocity + proper motion for solar
    reflex motion, using a distance that for now is assumed to 
    be known exactly. For real streams this would be based on 
    a dynamically-fit distance track and hopefully would also be 
    pretty well-constrained....

    data_dict should be one of the three choicdes of subdictionary (CoM, luminous, companions)
    so that it has position and velocity info. those should have astropy units. 

    returns a coordinate object, since all frames can be constructed in gala with great circles. 
    """
    pos, vel = data_dict['pos'], data_dict['vel']
    coords_ICRS = paf.galcen_to_ICRS(pos, vel) 
    coords_ICRS = reflex_correct(coords_ICRS) #<-- correct for solar reflex motion
    row = prog_tab[prog_tab['name']==orbit]
    

    if orbit=='gd1': #<-- pre-defined frame in gala. 
        # sc=coords_ICRS.transform_to(gc.GD1Koposov10)
        selected_streamframe = gc.GD1Koposov10
    if orbit=='pa5': #<-- pre-defined frame in gala. 
        # sc = coords_ICRS.transform_to(gc.Pal5PriceWhelan18)
        selected_streamframe = gc.Pal5PriceWhelan18

    if orbit=='c19': # defined by pole/origin in ibata+24, mohammed+26
        alpha_0 = 354.356*u.degree #<-- sets the phi1 zero point

        pole = SkyCoord( #<-- defines the pole
            ra= 81.45*u.degree,
            dec = -6.346*u.degree
        )
        prog_pos = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.kpc
        prog_vel = np.array([row['vx'][0], row['vy'][0], row['vz'][0]]) * u.km/u.s
        prog_coord = paf.galcen_to_ICRS(prog_pos, prog_vel) #<-- debug and print htis out to see if its between the endpoints...
        origin = SkyCoord(
            ra=prog_coord.ra, dec=prog_coord.dec #<-- choose an origin... in ra
        )

        C19Ibata24 = gc.GreatCircleICRSFrame.from_pole_ra0(
            pole=pole, ra0=alpha_0, origin_disambiguate=origin #[coordinate that the origin should be closest to?]
        )
        # sc = coords_ICRS.transform_to(C19Ibata24)
        selected_streamframe = C19Ibata24

    if orbit=='jet': #<-- this also has a pole/origin vibe: from Do+26
        pole = SkyCoord(
            ra=64.983*u.degree, 
            dec=34.747*u.degree
        )
        origin = SkyCoord(
            ra = 138.62*u.degree,
            dec = 22.10*u.degree
        )
        JetDo26 = gc.GreatCircleICRSFrame.from_pole_ra0(
            pole=pole, ra0=origin.ra, origin_disambiguate=origin
        )
        # sc = coords_ICRS.transform_to(JetDo26)
        selected_streamframe = JetDo26

    if orbit=='aau': #<-- endpoints defined for ATLAS in table 1 of Shipp+2018
        endpoints = SkyCoord(
            ra=[9.3, 30.7]*u.degree,
            dec=[-20.9, -33.2]*u.degree
        )
        # pole = SkyCoord( #<-- don't actually need in my definition. 
        #     ra=74.3*u.degree, dec=47.9*u.degree
        # )
        prog_pos = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.kpc
        prog_vel = np.array([row['vx'][0], row['vy'][0], row['vz'][0]]) * u.km/u.s
        prog_coord = paf.galcen_to_ICRS(prog_pos, prog_vel) #<-- debug and print htis out to see if its between the endpoints...
        origin = SkyCoord(
            ra=prog_coord.ra, dec=prog_coord.dec #<-- choose an origin... in ra
        )
        AAUShipp18 = gc.GreatCircleICRSFrame.from_endpoints( #<-- see table 1; for the ATLAS stream. 
            endpoints[0], endpoints[1], ra0=origin.ra#origin=origin, priority='pole'
        )
        # sc = coords_ICRS.transform_to(AAUShipp18)
        selected_streamframe = AAUShipp18


    if orbit=='m3': #<-- endpoints defined in Yang+23, sec 4.4
        # TODO: ics file gives something that is off from this frame; might want to redefine, and settle for "stream on an m3-like orbit." since this is mostly just a high e, low pericenter test. 
        endpoints = SkyCoord(
            ra=[186.45,197.20]*u.degree,
            dec=[19.06,27.76]*u.degree
        )
        # choose the origin based on the ICRS coordinate of the progenitor 
        # at present day. might be smart to read these in from data/FINAL_ics_nolmc.csv rather than hard coding here.
        prog_pos = np.array([row['x'][0], row['y'][0], row['z'][0]]) * u.kpc
        prog_vel = np.array([row['vx'][0], row['vy'][0], row['vz'][0]]) * u.km/u.s
        prog_coord = paf.galcen_to_ICRS(prog_pos, prog_vel) #<-- debug and print htis out to see if its between the endpoints...
        origin = SkyCoord(
            ra=prog_coord.ra, dec=prog_coord.dec #<-- choose an origin... in ra
        )
        M3Yang23 = gc.GreatCircleICRSFrame.from_endpoints(
            endpoints[0], endpoints[1], origin=origin, priority='origin' #<-- warns...
        )
        # sc = coords_ICRS.transform_to(M3Yang23) 
        selected_streamframe = M3Yang23
    
    # else: #<-- why do i get this behavior when orbit='gd1' ??
    #     raise ValueError('acceptable orbits are gd1, pa5, jet, aau, m3, c19. for circular orbits use nominal streamframe.')
    
    sc = coords_ICRS.transform_to(selected_streamframe)

    # coords_stream = {} #<-- a dictionary, like Jake', but will have v_gsr instead of vr, and distance instead of r (different keys. )
    # coords_stream['phi1'] = sc.phi1.to(u.degree).value
    # coords_stream['phi2'] = sc.phi2.to(u.degree).value
    # coords_stream['pm_phi1'] = sc.pm_phi1_cosphi2.to(u.mas/u.yr).value
    # coords_stream['pm_phi2'] = sc.pm_phi2.to(u.mas/u.yr).value
    # coords_stream['distance'] = sc.distance.to(u.kpc).value
    # coords_stream['v_gsr'] = sc.radial_velocity.to(u.km/u.s).value

    return sc, selected_streamframe #<-- returned as a coordinate object. 


def prog_orbit_track(w0, Dt):
    """
    assumes we're in mwp2014. w0 should be a gala phase space position. 
    """
    mwp = gp.BovyMWPotential2014(units=galactic)
    H=gp.Hamiltonian(mwp)
    orbit_forward = H.integrate_orbit(w0, dt=1, n_steps=Dt)
    orbit_backward = H.integrate_orbit(w0, dt=-1, n_steps=Dt)
    # some slick lambda usage from claude:
    P = lambda o: np.array([o.pos.x.to(u.kpc).value,
                            o.pos.y.to(u.kpc).value,
                            o.pos.z.to(u.kpc).value]).T
    V = lambda o: np.array([o.vel.d_x.to(u.km/u.s).value,
                            o.vel.d_y.to(u.km/u.s).value,
                            o.vel.d_z.to(u.km/u.s).value]).T
    orbit_pos = np.vstack([P(orbit_backward)[::-1][:-1], P(orbit_forward)])
    orbit_vel = np.vstack([V(orbit_backward)[::-1][:-1], V(orbit_forward)])
    # ^ note the first half of the orbit integration has the last element 
    # chopped off to avoid duplicating the progenitor position. 
    return orbit_pos*u.kpc, orbit_vel*u.km/u.s

def observed_orbit_track(orbit_pos, orbit_vel, obs_streamframe):
    """
    give orbit_pos and orbit_vel with astropy units. 
    """
    # orbit_pos, orbit_vel = orbit_w[:,:3] * u.kpc, orbit_w[:,3:] * u.kpc/u.Myr
    orbit_icrs = paf.galcen_to_ICRS(orbit_pos, orbit_vel)
    orbit_icrs = reflex_correct(orbit_icrs) #<-- reflex correct the orbit too, or it doesn't match the data.
    sc = orbit_icrs.transform_to(obs_streamframe)
    return sc

def chop_orbit_track(phi1, jump_threshold=45):
    """
    safe-guard against non-strictly-increasing phi1
    in the progenitor orbit interpolant. 
    """
    phi1_track = np.asarray(phi1, dtype=float)
    n = len(phi1_track)
    i_prog = n//2 #<-- the center, where the prog is
    dphi1 = np.diff(phi1_track)
    bad = np.abs(dphi1)>jump_threshold

    ## guard against when things turn around in phi1 without discontinuity. 
    s = np.sign(dphi1[i_prog]) if i_prog < n-1 else np.sign(dphi1[i_prog-1])
    bad = bad | (np.sign(dphi1) != s)


    lo=i_prog
    while lo>0 and not bad[lo-1]:
        lo-=1
    hi=i_prog
    while hi<n-1 and not bad[hi]:
        hi+=1
    idx = np.arange(lo, hi+1)
    if phi1_track[idx[0]] > phi1_track[idx[-1]]: #<-- make sure phi1 is increasing. 
        idx = idx[::-1]   # so it can go straight into np.interp
    return idx

def straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab, Dt=500):
    """
    should _this_ return a dictionary ??? 
    """
    coords_obs, obs_streamframe = streamframe_coords_observed(orbit, CMdict, prog_tab)
    row = prog_tab[prog_tab['name']==orbit]
    w0 = gd.PhaseSpacePosition(
        np.array([row['x'][0], row['y'][0], row['z'][0]])*u.kpc,
        np.array([row['vx'][0], row['vy'][0], row['vz'][0]])*u.km/u.s
        )
    orbit_pos, orbit_vel = prog_orbit_track(w0, Dt)
    orbit_sf = observed_orbit_track(orbit_pos, orbit_vel, obs_streamframe)
    idx = chop_orbit_track(orbit_sf.phi1.to(u.degree).value)



    scd = {} #<--"straight coord dict"
    scd['phi1'] = coords_obs.phi1.to(u.degree).value
    data_x = scd['phi1']
    data_y = [coords_obs.phi2.to(u.degree).value,
              coords_obs.pm_phi1_cosphi2.to(u.mas/u.yr).value,
              coords_obs.pm_phi2.to(u.mas/u.yr).value,
              coords_obs.radial_velocity.to(u.km/u.s).value,
              coords_obs.distance.to(u.kpc).value
              ]

    keys = ['phi2','pm_phi1','pm_phi2','v_gsr','distance']
    orbit_x = orbit_sf.phi1[idx].to(u.degree).value
    orbit_y = [orbit_sf.phi2[idx].to(u.degree).value,
               orbit_sf.pm_phi1_cosphi2[idx].to(u.mas/u.yr).value,
               orbit_sf.pm_phi2[idx].to(u.mas/u.yr).value,
               orbit_sf.radial_velocity[idx].to(u.km/u.s).value,
               orbit_sf.distance[idx].to(u.kpc).value
               ]

    # populate the coordinate dictionary:       
    for ii, key in enumerate(keys):
        scd[key] = data_y[ii] - np.interp(data_x, orbit_x, orbit_y[ii])
    return scd

def outlier_clip(vr, pmphi1, pmphi2):
    """
    i've decided that dvr should be clipped at 100 km/s
    and dpm should be clipped at 1.5 mas/yr
    this avoids biasing the cocoon dispersion 
    with one or two wack stars that really should
    have been taken out with inMW, trim
    """
    outlier_clip = (np.abs(pmphi1)<1.5) & (np.abs(pmphi2)<1.5) & (np.abs(vr)<100)
    return outlier_clip

def poly_straightening(coords, tc=None): #<-- Q; should i be doing in MW, trim first?
    """
    coords should be in dictionary form, 
    include phi1, and phi1 should be the first key
    also, when in dict form my convention is that 
    things don't have astropy units. 
    phi1/phi2 should be in degrees, 
    distance in kpc, 
    proper motions in mas/yr
    v_gsr in km/s
    """
    coords_straighter = {"phi1":coords['phi1']}
    if tc is None:
        tc = [np.ones(len(coords['phi1'])).astype(bool)]*2
    
    for k, key in enumerate(list(coords.keys())):
        if key=='phi1':
            continue
        y = coords[key]
        _, _, poly, fit = paf.straighten_stream_polynomial(coords['phi1'],y,
                                                            trim_criteria = tc,
                                                            degree=5, #<-- default, but choose explicitly here. 
                                                            return_poly_fn=True) #<-- all i care about
        y-=poly(coords_straighter['phi1'], *fit)
        coords_straighter[key] = y

    return coords_straighter

def clip_coords(coords, tc):
    inMW, trim = tc
    out = {}
    for key in coords.keys():
        out[key] = coords[key][inMW][trim]
    return out

def desi_RVerr(zmag, feh=-2.0):
    """
    get RV error for desi data model, which 
    depends on z magnitude and metallicity
    """
    log_err = -0.47 + 0.27*(zmag-16) - 0.23*feh
    return 10**log_err

### TODO: define Via RVerr(mag, metallicity) #<-- from viamock; will need to add to this env (?) check my machine.  

def add_noise(icrs_coords, survey='DESI'):
    """
    TODO: write this function lol. 
    i think in detail the transformation of proper motion errors (and phi1/phi2 errors if we have those)
    is non-trivial and Gala might have functions for transforming the covariance matrix or something. 
    will want to have "DESI" and "Via" options 
    for the survey velocity errors. 
    For Desi these will come from a Koposov paper, 
    for Via they will come from viamock. 
    """
    return


def get_cocoon_selection(coords, cuts):
    """
    assume that coords and cuts are in the same order,
    i.e., the convention of 
    phi2, pmphi1, pmphi2, vgsr
    """
    coords_to_cut = {"phi2":coords['phi2'],"pm_phi1":coords['pm_phi1'],"pm_phi2":coords['pm_phi2'], 'v_gsr':coords['v_gsr']}

    selections = []
    for ii, key in enumerate(coords_to_cut.keys()):
        selection = np.abs(coords_to_cut[key])>cuts[ii]
        selections.append(selection)
    cocoon_selection = np.logical_or.reduce(selections) #<-- this should be true if any of the selections are true. 

    return cocoon_selection
# %%
#### stuff i won't want to run when i import functions to other scripts below. 
# comment out the if __name__==... line and un-indent stuff if doing work in this notebook. 
# if __name__=='__main__': 
orbits = ['circ','gd1','aau','pa5','jet','m3','c19']
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
rvir_index=0
mass_index = 1


dicts = []
sf_coords_obs = []
straight_sf_coords_obs = []

copy_options = [0,1,2,3,4]
for ii, orbit in enumerate(tqdm(orbits)):
    print(orbit)

    (core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
        prepare_nbody_data_anycopy(
            orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
            include_photometry=False
        )

    dicts.append(data_dict)


    if orbit!='circ':
        coords_obs, sf = streamframe_coords_observed(orbit, CMdict, prog_tab)
        sf_coords_obs.append(coords_obs)

        scd = straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab)
        straight_sf_coords_obs.append(scd)

    if orbit=='circ':
        sf_coords_obs.append({})
        straight_sf_coords_obs.append({})

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']#,'distance']

# do cuts as phi2, pmphi1, pmphi2, vgsr
gd1_cuts = [0.75, 0.12,0.04, 1.5]
aau_cuts = [0.5, 0.05, 0.0125, 2]
pa5_cuts = [0.3, 0.02, 0.04, 4.5]
jet_cuts = [0.4, 0.013, 0.005, 2.]
m3_cuts = [1.8, 0.4, 0.15, 4.75]
c19_cuts = [0.32, 0.023, 0.018, 2.75]

orbit_cuts = [
    [], gd1_cuts, aau_cuts, pa5_cuts, jet_cuts, m3_cuts, c19_cuts
]




for ii, sc in enumerate(tqdm(straight_sf_coords_obs)):
    orbit=orbits[ii]
    if ii==0:
        continue
    fig, axs = plt.subplots(len(keys), 2, figsize=[10, 10], width_ratios = [4,1])

    plt.subplots_adjust(hspace=0.03, wspace=0.03)

    fig.suptitle(orbits[ii])
    
    cmdict = dicts[ii]['CoM']
    inMW, trim = cmdict['inMW'], cmdict['trim']

    unbound = ~cmdict['in_rtid']
    unbound = unbound[inMW][trim]

    trimmed_sc = clip_coords(sc, [inMW, trim])

    sc_straighter = poly_straightening(trimmed_sc)


    ol_clip = outlier_clip(
        sc_straighter['v_gsr'], sc_straighter['pm_phi1'], sc_straighter['pm_phi2']
    )

    cocoon_clips = orbit_cuts[ii]
    cocoon_selection = get_cocoon_selection(sc_straighter, cocoon_clips)

    key_labels = [
        r'$\phi_2~[\degree]$',
        r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$',
        r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$',
        r'$v_{\rm GSR}~[\rm km~s^{-1}]$'
    ]
    for jj, key in enumerate(keys):
        cut = cocoon_clips[jj]
        # ax_row = axs[jj]

        ax = axs[jj,0]
        ax.scatter(sc_straighter['phi1'][ol_clip & unbound], 
                   sc_straighter[key][ol_clip & unbound],
                    c='k', s=.1,
                    rasterized=True) 

        ax.scatter(sc_straighter['phi1'][ol_clip & unbound & cocoon_selection], 
                   sc_straighter[key][ol_clip & unbound & cocoon_selection],
                    c=cc[-1], s=20, edgecolor='k', lw=0.5,
                    rasterized=True) 

        ax.axhline(cut, c='k', lw=1)
        ax.axhline(-cut, c='k', lw=1)
        ax.set_ylim(-3*cut, 3*cut)
        ax.set_ylabel(key_labels[jj], fontsize=15)


        ax = axs[jj,1]
        bins = np.linspace(-3*cut, 3*cut, 50)
        ax.hist(sc_straighter[key][ol_clip & unbound & ~cocoon_selection], 
                alpha=0.2, density=True, color='k',orientation='horizontal',
                bins=bins)
        ax.hist(sc_straighter[key][ol_clip & unbound & cocoon_selection],
                histtype='step', density=True, lw=2, 
                color=cc[-1],orientation='horizontal',
                bins=bins)

        ax.set_yticklabels([])
        ax.set_xticks
        if jj<3:
            # print("REMOVING TICK LABLES>>>>>")
            axs[jj,0].set_xticklabels([])
            axs[jj,1].set_xticklabels([])


    axs[-1,0].set_xlabel(r'$\phi_1~[\degree]$')
    axs[-1,1].set_xlabel(r'density')

    # plt.savefig('/n/home02/amphillips/p27_nbody/plots/cocoon_separation/%s/%s_%.2f.pdf'%(masses[mass_index],orbit, rvirs[rvir_index]), 
    #             dpi=300, bbox_inches='tight')


# %%
