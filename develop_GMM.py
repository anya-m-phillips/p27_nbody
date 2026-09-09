#--------------------------------------------------------------#
#   in inspect_new_sims.py i defined a bunch of functions.     #
#   that will help me define stream frames and make good cuts  # 
#   on the data; here i will use those functions and 
#   develop cocoon separation using a GMM instead
#   of hard cuts. inspired by jarvis+2026, we will start
#   with a maximized likelihood of the GMM, then use 
#   emcee or similar to sample posteriors. hopefully this
#   will be more robust/useful for down-sampled "observed" 
#   data.                                                      #
#--------------------------------------------------------------#
# %%
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import petar
import numpy as np

from scipy.stats import binned_statistic, norm, multivariate_normal
from scipy.optimize import curve_fit, minimize
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
import inspect_new_sims as simspect # lol idk. but i need basically all the funcitons. 

# first import... very beefy...
from pygaia.errors.astrometric import parallax_uncertainty, proper_motion_uncertainty, total_proper_motion_uncertainty, total_position_uncertainty
# %%
#### this requires inputting a full covariance matrix. 
# def gmm_likelihood_multivariate(x_data, 
#                    f_1, Mu_1, Sigma_1, 
#                    Mu_2, Sigma_2
#                    ):
#     """
#     Mu_1, Mu_2 should be 4-vectors (phi2, transverse+radial velocities)
#     Sigma_1, Sigma_2 should I guess be covariance matrices. yikes. 
#     """
#     Q_1 = f_1
#     Q_2 = 1 - f_1 
#     component_1 = Q_1 * multivariate_normal.pdf(x_data, mean=Mu_1, cov=Sigma_1)
#     component_2 = Q_2 * multivariate_normal.pdf(x_data, mean=Mu_2, cov=Sigma_2)
#     li = component_1 + component_2

#     ln_li = np.log(li)
#     ln_L = np.sum(ln_li)

#     return ln_L # idk how to test. 
# except that idk what a covariance matrix is so let's take a different approach:
# note this "simple" version IS the multivariate one with a diagonal covariance
# matrix, i.e. sigma = np.diag(sigma_vec**2). no correlations between phi2 and the
# velocities within a single component. that's an assumption, not a hack, and it's
# the same one jarvis+26 make.


def component_likelihood(x_data, mu, sigma): #<-- this can stay the same...
    """
    log likelihood of each star under ONE gaussian component, with independent
    (uncorrelated) phase space dimensions.

    x_data : (N, K) array -- N stars, K phase space coords (phi2, pm_phi1, pm_phi2, v_gsr)
    mu     : (K,) array of means
    sigma  : (K,) array of standard deviations (NOT variances)

    returns : (N,) array of ln p(x_i | mu, sigma)

    the product over phase space dimensions,
        p(x_i) = prod_k N(x_ik | mu_k, sigma_k),
    becomes a SUM over k in log space -- this is the "inner" product, and it
    collapses to axis=1 here. broadcasting does the k loop for us.
    """
    mu = np.atleast_1d(mu)
    sigma = np.atleast_1d(sigma)

    if np.any(sigma <= 0) or not np.all(np.isfinite(sigma)):
        return np.full(x_data.shape[0], -np.inf) #<-- so an optimizer walking into
                                                 #    negative sigma gets rejected
                                                 #    instead of getting a nan.

    return norm.logpdf(x_data, loc=mu, scale=sigma).sum(axis=1)


### for the optimizer. 
def gmm_negative_loglikelihood(component_fractions, means, sigmas, data):
    """
    all inputs should be numpy arrays
    component fractions should have length (n_components - 1) -- the last weight
    is implicit, fixed by sum(Q) = 1, so only n-1 of them are free
    means, sigmas should have length n_components x four dimensions
    """
    ### number of components
    n_components = len(component_fractions)+1 #<-- provide the first n-1 components, last one is so that they sum to 1
    if (len(means)!= n_components) or (len(sigmas)!= n_components):
        raise ValueError("len(means)=%i, len(sigmas)=%i, ncomponents=%i"%(len(means), len(sigmas), n_components))

    # making sure we have valid component fractions
    if np.any(component_fractions <= 0) or np.sum(component_fractions) >= 1:
        return np.inf   # POSITIVE inf: we are MINIMIZING, so invalid must look BAD.

    # checking for positive sigmas.
    if np.any(sigmas <= 0) or not np.all(np.isfinite(sigmas)):
        return np.inf

    
    ln_ps = np.array([component_likelihood(data, means[j], sigmas[j])
                      for j in range(n_components)]) #<-- this is the likelihood over the four dimensions

    Qs = np.append(component_fractions, 1 - np.sum(component_fractions))


    ln_Li=logsumexp(np.log(Qs[:,None])+ln_ps, axis=0)

    return -np.sum(ln_Li)
    

#--- packing, so scipy.optimize.minimize can see the parameters ----------------#
# minimize() wants ONE flat 1-D array as the first argument to the objective.
# it cannot take (scalar, 4-vector, 4-vector, 4-vector, 4-vector) -- numpy turns
# that into a ragged object array and it dies. so we flatten to (17,) here.
#
# we also reparameterize while we're at it, so the fit is UNCONSTRAINED:
#   f_1   -> logit(f_1),  so any real number maps back into (0, 1)
#   sigma -> log(sigma),  so any real number maps back to sigma > 0
# this means the optimizer can never step into invalid territory and hit the
# +inf walls above, which is what wrecks a finite-difference gradient. it also
# fixes the conditioning: phi2 is ~0.1 deg while v_gsr is ~10 km/s, and in logs
# those are comparable steps.

# with n_components free, the (0,1) reparameterization of f_1 generalizes to a
# MULTINOMIAL LOGIT (softmax): the last component's alpha is pinned at 0 as the
# reference, so n-1 free reals map onto the interior of the simplex, i.e. all
# Q_j > 0 and sum(Q_j) = 1, automatically. for n_components=2 this reduces to
# exactly logit/expit, and the theta layout below is bit-identical to the old
# 2-component one, so an old result.x still unpacks correctly.
#
# layout, length (n-1) + 2nK:
#   [alpha_1 ... alpha_{n-1},  mu_1, ln sigma_1,  mu_2, ln sigma_2,  ...]
# the mu/ln-sigma pairs stay blocked BY COMPONENT (not all mus then all sigmas)
# to preserve that backward compatibility.

def pack_params(component_fractions, means, sigmas):
    """
    natural parameters -> flat ((n-1) + 2nK,) unconstrained vector.

    component_fractions : (n_components - 1,) -- the last weight is implicit
    means, sigmas       : (n_components, K)
    """
    fracs  = np.atleast_1d(np.asarray(component_fractions, dtype=float))
    means  = np.asarray(means,  dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)

    n_components = len(fracs) + 1
    if (len(means) != n_components) or (len(sigmas) != n_components):
        raise ValueError("len(means)=%i, len(sigmas)=%i, ncomponents=%i"
                         % (len(means), len(sigmas), n_components))
    if np.any(fracs <= 0) or np.sum(fracs) >= 1:
        raise ValueError("component_fractions must all be > 0 and sum to < 1; "
                         "got %s (sum %.6f)" % (fracs, np.sum(fracs)))
    if np.any(sigmas <= 0):
        raise ValueError("sigmas must be positive; got min %g" % np.min(sigmas))

    # alpha_j = ln(Q_j / Q_last). for n=2 this is ln(f_1/(1-f_1)) = logit(f_1).
    alpha = np.log(fracs) - np.log(1 - np.sum(fracs))
    blocks = np.hstack([means, np.log(sigmas)])          # (n_components, 2K)
    return np.concatenate([alpha, blocks.ravel()])


def unpack_params(theta, K=4):
    """
    flat unconstrained vector -> natural parameters. inverse of pack_params.

    n_components is INFERRED from the length: len(theta) = (n-1) + 2nK, so
    n = (len(theta) + 1) / (2K + 1). a length that isn't of that form is a
    K/n_components mismatch and raises rather than silently misreshaping.
    """
    theta = np.asarray(theta, dtype=float)

    n_components, remainder = divmod(len(theta) + 1, 2 * K + 1)
    if remainder != 0 or n_components < 2:
        raise ValueError("len(theta)=%i is not (n-1) + 2nK for any n >= 2 with "
                         "K=%i" % (len(theta), K))
    n_free = n_components - 1

    # softmax with the last component pinned at alpha=0. subtracting logsumexp
    # keeps it stable for large |alpha| -- no overflow, and the Qs sum to 1 by
    # construction rather than by luck.
    alpha_full = np.append(theta[:n_free], 0.0)
    Qs = np.exp(alpha_full - logsumexp(alpha_full))

    blocks = theta[n_free:].reshape(n_components, 2 * K)
    means  = blocks[:, :K]
    sigmas = np.exp(blocks[:, K:])

    # return only the n-1 free fractions, to match gmm_negative_loglikelihood's
    # signature -- the last one is always 1 - sum(the rest).
    return Qs[:-1], means, sigmas



def nll_flat(theta, x_data):
    """
    the objective actually handed to minimize().

    K is read off x_data.shape[1] rather than defaulted, so theta can never be
    unpacked against the wrong number of phase space dimensions.
    """
    return gmm_negative_loglikelihood(*unpack_params(theta, K=x_data.shape[1]),
                                      x_data)


def sort_components(component_fractions, means, sigmas, sort_dim=-1):
    """
    a mixture model has no idea which component you meant to call "thin" -- the
    likelihood is exactly invariant under relabelling the components (this is
    "label switching"), so we impose the convention ourselves: components come
    back ordered NARROWEST -> WIDEST, so component 0 is the thinnest part of the
    stream and the LAST one is the cocoon.

    sort_dim picks which phase space dimension sets the ordering, since a
    component can be the widest in one coordinate and not another. default -1 =
    v_gsr given keys = ['phi2','pm_phi1','pm_phi2','v_gsr'].

    component_fractions : (n_components - 1,) -- the last weight is implicit
    means, sigmas       : (n_components, K)

    NOTE the implicit weight is re-derived after sorting, so the returned
    component_fractions are the n-1 NARROWEST components and the cocoon fraction
    is 1 - sum(returned). that is the number you want to report.
    """
    fracs  = np.atleast_1d(np.asarray(component_fractions, dtype=float))
    means  = np.asarray(means,  dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)

    Qs = np.append(fracs, 1 - np.sum(fracs))       # reconstruct the full weights
    order = np.argsort(sigmas[:, sort_dim])        # narrowest -> widest
    return Qs[order][:-1], means[order], sigmas[order]


def membership_probability(x_data, component_fractions, means, sigmas, sort_dim=-1):
    """
    posterior probability that each star belongs to the THIN STREAM, defined as
    any of the n-1 narrower components -- i.e. everything except the final
    (widest) one, which is the cocoon. p_cocoon = 1 - this.

        p_thin(i) = sum_{j < n} Q_j p_j(x_i) / sum_{j} Q_j p_j(x_i)

    this is the thing that actually does the cocoon separation: a soft per-star
    weight rather than a boolean cut. with n_components > 2 the extra narrow
    components soak up the non-gaussian shape of the thin stream (epicyclic
    feathers etc.) instead of being mistaken for a cocoon.

    computed as a difference of logs so it never overflows.

    ORDER MATTERS: "the last component is the cocoon" is only true if the
    components are sorted, so this checks and refuses rather than quietly
    reporting the wrong population. run sort_components first.

    returns : (N,) array in [0, 1]
    """
    fracs  = np.atleast_1d(np.asarray(component_fractions, dtype=float))
    means  = np.asarray(means,  dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)

    n_components = len(fracs) + 1
    if n_components < 2:
        raise ValueError("need at least 2 components to separate a cocoon")

    widths = sigmas[:, sort_dim]
    if np.any(np.diff(widths) < 0):
        raise ValueError("components are not sorted narrowest -> widest in "
                         "sort_dim=%i (widths %s) -- call sort_components first"
                         % (sort_dim, widths))

    Qs = np.append(fracs, 1 - np.sum(fracs))

    # (n_components, N) -- ln Q_j + ln p_j(x_i)
    ln_w = np.array([np.log(Qs[j]) + component_likelihood(x_data, means[j], sigmas[j])
                     for j in range(n_components)])

    # numerator drops the last (widest) component; denominator keeps everything.
    return np.exp(logsumexp(ln_w[:-1], axis=0) - logsumexp(ln_w, axis=0))



# %%
#### START REWRITING FOR 3 COMPONENTS HERE: 
# define stuff about the grid:
grid_info = paf.extended_grid_info(scratch=True) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')


orbits = ['gd1','aau','pa5','jet','m3','c19']
init_displacements = [
    grid_info.gd1_init_displacement, 
    grid_info.aau_init_displacement,
    grid_info.pa5_init_displacement,
    grid_info.jet_init_displacement,
    grid_info.m3_init_displacement,
    grid_info.c19_init_displacement
]
masses = ['lm','hm']
rvirs = [0.75, 1.5, 3, 6]
copy_options = [0,1,2,3,4]
# copy_options = [4,3,2,1,0]

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']

ii = 2 # <--- 0 for gd1 i think. 
orbit = orbits[ii]
## do the orbit-wise check -- integrate prog orbit and find the pericenter. 
init_displacement = init_displacements[ii]
orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
peri = orbit_obj.pericenter()

mass_index = 1
rvir_index=0
(core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
    simspect.prepare_nbody_data_anycopy(
        orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
        include_photometry=False
    )

coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab) # observed frame
sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) # straight coords


# clip the straight coords before putting into polynomial straightening step
trimmed_sc = simspect.clip_coords(sc, [inMW, trim]) #<-- this applies inMW, trim to the coordinate dictionary
sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction


# select unbound + outlier clipping 
unbound = ~CMdict['in_rtid']
unbound = unbound[inMW][trim]
ol_clip = simspect.outlier_clip( #<-- avoid biasing the cocoon dispersion with a few crazy outliers. 
            sc_straighter['v_gsr'], sc_straighter['pm_phi1'], sc_straighter['pm_phi2'] 
        )

# %%
### for right now i am blindly putting these into the gmm, eventually will need 
# to convert proper motions -> transverse velocities and _then_ straighten before 
# doing the mixture modeling step. 
# (N, 4), column order set by `keys` = ['phi2','pm_phi1','pm_phi2','v_gsr'].
# mu/sigma have to be in this same order. `unbound` drops the bound progenitor so
# it doesn't get swallowed by the thin component.
x_data = np.column_stack([sc_straighter[k][unbound & ol_clip] for k in keys])

# --- initial guess -----------------------------------------------------------#
# NOT mu=0, sigma=1 for both components. two reasons:
#  1. identical components are an exact saddle point of the likelihood. if
#     p_1 == p_2 then every star has responsibility f_1 regardless of f_1, the
#     gradient wrt f_1 is exactly zero, and the two components can never split.
#     the initial guess HAS to break the symmetry.
#  2. sigma=1 is meaningless across mixed units -- 1 deg in phi2 is the whole
#     stream width, 1 km/s in v_gsr is nothing. scale off the data instead.
sd = x_data.std(axis=0)
mu_1, sigma_1 = np.zeros(4), 0.1 * sd  # thin: narrower than the data
mu_2, sigma_2 = np.zeros(4), 1.0 * sd  
mu_3, sigma_3 = np.zeros(4), 3. * sd   # cocoon: broader than the data
f_1 = 1/3                           # starting at 0.5 would "let the data decide." this is the thin stream fraction. 
f_2 = 1/3 # - 0.01

fracs_0 = np.array([f_1, f_2])
means_0 = np.array([mu_1, mu_2, mu_3])
sigmas_0 = np.array([sigma_1, sigma_2, sigma_3])
print(gmm_negative_loglikelihood(component_fractions = np.array([f_1, f_2]),
                                 means = means_0, 
                                 sigmas = sigmas_0,
                                 data=x_data
                                 )
)
# print(gmm_negative_loglikelihood(f_1, mu_1, sigma_1, mu_2, sigma_2, x_data)) # test function
# %%

# %%
# minimize() passes ONE flat array as the first argument and `args` as a TUPLE of
# everything after it -- `args=x_data` (no comma) gets iterated and splatted.
theta0 = pack_params(fracs_0, means_0, sigmas_0)

result = minimize(nll_flat, x0=theta0, args=(x_data,), method='Powell', # method='Nelder-Mead',
                  options={'maxiter': 100000, 'maxfev': 100000,
                           'fatol': 1e-6, 'xatol': 1e-6})

fracs_fit, means_fit, sigmas_fit = sort_components(*unpack_params(result.x, K=x_data.shape[1]))
p_thin = membership_probability(x_data, fracs_fit, means_fit, sigmas_fit)
cocoon_fraction = 1 - fracs_fit.sum()

print(result.success, result.message, '\nnll =', result.fun)
print("cocoon fraction:", cocoon_fraction)
print(fracs_fit)
print(sigmas_fit[:,-1])


# %%
for k in range(4):
    fig, axs = plt.subplots(1,2, figsize=[10,3], width_ratios=[3,1], sharey=True)
    x = sc_straighter['phi1'][unbound & ol_clip]
    y = x_data[:,k]

    order = np.argsort(1-p_thin)

    axs[0].scatter(x[order], y[order], c=p_thin[order], s=5, cmap='winter', rasterized=True)

    ts = p_thin>0.5
    cn = ~ts
    bins = np.linspace(min(y), max(y), 50)
    axs[1].hist(y[ts], bins=bins, histtype='step', color='magenta', orientation='horizontal', lw=3, density=True)
    axs[1].hist(y[cn], bins=bins, histtype='step', color='cyan', orientation='horizontal', lw=3, density=True)


# %%
#
#
#
#
#
#
#
#
#
############ RUN A LOOP OVER THE GRID. DO THINGS BEHAVE AS EXPECTED?? ############
grid_info = paf.extended_grid_info(scratch=True) 
lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')


# orbits = ["m3","pa5","aau","c19","gd1","jet"] #<-- ordered by pericenter. 
# orbits = ['gd1','pa5','m3']
# init_displacements = [
#     grid_info.gd1_init_displacement,
#     grid_info.pa5_init_displacement,
#     grid_info.m3_init_displacement
# ]


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
# copy_options = [0,1,2,3,4]
copy_options = [4,3,2,1,0]

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']

dicts = []
sf_coords_obs = [] # before straightening
straight_sf_coords_obs = [] 

f_cocoons = []
vgsr_dispersions = []
phi2_dispersions = []


pericenters_kpc = []


for ii, orbit in enumerate(tqdm(orbits)):

    ## do the orbit-wise check -- integrate prog orbit and find the pericenter. 
    init_displacement = init_displacements[ii]
    orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
    peri = orbit_obj.pericenter()
    pericenters_kpc.append(peri.to(u.kpc).value)

    ### eventually eventually another inner loop will go here for masses.
    f_cocoons_this_orbit = []
    vgsr_dispersions_this_orbit = []
    phi2_dispersions_this_orbit = []

    mass_index = 1 # <-- high mass stellar population... 
    for rvir_index in range(4):
        (core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
            simspect.prepare_nbody_data_anycopy(
                orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
                include_photometry=False
            )

        # dicts.append(data_dict)

        coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab)
        # sf_coords_obs.append(sf_coords_obs)

        # straightened coords
        sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab)

        unbound = ~CMdict['in_rtid']
        unbound = unbound[inMW][trim]

        trimmed_sc = simspect.clip_coords(sc, [inMW, trim]) #<-- this applies inMW, trim to the coordinate dictionary
        sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction
        ol_clip = simspect.outlier_clip( #<-- avoid biasing the cocoon dispersion with a few crazy outliers. 
            sc_straighter['v_gsr'], sc_straighter['pm_phi1'], sc_straighter['pm_phi2'] 
        )

        use = ol_clip & unbound

        #### assemble the data and perform the fit: 
        x_data = np.column_stack([sc_straighter[k][unbound & ol_clip] for k in keys])

        sd = x_data.std(axis=0)
        mu_1, sigma_1 = np.zeros(4), 0.1 * sd  # thin: narrower than the data
        mu_2, sigma_2 = np.zeros(4), 1.0 * sd  
        mu_3, sigma_3 = np.zeros(4), 10 * sd   # cocoon: broader than the data
        f_1 = 1/3                           # starting at even groups would "let the data decide." 
        f_2 = 1/3 # - 0.01

        fracs_0 = np.array([f_1, f_2])
        means_0 = np.array([mu_1, mu_2, mu_3])
        sigmas_0 = np.array([sigma_1, sigma_2, sigma_3])
        theta0 = pack_params(fracs_0, means_0, sigmas_0)

        result = minimize(nll_flat, x0=theta0, args=(x_data,), method='Powell', # method='Nelder-Mead',
                        options={'maxiter': 100000, 'maxfev': 100000,
                                'fatol': 1e-6, 'xatol': 1e-6})

        fracs_fit, means_fit, sigmas_fit = sort_components(*unpack_params(result.x, K=x_data.shape[1]))
        p_thin = membership_probability(x_data, fracs_fit, means_fit, sigmas_fit)
        p_cocoon = 1 - p_thin

        print(result.success, result.message, '\nnll =', result.fun) #<-- verbose? 



        # for now let's say I just care about the cocoon fraction... 
        # f_cocoon = len(p_thin[p_thin<0.5]) / len(p_thin)
        f_cocoon = 1-np.sum(fracs_fit)
        f_cocoons_this_orbit.append(f_cocoon)

        sigphi2, sigpmphi1, sigpmphi2, sigvgsr = sigmas_fit[-1] #<-- cocoon component. thin stream is components 0 and 1
        vgsr_dispersions_this_orbit.append(sigvgsr)
        phi2_dispersions_this_orbit.append(sigphi2)


    f_cocoons.append(f_cocoons_this_orbit)
    vgsr_dispersions.append(vgsr_dispersions_this_orbit)
    phi2_dispersions.append(phi2_dispersions_this_orbit)


# %%
pericenters_kpc = np.array(pericenters_kpc)
reordered = np.argsort(pericenters_kpc)
pericenters_kpc = pericenters_kpc[reordered]
f_cocoons = np.array(f_cocoons)[reordered]
vgsr_dispersions = np.array(vgsr_dispersions)[reordered]
phi2_dispersions = np.array(phi2_dispersions)[reordered]

orbits = np.array(orbits)[reordered]

ccc = cc[1:]
fig, axs = plt.subplots(1,3,figsize=[21,7], sharex=True)
for ii, orbit in enumerate(tqdm(orbits)):
    f_cocoons_this_orbit = f_cocoons[ii]
    phi2_dispersions_this_orbit = phi2_dispersions[ii]
    vgsr_dispersions_this_orbit = vgsr_dispersions[ii]

    x = rvirs
    axs[0].plot(x, f_cocoons_this_orbit, label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[ii],
                marker='o', color=ccc[ii], markersize=10)


    axs[1].plot(x, phi2_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10)
    axs[2].plot(x, vgsr_dispersions_this_orbit, marker='o', color=ccc[ii], markersize=10)
axs[0].legend(loc='upper right')
axs[0].legend(loc='upper right')
for ax in axs:
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$R_{\rm vir, 0}~[\rm pc]$')
axs[0].set_ylabel(r'$f_{\rm cocoon}$')
axs[1].set_ylabel(r'$\sigma_{\phi_2, \rm cocoon}~[\degree]$')
axs[2].set_ylabel(r'$\sigma_{v_{\rm GSR, cocoon}}~[\rm km~s^{-1}]$')

 