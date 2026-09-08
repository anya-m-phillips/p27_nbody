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
from scipy.special import expit, logit # inverse-logit / logit, for the f_1 reparameterization


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
# copy_options = [0,1,2,3,4]
copy_options = [4,3,2,1,0]

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']
# %%
def gmm_likelihood_multivariate(x_data, 
                   f_1, Mu_1, Sigma_1, 
                   Mu_2, Sigma_2
                   ):
    """
    Mu_1, Mu_2 should be 4-vectors (phi2, transverse+radial velocities)
    Sigma_1, Sigma_2 should I guess be covariance matrices. yikes. 
    """
    Q_1 = f_1
    Q_2 = 1 - f_1 
    component_1 = Q_1 * multivariate_normal.pdf(x_data, mean=Mu_1, cov=Sigma_1)
    component_2 = Q_2 * multivariate_normal.pdf(x_data, mean=Mu_2, cov=Sigma_2)
    li = component_1 + component_2

    ln_li = np.log(li)
    ln_L = np.sum(ln_li)

    return ln_L # idk how to test. 
# except that idk what a covariance matrix is so let's take a different approach:
# note this "simple" version IS the multivariate one with a diagonal covariance
# matrix, i.e. sigma = np.diag(sigma_vec**2). no correlations between phi2 and the
# velocities within a single component. that's an assumption, not a hack, and it's
# the same one jarvis+26 make.
def component_likelihood(x_data, mu, sigma):
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


def gmm_likelihood_simple(x_data,
                          f_1,
                          mu_1, sigma_1,
                          mu_2, sigma_2):
    """
    total log likelihood of a two-component (thin + cocoon) mixture.

    now mu's are still 4-vectors, but so are sigmas.

    the mixture weights go INSIDE the product over stars, not outside:
        L = prod_i [ Q_1 p_1(x_i) + Q_2 p_2(x_i) ]
    each star independently belongs to thin OR cocoon. the other ordering,
        L = Q_1 prod_i p_1(x_i) + Q_2 prod_i p_2(x_i),
    is a different (wrong) model: it says the WHOLE stream is thin or the whole
    stream is cocoon, and it will collapse onto whichever component wins overall.

    in log space the per-star mixing is the only place we'd have to exponentiate,
    and logaddexp does it stably (it factors out the larger term):
        ln L = sum_i logaddexp( ln Q_1 + ln p_1(x_i), ln Q_2 + ln p_2(x_i) )
    """
    if not (0 < f_1 < 1):
        return -np.inf
    Q_1 = f_1
    Q_2 = 1 - f_1

    ln_p1 = component_likelihood(x_data, mu_1, sigma_1) # (N,)
    ln_p2 = component_likelihood(x_data, mu_2, sigma_2) # (N,)

    # per-star mixture, still in logs:
    ln_li = np.logaddexp(np.log(Q_1) + ln_p1,
                         np.log(Q_2) + ln_p2) # (N,)

    # the "outer" product over stars -> a sum:
    return np.sum(ln_li) #<-- eventually will want to minimize the negative of this.

def gmm_negative_loglikelihood(f_1,
                                mu_1, sigma_1,
                                mu_2, sigma_2,
                                x_data):
    """
    same as gmm_likelihood_simple() but returns negative.
    """
    if not (0 < f_1 < 1):
        return np.inf #<-- POSITIVE inf here. this is what we're MINIMIZING, so an
                      #    invalid f_1 has to look BAD (+inf), not infinitely good.
    Q_1 = f_1
    Q_2 = 1 - f_1

    ln_p1 = component_likelihood(x_data, mu_1, sigma_1) # (N,)
    ln_p2 = component_likelihood(x_data, mu_2, sigma_2) # (N,)

    # per-star mixture, still in logs:
    ln_li = np.logaddexp(np.log(Q_1) + ln_p1,
                         np.log(Q_2) + ln_p2) # (N,)

    # the "outer" product over stars -> a sum:
    return -np.sum(ln_li)


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

def pack_params(f_1, mu_1, sigma_1, mu_2, sigma_2):
    """natural parameters -> flat (4K+1,) unconstrained vector."""
    return np.concatenate([[logit(f_1)],
                           mu_1, np.log(sigma_1),
                           mu_2, np.log(sigma_2)])


def unpack_params(theta, K=4):
    """flat unconstrained vector -> natural parameters. inverse of pack_params."""
    f_1     = expit(theta[0])
    mu_1    = theta[1       : 1 +   K]
    sigma_1 = np.exp(theta[1 +   K : 1 + 2*K])
    mu_2    = theta[1 + 2*K : 1 + 3*K]
    sigma_2 = np.exp(theta[1 + 3*K : 1 + 4*K])
    return f_1, mu_1, sigma_1, mu_2, sigma_2


def nll_flat(theta, x_data):
    """the objective actually handed to minimize()."""
    return gmm_negative_loglikelihood(*unpack_params(theta), x_data)


def sort_components(f_1, mu_1, sigma_1, mu_2, sigma_2, sort_dim=0):
    """
    a mixture model has no idea which component you meant to call "thin" -- the
    likelihood is identical if you swap 1<->2 and send f_1 -> 1-f_1 (this is
    "label switching"). so we impose the convention ourselves: component 1 is the
    NARROWER one in sort_dim (default 0 = phi2), i.e. the thin stream.
    """
    if sigma_1[sort_dim] <= sigma_2[sort_dim]:
        return f_1, mu_1, sigma_1, mu_2, sigma_2
    return 1 - f_1, mu_2, sigma_2, mu_1, sigma_1


def membership_probability(x_data, f_1, mu_1, sigma_1, mu_2, sigma_2):
    """
    posterior probability that each star belongs to component 1 (the thin stream),
    a.k.a. the "responsibility". this is the thing that actually does the cocoon
    separation -- p_cocoon = 1 - this.

        p_1(i) = Q_1 p_1(x_i) / [ Q_1 p_1(x_i) + Q_2 p_2(x_i) ]

    computed as a difference of logs so it never overflows.

    returns : (N,) array in [0, 1]
    """
    Q_1 = f_1
    Q_2 = 1 - f_1

    ln_w1 = np.log(Q_1) + component_likelihood(x_data, mu_1, sigma_1)
    ln_w2 = np.log(Q_2) + component_likelihood(x_data, mu_2, sigma_2)

    return np.exp(ln_w1 - np.logaddexp(ln_w1, ln_w2))

# %%
ii = 0 # <--- gd1 i think. 
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
mu_1, sigma_1 = np.zeros(4), 0.3 * sd  # thin: narrower than the data
mu_2, sigma_2 = np.zeros(4), 1.5 * sd  # cocoon: broader than the data
f_1 = 0.99                           # let the data decide, don't start at 0.999. this is the thin stream fraction. 

# print(gmm_negative_loglikelihood(f_1, mu_1, sigma_1, mu_2, sigma_2, x_data)) # test function

# %%
# minimize() passes ONE flat array as the first argument and `args` as a TUPLE of
# everything after it -- `args=x_data` (no comma) gets iterated and splatted.
theta0 = pack_params(f_1, mu_1, sigma_1, mu_2, sigma_2)

#### note this step takes a long time: 
result = minimize(nll_flat, x0=theta0, args=(x_data,), method='Nelder-Mead',
                  options={'maxiter': 100000, 'maxfev': 100000,
                           'fatol': 1e-6, 'xatol': 1e-6})

f_1_fit, mu_1_fit, sigma_1_fit, mu_2_fit, sigma_2_fit = \
    sort_components(*unpack_params(result.x))

print(result.success, result.message, '\nnll =', result.fun)
print('f_thin  =', np.round(f_1_fit, 4))
for k, m1, s1, m2, s2 in zip(keys, mu_1_fit, sigma_1_fit, mu_2_fit, sigma_2_fit):
    print(f'  {k:>8}   thin mu={m1:9.4f} sig={s1:8.4f} | cocoon mu={m2:9.4f} sig={s2:8.4f}')
# %%
# per-star cocoon membership -- the soft replacement for get_cocoon_selection()
p_thin = membership_probability(x_data, f_1_fit, mu_1_fit, sigma_1_fit,
                                mu_2_fit, sigma_2_fit)
p_cocoon = 1 - p_thin
# %%
for k in range(4):
    fig, axs = plt.subplots(1,2, figsize=[10,3], width_ratios=[3,1])
    x = sc_straighter['phi1'][unbound & ol_clip]
    y = x_data[:,k]

    order = np.argsort(p_cocoon)

    axs[0].scatter(x[order], y[order], c=p_thin[order], s=5, cmap='cool', rasterized=True)

    ts = p_thin>0.5
    cn = ~ts
    bins = np.linspace(min(y), max(y), 50)
    axs[1].hist(y[ts], bins=bins, histtype='step', color='magenta', orientation='horizontal', lw=3, density=True)
    axs[1].hist(y[cn], bins=bins, histtype='step', color='cyan', orientation='horizontal', lw=3, density=True)
# %%
