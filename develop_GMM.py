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


def unpack_params(theta, K=4, min_components=2):
    """
    flat unconstrained vector -> natural parameters. inverse of pack_params.

    n_components is INFERRED from the length: len(theta) = (n-1) + 2nK, so
    n = (len(theta) + 1) / (2K + 1). a length that isn't of that form is a
    K/n_components mismatch and raises rather than silently misreshaping.

    min_components stays at 2 by default because for a MIXTURE the n=1 case is
    almost always a wiring mistake (an empty fractions array that got flattened
    away somewhere) rather than something you meant. pass min_components=1 when
    you really do want the single-component model -- the AIC/BIC null -- which is
    what nll_flat_anyn does.
    """
    theta = np.asarray(theta, dtype=float)

    n_components, remainder = divmod(len(theta) + 1, 2 * K + 1)
    if remainder != 0 or n_components < min_components:
        raise ValueError("len(theta)=%i is not (n-1) + 2nK for any n >= %i with "
                         "K=%i" % (len(theta), min_components, K))
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



def nll_flat(theta, x_data, min_components=2):
    """
    the objective actually handed to minimize().

    K is read off x_data.shape[1] rather than defaulted, so theta can never be
    unpacked against the wrong number of phase space dimensions.

    min_components is passed straight through to unpack_params. leave it at 2 for
    a mixture, so an accidentally-emptied fractions array still fails loudly.
    pass min_components=1 for the single multivariate (diagonal) gaussian -- the
    null model for the AIC/BIC comparison -- where component_fractions really is
    empty. nothing downstream needs special-casing for that: pack_params turns
    (np.array([]), (1,K) means, (1,K) sigmas) into a length-2K theta with no
    weights in it, unpack_params sends it back to an empty fractions array, and
    gmm_negative_loglikelihood's validity guards are vacuously true on an empty
    array (np.any of nothing is False, np.sum of nothing is 0 < 1) so Qs = [1.0],
    ln Q = 0, and the single-row logsumexp is the identity -- what comes out is
    just the plain single-gaussian log likelihood.

    NOTE minimize() passes this via args, so the single-component call is
    args=(x_data, 1) -- the trailing comma rule applies as always.
    """
    return gmm_negative_loglikelihood(*unpack_params(theta, K=x_data.shape[1],
                                                     min_components=min_components),
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


def component_responsibilities(x_data, component_fractions, means, sigmas):
    """
    the FULL responsibility matrix: posterior probability that star i was drawn
    from component j, for every component.

        R[j, i] = Q_j p_j(x_i) / sum_k Q_k p_k(x_i)

    everything stays in logs and comes back through a single logsumexp, so this
    never overflows even when the individual ln p are ~ -3000. columns sum to 1
    by construction.

    component_fractions : (n_components - 1,) -- the last weight is implicit
    means, sigmas       : (n_components, K)

    returns : (n_components, N) array in [0, 1]

    NOTE this does NOT care whether the components are sorted -- it just labels
    them in whatever order you hand them over. if you want "component j" to mean
    something physical (j=0 thinnest, j=-1 cocoon), sort first.
    """
    fracs  = np.atleast_1d(np.asarray(component_fractions, dtype=float))
    means  = np.asarray(means,  dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)

    n_components = len(fracs) + 1
    if (len(means) != n_components) or (len(sigmas) != n_components):
        raise ValueError("len(means)=%i, len(sigmas)=%i, ncomponents=%i"
                         % (len(means), len(sigmas), n_components))

    Qs = np.append(fracs, 1 - np.sum(fracs))
    if np.any(Qs <= 0):
        raise ValueError("all component fractions must be > 0 (and sum to < 1); "
                         "got %s" % Qs)

    # (n_components, N) -- ln Q_j + ln p_j(x_i)
    ln_w = np.array([np.log(Qs[j]) + component_likelihood(x_data, means[j], sigmas[j])
                     for j in range(n_components)])

    return np.exp(ln_w - logsumexp(ln_w, axis=0))


def component_membership_probability(x_data, component_fractions, means, sigmas,
                                     component=0, sort_dim=None):
    """
    posterior probability that each star belongs to the component(s) you ask for.
    this is the indexable version of membership_probability.

    component : int, or a sequence/slice of ints. negative indices work the usual
                numpy way, so component=-1 is the widest (cocoon) component once
                things are sorted. a sequence sums the responsibilities, e.g.
                component=[0, 1] on a 3-component fit gives the whole thin
                stream (the two narrow components together), and
                component=slice(0, -1) gives "everything but the cocoon"
                regardless of n_components.

    sort_dim : None (default) to take the components exactly as given, or a
               phase space dimension index to CHECK that they are already sorted
               narrowest -> widest in that dimension before indexing. pass
               sort_dim=-1 when the index is meant to carry the physical
               narrow -> cocoon meaning, so a mislabelled fit raises instead of
               quietly reporting the wrong population.

    returns : (N,) array in [0, 1]
    """
    sigmas = np.asarray(sigmas, dtype=float)

    if sort_dim is not None:
        widths = sigmas[:, sort_dim]
        if np.any(np.diff(widths) < 0):
            raise ValueError("components are not sorted narrowest -> widest in "
                             "sort_dim=%i (widths %s) -- call sort_components "
                             "first" % (sort_dim, widths))

    R = component_responsibilities(x_data, component_fractions, means, sigmas)

    if isinstance(component, (int, np.integer)):
        return R[component]

    # a sequence or a slice: sum the selected components' responsibilities. the
    # atleast_2d is for the degenerate one-element case, where R[[j]] is already
    # 2-D but R[j:j+1] on a scalar-ish index would not be.
    return np.atleast_2d(R[component]).sum(axis=0)


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

    a thin wrapper on component_membership_probability with
    component=slice(0, -1) -- use that one directly if you want a single
    component rather than the thin/cocoon split.

    ORDER MATTERS: "the last component is the cocoon" is only true if the
    components are sorted, so this checks and refuses rather than quietly
    reporting the wrong population. run sort_components first.

    returns : (N,) array in [0, 1]
    """
    fracs = np.atleast_1d(np.asarray(component_fractions, dtype=float))
    if len(fracs) + 1 < 2:
        raise ValueError("need at least 2 components to separate a cocoon")

    return component_membership_probability(x_data, fracs, means, sigmas,
                                            component=slice(0, -1),
                                            sort_dim=sort_dim)


### for model choice...
def AIC(lnL0, k, N):
    t1 = -2 * lnL0
    t2 = 2*k
    t3 = (2*k*(k+1))/(N-k-1)

    return t1+t2+t3

def BIC(lnL0, k, N):
    t1 = -2*lnL0
    t2 = k*np.log(N)
    return t1+t2
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

ii = 2 # <--- pick orbit. 
orbit = orbits[ii]
## do the orbit-wise check -- integrate prog orbit and find the pericenter. 
init_displacement = init_displacements[ii]
orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
peri = orbit_obj.pericenter()

mass_index = 1 # <--- pick stellar population
rvir_index=-1 # <--- pick density
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

use = unbound & ol_clip
x_data = np.column_stack([sc_straighter[k][use] for k in keys])

sd = x_data.std(axis=0)
mu_1, sigma_1 = np.zeros(4), 0.1 * sd  # thin: narrower than the data
mu_2, sigma_2 = np.zeros(4), 10.0 * sd  
# mu_3, sigma_3 = np.zeros(4), 10 * sd   # cocoon: broader than the data
f_1 = 0.9                          # starting at even groups would "let the data decide." but for two components only, guessing 90% thin stream is sort of like a prior. 
# f_2 = 1/3 # - 0.01

fracs_0 = np.array([f_1])#, f_2])
means_0 = np.array([mu_1, mu_2])#, mu_3])
sigmas_0 = np.array([sigma_1, sigma_2])#, sigma_3])
theta0 = pack_params(fracs_0, means_0, sigmas_0)

ncomponents = len(fracs_0)+1

result = minimize(nll_flat, x0=theta0, args=(x_data,), method='Powell', # method='Nelder-Mead',
                options={'maxiter': 100000, 'maxfev': 100000})#,
                        # 'fatol': 1e-6, 'xatol': 1e-6}) #<-- those are for Nelder-Mead optimizer.
# %%
fracs_fit, means_fit, sigmas_fit = sort_components(*unpack_params(result.x, K=x_data.shape[1]))
p1, p2 = [component_membership_probability(x_data, fracs_fit, means_fit, sigmas_fit, component=ii) for ii in range(ncomponents)]

k = len(theta0) #<-- model complexity -- number of parameters here. 
N = len(x_data) #<-- number of data points

L0 = -gmm_negative_loglikelihood(fracs_fit, means_fit, sigmas_fit, 
                                        data=x_data)

aic = AIC(L0, k, N)
bic = BIC(L0, k, N)

print("AIC/BIC", aic, bic)

if len(fracs_fit)<ncomponents:
    fracs_fit = np.append(fracs_fit, 1-np.sum(fracs_fit))

cocoon_fraction = fracs_fit[-1]

print(result.success, result.message, '\nnll =', result.fun)
print("cocoon fraction:", cocoon_fraction)
print(fracs_fit)
print(sigmas_fit[:,-1])
# %%
#### NEXT: single-component model. 
fracs_1comp = np.array([]) #<-- empty by construction: with one component there are NO free weights, Q_1 = 1. pass min_components=1 so unpack_params allows it.
means_1comp = np.array([np.zeros(4)])
sigmas_1comp = np.array([np.ones(4) * sd])
theta0_1comp = pack_params(fracs_1comp, means_1comp, sigmas_1comp) #<-- length 2K = 8, no weights in it at all


result_1comp = minimize(nll_flat, x0=theta0_1comp, args=(x_data, 1), method='Powell', # <-- the 1 is min_components
                        options = {'maxiter':100000, "maxfev":100000})
fracs_fit_1comp, means_fit_1comp, sigmas_fit_1comp = \
    sort_components(*unpack_params(result_1comp.x, K=x_data.shape[1], min_components=1))
k_1comp = len(theta0_1comp)
L0_1comp = -gmm_negative_loglikelihood(fracs_fit_1comp, means_fit_1comp, sigmas_fit_1comp,
                                 data=x_data)
aic_1comp, bic_1comp = AIC(L0_1comp, k_1comp, N), BIC(L0_1comp, k_1comp, N)

print("AIC/BIC", aic_1comp, bic_1comp)

### ugH so officially two components is a better fit :((((


# %%
for k in range(4):
    fig, axs = plt.subplots(1,2, figsize=[10,3], width_ratios=[3,1], sharey=True)
    x = sc_straighter['phi1'][unbound & ol_clip]
    y = x_data[:,k]



    # ts = p_thin>0.5
    # ts = p3<(1-np.sum(fracs_fit))
    # if fracs_fit[1]+fracs_fit[2]<0.5: #<-- in this case, components 2 and 3 count as cocoon. 
    #     print("thin stream is only component 1")
    #     ts = p1>0.5
    #     p_thin = p1
    # if fracs_fit[1]+fracs_fit[2]>=0.5: #<-- in this case, only component 3 counts as cocoon
    #     print("thin stream is components one and two")
    #     ts = p3<0.5 #< ie both p1 and p2 count to thin stream. 
    #     p_thin = p1+p2

    order = np.argsort(1-p_thin)

    axs[0].scatter(x[order], y[order], c=p_thin[order], s=5, cmap='winter', rasterized=True)


    cn = ~ts
    bins = np.linspace(min(y), max(y), 50)

    axs[1].hist(y[ts], bins=bins, histtype='step', color='magenta', orientation='horizontal', lw=3, density=True)
    axs[1].hist(y[cn], bins=bins, histtype='step', color='cyan', orientation='horizontal', lw=3, density=True)


# %%













# %%
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
copy_options = [0,1,2,3,4]
# copy_options = [4,3,2,1,0]

keys = ['phi2','pm_phi1','pm_phi2','v_gsr']

dicts = []
sf_coords_obs = [] # before straightening
straight_sf_coords_obs = [] 

f_cocoons = []
vgsr_dispersions = []
phi2_dispersions = []


pericenters_kpc = []
apocenters_kpc = []
present_rs = []

for ii, orbit in enumerate(tqdm(orbits)):

    ## do the orbit-wise check -- integrate prog orbit and find the pericenter. 
    init_displacement = init_displacements[ii]
    orbit_obj = paf.integrate_prog_orbit(init_displacement, steps=100000, dt=1*u.Myr)
    peri = orbit_obj.pericenter()
    apo = orbit_obj.apocenter()
    pericenters_kpc.append(peri.to(u.kpc).value)
    apocenters_kpc.append(apo.to(u.kpc).value)

    x,y,z = init_displacement[:3]
    r = np.sqrt(x**2 + y**2 + z**2)
    present_rs.append(r) # kpc

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
        x_data = np.column_stack([sc_straighter[k][use] for k in keys])

        sd = x_data.std(axis=0)
        mu_1, sigma_1 = np.zeros(4), 0.1 * sd  # thin: narrower than the data
        mu_2, sigma_2 = np.zeros(4), 10.0 * sd  
        # mu_3, sigma_3 = np.zeros(4), 10 * sd   # cocoon: broader than the data
        f_1 = 0.9                          # starting at even groups would "let the data decide." but for two components only, guessing 90% thin stream is sort of like a prior. 
        # f_2 = 1/3 # - 0.01

        fracs_0 = np.array([f_1])#, f_2])
        means_0 = np.array([mu_1, mu_2])#, mu_3])
        sigmas_0 = np.array([sigma_1, sigma_2])#, sigma_3])
        theta0 = pack_params(fracs_0, means_0, sigmas_0)

        ncomponents = len(fracs_0)+1

        result = minimize(nll_flat, x0=theta0, args=(x_data,), method='Powell', # method='Nelder-Mead',
                        options={'maxiter': 100000, 'maxfev': 100000})#,
                                # 'fatol': 1e-6, 'xatol': 1e-6}) #<-- those are for Nelder-Mead optimizer.

        fracs_fit, means_fit, sigmas_fit = sort_components(*unpack_params(result.x, K=x_data.shape[1]))
        p1, p2 = [component_membership_probability(x_data, fracs_fit, means_fit, sigmas_fit, component=ii) for ii in range(ncomponents)]

        if len(fracs_fit)<ncomponents:
            fracs_fit = np.append(fracs_fit, 1-np.sum(fracs_fit))



        ##### stuff i was trying out letting either the cocoon or the thin stream have two components:
        # if fracs_fit[1]+fracs_fit[2]<0.5: #<-- in this case, components 2 and 3 count as cocoon. 
        #     print("thin stream is only component 1")
        #     ts = p1>0.5
        #     p_thin = p1
        #     f_cocoon = fracs_fit[1] + fracs_fit[2]
        # if fracs_fit[1]+fracs_fit[2]>=0.5: #<-- in this case, only component 3 counts as cocoon
        #     print("thin stream is components one and two")
        #     ts = p3<0.5 #< ie both p1 and p2 count to thin stream. 
        #     p_thin = p1+p2
        #     f_cocoon = fracs_fit[2]


        p_thin = p1
        ts = p1>0.5
        p_cocoon = 1-p_thin
        f_cocoon = fracs_fit[-1]

        print(result.success, result.message, '\nnll =', result.fun) #<-- verbose? 

        f_cocoons_this_orbit.append(f_cocoon)


        ######## THIS IS NOT RIGHT IF >2 components fit, THE COCOON COULD BE TWO OF THE THREE COMPONENTS IN SOME CASES. 
        sigphi2, sigpmphi1, sigpmphi2, sigvgsr = sigmas_fit[-1] #<-- cocoon component. thin stream is components 0 and 1
        vgsr_dispersions_this_orbit.append(sigvgsr)
        phi2_dispersions_this_orbit.append(sigphi2)

        ### TODO: add a step that plots everything and saves the folder so that I can visually inspect -- see inspect_new_sims.py for a nice plotting routine. 
        order = np.argsort(p_cocoon)
        fig, axs = plt.subplots(len(keys), 2, figsize=[10, 10], width_ratios = [4,1])

        plt.subplots_adjust(hspace=0.03, wspace=0.03)

        # fig.suptitle(orbits[ii])
        
        key_labels = [
            r'$\phi_2~[\degree]$',
            r'$\mu_{\phi_1}~[\rm mas~yr^{-1}]$',
            r'$\mu_{\phi_2}~[\rm mas~yr^{-1}]$',
            r'$v_{\rm GSR}~[\rm km~s^{-1}]$'
        ]
        for jj, key in enumerate(keys):
            # ax_row = axs[jj]
            cut = sigmas_fit[-1][jj]

            ax = axs[jj,0]

            ax.scatter(sc_straighter['phi1'][use][order],  # plot cocoon on top. 
                    sc_straighter[key][use][order], # plot cocoon on top. 
                    # x_data[:,ii],
                        c=p_thin[order], s=5, cmap='winter',
                        rasterized=True) 
 

            # ax.set_ylim(-3*cut, 3*cut)
            ax.set_ylabel(key_labels[jj], fontsize=15)


            ax = axs[jj,1]
            bins = np.linspace(-3*cut, 3*cut, 50)

            # tsd, _ = np.histogram(sc_straighter[key][ol_clip & unbound & ~cocoon_selection],
            #                       bins=bins, density=True)

            cocoon_selection = p_thin<0.5
            ax.hist(sc_straighter[key][ol_clip & unbound][~cocoon_selection], 
                    alpha=0.2, density=True, color='k',orientation='horizontal',
                    bins=bins)
            ax.hist(sc_straighter[key][ol_clip & unbound][cocoon_selection],
                    histtype='step', density=True, lw=2, 
                    color=cc[-1],orientation='horizontal',
                    bins=bins)

   

            # too annoying to get the limits to work out. being unrigorous for now...
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_yticklabels([])

            if jj<3:
                # print("REMOVING TICK LABLES>>>>>")
                axs[jj,0].set_xticklabels([])
                axs[jj,1].set_xticklabels([])


        axs[-1,0].set_xlabel(r'$\phi_1~[\degree]$')
        axs[-1,1].set_xlabel(r'density')

        plt.savefig("/n/home02/amphillips/p27_nbody/plots/cocoon_separation/gmm/%s_%.2f.pdf"%(orbit, rvirs[rvir_index]),
                    bbox_inches='tight')
        plt.close()

    f_cocoons.append(f_cocoons_this_orbit)
    vgsr_dispersions.append(vgsr_dispersions_this_orbit)
    phi2_dispersions.append(phi2_dispersions_this_orbit)


# %%
pericenters_kpc = np.array(pericenters_kpc)
# apocenters_kpc = np.array(apocenters_kpc)
# present_rs = np.array(present_rs)

f_cocoons = np.array(f_cocoons)
vgsr_dispersions = np.array(vgsr_dispersions)
phi2_dispersions = np.array(phi2_dispersions)

### roughly ~amount of the way through orbit
orbital_phases = (present_rs - pericenters_kpc) / (apocenters_kpc - pericenters_kpc)
orbits = np.array(orbits)

eccentricities = (apocenters_kpc - pericenters_kpc) / (apocenters_kpc + pericenters_kpc)

# %%
# reordered = np.argsort(pericenters_kpc)
reordered = np.argsort(orbital_phases)



ccc = cc[1:]
fig, axs = plt.subplots(1,3,figsize=[21,7], sharex=True)
for ii, orbit in enumerate(tqdm(orbits[reordered])):
    f_cocoons_this_orbit = f_cocoons[reordered][ii]
    phi2_dispersions_this_orbit = phi2_dispersions[reordered][ii]
    vgsr_dispersions_this_orbit = vgsr_dispersions[reordered][ii]

    x = rvirs
    axs[0].plot(x, f_cocoons_this_orbit, 
                # label=orbit+r"; $r_{\rm peri}=%.1f~\rm kpc$"%pericenters_kpc[reordered][ii],
                label = orbit+r'"orbital phase"=%.2f'%orbital_phases[reordered][ii],
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


# plt.savefig("plots/cocoon_separation/gmm/two_component_summary_orbPhase.pdf", dpi=300, bbox_inches='tight')
# %%
# idk man, that looks rly bad. let's get some MCMC going maybe... ... ... ... ... 

