#--------------------------------------------------------------#
#  clean functions to do the gaussian mixture modeling.        #
#   this script saves cocoon info + data_dicts at the end
#                                                              #
#--------------------------------------------------------------#
# %%
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
### helper functions here:
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

def _bound_pairs(b, shape, name):
    """
    normalize a group of bounds into a float array of shape + (2,).

    accepts None (the whole group unbounded), a single (lo, hi) pair (broadcast
    over the group), or anything reshapeable to (prod(shape), 2). None as an
    individual limit means unbounded on that side and becomes -+inf here --
    inf is easier to transform than None, and it gets turned back into None at
    the end of pack_bounds.
    """
    n = int(np.prod(shape))
    if b is None:
        flat = [(None, None)] * n
    else:
        flat = list(np.reshape(np.asarray(b, dtype=object), (-1, 2)))
        if len(flat) == 1 and n > 1:
            flat = flat * n                      # one pair, applied to the group
        if len(flat) != n:
            raise ValueError("%s: got %i (lo, hi) pairs, expected %i for shape %s"
                             % (name, len(flat), n, shape))

    out = np.empty((n, 2))
    for i, (lo, hi) in enumerate(flat):
        out[i, 0] = -np.inf if lo is None else float(lo)
        out[i, 1] = np.inf if hi is None else float(hi)
    bad = out[:, 0] > out[:, 1]
    if np.any(bad):
        raise ValueError("%s: lower bound above upper bound at index %s"
                         % (name, np.flatnonzero(bad)))
    return out.reshape(shape + (2,))

def pack_bounds(fracs_bounds=None, means_bounds=None, sigmas_bounds=None,
                n_components=None, K=4):
    """
    NATURAL-space bounds -> the ((n-1) + 2nK,) list of (lo, hi) pairs that
    minimize() wants, in the same THETA layout as pack_params.

    the whole point: theta is not the natural parameters, it is
    [alpha_1 ... alpha_{n-1}, mu_1, ln sigma_1, mu_2, ln sigma_2, ...], so a
    bound has to be pushed through the same monotonic reparameterization that
    pack_params applies to the value. both transforms are increasing, so
    (lo, hi) maps to (T(lo), T(hi)) with no reordering:
      f     -> alpha = logit(f)      [n_components == 2 only, see below]
      sigma -> ln sigma
    in particular sigma > 0 becomes ln sigma > -inf, i.e. UNBOUNDED -- positivity
    is already enforced by the parameterization, so those bounds are free and
    passing them costs nothing.

    fracs_bounds  : (n-1, 2) in Q space, or a single pair, or None.
    means_bounds  : (n, K, 2), or a single pair, or None.
    sigmas_bounds : (n, K, 2) in sigma (NOT ln sigma) space, or a pair, or None.
    n_components, K : only needed if they can't be inferred from the arguments.

    NOTE the n_components > 2 restriction on fracs_bounds. alpha_j = ln(Q_j/Q_last)
    depends on ALL the weights, so a box in Q space is not a box in alpha space
    and there is no honest elementwise bound to hand minimize(). for n=2 the map
    is just logit and it is exact -- which covers bounding the cocoon fraction,
    since f_cocoon = 1 - f_1 there. for n>2 use cocoon_fraction_constraint()
    with a constraint-aware method instead.
    """
    # ---- infer the grid shape from whatever was actually handed over --------
    if n_components is None:
        if fracs_bounds is not None and np.ndim(fracs_bounds) > 1:
            n_components = len(np.reshape(np.asarray(fracs_bounds, dtype=object),
                                          (-1, 2))) + 1
        elif means_bounds is not None and np.ndim(means_bounds) == 3:
            n_components = len(means_bounds)
        elif sigmas_bounds is not None and np.ndim(sigmas_bounds) == 3:
            n_components = len(sigmas_bounds)
        else:
            raise ValueError("can't infer n_components from these bounds; "
                             "pass n_components explicitly")
    if means_bounds is not None and np.ndim(means_bounds) == 3:
        K = np.shape(means_bounds)[1]
    elif sigmas_bounds is not None and np.ndim(sigmas_bounds) == 3:
        K = np.shape(sigmas_bounds)[1]
    n_free = n_components - 1

    f = _bound_pairs(fracs_bounds,  (n_free,),        'fracs_bounds')
    m = _bound_pairs(means_bounds,  (n_components, K), 'means_bounds')
    s = _bound_pairs(sigmas_bounds, (n_components, K), 'sigmas_bounds')

    # ---- fractions -> alpha -------------------------------------------------
    if np.isfinite(f).any():
        if n_components != 2:
            raise ValueError("finite fracs_bounds are only well defined for "
                             "n_components=2 (alpha_j couples all the weights); "
                             "got n_components=%i. use cocoon_fraction_constraint()"
                             % n_components)
        if np.any(f < 0) or np.any(f > 1):
            raise ValueError("fracs_bounds must lie in [0, 1]; got %s" % (f,))
        f_alpha = logit(f)          # logit(0) = -inf, logit(1) = +inf, monotonic
    else:
        f_alpha = f                 # already all +-inf

    # ---- sigmas -> ln sigma -------------------------------------------------
    if np.any(s < 0):
        raise ValueError("sigmas_bounds must be >= 0 (they are sigma, not ln sigma)")
    with np.errstate(divide='ignore'):
        s_ln = np.log(s)            # log(0) -> -inf, log(inf) -> inf

    # ---- assemble in pack_params' layout ------------------------------------
    # pack_params does hstack([means, log sigmas]) then ravel, so the mu/ln-sigma
    # pairs stay blocked BY COMPONENT. same thing here, with a trailing (2,).
    blocks = np.concatenate([m, s_ln], axis=1)                 # (n, 2K, 2)
    packed = np.concatenate([f_alpha, blocks.reshape(-1, 2)])  # ((n-1)+2nK, 2)

    expected = n_free + 2 * n_components * K
    if len(packed) != expected:
        raise ValueError("packed %i bounds, expected %i" % (len(packed), expected))

    # back to None for the unbounded sides: scipy takes either, but None is the
    # documented spelling and it reads better when you print the thing.
    return [(None if np.isneginf(lo) else lo,
             None if np.isposinf(hi) else hi) for lo, hi in packed]

def sigma_ratio_constraint(min_ratio, n_components=2, K=4, dims=None,
                           thin_component=0):
    """
    require the cocoon to be at least min_ratio times WIDER than the thin
    component, dimension by dimension -- with a DIFFERENT ratio per dimension
    if you want one.

    this is usually a better handle than bounding the cocoon fraction. bounding
    f_cocoon just clips the answer at whatever wall you put up -- the optimizer
    reports the bound, not a fit. what actually goes wrong for a puffy
    progenitor is that the two components stop being distinguishable, and the
    second one gets used to soak up the non-gaussian wings of the thin stream
    rather than a physically separate cocoon. constraining the SEPARATION says
    what you actually mean: "only call it a cocoon if it's much wider."

    it is exactly LINEAR in theta -- ln sigma is stored there directly, so
    sigma_cocoon / sigma_thin >= R is
        (ln sigma_cocoon)_k - (ln sigma_thin)_k >= ln R,
    one row per phase space dimension. no reparameterization needed. that row
    structure is also why a per-dimension ratio costs nothing: LinearConstraint
    takes a VECTOR lower bound, one entry per row, so ln R just stops being a
    scalar that broadcasts.

    dims : which phase space dimensions to apply it to (default all K). the
           keys order is ['phi2','pm_phi1','pm_phi2','v_gsr'], so dims=[0, 3]
           constrains phi2 and v_gsr only.

    min_ratio : either ONE ratio for every requested dimension, or a list with
           the same length as dims, ALIGNED WITH IT elementwise -- so
               dims=[0, 3], min_ratio=[10.0, 4.0]
           demands a cocoon 10x wider in phi2 and 4x wider in v_gsr. this is the
           physical case: from Jarvis+26 the GD-1 cocoon is ~10x wider in phi2
           but only ~3x in radial velocity, so a single ratio is either too weak
           in phi2 or too strong in v_gsr. the lengths are checked, since
           silently broadcasting a mismatched list would constrain the wrong
           dimension.

    needs 'SLSQP', 'trust-constr' or 'COBYLA' -- Powell ignores constraints.
    """
    from scipy.optimize import LinearConstraint

    dims = list(range(K)) if dims is None else list(dims)
    n_free = n_components - 1
    n_theta = n_free + 2 * n_components * K

    # one ratio per row of the constraint matrix. a scalar is expanded here
    # rather than left to broadcast, so len(ratios) == len(rows) always holds
    # and a wrong-length list is an error instead of a partial constraint.
    ratios = np.asarray(min_ratio, dtype=float)
    if ratios.ndim == 0:
        ratios = np.full(len(dims), float(ratios))
    elif ratios.shape != (len(dims),):
        raise ValueError("min_ratio has shape %s but %i dims were requested; "
                         "pass a single ratio or one per dim, in the same order "
                         "as dims" % (ratios.shape, len(dims)))
    if np.any(ratios <= 0):
        raise ValueError("min_ratio must be > 0 (it is a width ratio, and ln R "
                         "is the bound); got %s" % (ratios,))

    def lnsigma_index(j, k):
        # component j's block starts at n_free + j*2K; mus first, then ln sigmas
        return n_free + j * 2 * K + K + k

    rows = []
    for k in dims:
        row = np.zeros(n_theta)
        row[lnsigma_index(n_components - 1, k)] = 1.0    # cocoon = last
        row[lnsigma_index(thin_component, k)] = -1.0
        rows.append(row)

    return LinearConstraint(np.array(rows), np.log(ratios), np.inf)


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

# %%

# if __name__=="__main__":
make_plots=False
constrain_widths=False


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

def trim_obstream_percentile(sc, p=[1,99], 
                            trim_keys=['phi1','phi2','pm_phi1','pm_phi2','v_gsr']):
    criteria = []
    for key in trim_keys:
        key_low, key_high = np.percentile(sc[key], q=p)
        key_crit = (sc[key]<=key_high) & (sc[key]>=key_low)
        criteria.append(key_crit)

    trim_criteria = np.logical_and.reduce(criteria)
    return trim_criteria

for ii, orbit in enumerate(tqdm(orbits)): #<--- this i can do later i think. 


    mass_index = 1 # <-- LOW mass stellar population... should minimize cocoon contributions from stellar evolution-related kicks i think. 
    for rvir_index in range(4):

        (core, data_dict, CMdict, lumdict, inMW, trim), path, apo, age, init_displacement, copy = \
            simspect.prepare_nbody_data_anycopy(
                orbit, stellar_pop=masses[mass_index], rvir_index=rvir_index, copies=copy_options,
                include_photometry=False
            )

        # dicts.append(data_dict)

        # coords_obs, sf = simspect.streamframe_coords_observed(orbit, CMdict, prog_tab) #<-- i think i straight up never actually need these. 
        # sf_coords_obs.append(sf_coords_obs)

        # straightened coords
        sc = simspect.straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab) #<-- sc is returned as a DICTIONARY! 

        unbound = ~CMdict['in_rtid']
        # unbound = unbound[inMW][trim] # don't care about this. 


        ### NEW scheme for trimming the stream just dropped, no 'inMW' necessary now. 
        inMW_na = np.ones(len(sc['phi1']), dtype=bool) #<-- i don't actually want to do a "inMW" trim here. 
        trim_new = trim_obstream_percentile(sc)



        #### apply the trim to all of the coordinates:
        trimmed_sc = simspect.clip_coords(sc, [inMW_na, trim_new]) #<-- this applies inMW, trim to the coordinate dictionary
        sc_straighter = simspect.poly_straightening(trimmed_sc) #< subtract a polynomial on top of the orbit subtraction

        #### i am less worried about ol_clip now. deleting it since it's ~taken care of by trimming the 1-99th percentile in orbit-subtracted coordinates. 
        # ol_clip = simspect.outlier_clip( #<-- avoid biasing the cocoon dispersion with a few crazy outliers. 
        #     sc_straighter['v_gsr'], sc_straighter['pm_phi1'], sc_straighter['pm_phi2'] 
        # )

        # use = ol_clip & unbound
        unbound = unbound[trim_new]
        use = unbound

        ### confirming that esp for m3 this straightening is quite a bit better (subjectively)
        # fig, ax = plt.subplots(figsize=[9,3])
        # ax.scatter(sc_straighter['phi1'][use], 
        #            sc_straighter['phi2'][use], c='k', s=1)
        # ax.set_ylim(-1.5, 1.5)

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

        ### BOUNDS: mean box at +/- 1 sd of the DATA, per dimension. see the
        # interactive cell above for the full reasoning -- short version is that
        # a real cocoon shares the thin stream's mean (it is WIDER, not
        # displaced), so bounding mu near zero stops the second component from
        # wandering off-track and soaking up a diverging tail or the far end of
        # an epicyclic feather instead. has to be per dimension: sd is ~0.1 deg
        # in phi2 but ~10 km/s in v_gsr.
        # mu_halfwidth = 1.0 * sd

        ### for this i am not going to bound the mean, but i will save it and 
        #   use means far from zero to flag fits that might not have worked well. 
        # means_bound  = np.broadcast_to(np.column_stack([-mu_halfwidth, mu_halfwidth]),
        #                                (ncomponents, len(sd), 2))
        # fracs_bound  = None
        # sigmas_bound = (0, None)         # sigma > 0 -- already free in ln sigma

        # bounds = pack_bounds(fracs_bound, means_bound, sigmas_bound,
        #                      n_components=ncomponents, K=x_data.shape[1])
        bounds=None

        ### CONSTRAINTS: require the cocoon be at least min_ratio times WIDER
        # than the thin component, PER DIMENSION (aligned elementwise with
        # constraint_dims). keys order is ['phi2','pm_phi1','pm_phi2','v_gsr'].
        constraint_dims = [0, 2, 3]      # phi2, pm_phi2, v_gsr
        min_ratio = [10.0, 5.0, 5.0]     # same order as constraint_dims
        constraints = [sigma_ratio_constraint(min_ratio,
                                            n_components=ncomponents,
                                            K=x_data.shape[1],
                                            dims=constraint_dims)]


        # stage 1: Powell, to land in the right basin. it IGNORES constraints
        # (only warns) but DOES honour bounds, so pass them -- otherwise stage 1
        # can walk the mean outside the box and SLSQP silently CLIPS x0 back in,
        # throwing away the basin Powell was run to find.
        result_free = minimize(nll_flat, x0=theta0, args=(x_data,), method='Powell',
                            bounds=bounds,
                            options={'maxiter': 100000, 'maxfev': 100000})

        if constrain_widths==True:
            # stage 2: re-fit from there, WITH the constraint. SLSQP from a cold
            # start collapses this likelihood onto a single component, hence two stages.
            result = minimize(nll_flat, x0=result_free.x,
                            args=(x_data,),
                            method='SLSQP',
                            bounds=bounds,          # bounds and constraints coexist
                            constraints=constraints,
                            options={'maxiter': 5000})
        else:
            result=result_free

        # check result.fun, never result.success. the constrained nll is
        # necessarily >= the free one; the gap is how hard the data resist.
        # _demand = ", ".join("%gx %s" % (r, keys[d])
        #                     for d, r in zip(constraint_dims, np.atleast_1d(min_ratio)))
        # print("%s rvir=%.2f: free nll = %.2f | cnstr nll = %.2f "
        #     "(cost of demanding a wider cocoon [%s]: %.2f)"
        #     % (orbit, rvirs[rvir_index], result_free.fun, result.fun,
        #         _demand, result.fun - result_free.fun))

        # an ACTIVE mean bound is silent in scipy: a mu sitting exactly on the
        # wall is the optimizer reporting the bound, not a fitted mean. it means
        # the data wanted an OFFSET component, i.e. the straightening left
        # structure behind (see the m3/pa5 TODO) rather than a cocoon.
        # _fr, _mu, _sg = unpack_params(result.x, K=x_data.shape[1])
        # _on_wall = np.abs(np.abs(_mu) - mu_halfwidth) < 1e-6 * np.maximum(mu_halfwidth, 1)
        # if _on_wall.any():
        #     for j, d in zip(*np.nonzero(_on_wall)):
        #         print("  WARNING mu[comp %i, %s] = %+.4g is ON its +/-%.4g bound"
        #               % (j, keys[d], _mu[j, d], mu_halfwidth[d]))
        # else:
        #     print("  means all interior to the +/-1 sd box (max |mu|/sd = %.2f)"
        #           % np.max(np.abs(_mu) / mu_halfwidth))

        fracs_fit, means_fit, sigmas_fit = sort_components(*unpack_params(result.x, K=x_data.shape[1]))
        # NB: not `for ii in ...` -- that shadows the orbit-loop index.
        p1, p2 = [component_membership_probability(x_data, fracs_fit, means_fit, sigmas_fit, component=cc_i) for cc_i in range(ncomponents)]

        if len(fracs_fit)<ncomponents:
            fracs_fit = np.append(fracs_fit, 1-np.sum(fracs_fit))


        p_thin = p1
        ts = p1>0.5
        p_cocoon = 1-p_thin
        f_cocoon = fracs_fit[-1]


        ### modify the data dictionary with model information
        cocoon_info = {}
        cocoon_info['mu_thin'] = means_fit[0]
        cocoon_info['sigma_thin'] = sigmas_fit[0]
        cocoon_info['mu_cocoon'] = means_fit[-1]
        cocoon_info['sigma_cocoon'] = sigmas_fit[-1]
        cocoon_info['f_cocoon'] = f_cocoon

        ### things that are like per star
        cocoon_info['p_thin'] = p_thin
        cocoon_info['sc_straighter'] = sc_straighter #<-- already has [inMW][trim] applied, needs [ol_clip & unbound] applied.
        # cocoon_info['ol_clip'] = ol_clip # <-- with unbound is 'use'
        cocoon_info['unbound'] = unbound # <-- with ol_clip is 'use'
        cocoon_info['trim_new'] = trim_new #<-- but unbound is already trimmed to the trimmed length... lol

        data_dict['cocoon_info']= cocoon_info


        ### pickle the dictionary. 
        if constrain_widths==True:
            datapath = '/n/home02/amphillips/p27_nbody/data/data_dicts/constrained/'
        if constrain_widths==False:
            datapath = '/n/home02/amphillips/p27_nbody/data/data_dicts/unconstrained/'
        rvir = rvirs[rvir_index]
        print("dumping to pkl file...")
        with open(datapath+'%s_%.2f.pickle'%(orbit, rvir), 'wb') as handle:
            pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        if make_plots==True:
            order = np.argsort(p_cocoon)
            fig, axs = plt.subplots(len(keys), 2, figsize=[10, 10], width_ratios = [4,1])

            plt.subplots_adjust(hspace=0.03, wspace=0.03)

            
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
                ax.set_ylim(-cut,cut)
    

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
                ax.set_ylim(-cut,cut)

    

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


            if constrain_widths==False:
                plt.savefig("/n/home02/amphillips/p27_nbody/plots/cocoon_separation/gmm/%s_%.2f.pdf"%(orbit, rvirs[rvir_index]),
                            bbox_inches='tight')
            if constrain_widths==True:
                plt.savefig("/n/home02/amphillips/p27_nbody/plots/cocoon_separation/gmm_constrained/%s_%.2f.pdf"%(orbit, rvirs[rvir_index]),
                            bbox_inches='tight')

            plt.close()

# %%
