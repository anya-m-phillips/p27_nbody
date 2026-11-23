#--------------------------------------------------------------#
#  single-component dispersions of mock-observed GD-1 samples,  #
#  with sigma_binary from two-epoch RVs: fit or estimated.      #
#   (moved out of the last cell of gmm.py, 15 Nov)              #
#--------------------------------------------------------------#
# Assume you observe every star that passes the nominal selection of a survey (`survey` in the settings
# cell, = gmm.py's noise case: 'desi' = the N_jarvis = 679 brightest in catalog G, rverr < 10; 'via' /
# 'via10hr' = rverr < 5 with the Via 1 hr / 10 hr errors; good Gaia PMs either way; a binary is one source),
# and a random fraction f_repeat of them get a SECOND RV epoch, ~1 yr after the first. what intrinsic stream
# dispersion do you get out of a single-component fit, with sigma_binary either a free parameter
# (sb_mode = 'fit', Nibauer+26-like) or estimated from the repeats' epoch-to-epoch dv and added in
# quadrature (sb_mode = 'dv')? see the settings cell.
#    f_repeat = 0.2 is roughly another GD-1 measurement that accounts for binaries
#    (35 repeat observations of 160 unique sources; Nibauer+26)
#   **in this case I only care about fitting the single component model!!!**
#
# the mock RVs (load_mock_sample, multi_epoch=True): each binary gets a Keplerian RV curve from its
#   petar orbit (semi, ecc, masses; random inclination, periapsis argument and phase), about its
#   NOISELESS CoM v_gsr; singles sit at their CoM v_gsr. two epochs, Delta t ~ N(1 yr, 2 months);
#   each epoch gets a fresh noise draw with the star's rverr_<survey>.
#
# the likelihood (nll_multiepoch): star i's epochs v_ik are one draw from the stream (shared by its
#   epochs) plus, per epoch, binary motion and noise:
#       v_i ~ MVN(mu 1, C_i),   C_i = diag(err_ik^2 + sigma_binary^2) + sigma_stream^2 1 1^T
#   a star without a second epoch gets err = missing_epoch_err (1e4 km/s) for it. the within-star
#   scatter constrains sigma_binary, the between-star scatter sigma_stream^2 + sigma_binary^2.
#   (shared_stream=False gives C_i = diag(err_ik^2 + sigma_binary^2 + sigma_stream^2), every epoch
#   independent, in which sigma_binary and sigma_stream are exactly degenerate.)
#
# the scan: f_repeat = 0.2, 0.3, ..., 1.0. for each random draw the repeat stars are the first
#   round(f_repeat * N_obs) of ONE random permutation of the sample, so each step really is "add 10%
#   more stars to the repeat sample" (nested), and the step-to-step differences aren't drowned in
#   draw-to-draw scatter. the RV curves and noise are drawn once per rvir0 (in load_mock_sample).
#   (until 15 Nov the repeats replaced the binaries' RVs with their CoM RVs and an estimated
#    sigma_binary was added in quadrature for the rest -- table desi_noise_nominal_repeatfracscan.fits,
#    no longer written; see log.md.)
#
# outputs:
#   data/gmm_tables/<survey>_noise_nominal_multiepoch_scan_<tag>.fits  one row per (rvir, f_repeat);
#        <tag> = fit_shared / fit_indep / dv_deconvolved / dv_raw / ..., from sb_mode + shared_stream / sb_dv_method
#   plots/gd1_<survey>_nominal_multiepoch_scan_<tag>.pdf               sigma_v_gsr vs f_repeat
#   plots/jarvis_comparison_single_component[_with_errs].pdf       (moved here from results.py, 15 Nov)
#
# the cells below the experiment only read tables, so once the scan table exists they can be
# re-run on their own (run the imports/settings cell first).
# %%
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
sys.path.append(script_path)

import os
import pickle
import numpy as np
import petar
from scipy.optimize import minimize
import astropy.units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.style.use(script_path+'/vedant.mplstyle')

from tqdm import tqdm
import PETAR_ANALYSIS_FUNCTIONS as paf
from gmm import (measurement_errors, binary_excess_array,
                 pack_params, unpack_params, nll_flat, N_jarvis)

# %%
datapath = '/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
table_path = repo_path+"/data/gmm_tables/"
plot_path = repo_path+"/plots/"

rvirs = [0.75, 1.5, 3, 6] #<-- also the grid's rvir_index order (0-3)
f_repeat_nominal = 0.2 # ~ Nibauer+26 (35/160): fraction of sources with a second RV epoch
f_repeats = np.round(np.arange(f_repeat_nominal, 1.0 + 1e-9, 0.1), 2) # 0.2, 0.3, ..., 1.0

n_draws = 10        # random choices of the repeat stars, to see how much the choice matters

### which survey: RV errors + RV selection, exactly as gmm.py's noise case of the same name (load_mock_sample)
survey = 'via'         #<-- 'desi':    rverr_desi (Jarvis+26 V_ERRs by rank); the N_jarvis = 679 brightest in catalog G, rverr < 10.
                        #    'via':     rverr_via (viamock, 1 hr / 1 exposure); rverr < 5 (~1578 GD-1 stars, the whole catalog).
                        #    'via10hr': rverr_via10hr (10 hr / 10 exposures); rverr < 5. same stars as 'via' (see readme).
                        #    both: good Gaia PMs (pm_err < 0.5 mas/yr), catalog-matched, binaries case (primaries trim).
survey_label = {'desi': 'DESI', 'via': 'Via', 'via10hr': 'Via 10 hr'}[survey]

### the mock second epoch (load_mock_sample)
multi_epoch = True     #<-- True: two mock epochs on Keplerian RV curves (load_mock_sample). False: one observed epoch,
                        #    the sim's instantaneous primary v_gsr + the saved noise draw; only sb_mode='dv'
                        #    works then, with dv = observed epoch - true CoM v (estimate_sigma_binary_com).
dt_mean_yr = 5.0       # time between the two epochs ~ N(dt_mean_yr, dt_sd_yr) [yr]
dt_sd_yr = 2./12.
missing_epoch_err = 1e4 # km/s: the "error" on the epoch a non-repeat star doesn't have, so it drops out

### which sigma_binary measurement the scan makes
sb_mode = 'flag'         #<-- 'fit': Nibauer+26-like. sigma_binary is a FREE parameter of the two-epoch likelihood
                        #           (fit_multiepoch), fit jointly with sigma_stream. uses shared_stream.
                        #    'dv':  sigma_binary ESTIMATED from the repeats' epoch-to-epoch dv
                        #           (estimate_sigma_binary_dv), then the old way: one gaussian on epoch 0 with
                        #           that sigma_binary in quadrature on every star. uses sb_dv_method. 
                        #    'avg':  no sigma_binary. a repeat star's velocity is the MEAN of its two epochs (error rverr/sqrt 2),
                         #           a non-repeat star's is epoch 0; one gaussian on that. needs multi_epoch=True.
                        #    'flag': no sigma_binary. repeat stars flagged as RV variable by paf.get_detections (chi^2 of the
                         #           epochs about their mean, flagged if its p-value < 0.1, hardcoded there) are REMOVED; the
                         #           rest are fit with one gaussian (unflagged repeats: their two-epoch mean if flag_average,
                         #           else epoch 0; non-repeats: epoch 0, they can't be flagged). needs multi_epoch=True.

shared_stream = True    #<-- (sb_mode='fit') True: the stream term is shared by a star's epochs (C_i has
                        #    sigma_stream^2 off the diagonal). False: every epoch independent,
                        #    prod_k N(v_ik | mu, err_ik^2 + sigma_binary^2 + sigma_stream^2) (the Nibauer+26 form;
                        #    sigma_binary and sigma_stream are then degenerate).
sb_dv_method = 'deconvolved' #<-- (sb_mode='dv') 'deconvolved' or 'raw' (estimate_sigma_binary_dv / _com);
                             #    with multi_epoch=False also 'noiseless': std of the TRUE binary motion,
                             #    v_pri - v_com, over the repeats (no noise at all)
flag_average = True     #<-- (sb_mode='flag') unflagged repeat stars: True = two-epoch mean, False = epoch 0

if sb_mode == 'fit':
    if not multi_epoch:
        raise ValueError("sb_mode='fit' needs a second epoch to fit sigma_binary from: set multi_epoch=True")
    scan_tag = 'fit_shared' if shared_stream else 'fit_indep'
elif sb_mode == 'dv':
    if sb_dv_method == 'noiseless' and multi_epoch:
        raise ValueError("sb_dv_method='noiseless' is only set up for multi_epoch=False (dv vs the true CoM v)")
    scan_tag = ('dv_' if multi_epoch else 'dvcom_') + sb_dv_method #<-- dvcom: one epoch, dv against the true CoM v
elif sb_mode in ('avg', 'flag'):
    if not multi_epoch:
        raise ValueError("sb_mode=%r needs the second epoch: set multi_epoch=True" % (sb_mode,))
    scan_tag = 'avg' if sb_mode == 'avg' else 'flag_' + ('avg' if flag_average else 'ep0')
else:
    raise ValueError("sb_mode must be 'fit', 'dv', 'avg' or 'flag', got %r" % (sb_mode,))
# the epoch baseline goes in the file name (two-epoch modes only; with multi_epoch=False there's no second epoch,
#   so dt does nothing). e.g. dt_mean_yr = 1, dt_sd_yr = 2/12 -> '_dt1yr_sd0.17yr'; the exact values are in the
#   DTMEANYR / DTSDYR header keys.
dt_tag = ('_dt%gyr_sd%gyr' % (round(dt_mean_yr, 2), round(dt_sd_yr, 2))) if multi_epoch else ''
scan_name = survey + "_noise_nominal_multiepoch_scan_" + scan_tag + dt_tag #<-- one table per survey x measurement x baseline, so they don't overwrite each other
# one diagonal gaussian is separable, so fitting v_phi1/v_phi2 as well wouldn't change these two
exp_keys = ['phi2', 'v_gsr']
i_phi2, i_vgsr = exp_keys.index('phi2'), exp_keys.index('v_gsr')

# %%
def fit_one_component(x, xe, sb=None):
    """
    single gaussian, Powell on nll_flat (min_components=1), as in gmm.py's BIC_one fit. single epoch:
    phi2, the noiseless truth, and the "old way" dispersions with an estimated sigma_binary.
    x  : (N, K) data; xe : (N, K) 1 sd errors or None
    sb : optional fixed sigma_binary array, broadcastable to x's shape (gmm.binary_excess_array)
    returns mu, sigma (K-vectors; sigma intrinsic, i.e. deconvolved with xe and sb)
    """
    theta0 = pack_params(np.array([]), np.zeros((1, x.shape[1])), x.std(axis=0)[None, :])
    res = minimize(nll_flat, x0=theta0, args=(x, 1, xe, sb), method='Powell',
                    options={'maxiter': 100000, 'maxfev': 100000})
    _, mu, sig = unpack_params(res.x, K=x.shape[1], min_components=1)
    return mu[0], sig[0]


#---------------------------------------------------------------------#
#   multi-epoch likelihood: sigma_binary as a free parameter            #
#---------------------------------------------------------------------#
def multiepoch_nll_per_star(v, v_err, mu, sigma_stream, sigma_binary, shared_stream=True, xp=np):
    """
    -ln L for each star's epochs, (N,).
    v, v_err : (N, n_epochs) measured v_gsr and 1 sd errors [km/s]; a missing epoch is any finite v
               with v_err = missing_epoch_err (it then contributes ~nothing).
    xp : the array module, np (default) or jax.numpy (so the NUTS cell can differentiate it; needs x64,
         since missing_epoch_err^2 = 1e8).

    shared_stream=True (default):
        v_i ~ MVN(mu 1, C_i),  C_i = D_i + sigma_stream^2 1 1^T,  D_i = diag(v_err_ik^2 + sigma_binary^2)
        i.e. the star's systemic velocity is ONE draw from the stream, shared by its epochs, and the
        binary motion + noise are independent from epoch to epoch. with Sherman-Morrison
        (A = sum_k 1/d_k, B = sum_k r_k/d_k, r = v - mu):
            ln det C = sum_k ln d_k + ln(1 + sigma_stream^2 A)
            r C^-1 r = sum_k r_k^2/d_k - sigma_stream^2 B^2 / (1 + sigma_stream^2 A)
        one epoch reduces to N(v | mu, err^2 + sigma_binary^2 + sigma_stream^2).
    shared_stream=False:
        prod_k N(v_ik | mu, err_ik^2 + sigma_binary^2 + sigma_stream^2), every epoch independent.
        sigma_binary and sigma_stream only enter as a sum of squares, so they're degenerate.
    """
    v, v_err = xp.atleast_2d(v), xp.atleast_2d(v_err)
    s2 = sigma_stream**2
    d = v_err**2 + sigma_binary**2
    if not shared_stream:
        d, s2 = d + s2, 0.0
    r = v - mu
    A = xp.sum(1.0/d, axis=1)
    B = xp.sum(r/d, axis=1)
    denom = 1.0 + s2*A
    logdet = xp.sum(xp.log(d), axis=1) + xp.log(denom)
    chi2 = xp.sum(r**2/d, axis=1) - s2*B**2/denom
    return 0.5*(chi2 + logdet + v.shape[1]*np.log(2*np.pi))


def nll_multiepoch(theta, v, v_err, sigma_binary=None, shared_stream=True):
    """
    total -ln L = sum over stars of multiepoch_nll_per_star.
    theta = [mu, ln sigma_stream, ln sigma_binary]   if sigma_binary is None (fit it), or
            [mu, ln sigma_stream]                    if sigma_binary is a number (held fixed, e.g. 0).
    """
    mu, sigma_stream = theta[0], np.exp(theta[1])
    sb = np.exp(theta[2]) if sigma_binary is None else sigma_binary
    nll = np.sum(multiepoch_nll_per_star(v, v_err, mu, sigma_stream, sb, shared_stream=shared_stream))
    return nll if np.isfinite(nll) else np.inf


def fit_multiepoch(v, v_err, sigma_binary=None, shared_stream=True):
    """
    maximum likelihood of nll_multiepoch. sigma_binary=None fits it; a number holds it fixed.
    ln sigma is bounded to [-10, 6] (4.5e-5 to 400 km/s), so a width at ~0 = unresolved (see readme).
    NOT Powell: its first line search throws ln sigma_stream onto the lower bound, where the likelihood
    is flat (sigma^2 -> 0), and it stays there (on synthetic data: nll 8 worse than the truth). so
    L-BFGS-B and Nelder-Mead from the same start, keeping the lower nll (they agree on synthetic data).
    returns mu, sigma_stream, sigma_binary, result
    """
    v1 = v[:, 0]
    sd = np.std(v1)
    theta0 = [np.median(v1), np.log(sd/np.sqrt(2))]
    bounds = [(None, None), (-10., 6.)]
    if sigma_binary is None:
        theta0.append(np.log(sd/np.sqrt(2)))
        bounds.append((-10., 6.))
    fits = [minimize(nll_multiepoch, x0=np.array(theta0), args=(v, v_err, sigma_binary, shared_stream),
                     method=method, bounds=bounds, options={'maxiter': 100000})
            for method in ['L-BFGS-B', 'Nelder-Mead']]
    res = min(fits, key=lambda r: r.fun)
    sb = np.exp(res.x[2]) if sigma_binary is None else sigma_binary
    return res.x[0], np.exp(res.x[1]), sb, res


def epoch_arrays(v_epochs, rverr, repeat):
    """
    (N, 2) v and v_err for the likelihood: stars in `repeat` keep both epochs; the rest keep epoch 0
    and get missing_epoch_err on epoch 1 (its v set to epoch 0's, any finite value would do).
    """
    v = v_epochs.copy()
    v_err = np.repeat(rverr[:, None], v_epochs.shape[1], axis=1)
    v[~repeat, 1:] = v[~repeat, :1]
    v_err[~repeat, 1:] = missing_epoch_err
    return v, v_err


def estimate_sigma_binary_dv(v_epochs, rverr, repeat, method='deconvolved'):
    """
    sigma_binary from the repeat stars' epoch-to-epoch differences, dv = v_1 - v_0 (observed, so
    binaries: two noised points on the RV curve; singles: two noise draws).
    method = 'deconvolved' (default): sqrt( (var(dv) - mean(2 rverr^2)) / 2 ), clipped at 0.
                the noise (variance 2 rverr^2 per star, same error both epochs) is subtracted, and /2
                turns the variance of a difference of two independent epochs into a per-epoch one.
                the fit it goes into already deconvolves rverr, so this is the binary term alone.
                a binary that barely moves in ~1 yr (P >> 1 yr) contributes ~0 here.
             'raw': std(dv), as is. includes sqrt(2) x the noise, so it double counts the errors
                when added on top of x_err.
    """
    dv = v_epochs[repeat, 1] - v_epochs[repeat, 0]
    if method == 'raw':
        return np.std(dv)
    if method == 'deconvolved':
        return np.sqrt(max(0.0, (np.var(dv) - np.mean(2.*rverr[repeat]**2)) / 2.))
    raise ValueError("method must be 'deconvolved' or 'raw', got %r" % (method,))


def estimate_sigma_binary_com(v_obs, v_com, rverr, repeat, method='deconvolved', v_true=None):
    """
    single-epoch version (multi_epoch=False): sigma_binary from the repeat stars' dv = v_obs - v_com, the
    observed epoch (with binary motion and noise) minus the TRUE, noiseless CoM v_gsr, i.e. as if the
    repeats had pinned down each star's CoM velocity exactly. binaries: orbital motion + noise; singles: noise.
    method = 'deconvolved' (default): sqrt( var(dv) - mean(rverr^2) ), clipped at 0. one noise term
                (v_com is noiseless), so no factor 2 and no /2: the binary term alone.
             'raw': std(dv), as is (includes the noise).
             'noiseless': std(v_true - v_com) over the repeats, with v_true the epoch's velocity WITHOUT
                noise (the sim's primary v_gsr), i.e. the true binary motion (0 for singles). needs v_true.
    """
    if method == 'noiseless':
        if v_true is None:
            raise ValueError("method='noiseless' needs v_true (the observed epoch without noise)")
        return np.std(v_true[repeat] - v_com[repeat])
    dv = v_obs[repeat] - v_com[repeat]
    if method == 'raw':
        return np.std(dv)
    if method == 'deconvolved':
        return np.sqrt(max(0.0, np.var(dv) - np.mean(rverr[repeat]**2)))
    raise ValueError("method must be 'deconvolved', 'raw' or 'noiseless', got %r" % (method,))


def average_epochs(v_epochs, rverr, repeat):
    """
    (N,) velocity and error the 'avg' / 'flag' fits see: repeat stars get the mean of their epochs, with
    error rverr/sqrt(n_epochs) (same error every epoch); the rest keep epoch 0 and rverr.
    """
    n_ep = v_epochs.shape[1]
    v = np.where(repeat, v_epochs.mean(axis=1), v_epochs[:, 0])
    e = np.where(repeat, rverr/np.sqrt(n_ep), rverr)
    return v, e


def flag_rv_variable(v_epochs, rverr, repeat, v_com):
    """
    (N,) bool: repeat stars that paf.get_detections calls RV variable (chi^2 of the epochs about their mean
    with each star's rverr, flagged if its p-value < 0.1). non-repeat stars are never flagged.
    v_com only feeds get_detections' delta_vsys output, which isn't used here.
    """
    flagged = np.zeros(len(repeat), dtype=bool)
    undet, _ = paf.get_detections(rverr[repeat][:, None], v_epochs[repeat],
                                  v_com[repeat]*u.km/u.s, bool_arr='undet')
    flagged[repeat] = ~undet
    return flagged


#---------------------------------------------------------------------#
#   the mock sample                                                     #
#---------------------------------------------------------------------#
def load_present_day_binaries(dd, rvir, orbit='gd1', stellar_pop='hm', copies=range(5)):
    """
    the petar Binary object of the snapshot the pickle was built from, rows in the same order as the
    pickle's binaries (rows nsingles: of CoM/luminous). pickles built after 15 Nov 2026 carry it as
    dd['binaries']; for older ones, take the first grid copy whose present-day data.<age/10>.binary has
    the same binaries as the pickle (same number, same p1 ids in the same order). either way that's checked.
    """
    ns, nb = dd['nsingles'], dd['nbinaries']
    p1_ids = dd['IDs'][ns:ns+nb] #<-- top-level IDs are [singles, binary p1, binary p2]

    def matches(b):
        return len(b.mass) == nb and np.array_equal(b.p1.id, p1_ids)

    if 'binaries' in dd:
        if not matches(dd['binaries']):
            raise ValueError("dd['binaries'] doesn't line up with the pickle's binaries")
        return dd['binaries']

    grid = paf.extended_grid_info()
    for copy in copies:
        path, _, age, _ = grid.retrieve_sim_info(orbit, stellar_pop, rvirs.index(rvir), copy)
        binary_file = path + "data.%i.binary" % int(age/10)
        if not os.path.exists(binary_file):
            continue
        b = petar.Binary(member_particle_type=petar.Particle, G=petar.G_MSUN_PC_MYR,
                         interrupt_mode='bse', external_mode='galpy')
        b.loadtxt(binary_file)
        if matches(b):
            return b
    raise FileNotFoundError("no %s %s rvir %.2f copy has binaries matching the pickle" % (orbit, stellar_pop, rvir))


def load_mock_sample(rvir, survey=survey, multi_epoch=multi_epoch, rng=None,
                     dt_mean_yr=dt_mean_yr, dt_sd_yr=dt_sd_yr):
    """
    GD-1, `survey` noise ('desi', 'via' or 'via10hr'), binaries case (primary RVs): every system gmm.py's
    case with noise=survey would fit (include_binaries=True). returns a dict of everything the experiment
    needs, all length N_obs. always: v_com (true CoM v_gsr), v_pri (the sim's instantaneous primary v_gsr),
    v_noise (the saved noise draw for that survey), x_err, phi2, mG, is_binary. with multi_epoch=False the
    scan's one observed epoch is v_pri + v_noise.

    multi_epoch=True also builds two mock RV epochs per star (see the top of the file):
        t_obs         (N_obs, 2) [day]: 0 and Delta t ~ N(dt_mean_yr, dt_sd_yr) yr
        v_epochs_true (N_obs, 2) noiseless: singles at their CoM v_gsr; binaries on a Keplerian RV curve
                      about their CoM v_gsr (paf.get_orbital_params_from_binaries + paf.get_rvs)
        v_epochs_noise(N_obs, 2) fresh draws, N(0, rverr_<survey>) per epoch (NOT the saved draw)
        v_epochs_obs  = v_epochs_true + v_epochs_noise
        rverr         (N_obs,) rverr_<survey> [km/s]
    the random inclinations/phases, Delta t and noise all come from `rng`.
    (the selection is still the primaries one, i.e. made with the sim's instantaneous primary v_gsr.)
    """
    with open(datapath+"gd1_%.2f.pickle"%rvir, 'rb') as handle:
        dd = pickle.load(handle)

    sc_com, sc_pri = dd['sc_straighter'], dd['sc_straighter_primaries']
    distances = dd['coords_obs'].distance
    N_all = len(sc_com['phi1'])
    is_binary = np.arange(N_all) >= dd['nsingles'] #<-- CoM/luminous ordering is [singles, binaries]

    ### the binaries-case catalog products (primaries trim), as in gmm.py with include_binaries=True
    cut_p = dd['cut_for_catalog_photometry_primaries']
    matched_p = np.zeros(N_all, dtype=bool)
    matched_p[cut_p] = dd['matched_to_catalog_photometry_primaries']
    mG = dd['catalog_photometry_primaries']['mG']
    nd = dd['noise_catalog_photometry_primaries']

    ### the nominal selection, exactly as gmm.py's `use` for noise=survey:
    #   cut & matched & good_pm & good_RV & nonrem, with good_RV =
    #     'desi':           (N_jarvis brightest in catalog mG, NaNs last) & rverr_desi < 10
    #     'via', 'via10hr': rverr_<survey> < 5 (NaN rverr, i.e. NaN log_Teff, fails it)
    #   matched_p is inside cut_for_catalog_photometry_primaries (unbound & trim_new_primaries
    #   & nonrem & |v_gsr| < 30), so the cut/nonrem terms are already implied.
    if survey not in ('desi', 'via', 'via10hr'):
        raise ValueError("survey must be 'desi', 'via' or 'via10hr', got %r" % (survey,))
    rverr_all = nd['rverr_'+survey]
    good_pm = nd['pm_err_gaia']*(u.microarcsecond/u.yr).to(u.mas/u.yr) < 0.5
    if survey == 'desi':
        top_N_jarvis = np.zeros(N_all, dtype=bool)
        top_N_jarvis[np.argsort(mG)[:N_jarvis]] = True #<-- argsort puts NaNs last
        good_RV = top_N_jarvis & (rverr_all < 10.)
    else:
        good_RV = rverr_all < 5.
    use = cut_p & matched_p & good_pm & good_RV & dd['nonrem']
    idx = np.flatnonzero(use)
    idx = idx[np.argsort(mG[idx])] #<-- brightest first (order doesn't matter to the fits)
    x_err = measurement_errors(nd, survey, distances, exp_keys)[idx]
    if not np.all(np.isfinite(x_err)):
        raise ValueError("rvir %.2f: non-finite %s errors in the nominal %s sample" % (rvir, survey, survey))

    sample = {
        'mG': mG[idx],
        'is_binary': is_binary[idx],
        'x_err': x_err,
        'phi2_true': sc_com['phi2'][idx],
        'phi2_obs': sc_pri['phi2'][idx] + nd['phi2'][idx],
        'v_com': sc_com['v_gsr'][idx],    # true, no noise
        'v_pri': sc_pri['v_gsr'][idx],    # true, no noise, with binary motion (the sim's instantaneous one)
        'v_noise': nd['v_gsr_'+survey][idx], # the saved noise draw for this survey
    }
    if not multi_epoch:
        return sample

    if rng is None:
        rng = np.random.default_rng()
    ns = dd['nsingles']
    N_obs = len(idx)
    b = is_binary[idx]

    ### Keplerian parameters for every binary in the snapshot, v0 = its noiseless CoM v_gsr (residual about the track)
    binaries = load_present_day_binaries(dd, rvir)
    orb_params = paf.get_orbital_params_from_binaries(binaries, sc_com['v_gsr'][ns:], rng=rng)

    dt = rng.normal(dt_mean_yr, dt_sd_yr, size=N_obs) * 365.25 # days
    t_obs = np.column_stack([np.zeros(N_obs), dt])

    v_epochs_true = np.repeat(sample['v_com'][:, None], 2, axis=1)
    v_epochs_true[b] = paf.get_rvs(orb_params[idx[b] - ns], t_obs[b], verbose=False) #<-- binary row j of the pickle is binaries[j - nsingles]
    rverr = rverr_all[idx]
    v_epochs_noise = rng.normal(0., 1., size=v_epochs_true.shape) * rverr[:, None]

    sample.update({
        't_obs': t_obs,
        'v_epochs_true': v_epochs_true,
        'v_epochs_noise': v_epochs_noise,
        'v_epochs_obs': v_epochs_true + v_epochs_noise,
        'rverr': rverr,
        'orb_params': orb_params[idx[b] - ns], # [v0, K, w, phi0, e, P (day)] for the binaries in the sample, in sample order
    })
    return sample

# %%
run_scan = False #<-- False: skip the scan; the plot cells below just read the table of the current settings

if run_scan:
    rng_exp = np.random.default_rng(seed=1)
    scan_rows = [] #<-- one dict per (rvir, f_repeat)

    for rvir in rvirs:
        s = load_mock_sample(rvir, survey=survey, multi_epoch=multi_epoch, rng=rng_exp)
        bin_s, x_err = s['is_binary'], s['x_err']
        N_obs = len(bin_s) #<-- unique sources (a binary is one source); desi: 679 if every top-N_jarvis star passes
        all_epochs = np.ones(N_obs, dtype=bool)
        if multi_epoch:
            rverr = s['rverr']
            v_obs = s['v_epochs_obs']                          # (N, 2) mock epochs on the RV curves
            v_com_obs0 = s['v_com'] + s['v_epochs_noise'][:, 0] # CoM v + epoch 0's noise
            dv_bin = s['v_epochs_true'] - s['v_com'][:, None]   #<-- binary orbital motion per epoch (0 for singles)
        else:
            rverr = x_err[:, i_vgsr]                            # rverr_<survey>
            v_obs = (s['v_pri'] + s['v_noise'])[:, None]        # (N, 1): the one observed epoch, sim primary v + saved noise
            v_com_obs0 = s['v_com'] + s['v_noise']              # CoM v + the same noise
            dv_bin = (s['v_pri'] - s['v_com'])[:, None]         #<-- the sim's instantaneous binary motion (0 for singles)

        ### references
        #   truth: noiseless CoM v_gsr of these stars, one gaussian.
        #   none:  epoch 0 only, binary motion in, sigma_binary fixed at 0 (one epoch can't separate it).
        #   all:   the noise-limited best case, binary motion removed (CoM v + the same noise draws), sigma_binary = 0.
        #          sb_mode 'fit':         both epochs of every star, the fit's likelihood (shared_stream);
        #          sb_mode 'dv':          epoch 0 of every star, one gaussian (the data the dv-mode fit sees);
        #          sb_mode 'avg'/'flag':  the two-epoch mean of every star (error rverr/sqrt 2), one gaussian.
        #   phi2:  one gaussian on the observed phi2 (no position errors), the same for every case.
        _, sig_true = fit_one_component(np.column_stack([s['phi2_true'], s['v_com']]), None)
        _, sig_phi2 = fit_one_component(s['phi2_obs'][:, None], None)
        sig_phi2 = sig_phi2[0]
        _, ss_none, _, _ = fit_multiepoch(v_obs[:, :1], rverr[:, None], sigma_binary=0.)
        if sb_mode == 'fit':
            v_com_epochs = s['v_com'][:, None] + s['v_epochs_noise']
            _, ss_all, _, _ = fit_multiepoch(*epoch_arrays(v_com_epochs, rverr, all_epochs), sigma_binary=0.,
                                             shared_stream=shared_stream)
        elif sb_mode in ('avg', 'flag'):
            v_com_avg, e_avg = average_epochs(s['v_com'][:, None] + s['v_epochs_noise'], rverr, all_epochs)
            _, ss_all, _, _ = fit_multiepoch(v_com_avg[:, None], e_avg[:, None], sigma_binary=0.)
        else:
            _, ss_all, _, _ = fit_multiepoch(v_com_obs0[:, None], rverr[:, None], sigma_binary=0.)
        sig_none = np.array([sig_phi2, ss_none])
        sig_all = np.array([sig_phi2, ss_all])
        ### what a gaussian sigma_binary "should" be for this sample: rms binary motion over every star and
        #   epoch (singles count as 0, as in the model); and the part that changes between the two epochs,
        #   rms(dv_1 - dv_0)/sqrt(2), which is all a second epoch can see (slow binaries look like stream).
        #   (one epoch: NaN.)
        sb_rms = np.sqrt(np.mean(dv_bin**2))
        sb_rms_diff = np.sqrt(np.mean((dv_bin[:, 1] - dv_bin[:, 0])**2) / 2.) if multi_epoch else np.nan

        ### the "old way": one gaussian in (phi2, v_gsr) on EPOCH 0 for every star, with x_err, and a fixed
        #   sigma_binary in quadrature on every star. oracle = sigma_binary set to the true rms binary motion
        #   at epoch 0 over the sample, i.e. the best a single gaussian sigma_binary can do.
        x_ep0 = np.column_stack([s['phi2_obs'], v_obs[:, 0]])
        x_err_ep0 = np.column_stack([np.zeros(N_obs), rverr])
        sb_oracle = np.sqrt(np.mean(dv_bin[:, 0]**2))
        _, sig_oracle = fit_one_component(x_ep0, x_err_ep0, binary_excess_array(sb_oracle, exp_keys, N_obs))

        ### one permutation per draw, shared by every f_repeat -> nested repeat samples
        perms = np.array([rng_exp.permutation(N_obs) for _ in range(n_draws)])

        print("\nGD-1 rvir0 = %.2f pc: %i in the nominal %s sample (mG %.1f-%.1f), %i binaries (f_bin %.2f), median %s rverr %.2f km/s"
                % (rvir, N_obs, survey_label, s['mG'].min(), s['mG'].max(), bin_s.sum(), bin_s.mean(),
                   survey_label, np.median(rverr)))
        print("  sigma_phi2 = %.3f deg  (noiseless CoM: %.3f)" % (sig_phi2, sig_true[i_phi2]))
        print("  sigma_v_gsr [km/s]: truth (noiseless CoM) %.2f | %s 1 epoch, binaries in %.2f | %s, CoM (%s) %.2f"
                % (sig_true[i_vgsr], survey_label, ss_none, survey_label,
                   {'fit': '2 epochs', 'dv': 'epoch 0'}.get(sb_mode, '2-epoch mean'), ss_all))
        print("  rms binary motion per epoch %.2f km/s; epoch-to-epoch part %.2f km/s" % (sb_rms, sb_rms_diff))
        print("  old way, epoch 0, sigma_binary = true rms at epoch 0 (%.2f): sigma_v_gsr %.2f" % (sb_oracle, sig_oracle[i_vgsr]))
        print("  sb_mode = %r (%s)" % (sb_mode, scan_tag))
        print("  %8s %8s %8s | %-22s | %-22s%s" % ('f_repeat', 'N_repeat', 'n_bin', 'sigma_v_gsr (stream)',
              {'fit': 'sigma_binary (fit)', 'dv': 'sigma_binary (from dv)'}.get(sb_mode, '(no sigma_binary)'),
              {'fit': ' | %-22s' % 'sigma_stream (sb = 0)',
               'flag': ' | flagged: all / binaries / singles, N fit'}.get(sb_mode, '')))

        for f_repeat in f_repeats:
            N_repeat = int(round(f_repeat * N_obs))
            ### the experiment: N_repeat random sources have both epochs, the rest only epoch 0. then
            #   sb_mode 'fit':  fit (mu, sigma_stream, sigma_binary) with nll_multiepoch; and, for comparison,
            #                   the same with sigma_binary held at 0.
            #   sb_mode 'dv':   estimate sigma_binary from the repeats' dv, then fit one gaussian to epoch 0 of
            #                   every star with it in quadrature on every star.
            #   sb_mode 'avg':  one gaussian on the repeats' two-epoch means + the rest's epoch 0.
            #   sb_mode 'flag': drop the repeats paf.get_detections flags as RV variable, one gaussian on the rest.
            ss_draws, sb_draws, ss_nosb_draws, n_bin_repeat = [], [], [], []
            n_flag, n_flag_bin, n_fit = [], [], []
            for perm in perms:
                repeat = np.zeros(N_obs, dtype=bool)
                repeat[perm[:N_repeat]] = True

                if sb_mode == 'fit':
                    v, v_err = epoch_arrays(v_obs, rverr, repeat)
                    _, ss, sb, _ = fit_multiepoch(v, v_err, shared_stream=shared_stream)
                    _, ss_nosb, _, _ = fit_multiepoch(v, v_err, sigma_binary=0., shared_stream=shared_stream)
                    ss_nosb_draws.append(ss_nosb)
                elif sb_mode == 'dv':
                    if multi_epoch:  # dv = epoch 1 - epoch 0 of the repeats
                        sb = estimate_sigma_binary_dv(v_obs, rverr, repeat, method=sb_dv_method)
                    else:            # dv = the observed epoch - the true CoM v of the repeats
                        sb = estimate_sigma_binary_com(v_obs[:, 0], s['v_com'], rverr, repeat, method=sb_dv_method,
                                                       v_true=s['v_pri'])
                    _, sig_dv = fit_one_component(x_ep0, x_err_ep0, binary_excess_array(sb, exp_keys, N_obs))
                    ss = sig_dv[i_vgsr]
                else:  # 'avg' / 'flag': no sigma_binary
                    if sb_mode == 'avg' or flag_average:
                        v_use, e_use = average_epochs(v_obs, rverr, repeat)
                    else:
                        v_use, e_use = v_obs[:, 0], rverr
                    keep = np.ones(N_obs, dtype=bool)
                    if sb_mode == 'flag':
                        flagged = flag_rv_variable(v_obs, rverr, repeat, s['v_com'])
                        keep = ~flagged
                        n_flag.append(flagged.sum())
                        n_flag_bin.append((flagged & bin_s).sum())
                    _, sig_k = fit_one_component(np.column_stack([s['phi2_obs'][keep], v_use[keep]]),
                                                 np.column_stack([np.zeros(keep.sum()), e_use[keep]]))
                    ss, sb = sig_k[i_vgsr], np.nan
                    n_fit.append(keep.sum())

                ss_draws.append(ss)
                sb_draws.append(sb)
                n_bin_repeat.append((repeat & bin_s).sum())

            ss16, ss50, ss84 = np.percentile(ss_draws, [16, 50, 84])
            sb16, sb50, sb84 = np.percentile(sb_draws, [16, 50, 84])
            line = ("  %8.1f %8i %8.1f | %5.2f [%5.2f-%5.2f]    | %5.2f [%5.2f-%5.2f]   "
                    % (f_repeat, N_repeat, np.mean(n_bin_repeat), ss50, ss16, ss84, sb50, sb16, sb84))
            mode_cols = {}
            if sb_mode == 'fit':
                ns16, ns50, ns84 = np.percentile(ss_nosb_draws, [16, 50, 84])
                line += " | %5.2f [%5.2f-%5.2f]" % (ns50, ns16, ns84)
                mode_cols = {'S_vgsr_nosb': ns50,            # two-epoch fit, sigma_binary held at 0
                             'S_vgsr_nosb_p16': ns16,
                             'S_vgsr_nosb_p84': ns84}
            if sb_mode in ('avg', 'flag'):
                mode_cols['N_fit_mean'] = np.mean(n_fit)     # stars in the fit, on average
            if sb_mode == 'flag':
                line += " | %6.1f / %6.1f / %6.1f, %6.1f" % (np.mean(n_flag), np.mean(n_flag_bin),
                                                            np.mean(n_flag) - np.mean(n_flag_bin), np.mean(n_fit))
                mode_cols.update({
                    'n_flagged_mean': np.mean(n_flag),                               # repeats flagged RV variable, on average
                    'n_flagged_bin_mean': np.mean(n_flag_bin),                       # ... of which true binaries
                    'n_flagged_single_mean': np.mean(n_flag) - np.mean(n_flag_bin),  # ... of which singles (false positives)
                })
            print(line)

            ### sigmas are INTRINSIC (deconvolved with the RV errors), except S_true (noiseless).
            #   S_true / S_none / S_all are 2-vectors in exp_keys order: (phi2 [deg], v_gsr [km/s]).
            #   S_vgsr / sigma_binary are the sb_mode measurement (same column names in every mode;
            #   sigma_binary is NaN for 'avg' / 'flag', which have none).
            scan_rows.append({
                'sim': 'gd1', 'orbit': 'gd1', 'Rvir0': rvir,
                'N_obs': N_obs, 'f_repeat': f_repeat, 'N_repeat': N_repeat, 'n_draws': n_draws,
                'f_bin': bin_s.mean(),                       # binary fraction among the N_obs
                'n_bin_repeat_mean': np.mean(n_bin_repeat),  # binaries landing in the repeat sample, on average
                'S_true': sig_true,                          # noiseless, CoM RVs
                'S_none': sig_none,                          # survey noise, epoch 0 only, binary motion in, sigma_binary = 0
                'S_all': sig_all,                            # survey noise, CoM v, sigma_binary = 0 (see the references comment)
                'S_phi2': sig_phi2,                          # same in every sb_mode for a given survey (only v_gsr changes)
                'S_vgsr': ss50,                              # the stream sigma_v_gsr: median over draws
                'S_vgsr_p16': ss16,
                'S_vgsr_p84': ss84,
                'sigma_binary': sb50,                        # fitted ('fit') or estimated from dv ('dv'): median over draws
                'sigma_binary_p16': sb16,
                'sigma_binary_p84': sb84,
                **mode_cols,
                'sigma_binary_rms': sb_rms,                  # rms binary motion over all stars and epochs (reference)
                'sigma_binary_rms_diff': sb_rms_diff,        # rms(dv_1 - dv_0)/sqrt(2) (reference)
                'sigma_binary_oracle': sb_oracle,            # true rms binary motion at epoch 0 over the sample
                'S_vgsr_sboracle': sig_oracle[i_vgsr],       # old way with sigma_binary = sigma_binary_oracle
            })

    t_scan = Table(rows=scan_rows)
    t_scan.meta['KEYS'] = ','.join(exp_keys) #<-- order of the S_true/S_none/S_all columns
    t_scan.meta['NESTED'] = True #<-- repeat samples nested across f_repeat within a draw
    t_scan.meta['SURVEY'] = survey
    t_scan.meta['SAMPLE'] = '%s nominal (gmm.py %s, binaries case)' % (survey, survey)
    t_scan.meta['DTMEANYR'] = dt_mean_yr
    t_scan.meta['DTSDYR'] = dt_sd_yr
    t_scan.meta['MISSERR'] = missing_epoch_err
    t_scan.meta['SBMODE'] = sb_mode
    t_scan.meta['MULTIEP'] = multi_epoch #<-- False: one observed epoch (sim primary v + saved noise), dv vs the true CoM v
    if sb_mode == 'fit':
        t_scan.meta['LIKELIHD'] = ('multiepoch, stream term shared by a star\'s epochs' if shared_stream
                                   else 'multiepoch, every epoch independent (Nibauer+26 form)')
    elif sb_mode == 'dv':
        t_scan.meta['SBDVMETH'] = sb_dv_method #<-- epoch 0, sigma_binary from the repeats' dv, on every star
    elif sb_mode == 'flag':
        t_scan.meta['FLAGTEST'] = 'paf.get_detections, p < 0.1'
        t_scan.meta['FLAGAVG'] = flag_average #<-- unflagged repeats: two-epoch mean (True) or epoch 0
    t_scan.write(table_path+scan_name+".fits", overwrite=True)
    print("wrote", table_path+scan_name+".fits")

# %%
### sigma_v_gsr vs f_repeat, one square panel per rvir0, for the table of the current sb_mode settings.
#   solid + band = the sb_mode measurement (median, 16-84% over draws); 'fit' also: dashed = the same
#   fit with sigma_binary = 0. horizontal lines = truth (noiseless CoM), `survey` one epoch with binary
#   motion, `survey` at CoM v, the old way with the true rms sigma_binary.
t_scan = Table.read(table_path+scan_name+".fits", format='fits')
mode_label = {'fit_shared': '2 epochs, $\\sigma_{\\rm binary}$ free',
              'fit_indep': '2 epochs (independent), $\\sigma_{\\rm binary}$ free',
              'dv_deconvolved': '1 epoch + $\\sigma_{\\rm binary}$ from repeats',
              'dv_raw': '1 epoch + raw std($\\Delta v$) as $\\sigma_{\\rm binary}$',
              'dvcom_deconvolved': '1 epoch + $\\sigma_{\\rm binary}$ from repeats (vs CoM)',
              'dvcom_raw': '1 epoch + raw std($v - v_{\\rm CoM}$) as $\\sigma_{\\rm binary}$',
              'dvcom_noiseless': '1 epoch + noiseless std($v - v_{\\rm CoM}$) as $\\sigma_{\\rm binary}$',
              'avg': 'repeats: 2-epoch mean, no $\\sigma_{\\rm binary}$',
              'flag_avg': 'RV-variable repeats removed (rest: 2-epoch mean)',
              'flag_ep0': 'RV-variable repeats removed (rest: epoch 0)'}[scan_tag]

fig, axs = plt.subplots(2, 2, figsize=(14, 14), sharex=True)
for ax, rvir in zip(axs.flat, rvirs):
    t = t_scan[t_scan['Rvir0']==rvir]
    t = t[np.argsort(t['f_repeat'])]
    ax.fill_between(t['f_repeat'], t['S_vgsr_p16'], t['S_vgsr_p84'], color='C0', alpha=0.25, lw=0)
    ax.plot(t['f_repeat'], t['S_vgsr'], marker='o', color='C0', label=mode_label)
    ymax = t['S_none'][0][i_vgsr]
    if 'S_vgsr_nosb' in t.colnames:
        ax.plot(t['f_repeat'], t['S_vgsr_nosb'], marker='s', color='C1', ls='--', label='2 epochs, $\\sigma_{\\rm binary}$ = 0')
        ymax = max(ymax, t['S_vgsr_nosb'].max())
    ax.axhline(t['S_vgsr_sboracle'][0], color='C2', lw=1, ls=':', label='1 epoch + true rms $\\sigma_{\\rm binary}$')
    ax.axhline(t['S_true'][0][i_vgsr], color='k', lw=1.5, label='truth (noiseless CoM)')
    ax.axhline(t['S_all'][0][i_vgsr], color='k', lw=1, ls=':',
               label='%s, CoM (%s)' % (survey_label, '2 epochs' if t_scan.meta['SBMODE'] == 'fit' else 'epoch 0'))
    ax.axhline(t['S_none'][0][i_vgsr], color='0.5', lw=1, ls='-.', label='%s, 1 epoch, binaries in' % survey_label)
    ax.set_title("$r_{\\rm vir,0}$ = %.2f pc ($N$ = %i)" % (rvir, t['N_obs'][0]), fontsize=18)
    ax.set_ylim(0, 1.35*ymax) #<-- headroom for the legend
for ax in axs[-1]:
    ax.set_xlabel("$f_{\\rm repeat}$")
for ax in axs[:, 0]:
    ax.set_ylabel("single-component $\\sigma_{v_{\\rm gsr}}$ [km/s]")
axs[0, 0].legend(fontsize=12, loc='upper right')
# fig.savefig(plot_path+"gd1_%s_nominal_multiepoch_scan_%s%s.pdf" % (survey, scan_tag, dt_tag))

# %%
### sigma_binary vs f_repeat (median, 16-84% over draws), fitted ('fit') or estimated from dv ('dv'),
#   against the rms binary motion (dotted: per epoch, over every star; dashed: the epoch-to-epoch part).
t_scan = Table.read(table_path+scan_name+".fits", format='fits')
fig, ax = plt.subplots()
for k, rvir in enumerate(rvirs):
    t = t_scan[t_scan['Rvir0']==rvir]
    t = t[np.argsort(t['f_repeat'])]
    ax.fill_between(t['f_repeat'], t['sigma_binary_p16'], t['sigma_binary_p84'], color='C%i'%k, alpha=0.2, lw=0)
    ax.plot(t['f_repeat'], t['sigma_binary'], marker='o', color='C%i'%k, label='%.2f pc' % rvir)
    ax.axhline(t['sigma_binary_rms'][0], color='C%i'%k, ls=':', lw=1)
    ax.axhline(t['sigma_binary_rms_diff'][0], color='C%i'%k, ls='--', lw=1)
ax.set_xlabel("$f_{\\rm repeat}$")
ax.set_ylabel("$\\sigma_{\\rm binary}$ (%s) [km/s]" % ('fitted' if t_scan.meta['SBMODE'] == 'fit' else 'from $\\Delta v$'))
ax.legend(fontsize=12)
# %%
#-------------------------------------------------------#
#   A FINAL FIG THAT JUST HAS THE JARVIS BOX FOR THE THIN STREAM AND THE SINGLE COMPONENT FITS WITH `survey` NOISE
#   (moved here from results.py, 15 Nov). GD-1, <survey>_noise_binaries table, S_1comp = one gaussian fit
#   to the survey-selected stars (desi: the 679 brightest), with binary orbital motion in the RVs. joined in
#   rvir0 order, marker size = rvir0. (there's no via10hr gmm.py table, so survey='via10hr' fails here.)
show_CoM = True #<-- also draw <survey>_noise_CoM's S_1comp (same selection idea, CoM RVs, i.e. no binary motion)
add_errbars = False
show_noiseless_thin = True #<-- also draw the THIN component (S_ts) of the noiseless two-component fit, for comparison
noiseless_thin_case = "via_noise_CoM" #<-- or "noiseless_binaries" (thin sigma_vgsr ~0.05 km/s higher)
f_repeat_plot = 1.0#f_repeat_nominal #<-- which f_repeat of the repeat scan to draw as "analogous to Nibauer+26"
show_nibauer = True #<-- also draw the NUTS sigma_stream with Nibauer+26's Eq. 1 (independent epochs, shared=False),
                    #    pooled over random draws of the repeat stars (the "NUTS over random draws" cell below)
nibauer_table = "jake_nuts_N160_rep35_indep_capstd_ndraw10_pooled" #<-- the pooled table (one row per survey x rvir0 x baseline)
nibauer_survey, nibauer_baseline = 'desi', '5yr' #<-- the only ones in that table so far
nibauer_err = 'pooled' #<-- error bar: 'pooled' = 16-84% of the posterior pooled over the draws;
                       #    'draws' = 16-84% of the per-draw posterior medians (the spread from the draws only)
nibauer_focus = True   #<-- True: draw ONLY the rvir0 = nibauer_focus_rvir point, from nibauer_focus_table (50 draws of
                       #    both the 160 stars, from the full catalog, and the 35 repeats); False: every rvir0 from
                       #    nibauer_table (10 draws of the repeats, 160 brightest)
nibauer_focus_rvir = 3.
nibauer_focus_table = "jake_nuts_N160_rep35_indep_capstd_ndraw50_varystars_pooled"

### jarvis+26 published thin stream (as hardcoded at the top of results.py's summary figure)
jarvis_svgsr_ts = 2.49
jarvis_svgsr_ts_err = 0.28
jarvis_sphi2_ts = 0.23
jarvis_sphi2_ts_err = 0.01

def size_from_rvir(r):
    return 60*(np.asarray(r, dtype=float)/0.75)**0.7 # 60 -> ~250 pt^2 from 0.75 -> 6 pc

def gd1_color_from_results():
    """
    results.py's gd1_color: its palette (ccc) is handed out to the orbits in order of orbital
    phase, (r_now - r_peri)/(r_apo - r_peri). same orbit integration, so the colour matches.
    """
    lm_colors, hm_colors, _ = paf.define_simcolors()
    ccc = (hm_colors[:-1] + lm_colors[::-1])[:-1][1:]
    prog_tab = Table.read(repo_path+'/data/FINAL_ics_nolmc.csv')
    orbits = np.array(['gd1','aau','pa5','jet','c19'])
    phases = []
    for orbit in orbits:
        pr = prog_tab[prog_tab['name']==orbit]
        disp = [pr['x'][0], pr['y'][0], pr['z'][0], pr['vx'][0], pr['vy'][0], pr['vz'][0]]
        orbit_obj = paf.integrate_prog_orbit(disp, steps=10000, dt=1*u.Myr)
        peri = orbit_obj.pericenter().to(u.kpc).value
        apo = orbit_obj.apocenter().to(u.kpc).value
        phases.append((np.sqrt(np.sum(np.square(disp[:3]))) - peri) / (apo - peri))
    reordered = np.argsort(phases)
    return ccc[list(orbits[reordered]).index('gd1')]

gd1_color = gd1_color_from_results()

t_1c = Table.read(table_path+"%s_noise_binaries.fits" % survey, format='fits')
t_1c = t_1c[t_1c['orbit']=='gd1']
t_1c = t_1c[np.argsort(t_1c['Rvir0'])]

fig, ax = plt.subplots()

# jarvis+26 published thin stream, and this model fit to their data (as in the panel above)
ax.add_patch(Rectangle((jarvis_sphi2_ts-jarvis_sphi2_ts_err, jarvis_svgsr_ts-jarvis_svgsr_ts_err),
                       2*jarvis_sphi2_ts_err, 2*jarvis_svgsr_ts_err,
                       color='k', alpha=0.1, lw=0, zorder=0, label='GD-1 thin stream\nfrom Jarvis+26'))
# ax.plot(t_jarvis_fit['S_ts'][0], t_jarvis_fit['S_ts'][-1], marker='x', ls='--', color='k', markersize=10,
#         lw=1, zorder=6, label='This model fit to Jarvis+26 data')

# DESI-like noise, with binaries, one component: open grey circles, as they were on the panel above
ax.plot(t_1c['S_1comp'][:,0], t_1c['S_1comp'][:,-1], color='k', ls=':', lw=1, zorder=3, alpha=0.5)
ax.scatter(t_1c['S_1comp'][:,0], t_1c['S_1comp'][:,-1], marker='o',
           s=size_from_rvir(t_1c['Rvir0']),
           color='none', edgecolor='k', lw=1.5, zorder=4, alpha=0.5)#,
        #    label='DESI-like noise, with binaries;\nsingle component')

# the multi-epoch experiment (above): `survey` noise, its nominal sample, a fraction f_repeat_plot of
#   them at random with a second RV epoch; S_vgsr from the table of the current sb_mode settings (scan_name).
#   y = median over the random choices of the repeats, error bar = 16-84%.
rep_series = 'S_vgsr'
t_rep = Table.read(table_path+scan_name+".fits", format='fits')
t_rep = t_rep[np.isclose(t_rep['f_repeat'], f_repeat_plot)]
t_rep = t_rep[np.argsort(t_rep['Rvir0'])]
if add_errbars==True:
    ax.errorbar(t_rep['S_phi2'], t_rep[rep_series],
                yerr=[t_rep[rep_series] - t_rep[rep_series+'_p16'],
                t_rep[rep_series+'_p84'] - t_rep[rep_series]],
                fmt='none', ecolor=gd1_color, elinewidth=1.5, capsize=3, zorder=3)
ax.plot(t_rep['S_phi2'], t_rep[rep_series], color=gd1_color, ls='-', lw=1, zorder=3)
ax.scatter(t_rep['S_phi2'], t_rep[rep_series], marker='o',
           s=size_from_rvir(t_rep['Rvir0']),
           color=gd1_color, edgecolor='k', lw=1.5, zorder=4)#,
        #    label="Analogous measurement to Nibauer+26")#'DESI-like noise, %i brightest,\n%i repeats + $\\sigma_{\\rm binary}$;\nsingle component'
                #  % (t_rep['N_obs'][0], t_rep['N_repeat'][0]))

if show_CoM:
    t_1c_com = Table.read(table_path+"%s_noise_CoM.fits" % survey, format='fits')
    t_1c_com = t_1c_com[t_1c_com['orbit']=='gd1']
    t_1c_com = t_1c_com[np.argsort(t_1c_com['Rvir0'])]
    ax.plot(t_1c_com['S_1comp'][:,0], t_1c_com['S_1comp'][:,-1], color='k', ls=':', lw=1, zorder=3, alpha=0.5)
    ax.scatter(t_1c_com['S_1comp'][:,0], t_1c_com['S_1comp'][:,-1], marker='o',
               s=size_from_rvir(t_1c_com['Rvir0']),
               color='k', edgecolor='w', lw=1.5, zorder=4, alpha=0.5)#,
            #    label='DESI-like noise, CMs;\nsingle component')

# noiseless, every unbound star (trim_new & unbound), two-component fit: just the thin component.
#   i.e. what's left once the mixture model has taken the strong ejections out into the wide component.
if show_noiseless_thin:
    t_ts = Table.read(table_path+noiseless_thin_case+".fits", format='fits')
    t_ts = t_ts[t_ts['orbit']=='gd1']
    t_ts = t_ts[np.argsort(t_ts['Rvir0'])]
    ts_is_bin = 'binaries' in noiseless_thin_case
    ax.plot(t_ts['S_ts'][:,0], t_ts['S_ts'][:,-1], color=gd1_color, ls='--', lw=1, zorder=3)
    ax.scatter(t_ts['S_ts'][:,0], t_ts['S_ts'][:,-1], marker='o',
               s=size_from_rvir(t_ts['Rvir0']),
               color='none', edgecolor=gd1_color, lw=1.5, zorder=4)

# Nibauer+26's measurement (NUTS, independent epochs, 160 stars, 35 random repeats, pooled over draws). it has no
#   stream width, so x = the binary-CoM single-component sigma_phi2 (<survey>_noise_CoM's S_1comp) at the same
#   rvir0. y = pooled posterior median, error bar per nibauer_err. nibauer_focus: just the one rvir0, no line.
if show_nibauer:
    t_nib = Table.read(table_path+(nibauer_focus_table if nibauer_focus else nibauer_table)+".fits", format='fits')
    t_nib = t_nib[(t_nib['survey'] == nibauer_survey) & (t_nib['baseline'] == nibauer_baseline)]
    if nibauer_focus:
        t_nib = t_nib[np.isclose(t_nib['Rvir0'], nibauer_focus_rvir)]
        if len(t_nib) != 1:
            raise ValueError("expected one rvir0 = %g row in %s, got %i" % (nibauer_focus_rvir, nibauer_focus_table, len(t_nib)))
    t_nib = t_nib[np.argsort(t_nib['Rvir0'])]
    t_nib_x = Table.read(table_path+"%s_noise_CoM.fits" % survey, format='fits')
    t_nib_x = t_nib_x[t_nib_x['orbit']=='gd1']
    nib_x = np.array([t_nib_x['S_1comp'][np.isclose(t_nib_x['Rvir0'], r), 0][0] for r in t_nib['Rvir0']])
    nib_y = np.asarray(t_nib['sigma_stream'])
    if nibauer_err == 'pooled':
        nib_lo, nib_hi = np.asarray(t_nib['sigma_stream_p16']), np.asarray(t_nib['sigma_stream_p84'])
    elif nibauer_err == 'draws':
        nib_lo, nib_hi = np.asarray(t_nib['sigma_stream_draw_p16']), np.asarray(t_nib['sigma_stream_draw_p84'])
    else:
        raise ValueError("nibauer_err must be 'pooled' or 'draws', got %r" % (nibauer_err,))
    ax.errorbar(nib_x, nib_y, yerr=[nib_y - nib_lo, nib_hi - nib_y],
                fmt='none', ecolor=gd1_color, elinewidth=1.5, capsize=3, zorder=3)
    if not nibauer_focus:
        ax.plot(nib_x, nib_y, color=gd1_color, ls='-.', lw=1, zorder=3)
    ax.scatter(nib_x, nib_y, marker='D', s=0.8*size_from_rvir(t_nib['Rvir0']),
               color=gd1_color, edgecolor='k', lw=1.5, zorder=4)

# annotations: text at the left edge of the panel, arrow to the leftmost (smallest sigma_phi2) point of each series.
#   text_dy = vertical nudge of the text relative to that point [km/s], to move labels off each other / the data.
#   rank = which point to point at, counting up in sigma_phi2 (0 = leftmost = most diffuse progenitor, 1 = second, ...).
#   offset_pts = (dx, dy) [points]: if given, put the text at that offset from the point instead (left-aligned,
#   so (+dx, -dy) = below and to the right), and ignore text_x / text_dy.
def annotate_leftmost(xs, ys, text, text_dy=0.0, text_x=0.3, color='k', rank=0, offset_pts=None):
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    k = np.argsort(xs)[rank]
    if offset_pts is None:
        xytext, textcoords, ha = (text_x, ys[k]+text_dy), ('axes fraction', 'data'), 'right'
    else:
        xytext, textcoords, ha = offset_pts, 'offset points', 'left'
    ax.annotate(text, xy=(xs[k], ys[k]), xycoords='data',
                xytext=xytext, textcoords=textcoords,
                ha=ha, va='center', fontsize=13, zorder=7,
                arrowprops=dict(arrowstyle='->', color=color,
                                 lw=1, shrinkA=2, shrinkB=6))

annotate_leftmost(t_1c['S_1comp'][:,0], t_1c['S_1comp'][:,-1],
                  'including all\nbinary orbits', text_dy=-0.4, text_x=0.3, color='0.5')


annotate_leftmost(t_rep['S_phi2'], t_rep[rep_series],
                #   'analogous binary\naccounting\nto Nibauer+26',
                'Rejecting\nvariables',
                  text_dy=0.4,
                  color=gd1_color)
if show_CoM:
    annotate_leftmost(t_1c_com['S_1comp'][:,0],
                      t_1c_com['S_1comp'][:,-1],
                      'binary CMs',
                      rank=1, offset_pts=(25, -30)) #<-- second most diffuse point; text below + right of it
if show_noiseless_thin:
    annotate_leftmost(t_ts['S_ts'][:,0], t_ts['S_ts'][:,-1],
                      'thin component only' + (' (binaries)' if ts_is_bin else ''),
                      text_dy=0.25, text_x=0.36, color=gd1_color) #<-- text_x was 0.3 (the default)
if show_nibauer:
    if nibauer_focus:
        annotate_leftmost(nib_x, nib_y, 'Nibauer+26 method\n(NUTS, indep. epochs)',
                          offset_pts=(50, 100), color=gd1_color) #<-- hand-tuned: well above + right of the one point
    else:
        annotate_leftmost(nib_x, nib_y, 'Nibauer+26 method\n(NUTS, indep. epochs)',
                          rank=2, offset_pts=(25, 5), color=gd1_color) #<-- hand-tuned: right of the 3rd (highest) point

# marker-size key
# for r in [0.75, 1.5, 3., 6.]:
#     ax.scatter([None], [None],
#                marker='o',
#                s=size_from_rvir(r),
#                color='grey',
#                label=r'$R_{\rm vir,0}=%g$ pc'%r
#                )

# ax.set_title('GD-1, DESI-like noise, %i brightest; single component, with binaries' % n_obs, fontsize=14)
ax.set_xlabel(r'$\sigma_{\phi_2}~[\degree]$')
ax.set_ylabel(r'$\sigma_{v_{\rm GSR}}~[\rm km~s^{-1}]$')
ax.set_xlim(left=0)#, 2*max(np.max(t_1c['S_1comp'][:,0]), np.max(t_rep['S_phi2']),
                        # jarvis_sphi2_ts+jarvis_sphi2_ts_err))
ax.set_ylim(bottom=0)
# ax.legend(loc='lower left',
#           bbox_to_anchor=[0.45, 0.75],#[1.02, 1],
#            fontsize=11,
#            labelspacing=1.2)
# (legend replaced by an annotation on the jarvis box, arrow to its bottom right corner)
ax.annotate('GD-1 thin stream\nfrom Jarvis+26',
            xy=(jarvis_sphi2_ts+jarvis_sphi2_ts_err, jarvis_svgsr_ts-jarvis_svgsr_ts_err), xycoords='data',
            xytext=(30, -30), textcoords='offset points', #<-- text position relative to the corner [points]
            ha='left', va='top', fontsize=13, zorder=7,
            arrowprops=dict(arrowstyle='->', color='0.7', lw=1, shrinkA=2, shrinkB=0))

fname = plot_path+"jarvis_comparison_single_component"
if add_errbars==True:
    fname+="_with_errs"
plt.savefig(fname+".pdf",
            dpi=300, bbox_inches='tight')

# %%
### per-draw NUTS medians: sigma_stream vs sigma_binary, one point per draw of the 160 stars + 35 repeats ###
# the raw posterior samples weren't saved, so this is NOT the joint posterior: each point is one draw's marginal
#   posterior medians, with that draw's marginal 16-84% as error bars in x and y (the within-draw correlation isn't
#   known; with independent epochs each draw's posterior is a thin arc sigma_stream^2 + sigma_binary^2 ~ const).
# big square: the pooled (all draws' samples concatenated) median and 16-84% in each sigma.
# reads the per-draw and pooled tables of the Jarvis figure's nibauer_focus_table (run that cell first), for
#   nibauer_survey / nibauer_baseline / nibauer_focus_rvir. points coloured by the number of binaries among the
#   draw's 35 repeats. truth: the noiseless CoM dispersion of the draw's 160 stars (horizontal) and the rms true
#   binary motion over them at the two epochs used (vertical); line = median over draws, band = their 16-84%.
med_tab = Table.read(table_path+nibauer_focus_table.replace('_pooled', '')+".fits", format='fits')
med_pool = Table.read(table_path+nibauer_focus_table+".fits", format='fits')
sel_ = lambda t: t[(t['survey'] == nibauer_survey) & (t['baseline'] == nibauer_baseline)
                   & np.isclose(t['Rvir0'], nibauer_focus_rvir)]
med_tab, med_pool = sel_(med_tab), sel_(med_pool)[0]

fig, ax = plt.subplots(figsize=(8, 8))
sb, ss = np.asarray(med_tab['sigma_binary']), np.asarray(med_tab['sigma_stream'])
ax.errorbar(sb, ss, xerr=[sb - med_tab['sigma_binary_p16'], med_tab['sigma_binary_p84'] - sb],
            yerr=[ss - med_tab['sigma_stream_p16'], med_tab['sigma_stream_p84'] - ss],
            fmt='none', ecolor='0.6', elinewidth=0.8, alpha=0.5, zorder=2)
sc_ = ax.scatter(sb, ss, c=med_tab['n_bin_repeat'], cmap='viridis', s=45, edgecolor='k', lw=0.6, zorder=3)
fig.colorbar(sc_, ax=ax, label='binaries among the %i repeats' % med_pool['N_repeat'], fraction=0.046, pad=0.03)

# pooled
ax.errorbar(med_pool['sigma_binary'], med_pool['sigma_stream'],
            xerr=[[med_pool['sigma_binary'] - med_pool['sigma_binary_p16']],
                  [med_pool['sigma_binary_p84'] - med_pool['sigma_binary']]],
            yerr=[[med_pool['sigma_stream'] - med_pool['sigma_stream_p16']],
                  [med_pool['sigma_stream_p84'] - med_pool['sigma_stream']]],
            fmt='s', color=gd1_color, mec='k', ms=12, elinewidth=2.5, capsize=0, zorder=5,
            label='pooled over %i draws' % len(med_tab))

# truths: median over draws (line) and 16-84% over draws (band)
for col, axline, axspan in [('S_true', ax.axhline, ax.axhspan), ('sigma_binary_rms', ax.axvline, ax.axvspan)]:
    lo, mid, hi = np.percentile(med_tab[col], [16, 50, 84])
    axspan(lo, hi, color='k', alpha=0.08, lw=0, zorder=0)
    axline(mid, color='k', ls='--', lw=1.2, zorder=1,
           label='truth (noiseless CoM)' if col == 'S_true' else 'true rms binary motion')

lim = 1.05*max(np.max(med_tab['sigma_binary_p84']), np.max(med_tab['sigma_stream_p84']))
ax.plot([0, lim], [0, lim], color='k', ls=':', lw=1, zorder=1, label='1:1')
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.set_xlabel(r'$\sigma_{\rm binary}$ [km/s] (per-draw posterior median)')
ax.set_ylabel(r'$\sigma_{v_r,\,\rm stream}$ [km/s] (per-draw posterior median)')
ax.set_title('%s-like, $r_{\\rm vir,0}$ = %g pc, %s baseline, %s epochs'
             % ({'desi': 'DESI', 'via': 'Via'}[nibauer_survey], nibauer_focus_rvir,
                nibauer_baseline.replace('yr', ' yr'), med_tab.meta['LIKELIHD']), fontsize=14)
ax.legend(fontsize=11, loc='upper left')
# plt.savefig(plot_path+"gd1_%s_%s_rvir%g_%s_medians.pdf" % (nibauer_survey, nibauer_focus_table.replace('_pooled', ''),
#                                                           nibauer_focus_rvir, nibauer_baseline),
#             dpi=300, bbox_inches='tight')

# %%
### joint and marginal posteriors of sigma_binary and sigma_stream, from the saved NUTS samples ###
# reads data/nuts_samples/<per-draw table name>_<survey>_rvir<X>_<baseline>.npz, written by the "NUTS over random
#   draws" cell (same nibauer_* settings as the two cells above; run the Jarvis cell first).
# main panel: the POOLED joint posterior (every draw's samples, equal weight) as a 2D histogram, with its 68% and
#   95% highest-density contours; optionally (show_draw_contours) each draw's own 68% contour, faint.
# side panels: the marginals, all as unit-area densities (so a single draw's is taller than the pooled one).
#   pooled = filled histogram with its median (solid) and 16-84% (dotted); each draw's own = faint step histograms.
# truth (black dashed): the noiseless CoM dispersion of the 160 stars and the rms true binary motion over them,
#   median over draws. 1:1 dotted. note sigma_binary's prior is truncated at each draw's own sb_cap.
# the contours are drawn on 2D histograms smoothed by a gaussian of contour_smooth bins (the grey image isn't).
show_draw_contours = True #<-- each draw's 68% contour in the joint panel
n_bins_joint = 80
contour_smooth = 2.5      #<-- [bins] gaussian smoothing before contouring
marg_ymax = 2.5           #<-- marginal panels' density axis limit, in units of the pooled histogram's peak
from scipy.ndimage import gaussian_filter

samp_file = repo_path+"/data/nuts_samples/%s_%s_rvir%g_%s.npz" % (nibauer_focus_table.replace('_pooled', ''),
                                                                  nibauer_survey, nibauer_focus_rvir, nibauer_baseline)
samp = np.load(samp_file)
ss_d, sb_d = samp['sigma_stream'], samp['sigma_binary'] #<-- (n_draws, n_samples)
ss_all, sb_all = ss_d.ravel(), sb_d.ravel()
n_draws_s = ss_d.shape[0]

def hpd_levels(H, fracs):
    """density thresholds of a 2D histogram enclosing the given fractions of its mass, highest density first."""
    h = np.sort(H.ravel())[::-1]
    c = np.cumsum(h) / h.sum()
    return [h[min(np.searchsorted(c, f), len(h) - 1)] for f in fracs]

lim_j = np.percentile(np.concatenate([ss_all, sb_all]), 99.5)
edges = np.linspace(0, lim_j, n_bins_joint + 1)
centers = 0.5*(edges[1:] + edges[:-1])

fig = plt.figure(figsize=(9, 9))
gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], wspace=0.04, hspace=0.04)
ax_j = fig.add_subplot(gs[1, 0])
ax_x = fig.add_subplot(gs[0, 0], sharex=ax_j)
ax_y = fig.add_subplot(gs[1, 1], sharey=ax_j)

# joint: pooled 2D histogram + 68/95% contours (+ each draw's 68%)
H, _, _ = np.histogram2d(sb_all, ss_all, bins=[edges, edges])
ax_j.pcolormesh(edges, edges, H.T, cmap='Greys', shading='flat', rasterized=True, zorder=0)
Hs = gaussian_filter(H, contour_smooth)
l95, l68 = hpd_levels(Hs, [0.95, 0.68])
for lev, a in [(l95, 0.15), (l68, 0.25)]: #<-- shade inside the 95% and 68% regions
    ax_j.contourf(centers, centers, Hs.T, levels=[lev, 1.01*Hs.max()], colors=[gd1_color], alpha=a, zorder=1)
ax_j.contour(centers, centers, Hs.T, levels=[l95, l68], colors=[gd1_color], linewidths=[1.5, 2.5], zorder=3)
if show_draw_contours:
    for k in range(n_draws_s):
        Hk, _, _ = np.histogram2d(sb_d[k], ss_d[k], bins=[edges, edges])
        Hk = gaussian_filter(Hk, contour_smooth)
        ax_j.contour(centers, centers, Hk.T, levels=hpd_levels(Hk, [0.68]), colors=['C0'], linewidths=0.6,
                     alpha=0.35, zorder=2)
ax_j.plot([0, lim_j], [0, lim_j], color='k', ls=':', lw=1, zorder=1)
S_true_med, sb_rms_med = np.median(samp['S_true']), np.median(samp['sigma_binary_rms'])
ax_j.axhline(S_true_med, color='k', ls='--', lw=1.2, zorder=4)
ax_j.axvline(sb_rms_med, color='k', ls='--', lw=1.2, zorder=4)
ax_j.set_xlim(0, lim_j)
ax_j.set_ylim(0, lim_j)
ax_j.set_xlabel(r'$\sigma_{\rm binary}$ [km/s]')
ax_j.set_ylabel(r'$\sigma_{v_r,\,\rm stream}$ [km/s]')

# marginals
for ax_m, d_draws, d_all, truth_val, orient in [(ax_x, sb_d, sb_all, sb_rms_med, 'vertical'),
                                                 (ax_y, ss_d, ss_all, S_true_med, 'horizontal')]:
    for k in range(n_draws_s):
        ax_m.hist(d_draws[k], bins=edges, density=True, orientation=orient, histtype='step', color='C0', lw=0.6,
                  alpha=0.35, zorder=1)
    n_pool, _, _ = ax_m.hist(d_all, bins=edges, density=True, orientation=orient, color=gd1_color, alpha=0.35,
                             lw=0, zorder=2)
    ax_m.hist(d_all, bins=edges, density=True, orientation=orient, histtype='step', color=gd1_color, lw=2, zorder=3)
    (ax_m.set_ylim if orient == 'vertical' else ax_m.set_xlim)(0, marg_ymax*np.max(n_pool))
    p16, p50, p84 = np.percentile(d_all, [16, 50, 84])
    line = ax_m.axvline if orient == 'vertical' else ax_m.axhline
    line(p50, color=gd1_color, ls='-', lw=1.5)
    for q in (p16, p84):
        line(q, color=gd1_color, ls=':', lw=1.5)
    line(truth_val, color='k', ls='--', lw=1.2)
ax_x.tick_params(labelbottom=False, labelleft=False)
ax_y.tick_params(labelleft=False, labelbottom=False)
ss16, ss50, ss84 = np.percentile(ss_all, [16, 50, 84])
sb16, sb50, sb84 = np.percentile(sb_all, [16, 50, 84])
ax_x.set_title('%s-like, $r_{\\rm vir,0}$ = %g pc, %s, %i draws\n'
               % ({'desi': 'DESI', 'via': 'Via'}[nibauer_survey], nibauer_focus_rvir,
                  nibauer_baseline.replace('yr', ' yr'), n_draws_s)
               + r'$\sigma_{\rm stream} = %.2f^{+%.2f}_{-%.2f}$, $\sigma_{\rm binary} = %.2f^{+%.2f}_{-%.2f}$ km/s'
               % (ss50, ss84 - ss50, ss50 - ss16, sb50, sb84 - sb50, sb50 - sb16), fontsize=13)
ax_c = fig.add_subplot(gs[0, 1]) # legend in the empty corner
ax_c.axis('off')
ax_c.legend(handles=[plt.Line2D([], [], color=gd1_color, lw=2.5, label='pooled\n(68, 95%)'),
                     plt.Line2D([], [], color='C0', lw=1, alpha=0.6, label='single draws'),
                     plt.Line2D([], [], color='k', ls='--', lw=1.2, label='truth')],
            loc='center', fontsize=10, frameon=False)
# plt.savefig(plot_path+"gd1_%s_%s_rvir%g_%s_joint.pdf" % (nibauer_survey, nibauer_focus_table.replace('_pooled', ''),
#                                                         nibauer_focus_rvir, nibauer_baseline),
#             dpi=300, bbox_inches='tight')

# %%
### joint and marginal posteriors for ONE draw of 160 stars / 35 repeats ###
# the same plot as above for a single NUTS fit, i.e. what one Nibauer+26-like measurement gives. no new sampling:
#   each draw in the .npz above is already its own fit, so this reads draw single_draw from samp (run the cell
#   above first; it also defines hpd_levels, gaussian_filter, n_bins_joint, contour_smooth).
# single_draw: an int (0 .. n_draws-1), or 'median' = the draw whose median total width
#   sqrt(sigma_stream^2 + sigma_binary^2) is the median over draws (a typical draw; the totals range ~1-6 km/s).
# truth (black dashed): THIS draw's noiseless CoM dispersion of its 160 stars and rms true binary motion; grey
#   dashed arc: the true total width sqrt(S_true^2 + rms^2). red dash-dot in the top marginal: this draw's sb_cap
#   (the truncation of sigma_binary's prior).
single_draw = 16# 'median'

r_draws = np.median(np.hypot(samp['sigma_stream'], samp['sigma_binary']), axis=1)
k1 = int(np.argsort(r_draws)[len(r_draws)//2]) if single_draw == 'median' else int(single_draw)
ss1, sb1 = samp['sigma_stream'][k1], samp['sigma_binary'][k1]
S_true1, sb_rms1, sb_cap1 = samp['S_true'][k1], samp['sigma_binary_rms'][k1], samp['sb_cap'][k1]

lim_1 = 1.05*max(np.percentile(np.concatenate([ss1, sb1]), 99.5), np.hypot(S_true1, sb_rms1))
edges1 = np.linspace(0, lim_1, n_bins_joint + 1)
centers1 = 0.5*(edges1[1:] + edges1[:-1])

fig = plt.figure(figsize=(9, 9))
gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], wspace=0.04, hspace=0.04)
ax_j = fig.add_subplot(gs[1, 0])
ax_x = fig.add_subplot(gs[0, 0], sharex=ax_j)
ax_y = fig.add_subplot(gs[1, 1], sharey=ax_j)

H1, _, _ = np.histogram2d(sb1, ss1, bins=[edges1, edges1])
ax_j.pcolormesh(edges1, edges1, H1.T, cmap='Greys', shading='flat', rasterized=True, zorder=0)
H1s = gaussian_filter(H1, contour_smooth)
l95, l68 = hpd_levels(H1s, [0.95, 0.68])
for lev, a in [(l95, 0.15), (l68, 0.25)]:
    ax_j.contourf(centers1, centers1, H1s.T, levels=[lev, 1.01*H1s.max()], colors=[gd1_color], alpha=a, zorder=1)
ax_j.contour(centers1, centers1, H1s.T, levels=[l95, l68], colors=[gd1_color], linewidths=[1.5, 2.5], zorder=3)
ax_j.plot([0, lim_1], [0, lim_1], color='k', ls=':', lw=1, zorder=1)
th = np.linspace(0, np.pi/2, 200)
# ax_j.plot(np.hypot(S_true1, sb_rms1)*np.cos(th), np.hypot(S_true1, sb_rms1)*np.sin(th), color='0.4', ls='--', lw=1,
#           zorder=4)
# ax_j.axhline(S_true1, color='k', ls='--', lw=1.2, zorder=4)
# ax_j.axvline(sb_rms1, color='k', ls='--', lw=1.2, zorder=4)
ax_j.set_xlim(0, lim_1)
ax_j.set_ylim(0, lim_1)
ax_j.set_xlabel(r'$\sigma_{\rm binary}$ [km/s]')
ax_j.set_ylabel(r'$\sigma_{v_r,\,\rm stream}$ [km/s]')

for ax_m, d1, truth_val, orient in [(ax_x, sb1, sb_rms1, 'vertical'), (ax_y, ss1, S_true1, 'horizontal')]:
    ax_m.hist(d1, bins=edges1, density=True, orientation=orient, color=gd1_color, alpha=0.35, lw=0)
    ax_m.hist(d1, bins=edges1, density=True, orientation=orient, histtype='step', color=gd1_color, lw=2)
    p16, p50, p84 = np.percentile(d1, [16, 50, 84])
    line = ax_m.axvline if orient == 'vertical' else ax_m.axhline
    line(p50, color=gd1_color, ls='-', lw=1.5)
    for q in (p16, p84):
        line(q, color=gd1_color, ls=':', lw=1.5)
    # line(truth_val, color='k', ls='--', lw=1.2)
# if sb_cap1 < lim_1:
#     ax_x.axvline(sb_cap1, color=gd1_color, ls='-.', lw=1.2)
ax_x.set_ylim(bottom=0)
ax_y.set_xlim(left=0)
ax_x.tick_params(labelbottom=False, labelleft=False)
ax_y.tick_params(labelleft=False, labelbottom=False)
ss16, ss50, ss84 = np.percentile(ss1, [16, 50, 84])
sb16, sb50, sb84 = np.percentile(sb1, [16, 50, 84])
ax_x.set_title('%s-like, $r_{\\rm vir,0}$ = %g pc, %s, draw %i (%i binaries, %i in repeats)\n'
               % ({'desi': 'DESI', 'via': 'Via'}[nibauer_survey], nibauer_focus_rvir,
                  nibauer_baseline.replace('yr', ' yr'), k1, samp['n_bin'][k1], samp['n_bin_repeat'][k1])
               + r'$\sigma_{\rm stream} = %.2f^{+%.2f}_{-%.2f}$, $\sigma_{\rm binary} = %.2f^{+%.2f}_{-%.2f}$ km/s'
               % (ss50, ss84 - ss50, ss50 - ss16, sb50, sb84 - sb50, sb50 - sb16), fontsize=13)
ax_c = fig.add_subplot(gs[0, 1])
ax_c.axis('off')
# ax_c.legend(handles=[plt.Line2D([], [], color=gd1_color, lw=2.5, label='posterior\n(68, 95%)'),
#                      plt.Line2D([], [], color='k', ls='--', lw=1.2, label='truth'),
#                      plt.Line2D([], [], color='0.4', ls='--', lw=1, label='true total'),
#                      plt.Line2D([], [], color=gd1_color, ls='-.', lw=1.2, label=r'$\sigma_{\rm binary}$ cap')],
#             loc='center', fontsize=10, frameon=False)
# plt.savefig(plot_path+"gd1_%s_%s_rvir%g_%s_draw%i_joint.pdf" % (nibauer_survey,
#             nibauer_focus_table.replace('_pooled', ''), nibauer_focus_rvir, nibauer_baseline, k1),
#             dpi=300, bbox_inches='tight')

 # %%
### assemble catalogs for jake ###
# mock GD-1 RV catalogs for Jake Nibauer to run his own dispersion measurement on: one file per
#   rvir0 x survey, for_jake/rvir0<X>_<desi|via>.fits. needs only the imports/settings cell and the
#   function cell above (not the scan or the plot cells).
#
# which stars: the catalog-matched stars of the binaries case (cut_for_catalog_photometry_primaries &
#   matched_to_catalog_photometry_primaries), i.e. the progenitor (2 rtid) and remnants removed, and the
#   1-99 percentile trim + |v_gsr,primary| < 30 km/s outlier cuts applied. NO error-based quality cuts
#   (no pm_err, no rverr cut; Jake does those):
#     'via':  every matched star = the length of the Bonaca & Price-Whelan catalog (1578 for GD-1). the
#             stars whose Via error is NaN (log_Teff outside the colour-Teff fit, 1 for GD-1) are kept,
#             with NaN vr_k and vr_k_err.
#     'desi': the N_jarvis = 679 brightest of those in catalog G, the stars that get the Jarvis+26 V_ERRs.
#             so the DESI stars are a subset of the Via stars.
# three RV epochs per star at t_0 = 0, t_1 ~ N(1 yr, 2 months), t_2 ~ N(5 yr, 6 months). singles sit at their
#   CoM v_gsr; binaries move on a Keplerian RV curve about it (paf.get_orbital_params_from_binaries +
#   paf.get_rvs, as in load_mock_sample). every epoch gets a fresh noise draw N(0, rverr_<survey>); the error
#   is the same for all three epochs. the binary orbits (inclination, periapsis argument, phase) are drawn once
#   per rvir0 and shared by the two surveys; epoch times and noise are drawn per survey.
# positions: phi1 and the luminous component's phi2, plus the saved Gaia DR3 position noise draw
#   (noise_catalog_photometry_primaries), with err = pygaia's total position uncertainty / sqrt 2, as in
#   final_datasets.survey_noise.
# phi2 and vr_k are RESIDUALS about our stream track (orbit + degree-5 polynomial straightening, as fit by
#   gmm.py); phi2_track / vr_track are that track at each star's phi1, so phi2 + phi2_track and
#   vr_k + vr_track are the absolute GD1Koposov10 phi2 and v_gsr (reflex-corrected) of the system's CoM track.
#   phi1 is absolute.
# the truth (is_binary, noiseless CoM and per-epoch v_gsr, orbital params) goes to data/for_jake_truth/, which
#   is not for sending.
from pygaia.errors.astrometric import total_position_uncertainty

jake_path = repo_path+"/for_jake/"
jake_truth_path = repo_path+"/data/for_jake_truth/"
jake_surveys = ['desi', 'via']
jake_epochs_yr = [(0., 0.), (1., 2./12.), (5., 6./12.)] #<-- (mean, sd) of each epoch's time after epoch 0 [yr]
jake_seed = 2026

def jake_catalog(dd, survey, orb_params_all, rng):
    """
    one mock catalog (see the cell comment). dd = the gd1 pickle; orb_params_all = Keplerian params for
    every binary in the snapshot, row j = binary j (paf.get_orbital_params_from_binaries).
    returns (catalog Table, truth Table), both ordered by phi1.
    """
    sc_com, sc_pri = dd['sc_straighter'], dd['sc_straighter_primaries']
    coords_obs = dd['coords_obs']
    N_all, ns = len(sc_com['phi1']), dd['nsingles']
    is_binary = np.arange(N_all) >= ns

    cut_p = dd['cut_for_catalog_photometry_primaries']
    matched_p = np.zeros(N_all, dtype=bool)
    matched_p[cut_p] = dd['matched_to_catalog_photometry_primaries']
    mG = dd['catalog_photometry_primaries']['mG']
    nd = dd['noise_catalog_photometry_primaries']

    idx = np.flatnonzero(matched_p)
    if survey == 'desi':
        idx = idx[np.argsort(mG[idx])[:N_jarvis]] #<-- the N_jarvis brightest = the ones with Jarvis V_ERRs
        if not np.all(np.isfinite(nd['rverr_desi'][idx])):
            raise ValueError("non-finite DESI errors among the %i brightest matched stars" % N_jarvis)
    elif survey != 'via':
        raise ValueError("survey must be 'desi' or 'via', got %r" % (survey,))
    idx = idx[np.argsort(sc_com['phi1'][idx])]
    N = len(idx)
    b = is_binary[idx]

    rverr = nd['rverr_'+survey][idx]
    pos_err = ((total_position_uncertainty(mG[idx], 'dr3') / np.sqrt(2)) * u.microarcsecond).to(u.deg).value

    # epoch times [yr] and noiseless velocities (residuals about the track)
    t_yr = np.column_stack([rng.normal(m, s, size=N) if s > 0 else np.full(N, m) for m, s in jake_epochs_yr])
    v_com = sc_com['v_gsr'][idx]
    v_true = np.repeat(v_com[:, None], len(jake_epochs_yr), axis=1)
    v_true[b] = paf.get_rvs(orb_params_all[idx[b] - ns], t_yr[b]*365.25, verbose=False) #<-- days, as P is
    v_obs = v_true + rng.normal(0., 1., size=v_true.shape) * rverr[:, None]

    # track at each star's phi1 = absolute (CoM) - residual
    phi2_track = coords_obs.phi2.deg[idx] - sc_com['phi2'][idx]
    vr_track = coords_obs.radial_velocity.to(u.km/u.s).value[idx] - sc_com['v_gsr'][idx]

    cat = Table()
    cat['phi1'] = (sc_com['phi1'][idx] + nd['phi1'][idx]) * u.deg
    cat['phi1_err'] = pos_err * u.deg
    cat['phi2'] = (sc_pri['phi2'][idx] + nd['phi2'][idx]) * u.deg
    cat['phi2_err'] = pos_err * u.deg
    for k in range(len(jake_epochs_yr)):
        cat['vr_%i' % k] = v_obs[:, k] * u.km/u.s
        cat['vr_%i_err' % k] = rverr * u.km/u.s
    for k in range(len(jake_epochs_yr)):
        cat['t_%i' % k] = t_yr[:, k] * u.yr
    cat['phi2_track'] = phi2_track * u.deg
    cat['vr_track'] = vr_track * u.km/u.s

    truth = Table()
    truth['phi1'] = cat['phi1']
    truth['is_binary'] = b
    truth['phi2_true'] = sc_pri['phi2'][idx] * u.deg
    truth['v_com'] = v_com * u.km/u.s
    for k in range(len(jake_epochs_yr)):
        truth['v_true_%i' % k] = v_true[:, k] * u.km/u.s
    op = np.full((N, 6), np.nan)
    op[b] = orb_params_all[idx[b] - ns]
    truth['orb_params'] = op #<-- [v0, K, w, phi0, e, P (day)], NaN for singles
    truth['mG'] = mG[idx]
    return cat, truth

os.makedirs(jake_truth_path, exist_ok=True)
for i_r, rvir in enumerate(rvirs):
    with open(datapath+"gd1_%.2f.pickle"%rvir, 'rb') as handle:
        dd = pickle.load(handle)
    rng_jake = np.random.default_rng(seed=[jake_seed, i_r]) #<-- one stream per rvir0, so files don't depend on run order
    ns = dd['nsingles']
    orb_params_all = paf.get_orbital_params_from_binaries(load_present_day_binaries(dd, rvir),
                                                          dd['sc_straighter']['v_gsr'][ns:], rng=rng_jake)
    for survey_j in jake_surveys:
        cat, truth = jake_catalog(dd, survey_j, orb_params_all, rng_jake)
        for t in (cat, truth):
            t.meta.update({'STREAM': 'GD-1 (mock)', 'RVIR0': rvir, 'SURVEY': survey_j, 'NEPOCH': len(jake_epochs_yr),
                           'SEED': jake_seed})
        for k, (m, s) in enumerate(jake_epochs_yr):
            cat.meta['T%iMEAN' % k], cat.meta['T%iSD' % k] = m, s
        cat.meta['COMMENTS'] = ['phi2, vr_k: residuals about the stream track; add phi2_track, vr_track for absolute',
                                'phi1, phi2 in the GD1Koposov10 frame; vr = reflex-corrected v_gsr',
                                'noise included; no error-based quality cuts']
        fname = "rvir0%g_%s" % (rvir, survey_j)
        cat.write(jake_path+fname+".fits", overwrite=True)
        truth.write(jake_truth_path+fname+"_truth.fits", overwrite=True)
        print("%s: N = %i, binaries %i, NaN vr_0 %i" % (fname, len(cat), truth['is_binary'].sum(),
                                                         np.sum(~np.isfinite(cat['vr_0']))))

# %%
### check that input table looks good.
# t = Table.read('for_jake/rvir03_via.fits', format='fits')
# %%
### NUTS on Jake's catalogs: sigma_binary and sigma_stream (v_gsr), Nibauer+26 priors ###
# posteriors for the stream's v_gsr dispersion and a global sigma_binary from the mock catalogs in for_jake/
#   (the cell above), sampled with numpyro's NUTS as Nibauer+26 do. no phi2 / stream width.
# the sample: per rvir0 x survey, the nuts_n_obs = 160 brightest stars in catalog G (mG from the truth file;
#   Nibauer+26: 160 GD-1 stars) with finite errors, of which nuts_n_repeat = 35 random ones have a second epoch
#   (Nibauer+26: 35). the same 35 stars for both baselines and both surveys (one rng per rvir0).
# the baselines: every star has vr_0; a repeat star's second epoch is vr_1 (~1 yr later) or vr_2 (~5 yr later).
#   one table row each.
# the likelihood: multiepoch_nll_per_star (with xp = jax.numpy), on the (N, 2) arrays from epoch_arrays,
#   i.e. exactly the ML fits' likelihood. nuts_shared_stream = True is the shared stream term; False is the
#   independent-epoch form of Nibauer+26's Eq. 1 (sigma_binary and sigma_stream degenerate, see the readme).
# the priors: sigma_stream ~ HalfNormal(10 km/s), sigma_binary ~ HalfNormal(10 km/s) TRUNCATED at sb_cap, the
#   repeats' epoch-to-epoch scatter: 'std' = std(dv)/sqrt 2 (dv = second epoch - vr_0; the per-epoch scatter a
#   difference of two epochs implies, noise included, so an upper bound on the binary term), or 'max' =
#   max |dv|/sqrt 2 (the largest single star's). mu ~ N(0, 10 km/s) (velocities are residuals about the track;
#   Nibauer+26's zero-point prior).
# the truth columns (from data/for_jake_truth/): S_true = std of the noiseless CoM v_gsr over the sample;
#   sigma_binary_rms = rms of the true binary motion (v_true_k - v_com) over the sample at both epochs used.
# outputs: data/gmm_tables/jake_nuts_N<n_obs>_rep<n_repeat>_<shared|indep>_cap<std|max>.fits (one row per
#   survey x rvir0 x baseline) and plots/gd1_desi_jake_nuts_<same tag>.pdf. needs the imports/settings and
#   function cells and the for_jake/ files (not the scan); ~1 min.
import jax
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
from numpyro.diagnostics import summary as numpyro_summary
numpyro.enable_x64() #<-- missing epochs have err^2 = missing_epoch_err^2 = 1e8: float32 isn't enough
jnp = jax.numpy

jake_path = repo_path+"/for_jake/"
jake_truth_path = repo_path+"/data/for_jake_truth/"

nuts_n_obs = 160          # the brightest N stars (Nibauer+26: 160)
nuts_n_repeat = 35        # of which this many (random) have a second epoch (Nibauer+26: 35)
nuts_baselines = {'1yr': 1, '5yr': 2} #<-- name -> which of Jake's epochs is the repeats' second one
nuts_surveys = ['desi', 'via']
nuts_prior_scale = 10.    # km/s, half-normal scale on sigma_stream and sigma_binary (Nibauer+26)
nuts_mu_prior_scale = 10. # km/s, N(0, .) on mu
nuts_shared_stream = False # True: shared stream term; False: Nibauer+26's Eq. 1 (independent epochs)
nuts_sb_cap = 'std'       # 'std': sigma_binary <= std(dv)/sqrt 2 over the repeats; 'max': <= max |dv|/sqrt 2
nuts_seed = 7
nuts_mcmc_kwargs = dict(num_warmup=1000, num_samples=2000, num_chains=4, chain_method='sequential',
                        progress_bar=False)
nuts_tag = "N%i_rep%i_%s_cap%s" % (nuts_n_obs, nuts_n_repeat, 'shared' if nuts_shared_stream else 'indep',
                                   nuts_sb_cap)

def jake_nuts_sample(rvir, survey, n_obs, n_repeat, rng):
    """
    the brightest n_obs stars (finite errors and velocities) of Jake's catalog for (rvir, survey), and a random
    n_repeat of them flagged as repeats. returns (cat, truth, repeat), rows in catalog order.
    """
    fname = "rvir0%g_%s" % (rvir, survey)
    cat = Table.read(jake_path+fname+".fits")
    truth = Table.read(jake_truth_path+fname+"_truth.fits")
    if not np.array_equal(np.asarray(cat['phi1']), np.asarray(truth['phi1'])):
        raise ValueError(fname+": catalog and truth rows don't line up")
    v3 = np.column_stack([np.asarray(cat['vr_%i' % k]) for k in range(cat.meta['NEPOCH'])])
    ok = np.isfinite(np.asarray(cat['vr_0_err'])) & np.all(np.isfinite(v3), axis=1)
    idx = np.flatnonzero(ok)
    idx = np.sort(idx[np.argsort(np.asarray(truth['mG'])[idx])[:n_obs]])
    repeat = np.zeros(len(idx), dtype=bool)
    repeat[rng.choice(len(idx), size=n_repeat, replace=False)] = True
    return cat[idx], truth[idx], repeat

def sigma_binary_cap(v_epochs, repeat, how='std'):
    """upper bound on sigma_binary from the repeats' dv = v_1 - v_0 (observed): std(dv)/sqrt 2 or max|dv|/sqrt 2."""
    dv = v_epochs[repeat, 1] - v_epochs[repeat, 0]
    if how == 'std':
        return np.std(dv, ddof=1)/np.sqrt(2)
    if how == 'max':
        return np.max(np.abs(dv))/np.sqrt(2)
    raise ValueError("how must be 'std' or 'max', got %r" % (how,))

def nuts_model(v, v_err, sb_cap, shared_stream=True):
    """numpyro model: Nibauer+26's half-normal priors, sigma_binary truncated at sb_cap; the multi-epoch likelihood."""
    mu = numpyro.sample('mu', dist.Normal(0., nuts_mu_prior_scale))
    sigma_stream = numpyro.sample('sigma_stream', dist.HalfNormal(nuts_prior_scale))
    sigma_binary = numpyro.sample('sigma_binary', dist.TruncatedNormal(0., nuts_prior_scale, low=0., high=sb_cap))
    numpyro.factor('lnL', -jnp.sum(multiepoch_nll_per_star(v, v_err, mu, sigma_stream, sigma_binary,
                                                            shared_stream=shared_stream, xp=jnp)))

print("running NUTS sampler")
nuts_rows, nuts_samples = [], {}
for i_r, rvir in tqdm(enumerate(rvirs)):
    for survey_n in nuts_surveys:
        rng_nuts = np.random.default_rng(seed=[nuts_seed, i_r]) #<-- same repeats for both surveys and baselines
        cat, truth, repeat = jake_nuts_sample(rvir, survey_n, nuts_n_obs, nuts_n_repeat, rng_nuts)
        rverr = np.asarray(cat['vr_0_err'])
        is_bin = np.asarray(truth['is_binary'])
        v_com = np.asarray(truth['v_com'])
        for bl, k2 in nuts_baselines.items():
            v_epochs = np.column_stack([np.asarray(cat['vr_0']), np.asarray(cat['vr_%i' % k2])])
            v, v_err = epoch_arrays(v_epochs, rverr, repeat)
            sb_cap = sigma_binary_cap(v_epochs, repeat, nuts_sb_cap)
            mcmc = MCMC(NUTS(nuts_model), **nuts_mcmc_kwargs)
            mcmc.run(jax.random.PRNGKey(nuts_seed), #<-- FITS columns are big-endian (>f8), which jax refuses
                     v=jnp.asarray(v.astype(np.float64)), v_err=jnp.asarray(v_err.astype(np.float64)), sb_cap=sb_cap,
                     shared_stream=nuts_shared_stream)
            s = {key: np.asarray(val) for key, val in mcmc.get_samples().items()}
            nuts_samples[(survey_n, rvir, bl)] = s
            diag = numpyro_summary(mcmc.get_samples(group_by_chain=True))
            n_div = int(np.sum(mcmc.get_extra_fields()['diverging']))
            # truth: noiseless CoM dispersion, and the rms true binary motion at the two epochs used
            v_true2 = np.column_stack([np.asarray(truth['v_true_0']), np.asarray(truth['v_true_%i' % k2])])
            row = {'survey': survey_n, 'Rvir0': rvir, 'baseline': bl, 'epoch2': k2,
                   'dt_mean_yr': cat.meta['T%iMEAN' % k2], 'dt_sd_yr': cat.meta['T%iSD' % k2],
                   'N_obs': len(cat), 'N_repeat': int(repeat.sum()), 'n_bin': int(is_bin.sum()),
                   'n_bin_repeat': int((is_bin & repeat).sum()), 'sb_cap': sb_cap,
                   'S_true': np.std(v_com), 'sigma_binary_rms': np.sqrt(np.mean((v_true2 - v_com[:, None])**2))}
            for p in ['sigma_stream', 'sigma_binary', 'mu']:
                p16, p50, p84 = np.percentile(s[p], [16, 50, 84])
                row.update({p: p50, p+'_p16': p16, p+'_p84': p84, p+'_err': 0.5*(p84 - p16),
                            p+'_mean': np.mean(s[p]), p+'_std': np.std(s[p])})
            row['rhat_max'] = max(float(np.max(d['r_hat'])) for d in diag.values())
            row['n_divergent'] = n_div
            nuts_rows.append(row)
            print("%-4s rvir0 %4g %s: sigma_stream %.2f [%.2f, %.2f] (true %.2f) | sigma_binary %.2f [%.2f, %.2f]"
                  " (cap %.2f, true rms %.2f) | r_hat %.3f, %i divergent"
                  % (survey_n, rvir, bl, row['sigma_stream'], row['sigma_stream_p16'], row['sigma_stream_p84'],
                     row['S_true'], row['sigma_binary'], row['sigma_binary_p16'], row['sigma_binary_p84'],
                     sb_cap, row['sigma_binary_rms'], row['rhat_max'], n_div))

nuts_table = Table(rows=nuts_rows)
nuts_table.meta.update({'STREAM': 'GD-1 (mock)', 'SOURCE': 'for_jake/ catalogs', 'NOBS': nuts_n_obs,
                        'NREPEAT': nuts_n_repeat, 'PRIOR': nuts_prior_scale, 'MUPRIOR': nuts_mu_prior_scale,
                        'SBCAP': nuts_sb_cap, 'LIKELIHD': 'shared' if nuts_shared_stream else 'indep',
                        'SEED': nuts_seed, 'NWARMUP': nuts_mcmc_kwargs['num_warmup'],
                        'NSAMPLE': nuts_mcmc_kwargs['num_samples'], 'NCHAIN': nuts_mcmc_kwargs['num_chains']})
nuts_table.write(table_path+"jake_nuts_"+nuts_tag+".fits", overwrite=True)

# %%
### summary plot: DESI-like, sigma_binary and sigma_stream vs rvir0, 1 yr vs 5 yr (median, 16-84%) ###

nuts_table = Table.read(table_path+"jake_nuts_"+nuts_tag+".fits", format='fits')
nuts_plot_survey = 'desi'
nuts_style = {'1yr': dict(color='#E84855', marker='o', dx=0.93), '5yr': dict(color='#2B59C3', marker='s', dx=1.07)}
tp = nuts_table[nuts_table['survey'] == nuts_plot_survey]
fig, axs = plt.subplots(1, 2, figsize=(15, 7.5))
for ax, p, truth_col, ylabel in [(axs[0], 'sigma_binary', 'sigma_binary_rms', r'$\sigma_{\rm binary}$ [km/s]'),
                                 (axs[1], 'sigma_stream', 'S_true', r'$\sigma_{v_r,\,\rm stream}$ [km/s]')]:
    for bl, st in nuts_style.items():
        tb = tp[tp['baseline'] == bl]
        x = np.asarray(tb['Rvir0'])*st['dx']
        y = np.asarray(tb[p])
        ax.errorbar(x, y, yerr=[y - tb[p+'_p16'], tb[p+'_p84'] - y], color=st['color'], marker=st['marker'],
                    ms=9, lw=2, capsize=0, ls='none',
                    label='NUTS, %s baseline' % bl.replace('yr', ' yr'))
        if p == 'sigma_binary':  # the true binary rms depends on which epochs are used
            ax.plot(x, tb[truth_col], color=st['color'], marker='x', ms=10, mew=2, ls='none',
                    label='true rms binary motion (%s)' % bl.replace('yr', ' yr'))
    if p == 'sigma_stream':
        tb = tp[tp['baseline'] == '1yr']
        ax.plot(tb['Rvir0'], tb[truth_col], color='k', marker='x', ms=10, mew=2, ls='none',
                label='truth (noiseless CoM)')
    ax.set_xscale('log')
    ax.set_xticks(rvirs)
    ax.set_xticklabels(['%g' % r for r in rvirs])
    ax.minorticks_off()
    ax.set_xlabel(r'$r_{\rm vir,0}$ [pc]')
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=13, loc='upper left' if p == 'sigma_binary' else 'upper right')
axs[0].set_title('%s-like, %i stars, %i repeats' % ({'desi': 'DESI', 'via': 'Via'}[nuts_plot_survey],
                                                     nuts_n_obs, nuts_n_repeat), fontsize=15)
plt.tight_layout()
plt.savefig(plot_path+"gd1_%s_jake_nuts_%s.pdf" % (nuts_plot_survey, nuts_tag), dpi=300, bbox_inches='tight')
# %%
### NUTS over random draws of the 35 repeat stars (and, optionally, of the 160 stars) ###
# the single-draw posteriors above are tight, but which binaries land among the 35 repeats moves them by more
#   than their width. so: the same NUTS fit (same priors, cap, likelihood; always nuts_n_obs = 160 stars of which
#   nuts_n_repeat = 35 are repeats) for nuts_n_draws random draws, per rvir0 in nuts_draw_rvirs, for
#   nuts_draw_surveys x nuts_draw_baselines only.
# nuts_draw_vary_stars = False: every draw has the same 160 BRIGHTEST stars and only the 35 repeats change. draws
#   come from the same per-rvir0 rng as the single-draw cell, so draw 0 is that cell's choice of repeats (and with
#   jax key nuts_seed + 0 it reproduces that cell's row exactly).
# nuts_draw_vary_stars = True: each draw picks its 160 stars at random (uniformly, without replacement) from the
#   FULL catalog of the survey (jake_nuts_pool: every star in Jake's catalog with finite errors, i.e. every
#   star matched to a Bonaca & Price-Whelan catalog row that has an RV error for the survey: desi 679, the
#   Jarvis-error stars; via 1577 of the 1578), then 35 of those 160 at random as the repeats.
# draw d uses jax key nuts_seed + d.
# pooled = all draws' posterior samples concatenated with equal weight, i.e. the posterior marginalized over the
#   draws (of the repeats, and of the 160 stars if they vary).
# outputs: data/gmm_tables/jake_nuts_<nuts_draw_tag>.fits (one row per survey x rvir0 x baseline x draw) and
#   ..._pooled.fits (one row per survey x rvir0 x baseline), nuts_draw_tag = <nuts_tag>_ndraw<n>[_varystars].
#   needs the single-draw cell's definitions (run it with its loop commented out); ~1.5 min per 10 fits.
nuts_n_draws = 50
nuts_draw_rvirs = [3] #<-- rvir0 = 3 pc only: its one-component width best matches the GD-1 thin stream
nuts_draw_surveys = ['desi']
nuts_draw_baselines = {'5yr': 2}
nuts_draw_vary_stars = True #<-- True: redraw the 160 stars too, from the full catalog; False: always the 160 brightest
nuts_draw_tag = nuts_tag + "_ndraw%i" % nuts_n_draws + ("_varystars" if nuts_draw_vary_stars else "")
nuts_samples_path = repo_path+"/data/nuts_samples/" #<-- raw samples, one .npz per survey x rvir0 x baseline (gitignored)

def jake_nuts_pool(rvir, survey):
    """
    every star of Jake's catalog for (rvir, survey) with finite errors and velocities (the pool the 160 are drawn
    from when nuts_draw_vary_stars). returns (cat, truth), rows in catalog order.
    """
    fname = "rvir0%g_%s" % (rvir, survey)
    cat = Table.read(jake_path+fname+".fits")
    truth = Table.read(jake_truth_path+fname+"_truth.fits")
    if not np.array_equal(np.asarray(cat['phi1']), np.asarray(truth['phi1'])):
        raise ValueError(fname+": catalog and truth rows don't line up")
    v3 = np.column_stack([np.asarray(cat['vr_%i' % k]) for k in range(cat.meta['NEPOCH'])])
    ok = np.isfinite(np.asarray(cat['vr_0_err'])) & np.all(np.isfinite(v3), axis=1)
    return cat[ok], truth[ok]

def nuts_fit_row(cat, truth, repeat, k2, key_seed):
    """one NUTS fit (nuts_model) with the repeats' second epoch = vr_<k2>. returns (row dict, samples dict)."""
    rverr = np.asarray(cat['vr_0_err'])
    is_bin = np.asarray(truth['is_binary'])
    v_com = np.asarray(truth['v_com'])
    v_epochs = np.column_stack([np.asarray(cat['vr_0']), np.asarray(cat['vr_%i' % k2])])
    v, v_err = epoch_arrays(v_epochs, rverr, repeat)
    sb_cap = sigma_binary_cap(v_epochs, repeat, nuts_sb_cap)
    mcmc = MCMC(NUTS(nuts_model), **nuts_mcmc_kwargs)
    mcmc.run(jax.random.PRNGKey(key_seed), #<-- FITS columns are big-endian (>f8), which jax refuses
             v=jnp.asarray(v.astype(np.float64)), v_err=jnp.asarray(v_err.astype(np.float64)), sb_cap=sb_cap,
             shared_stream=nuts_shared_stream)
    s = {key: np.asarray(val) for key, val in mcmc.get_samples().items()}
    diag = numpyro_summary(mcmc.get_samples(group_by_chain=True))
    v_true2 = np.column_stack([np.asarray(truth['v_true_0']), np.asarray(truth['v_true_%i' % k2])])
    row = {'epoch2': k2, 'dt_mean_yr': cat.meta['T%iMEAN' % k2], 'dt_sd_yr': cat.meta['T%iSD' % k2],
           'N_obs': len(cat), 'N_repeat': int(repeat.sum()), 'n_bin': int(is_bin.sum()),
           'n_bin_repeat': int((is_bin & repeat).sum()), 'sb_cap': sb_cap,
           'S_true': np.std(v_com), 'sigma_binary_rms': np.sqrt(np.mean((v_true2 - v_com[:, None])**2))}
    for p in ['sigma_stream', 'sigma_binary', 'mu']:
        p16, p50, p84 = np.percentile(s[p], [16, 50, 84])
        row.update({p: p50, p+'_p16': p16, p+'_p84': p84, p+'_err': 0.5*(p84 - p16),
                    p+'_mean': np.mean(s[p]), p+'_std': np.std(s[p])})
    row['rhat_max'] = max(float(np.max(d['r_hat'])) for d in diag.values())
    row['n_divergent'] = int(np.sum(mcmc.get_extra_fields()['diverging']))
    return row, s

nuts_draw_rows, nuts_pooled_rows = [], []
for i_r, rvir in enumerate(rvirs):
    if rvir not in nuts_draw_rvirs:
        continue
    for survey_n in nuts_draw_surveys:
        rng_nuts = np.random.default_rng(seed=[nuts_seed, i_r]) #<-- as in the single-draw cell: draw 0 = its repeats
        # draws = list of (cat, truth, repeat), one per draw
        if nuts_draw_vary_stars:
            cat_pool, truth_pool = jake_nuts_pool(rvir, survey_n)
            if len(cat_pool) < nuts_n_obs:
                raise ValueError("pool of %i stars < nuts_n_obs = %i" % (len(cat_pool), nuts_n_obs))
            S_true_pool = np.std(np.asarray(truth_pool['v_com']))
            draws = []
            for _ in range(nuts_n_draws):
                sel = np.sort(rng_nuts.choice(len(cat_pool), size=nuts_n_obs, replace=False))
                r = np.zeros(nuts_n_obs, dtype=bool)
                r[rng_nuts.choice(nuts_n_obs, size=nuts_n_repeat, replace=False)] = True
                draws.append((cat_pool[sel], truth_pool[sel], r))
        else:
            cat, truth, repeat0 = jake_nuts_sample(rvir, survey_n, nuts_n_obs, nuts_n_repeat, rng_nuts)
            S_true_pool = np.std(np.asarray(truth['v_com']))
            draws = [(cat, truth, repeat0)]
            for _ in range(nuts_n_draws - 1):
                r = np.zeros(len(cat), dtype=bool)
                r[rng_nuts.choice(len(cat), size=nuts_n_repeat, replace=False)] = True
                draws.append((cat, truth, r))
        for bl, k2 in nuts_draw_baselines.items():
            draw_samples = []
            for i_d, (cat, truth, repeat) in enumerate(draws):
                row, s = nuts_fit_row(cat, truth, repeat, k2, nuts_seed + i_d)
                row = {'survey': survey_n, 'Rvir0': rvir, 'baseline': bl, 'draw': i_d, **row}
                nuts_draw_rows.append(row)
                draw_samples.append(s)
                print("%-4s rvir0 %4g %s draw %i: sigma_stream %.2f [%.2f, %.2f] (true %.2f) | sigma_binary %.2f"
                      " [%.2f, %.2f] (cap %.2f, true rms %.2f) | %2i binaries in repeats | r_hat %.3f, %i div"
                      % (survey_n, rvir, bl, i_d, row['sigma_stream'], row['sigma_stream_p16'],
                         row['sigma_stream_p84'], row['S_true'], row['sigma_binary'], row['sigma_binary_p16'],
                         row['sigma_binary_p84'], row['sb_cap'], row['sigma_binary_rms'], row['n_bin_repeat'],
                         row['rhat_max'], row['n_divergent']))
            # pooled over draws, and the spread of the per-draw medians
            rows_d = nuts_draw_rows[-nuts_n_draws:]
            # S_true, sigma_binary_rms: median over draws (constant over draws unless the stars vary);
            #   S_true_pool: the noiseless CoM dispersion of the whole pool (= S_true if the stars don't vary)
            prow = {'survey': survey_n, 'Rvir0': rvir, 'baseline': bl, 'epoch2': k2, 'n_draws': nuts_n_draws,
                    'vary_stars': nuts_draw_vary_stars, 'N_pool': len(cat_pool) if nuts_draw_vary_stars else len(cat),
                    'N_obs': len(cat), 'N_repeat': nuts_n_repeat,
                    'n_bin_min': min(r['n_bin'] for r in rows_d), 'n_bin_max': max(r['n_bin'] for r in rows_d),
                    'n_bin_repeat_min': min(r['n_bin_repeat'] for r in rows_d),
                    'n_bin_repeat_max': max(r['n_bin_repeat'] for r in rows_d),
                    'S_true': np.median([r['S_true'] for r in rows_d]), 'S_true_pool': S_true_pool,
                    'sigma_binary_rms': np.median([r['sigma_binary_rms'] for r in rows_d])}
            for p in ['sigma_stream', 'sigma_binary', 'mu']:
                pooled = np.concatenate([s[p] for s in draw_samples])
                p16, p50, p84 = np.percentile(pooled, [16, 50, 84])
                m16, m50, m84 = np.percentile([r[p] for r in rows_d], [16, 50, 84])
                prow.update({p: p50, p+'_p16': p16, p+'_p84': p84, p+'_err': 0.5*(p84 - p16),
                             p+'_draw_med': m50, p+'_draw_p16': m16, p+'_draw_p84': m84})
            prow['rhat_max'] = max(r['rhat_max'] for r in rows_d)
            prow['n_divergent'] = sum(r['n_divergent'] for r in rows_d)
            nuts_pooled_rows.append(prow)
            # the raw posterior samples, (n_draws, n_samples) per parameter (float32), plus per-draw truths/counts.
            #   data/nuts_samples/ is gitignored (~5 MB per file).
            os.makedirs(nuts_samples_path, exist_ok=True)
            np.savez_compressed(nuts_samples_path+"jake_nuts_%s_%s_rvir%g_%s.npz" % (nuts_draw_tag, survey_n, rvir, bl),
                                **{p: np.stack([s[p] for s in draw_samples]).astype(np.float32)
                                   for p in ['sigma_stream', 'sigma_binary', 'mu']},
                                **{c: np.array([r[c] for r in rows_d])
                                   for c in ['S_true', 'sigma_binary_rms', 'sb_cap', 'n_bin', 'n_bin_repeat']})
            print("  POOLED: sigma_stream %.2f [%.2f, %.2f] | sigma_binary %.2f [%.2f, %.2f]"
                  % (prow['sigma_stream'], prow['sigma_stream_p16'], prow['sigma_stream_p84'],
                     prow['sigma_binary'], prow['sigma_binary_p16'], prow['sigma_binary_p84']))

nuts_draw_meta = {'STREAM': 'GD-1 (mock)', 'SOURCE': 'for_jake/ catalogs', 'NOBS': nuts_n_obs,
                  'NREPEAT': nuts_n_repeat, 'NDRAWS': nuts_n_draws, 'PRIOR': nuts_prior_scale,
                  'MUPRIOR': nuts_mu_prior_scale, 'SBCAP': nuts_sb_cap,
                  'LIKELIHD': 'shared' if nuts_shared_stream else 'indep', 'SEED': nuts_seed,
                  'NWARMUP': nuts_mcmc_kwargs['num_warmup'], 'NSAMPLE': nuts_mcmc_kwargs['num_samples'],
                  'NCHAIN': nuts_mcmc_kwargs['num_chains']}
for rows_, sfx_ in [(nuts_draw_rows, ''), (nuts_pooled_rows, '_pooled')]:
    t_ = Table(rows=rows_)
    t_.meta.update(nuts_draw_meta)
    t_.write(table_path+"jake_nuts_"+nuts_draw_tag+sfx_+".fits", overwrite=True)

# %%
### summary plot over the draws: sigma_binary and sigma_stream vs rvir0 ###
# small faint points: each draw's posterior median and 16-84%, spread out in x to the left of the tick;
#   big square: pooled over draws (median, 16-84%); cross: the truth, at the pooled point's x.
nuts_draw_table = Table.read(table_path+"jake_nuts_"+nuts_draw_tag+".fits", format='fits')
nuts_pooled_table = Table.read(table_path+"jake_nuts_"+nuts_draw_tag+"_pooled.fits", format='fits')
nuts_draw_plot_survey, nuts_draw_plot_baseline = 'desi', '5yr'
nuts_draw_color = '#2B59C3'
td = nuts_draw_table[(nuts_draw_table['survey'] == nuts_draw_plot_survey)
                     & (nuts_draw_table['baseline'] == nuts_draw_plot_baseline)]
tpool = nuts_pooled_table[(nuts_pooled_table['survey'] == nuts_draw_plot_survey)
                          & (nuts_pooled_table['baseline'] == nuts_draw_plot_baseline)]
draw_dlogx = np.linspace(-0.11, -0.03, nuts_n_draws) #<-- log10 x offsets of the draws
pool_dlogx = 0.04
fig, axs = plt.subplots(1, 2, figsize=(15, 7.5))
for ax, p, truth_col, ylabel, truth_label in [
        (axs[0], 'sigma_binary', 'sigma_binary_rms', r'$\sigma_{\rm binary}$ [km/s]', 'true rms binary motion'),
        (axs[1], 'sigma_stream', 'S_true', r'$\sigma_{v_r,\,\rm stream}$ [km/s]', 'truth (noiseless CoM)')]:
    for i_r, rvir in enumerate(np.unique(td['Rvir0'])):
        tr = td[td['Rvir0'] == rvir]
        tr = tr[np.argsort(tr['draw'])]
        x = rvir*10**draw_dlogx[:len(tr)]
        y = np.asarray(tr[p])
        ax.errorbar(x, y, yerr=[y - tr[p+'_p16'], tr[p+'_p84'] - y], color=nuts_draw_color,
                    alpha=0.45 if nuts_n_draws <= 10 else 0.25, marker='o', ms=5, lw=1.2, capsize=0, ls='none',
                    label=('single draw of %s' % ('stars + repeats' if nuts_draw_vary_stars else 'repeats'))
                    if i_r == 0 else None)
    xp_ = np.asarray(tpool['Rvir0'])*10**pool_dlogx
    y = np.asarray(tpool[p])
    ax.errorbar(xp_, y, yerr=[y - tpool[p+'_p16'], tpool[p+'_p84'] - y], color=nuts_draw_color, marker='s',
                ms=10, lw=2.5, capsize=0, ls='none', label='pooled over %i draws' % nuts_n_draws)
    ax.plot(xp_, tpool[truth_col], color='k', marker='x', ms=11, mew=2.5, ls='none', label=truth_label)
    ax.set_xscale('log')
    ax.set_xticks(rvirs)
    ax.set_xticklabels(['%g' % r for r in rvirs])
    ax.minorticks_off()
    ax.set_xlabel(r'$r_{\rm vir,0}$ [pc]')
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=13, loc='upper right')
axs[0].set_title('%s-like, %i stars, %i repeats (%s baseline)'
                 % ({'desi': 'DESI', 'via': 'Via'}[nuts_draw_plot_survey], nuts_n_obs, nuts_n_repeat,
                    nuts_draw_plot_baseline.replace('yr', ' yr')), fontsize=15)
plt.tight_layout()
plt.savefig(plot_path+"gd1_%s_jake_nuts_%s_%s.pdf" % (nuts_draw_plot_survey, nuts_draw_tag, nuts_draw_plot_baseline),
            dpi=300, bbox_inches='tight')
# %%
