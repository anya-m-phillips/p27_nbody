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
def multiepoch_nll_per_star(v, v_err, mu, sigma_stream, sigma_binary, shared_stream=True):
    """
    -ln L for each star's epochs, (N,).
    v, v_err : (N, n_epochs) measured v_gsr and 1 sd errors [km/s]; a missing epoch is any finite v
               with v_err = missing_epoch_err (it then contributes ~nothing).

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
    v, v_err = np.atleast_2d(v), np.atleast_2d(v_err)
    s2 = sigma_stream**2
    d = v_err**2 + sigma_binary**2
    if not shared_stream:
        d, s2 = d + s2, 0.0
    r = v - mu
    A = np.sum(1.0/d, axis=1)
    B = np.sum(r/d, axis=1)
    denom = 1.0 + s2*A
    logdet = np.sum(np.log(d), axis=1) + np.log(denom)
    chi2 = np.sum(r**2/d, axis=1) - s2*B**2/denom
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
                'observational\nbinary treatment',
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
t = Table.read('for_jake/rvir03_via.fits', format='fits')
# %%
# fig, ax = plt.subplots(figsize=[8,3])
# ax.scatter(t['phi1'], t['vr_2'])