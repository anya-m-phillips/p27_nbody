#--------------------------------------------------------------#
#  single-component dispersions of mock-observed GD-1 samples,  #
#  and how well repeat RVs + sigma_binary undo binary motion.   #
#   (moved out of the last cell of gmm.py, 15 Nov)              #
#--------------------------------------------------------------#
# Assume you observe every star that passes the nominal DESI selection (gmm.py's desi case: the
# N_jarvis = 679 brightest in catalog G, rverr < 10, good Gaia PMs; just count binaries as a single
# star here), then ask if I randomly give a fraction f_repeat of them repeat RVs and, if it is a binary,
# replace the binary motions with the COM motion, what dispersion would you get out?
#    f_repeat = 0.2 is roughly another GD-1 measurement that allegedly accounts for binaries
#    (35 repeat observations of 160 unique sources; Nibauer+26)
#   **in this case I only care about fitting the single component model!!!**
#   (until 15 Nov this used the 160 brightest and scanned N_repeat = 35, 60, ..., 160 -- table
#    desi_noise_brightest160_repeatscan.fits, no longer written.)
#
# the scan: f_repeat = 0.2, 0.3, ..., 1.0. for each random draw the repeat stars are the first
#   round(f_repeat * N_obs) of ONE random permutation of the sample, so each step really is "add 10%
#   more stars to the repeat sample" (nested), and the step-to-step differences aren't drowned in
#   draw-to-draw scatter.
#
# outputs:
#   data/gmm_tables/desi_noise_nominal_repeatfracscan.fits         one row per (rvir, f_repeat)
#   plots/gd1_desi_nominal_repeatfracscan.pdf                      sigma_v_gsr vs f_repeat
#   plots/jarvis_comparison_single_component[_with_errs].pdf       (moved here from results.py, 15 Nov)
#
# the cells below the experiment only read tables, so once the scan table exists they can be
# re-run on their own (run the imports/settings cell first).
# %%
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
sys.path.append(script_path)

import pickle
import numpy as np
from scipy.optimize import minimize
import astropy.units as u
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.style.use(script_path+'/vedant.mplstyle')

import PETAR_ANALYSIS_FUNCTIONS as paf
from gmm import (measurement_errors, binary_excess_array,
                 pack_params, unpack_params, nll_flat, N_jarvis)

# %%
datapath = '/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
table_path = repo_path+"/data/gmm_tables/"
plot_path = repo_path+"/plots/"

rvirs = [0.75, 1.5, 3, 6]
f_repeat_nominal = 0.2 # ~ Nibauer+26 (35/160): fraction of sources with repeat RVs -> binary motion removed
f_repeats = np.round(np.arange(f_repeat_nominal, 1.0 + 1e-9, 0.1), 2) # 0.2, 0.3, ..., 1.0

n_draws = 200        # random choices of the repeat stars, to see how much the choice matters
sb_method = 'rms_dv' # sigma_binary estimator, see estimate_sigma_binary
scan_name = "desi_noise_nominal_repeatfracscan"

# one diagonal gaussian is separable, so fitting v_phi1/v_phi2 as well wouldn't change these two
exp_keys = ['phi2', 'v_gsr']
i_phi2, i_vgsr = exp_keys.index('phi2'), exp_keys.index('v_gsr')

# %%
def fit_one_component(x, xe, sb=None):
    """
    single gaussian, Powell on nll_flat (min_components=1), as in gmm.py's BIC_one fit.
    x  : (N, K) data; xe : (N, K) 1 sd errors or None
    sb : optional sigma_binary array, broadcastable to x's shape (see gmm.binary_excess_array)
    returns mu, sigma (K-vectors; sigma intrinsic, i.e. deconvolved with xe and sb)
    """
    theta0 = pack_params(np.array([]), np.zeros((1, x.shape[1])), x.std(axis=0)[None, :])
    res = minimize(nll_flat, x0=theta0, args=(x, 1, xe, sb), method='Powell',
                    options={'maxiter': 100000, 'maxfev': 100000})
    _, mu, sig = unpack_params(res.x, K=x.shape[1], min_components=1)
    return mu[0], sig[0]


def estimate_sigma_binary(v_pri, v_com, sel, method=sb_method):
    """
    v_pri : observed v_gsr WITH binary motion; v_com : the same with every binary at its CoM v.
    both carry the SAME saved DESI noise draw, and singles have identical values in both,
    so dv = v_pri - v_com is exactly the binary orbital motion (0 for singles).

    method = 'rms_dv' (default): sqrt(mean(dv^2)) over `sel`, the rms orbital-motion offset.
                always >= 0, and 0 only if `sel` has no binaries. no cross term.
             'var_diff': sqrt(var(v_pri) - var(v_com)) over `sel`, clipped at 0. the
                difference is var(dv) + 2 cov(v_CoM + noise, dv), and the covariance can
                drive it negative when `sel` holds only a few, slow binaries.
    """
    if method == 'rms_dv':
        dv = v_pri[sel] - v_com[sel]
        return np.sqrt(np.mean(dv**2))
    if method == 'var_diff':
        return np.sqrt(max(0.0, np.var(v_pri[sel]) - np.var(v_com[sel])))
    raise ValueError("method must be 'rms_dv' or 'var_diff', got %r" % (method,))


def load_mock_sample(rvir):
    """
    GD-1, DESI noise, binaries case (primary RVs): every system gmm.py's desi case would fit
    (include_binaries=True). returns a dict of everything the experiment needs, all length N_obs.
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

    ### the nominal DESI selection, exactly as gmm.py's `use` for noise='desi':
    #   cut & matched & good_pm & (N_jarvis brightest in catalog mG, NaNs last) & rverr < 10 & nonrem.
    #   matched_p is inside cut_for_catalog_photometry_primaries (unbound & trim_new_primaries
    #   & nonrem & |v_gsr| < 30), so the cut/nonrem terms are already implied.
    good_pm = nd['pm_err_gaia']*(u.microarcsecond/u.yr).to(u.mas/u.yr) < 0.5
    top_N_jarvis = np.zeros(N_all, dtype=bool)
    top_N_jarvis[np.argsort(mG)[:N_jarvis]] = True #<-- argsort puts NaNs last
    good_RV = top_N_jarvis & (nd['rverr_desi'] < 10.)
    use = cut_p & matched_p & good_pm & good_RV & dd['nonrem']
    idx = np.flatnonzero(use)
    idx = idx[np.argsort(mG[idx])] #<-- brightest first (order doesn't matter to the fits)
    x_err = measurement_errors(nd, 'desi', distances, exp_keys)[idx]
    if not np.all(np.isfinite(x_err)):
        raise ValueError("rvir %.2f: non-finite DESI errors in the nominal DESI sample" % rvir)

    return {
        'mG': mG[idx],
        'is_binary': is_binary[idx],
        'x_err': x_err,
        'phi2_true': sc_com['phi2'][idx],
        'phi2_obs': sc_pri['phi2'][idx] + nd['phi2'][idx],
        'v_com': sc_com['v_gsr'][idx],    # true, no noise
        'v_pri': sc_pri['v_gsr'][idx],    # true, no noise, with binary motion
        'v_noise': nd['v_gsr_desi'][idx], # the saved DESI noise draw
    }

# %%
rng_exp = np.random.default_rng(seed=1)
scan_rows = [] #<-- one dict per (rvir, f_repeat)

for rvir in rvirs:
    s = load_mock_sample(rvir)
    bin_s, x_err = s['is_binary'], s['x_err']
    N_obs = len(bin_s) #<-- unique sources (a binary is one source); 679 if every top-N_jarvis star passes

    def observed(v_true):
        return np.column_stack([s['phi2_obs'], v_true + s['v_noise']])

    v_com_all = np.where(bin_s, s['v_com'], s['v_pri']) #<-- == v_com (singles: pri == com); kept explicit
    v_obs_pri = observed(s['v_pri'])[:, i_vgsr]
    v_obs_com_all = observed(v_com_all)[:, i_vgsr]

    ### references: true CoM RVs with no noise (truth for these stars); every RV with
    #   binary motion; every binary corrected (the noise-limited best case).
    _, sig_true = fit_one_component(np.column_stack([s['phi2_true'], s['v_com']]), None)
    _, sig_none = fit_one_component(observed(s['v_pri']), x_err)
    _, sig_all = fit_one_component(observed(v_com_all), x_err)
    sb_allN = estimate_sigma_binary(v_obs_pri, v_obs_com_all, np.ones(N_obs, dtype=bool))

    ### one permutation per draw, shared by every f_repeat -> nested repeat samples
    perms = np.array([rng_exp.permutation(N_obs) for _ in range(n_draws)])

    print("\nGD-1 rvir0 = %.2f pc: %i in the nominal DESI sample (mG %.1f-%.1f), %i binaries (f_bin %.2f), median DESI rverr %.2f km/s"
            % (rvir, N_obs, s['mG'].min(), s['mG'].max(), bin_s.sum(), bin_s.mean(), np.median(x_err[:, i_vgsr])))
    print("  sigma_phi2 = %.3f deg  (noiseless CoM: %.3f)" % (sig_none[i_phi2], sig_true[i_phi2]))
    print("  sigma_v_gsr [km/s]: truth (noiseless CoM) %.2f | DESI none corrected %.2f | DESI all corrected %.2f"
            % (sig_true[i_vgsr], sig_none[i_vgsr], sig_all[i_vgsr]))
    print("  %8s %8s %8s | %-22s | %-22s | %-22s" % ('f_repeat', 'N_repeat', 'n_bin', 'repeats only', 'sigma_binary', 'repeats + sigma_binary'))

    for f_repeat in f_repeats:
        N_repeat = int(round(f_repeat * N_obs))
        ### the experiment: N_repeat random sources; the binaries among them get CoM RVs.
        #   then (a) fit as before, and (b) estimate sigma_binary from the N_repeat and
        #   add it in quadrature (v_gsr only) for the stars WITHOUT repeats, whose RVs
        #   still include binary motion. the repeat stars get 0 (their CoM v is known).
        sig_vgsr_draws, sig_vgsr_sb_draws, sb_draws, n_bin_fixed = [], [], [], []
        for perm in perms:
            repeat = np.zeros(N_obs, dtype=bool)
            repeat[perm[:N_repeat]] = True
            fix = repeat & bin_s
            x_obs = observed(np.where(fix, s['v_com'], s['v_pri']))

            _, sig = fit_one_component(x_obs, x_err)

            sigma_binary = estimate_sigma_binary(v_obs_pri, v_obs_com_all, repeat)
            sb_arr = binary_excess_array(sigma_binary, exp_keys, N_obs, apply_to=~repeat)
            _, sig_sb = fit_one_component(x_obs, x_err, sb_arr)

            sig_vgsr_draws.append(sig[i_vgsr])
            sig_vgsr_sb_draws.append(sig_sb[i_vgsr])
            sb_draws.append(sigma_binary)
            n_bin_fixed.append(fix.sum())

        sig_vgsr_draws = np.array(sig_vgsr_draws)
        sig_vgsr_sb_draws = np.array(sig_vgsr_sb_draws)
        sb_draws = np.array(sb_draws)
        p16, p50, p84 = np.percentile(sig_vgsr_draws, [16, 50, 84])
        p16_sb, p50_sb, p84_sb = np.percentile(sig_vgsr_sb_draws, [16, 50, 84])
        sb16, sb50, sb84 = np.percentile(sb_draws, [16, 50, 84])

        print("  %8.1f %8i %8.1f | %5.2f [%5.2f-%5.2f]    | %5.2f [%5.2f-%5.2f]    | %5.2f [%5.2f-%5.2f]  (/all %.2f)"
                % (f_repeat, N_repeat, np.mean(n_bin_fixed), p50, p16, p84, sb50, sb16, sb84,
                   p50_sb, p16_sb, p84_sb, p50_sb/sig_all[i_vgsr]))

        ### single-component sigmas are INTRINSIC (deconvolved with x_err), except S_true (noiseless).
        #   S_* columns are 2-vectors in exp_keys order: (phi2 [deg], v_gsr [km/s]).
        scan_rows.append({
            'sim': 'gd1', 'orbit': 'gd1', 'Rvir0': rvir,
            'N_obs': N_obs, 'f_repeat': f_repeat, 'N_repeat': N_repeat, 'n_draws': n_draws,
            'f_bin': bin_s.mean(),                       # binary fraction among the N_obs
            'n_bin_fixed_mean': np.mean(n_bin_fixed),    # binaries landing in the repeat sample, on average
            'S_true': sig_true,                          # noiseless, CoM RVs
            'S_none': sig_none,                          # DESI, no binaries corrected
            'S_all': sig_all,                            # DESI, every binary corrected
            'S_phi2': sig_none[i_phi2],                  # same in every DESI case (only v_gsr changes)
            'S_vgsr_repeat': p50,                        # DESI, N_repeat random corrected: median over draws
            'S_vgsr_repeat_p16': p16,
            'S_vgsr_repeat_p84': p84,
            'sigma_binary': sb50,                        # sb_method estimator over the N_repeat: median over draws
            'sigma_binary_p16': sb16,
            'sigma_binary_p84': sb84,
            'sigma_binary_allN': sb_allN,                # same estimator over all N_obs (reference)
            'S_vgsr_repeat_sb': p50_sb,                  # as S_vgsr_repeat, + sigma_binary in quadrature for the non-repeat stars
            'S_vgsr_repeat_sb_p16': p16_sb,
            'S_vgsr_repeat_sb_p84': p84_sb,
        })

t_scan = Table(rows=scan_rows)
t_scan.meta['KEYS'] = ','.join(exp_keys) #<-- order of the S_true/S_none/S_all columns
t_scan.meta['SBMETHOD'] = sb_method
t_scan.meta['NESTED'] = True #<-- repeat samples nested across f_repeat within a draw
t_scan.meta['SAMPLE'] = 'desi nominal (gmm.py desi, binaries case)'
t_scan.write(table_path+scan_name+".fits", overwrite=True)
print("wrote", table_path+scan_name+".fits")

# %%
### sigma_v_gsr vs f_repeat, one square panel per rvir0. solid + band = repeats + sigma_binary
#   (median, 16-84% over draws); dashed = repeats only; horizontal lines = truth (noiseless CoM),
#   DESI with no binaries corrected, DESI with all of them corrected.
t_scan = Table.read(table_path+scan_name+".fits", format='fits')

fig, axs = plt.subplots(2, 2, figsize=(14, 14), sharex=True)
for ax, rvir in zip(axs.flat, rvirs):
    t = t_scan[t_scan['Rvir0']==rvir]
    t = t[np.argsort(t['f_repeat'])]
    ax.fill_between(t['f_repeat'], t['S_vgsr_repeat_sb_p16'], t['S_vgsr_repeat_sb_p84'],
                    color='C0', alpha=0.25, lw=0)
    ax.plot(t['f_repeat'], t['S_vgsr_repeat_sb'], marker='o', color='C0', label='repeats + $\\sigma_{\\rm binary}$')
    ax.plot(t['f_repeat'], t['S_vgsr_repeat'], marker='s', color='C1', ls='--', label='repeats only')
    ax.axhline(t['S_true'][0][i_vgsr], color='k', lw=1.5, label='truth (noiseless CoM)')
    ax.axhline(t['S_all'][0][i_vgsr], color='k', lw=1, ls=':', label='DESI, all binaries corrected')
    ax.axhline(t['S_none'][0][i_vgsr], color='0.5', lw=1, ls='-.', label='DESI, none corrected')
    ax.set_title("$r_{\\rm vir,0}$ = %.2f pc ($N$ = %i)" % (rvir, t['N_obs'][0]), fontsize=18)
    ax.set_ylim(0, 1.35*t['S_none'][0][i_vgsr]) #<-- headroom for the legend
for ax in axs[-1]:
    ax.set_xlabel("$f_{\\rm repeat}$")
for ax in axs[:, 0]:
    ax.set_ylabel("single-component $\\sigma_{v_{\\rm gsr}}$ [km/s]")
axs[0, 0].legend(fontsize=12, loc='upper right')
# fig.savefig(plot_path+"gd1_desi_nominal_repeatfracscan.pdf")
# %%
t_scan = Table.read(table_path+scan_name+".fits", format='fits')
t_scan

# %%
#-------------------------------------------------------#
#   A FINAL FIG THAT JUST HAS THE JARVIS BOX FOR THE THIN STREAM AND THE SINGLE COMPONENT FITS WITH DESI NOISE
#   (moved here from results.py, 15 Nov). GD-1, desi_noise_binaries table, S_1comp = one gaussian fit
#   to the DESI-selected stars (the 679 brightest), with binary orbital motion in the RVs. joined in
#   rvir0 order, marker size = rvir0.
show_CoM = True #<-- also draw desi_noise_CoM's S_1comp (same selection idea, CoM RVs, i.e. no binary motion)
add_errbars = False
show_noiseless_thin = True #<-- also draw the THIN component (S_ts) of the noiseless two-component fit, for comparison
noiseless_thin_case = "noiseless_CoM" #<-- or "noiseless_binaries" (thin sigma_vgsr ~0.05 km/s higher)
f_repeat_plot = 0.7#f_repeat_nominal #<-- which f_repeat of the repeat scan to draw as "analogous to Nibauer+26"

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

t_1c = Table.read(table_path+"desi_noise_binaries.fits", format='fits')
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

# the repeat-RV experiment (above): DESI noise, the nominal DESI sample, a fraction f_repeat_plot of
#   them at random with binary motion removed, + sigma_binary (estimated from those) in quadrature for
#   the rest; one gaussian. y = median over the random choices of the repeats, error bar = 16-84%.
t_rep = Table.read(table_path+scan_name+".fits", format='fits')
t_rep = t_rep[np.isclose(t_rep['f_repeat'], f_repeat_plot)]
t_rep = t_rep[np.argsort(t_rep['Rvir0'])]
if add_errbars==True:
    ax.errorbar(t_rep['S_phi2'], t_rep['S_vgsr_repeat_sb'],
                yerr=[t_rep['S_vgsr_repeat_sb'] - t_rep['S_vgsr_repeat_sb_p16'],
                t_rep['S_vgsr_repeat_sb_p84'] - t_rep['S_vgsr_repeat_sb']],
                fmt='none', ecolor=gd1_color, elinewidth=1.5, capsize=3, zorder=3)
ax.plot(t_rep['S_phi2'], t_rep['S_vgsr_repeat_sb'], color=gd1_color, ls='-', lw=1, zorder=3)
ax.scatter(t_rep['S_phi2'], t_rep['S_vgsr_repeat_sb'], marker='o',
           s=size_from_rvir(t_rep['Rvir0']),
           color=gd1_color, edgecolor='k', lw=1.5, zorder=4)#,
        #    label="Analogous measurement to Nibauer+26")#'DESI-like noise, %i brightest,\n%i repeats + $\\sigma_{\\rm binary}$;\nsingle component'
                #  % (t_rep['N_obs'][0], t_rep['N_repeat'][0]))

if show_CoM:
    t_1c_com = Table.read(table_path+"desi_noise_CoM.fits", format='fits')
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
annotate_leftmost(t_rep['S_phi2'], t_rep['S_vgsr_repeat_sb'],
                  'analogous binary\naccounting\nto Nibauer+26',
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
# plt.savefig(fname+".pdf",
#             dpi=300, bbox_inches='tight')
# %%
