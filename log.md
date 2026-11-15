# log

dated record of changes to the repository and of results / checks as they were made. newest first. the readme describes only the current state; history, superseded results and "used to be" notes live here.

entries before 15 Nov 2026 were moved out of the readme when the log was started, so they're only as complete as the readme was.

---

## 2026-11-15
- **two more sb_modes in `single_component_comparisons_mockobs.py`: `'avg'` and `'flag'`** (both need `multi_epoch = True`, neither has a sigma_binary). `'avg'`: repeat stars' velocity = the mean of their two epochs, error rverr/sqrt 2 (`average_epochs`). `'flag'`: repeat stars flagged RV variable by `paf.get_detections` (p < 0.1, hardcoded there; used at the user's request instead of a new chi^2 test, so the threshold isn't a setting) are removed (`flag_rv_variable`), with `flag_average` choosing two-epoch mean or epoch 0 for the unflagged repeats. tags `avg`, `flag_avg`, `flag_ep0`; extra columns `N_fit_mean`, `n_flagged_mean`, `n_flagged_bin_mean`, `n_flagged_single_mean`; headers `FLAGTEST`, `FLAGAVG`. the scan cell had been commented out by the user; it's restored as live code behind a new `run_scan` switch (False skips it). ran all three. results in the readme: `avg` = `fit_shared`'s `S_vgsr_nosb` exactly (as it should); `flag_avg` at f_repeat = 1 gives 1.02 / 0.88 / 1.32 / 1.09 km/s (truths 0.84 / 0.67 / 0.52 / 0.39), the first method that doesn't collapse to 0, but only 18-24% of the binaries get flagged; `flag_ep0` is lower but biased low by the selection (epoch-0 noise correlates with the epoch difference the flag cuts on).
- **`sb_dv_method = 'noiseless'`** for the one-epoch `'dv'` case: `estimate_sigma_binary_com(..., v_true=None)` with `method='noiseless'` returns std(`v_true - v_com`) over the repeats, the true binary motion with no noise; the scan passes `v_true = s['v_pri']`. only with `multi_epoch = False` (the settings cell raises otherwise). tag `dvcom_noiseless`, table `desi_noise_nominal_multiepoch_scan_dvcom_noiseless.fits`. check: at f_repeat = 1 sigma_binary equals the true rms and S_vgsr the oracle, exactly. numbers in the readme.
- **`sb_mode = 'dv'` now works with `multi_epoch = False`** (`single_component_comparisons_mockobs.py`; it used to stop at `assert multi_epoch`). one observed epoch = `v_pri + v_noise` (sim primary v_gsr + the saved DESI draw); dv for the repeats = observed epoch - true noiseless CoM v, via the new `estimate_sigma_binary_com` (`'deconvolved'` = sqrt(var(dv) - <rverr^2>), `'raw'` = std(dv)). the scan's per-rvir setup builds `v_obs`, `rverr`, `v_com_obs0`, `dv_bin` from either source; `sigma_binary_rms_diff` is NaN with one epoch; tag `dvcom_<method>`; header `MULTIEP`; `'fit'` + `multi_epoch = False` raises. ran `dvcom_deconvolved` (the file's settings), `dvcom_raw`, and `dv_deconvolved` again (identical to before). results in the readme: still ~0 except at 1.5 pc. the deconvolved estimate is 0.3-0.7 km/s above the true rms at f_repeat = 1; checked that dv = b + noise exactly and that the excess is the saved noise draw's realised variance (0.8-4.4 km^2/s^2 above <rverr^2>) plus the b x noise cross term. the user's settings at this point: `multi_epoch = False`, `sb_mode = 'dv'`, `shared_stream = True`.
- **sigma_binary measurement made a toggle** in `single_component_comparisons_mockobs.py`: `sb_mode` = `'fit'` (free parameter, uses `shared_stream`) or `'dv'` (estimated from the repeats' dv, uses `sb_dv_method`, which replaces the `sb_dv_methods` list: one estimator per run). the scan only runs the chosen measurement and writes it to generic columns `S_vgsr[_p16,_p84]` / `sigma_binary[_p16,_p84]` (plus `S_vgsr_nosb*` for `'fit'`; oracle columns in every mode), to `desi_noise_nominal_multiepoch_scan_<tag>.fits`, `<tag>` = `fit_shared` / `fit_indep` / `dv_deconvolved` / `dv_raw`; headers `SBMODE` + `LIKELIHD` or `SBDVMETH`. `S_all` is now mode-dependent (two epochs for `'fit'`, epoch 0 for `'dv'`). the plot cells and the Jarvis figure (`rep_series = 'S_vgsr'`) read the current settings' table; plot filename gets the tag. all four combinations run (settings swapped in memory by a scratch runner; the file stays at `'fit'`, `True`, `'deconvolved'`) and reproduce the earlier numbers exactly. the two earlier tables (`desi_noise_nominal_multiepoch_scan.fits`, `_indep.fits`) are superseded and left on disk.
- **"old way" with sigma_binary estimated from the repeats' epoch-to-epoch dv** added to the scan (`single_component_comparisons_mockobs.py`): `estimate_sigma_binary_dv` (`'deconvolved'` = sqrt((var(dv) - <2 rverr^2>)/2), default; `'raw'` = std(dv), as first proposed: kept as an option, but it double counts the noise), `fit_one_component` takes `sb` again, `binary_excess_array` imported again, `sb_dv_methods` setting, an oracle fit with the true rms sigma_binary, new columns `sigma_binary_dv*`, `S_vgsr_sbdv*`, `sigma_binary_dvraw*`, `S_vgsr_sbdvraw*`, `sigma_binary_oracle`, `S_vgsr_sboracle` (header `SBDV`). plots: the panel figure and the sigma_binary figure show it; the Jarvis figure draws `rep_series` (default `'S_vgsr_sbdv'`). `shared_stream` set back to `True`; the scan rewrote `desi_noise_nominal_multiepoch_scan.fits` (its multi-epoch columns are unchanged: same seed, same draws). the file's extra empty `# %%` cells were the user's.
  - result: the deconvolved sigma_binary is about right, but the dispersion still comes out ~0; even the oracle gives 0.00 at 0.75 / 6 pc and 2.4-2.6 km/s at 1.5 / 3 pc (truths 0.39-0.84). the stream variance is smaller than the random cross terms between noise and binary motion in a 679-star sample, so subtracting binaries in quadrature can't recover it. numbers in the readme.
  - context from the user: Nibauer+26 add sigma_binary to every star, with no systemic-velocity parameter and no binary population model, and put a half-normal (scale 10 km/s) prior on it, sampled with NUTS; a 30 km/s scale doesn't change their results.
- **independent-epoch (Nibauer+26 as written) likelihood run, on request.** new module-level `shared_stream` switch in `single_component_comparisons_mockobs.py`, passed to every `fit_multiepoch` call in the scan; `False` appends `_indep` to `scan_name` and the `LIKELIHD` header says which. set to `False` and run: `desi_noise_nominal_multiepoch_scan_indep.fits`. result: degenerate, as expected. on the 3 pc sample (f_repeat 0.5) the nll is identical to 1e-6 along sigma_stream^2 + sigma_binary^2 = 5.67^2, from all-binary to all-stream. the scan's sigma_stream / sigma_binary are arbitrary splits (e.g. 0.75 pc: sigma_binary ~0, sigma_stream 1.9-2.6; 1.5 pc: ~3.3 / ~3.5), and their quadrature sum reproduces `S_vgsr_nosb` (1.5 pc, f_repeat 1: hypot(2.30, 4.05) = 4.66). the shared-term table `desi_noise_nominal_multiepoch_scan.fits` from the earlier run is untouched.
- **repeat-RV experiment rebuilt around two-epoch mock RVs and a fitted sigma_binary** (`single_component_comparisons_mockobs.py`). details are in the readme section; what changed:
  - `load_mock_sample(rvir, multi_epoch=True, rng, dt_mean_yr=1, dt_sd_yr=2/12)`: binaries get Keplerian RV curves about their noiseless CoM v_gsr (`paf.get_orbital_params_from_binaries` + `paf.get_rvs`); two epochs, Delta t ~ N(1 yr, 2 months); fresh DESI noise per epoch from `rverr_desi`. `multi_epoch=False` gives the old return dict.
  - new likelihood in the script (not `gmm.py`): `multiepoch_nll_per_star` / `nll_multiepoch` / `fit_multiepoch` / `epoch_arrays`. per star, MVN with C = diag(err^2 + sigma_binary^2) + sigma_stream^2 11^T; missing epochs get err = 1e4 km/s. **deviation from the formula as first specified** (a product over epochs of N(v | mu, err^2 + sigma_binary^2 + sigma_stream^2)): with independent epochs the two sigmas only enter as a sum and can't be separated, so the stream term is shared by a star's epochs. the independent version is kept as `shared_stream=False`; on synthetic data it returns the same total width with an arbitrary split.
  - optimizer: Powell sticks at sigma_stream on its lower bound (nll 8 worse than the truth on synthetic data), so `fit_multiepoch` runs L-BFGS-B and Nelder-Mead and keeps the lower nll.
  - removed `estimate_sigma_binary`, the `sb` argument of `fit_one_component`, `sb_method`, and the `binary_excess_array` import. new table `desi_noise_nominal_multiepoch_scan.fits` (the old `desi_noise_nominal_repeatfracscan.fits` stays on disk, unread). the plot cells and the Jarvis comparison figure read the new table (`S_vgsr_stream` in place of `S_vgsr_repeat_sb`); the sigma_binary-vs-f_repeat cell now shows every rvir0 against the rms binary motion.
  - first run (10 draws, ~2 min): **sigma_stream = 0 at every rvir0 and f_repeat**, sigma_binary 2.2-5.0 km/s takes everything. heavy-tailed dv; same mis-specification as before, worse. a scratch per-star mixture (binary w.p. f_b) gives sigma_stream 0.6-0.95 km/s against truths 0.39-0.84 (f_b ~ 0.02-0.06, sigma_b ~ 15-21 km/s). numbers in the readme; TODO updated.
- **`paf.radial_velocity` phase fix:** `phi0` is drawn uniform in [0, 1) but was subtracted in radians (`M = 2 pi t/P - phi0`), so every binary started within 1 rad of mean anomaly 0. now `M = 2 pi (t/P - phi0)`. affects anything that used `get_rvs` / `radial_velocity` before (the old detection-fraction code).
- **`paf.intrinsic_stream_data_v3` now saves the raw petar `singles` and `binaries` objects** in its dict, so pickles rebuilt by `final_datasets.py` carry them (~1 MB per 2000 binaries; checked that a `Binary` round-trips through pickle). nothing else in the pickle changes. none of the pickles on disk have been rebuilt; until then `load_present_day_binaries` reads the present-day `data.<k>.binary` and checks it against the pickle by p1 id (copy 0 matches for all four GD-1 rvir0).
- readme: new section "what's in a data_dict (the pickles)" listing every key and subkey (from `gd1_0.75.pickle`), replacing the partial table under `final_datasets.py`; new paf subsection on the mock binary RV functions. confirmed while listing: `companions['in_rtid']` has length N, not N2 (known bug).
- found by hand in the working tree (not committed when this entry was written): `paf.get_orbital_params_from_binaries(binaries, v_sys_vals, rng)` added; `paf.get_obstimes(N_obsm, N, rng)` -> `get_obstimes(N, rng)`; in the comparison script `n_draws` 200 -> 10, a `tqdm` import and progress bar on the f_repeat loop, the per-f_repeat print commented out, a scratch cell plotting sigma_binary vs f_repeat at 3 pc. `desi_noise_nominal_repeatfracscan.fits` was rewritten, presumably with 10 draws, so it no longer matches the 200-draw numbers logged below.
- **docs reorganized:** `log.md` started; everything dated or historical moved here out of `readme.md`, which now describes only the current state. new working preference in the readme: update `readme.md` and `log.md` every time the code is touched.
- `N_repeat` for the repeat-RV experiment changed 34 -> 35 (by hand, in `gmm.py`), to match Nibauer+26 (35 repeat observations of 160 unique sources).
- **new script `single_component_comparisons_mockobs.py`.** moved into it: the repeat-RV / sigma_binary experiment (formerly the last cell of `gmm.py`, from line 1198 down; removed from `gmm.py` and replaced with a pointer comment), and the single-component Jarvis comparison figure cell (formerly in `results.py`, `plots/jarvis_comparison_single_component[_with_errs].pdf`; replaced with a pointer comment). `fit_one_component` / `estimate_sigma_binary` / `load_mock_sample` are now module-level functions; `estimate_sigma_binary(v_pri, v_com, sel, method)` takes the velocities as arguments instead of closing over them. `gd1_color` is recomputed with `results.py`'s orbital-phase ordering (`#E84855`).
  - side effect: `import gmm` no longer runs the experiment. (before this, `gmm.py`'s last cell ran on import because its `__main__` guard was commented out, so every `import gmm`, including `results.py`'s, refit it and rewrote `desi_noise_brightest160_repeat34.fits`. checked 5 Nov that the rewrite was byte-identical to the committed file.)
- **first version of the scan: the 160 brightest, N_repeat = 35, 60, 85, 110, 135, 160** (200 draws, nested repeat samples). table `data/gmm_tables/desi_noise_brightest160_repeatscan.fits`, plot `plots/gd1_desi_brightest160_repeatscan.pdf` (both superseded the same day, still on disk). ~5 min. median DESI rverr 1.24 km/s. median sigma_v_gsr [km/s], repeats + sigma_binary, [16-84%] over draws; "all" = every binary corrected:

  | rvir0 | truth (noiseless CoM) | all | none corrected | 35 | 60 | 85 | 110 | 135 |
  |---|---|---|---|---|---|---|---|---|
  | 0.75 | 0.81 | 0.43 | 1.24 | 0.59 [0.00-1.09] | 0.43 [0.00-0.94] | 0.39 [0.00-0.74] | 0.38 [0.14-0.62] | 0.41 [0.29-0.50] |
  | 1.5 | 0.65 | 0.59 | 3.50 | 1.24 [0.38-3.28] | 0.77 [0.42-2.71] | 0.64 [0.47-1.73] | 0.61 [0.50-0.92] | 0.60 [0.53-0.67] |
  | 3 | 0.43 | 0.48 | 2.75 | 1.27 [0.32-2.56] | 0.59 [0.29-1.95] | 0.52 [0.34-1.38] | 0.50 [0.38-0.70] | 0.48 [0.42-0.59] |
  | 6 | 0.44 | 0.52 | 2.95 | 0.87 [0.18-2.75] | 0.65 [0.37-1.71] | 0.57 [0.40-1.17] | 0.54 [0.43-0.69] | 0.52 [0.47-0.57] |

  - at 35 repeats the median was 0.2-0.8 km/s above "all", and the 16-84% range ran out to near the uncorrected value: whether the ~8-11 repeat binaries include the fast ones decides the answer. by 60-85 the median was within ~0.1-0.2 km/s; by 110 the 16-84% range was within a few tenths of a km/s.
  - repeats only (no sigma_binary) fell roughly linearly in N_repeat.
  - at 0.75 pc "all" (0.43) was well below the noiseless truth (0.81): the DESI noise draw on those 160 stars.
  - sigma_binary itself was high-variance at 35 (e.g. 2.76 [1.06-4.12] at 1.5 pc).
- **the scan switched to the nominal DESI sample and a repeat fraction**: every star in `gmm.py`'s `desi` selection (679 at every rvir0), f_repeat = 0.2, 0.3, ..., 1.0. table `desi_noise_nominal_repeatfracscan.fits`, plot `plots/gd1_desi_nominal_repeatfracscan.pdf`. ~7 min. results are in the readme. new finding: at 0.75 pc sigma_binary over-corrects and the f_repeat = 0.2 fit collapses to 0 (heavy-tailed dv; TODO added).
- by hand in `single_component_comparisons_mockobs.py`: `f_repeat_plot` set to 0.7, and both `savefig` calls commented out.
- `desi_noise_brightest160_repeat35.fits` was written once by the first run and deleted (nothing read it).
- readme's bug list: dropped the "importing gmm runs the repeat-RV experiment" item (fixed by the move).
- on-disk mtimes, for the record: grid pickles 4 Nov; supplemental pickles 11-29 Oct (so the 4 Nov `--which supplemental` regeneration didn't happen, and they still have nearest-point Teffs); `via_noise_binaries.fits` 4 Nov; `desi_noise_{CoM,binaries}`, `via_noise_CoM`, `noiseless_{CoM,binaries}`, `jarvis26` 19 Oct; `*_supplemental` tables 29 Oct. (git checkouts also set mtimes, so the 19 Oct ones may be older.) the readme's refit TODO and the netscratch purge dates were updated from this.
- the readme's seed TODO corrected: `desi_noise_*_supplemental` tables exist, fit on the supplemental loop's seed-42 draws (the readme had said no noisy supplemental case had been fit).

## 2026-11-06
- `results.py` populations cells: `mass_stat` toggle added, default `'mean_log'` (<log10 M_sys>, same treatment as <G>); `'log_mean'` (log10 <M_sys>) was the only option before, so the older `plots/gd1_combined_msys_fbin_meanG_<mass_sample>.pdf` files are `log_mean`.
- tried medians with 5-95% intervals instead of means in those cells. looked bad, reverted.

## 2026-11-05
- `results.py`: thin stream vs cocoon population cells (f_bin, <G>, <M_sys>, combined) added for GD-1. result: cocoon systems are more often binaries, more massive, and brighter (Weatherford+2025-like selection). numbers in the readme.
- checked: the repeat-RV table rewritten on `import gmm` was byte-identical to the committed file.

## 2026-11-04
- **catalog `log_Teff` switched to `noise.Teff_from_gaia_color_temp_relation`** (degree-2 fit of log Teff vs BP-RP to the isochrone track), from `nearest_isochrone_point` (CMD nearest point, which needed absolute G and so depended on the N-body distance each catalog star was painted onto). `nearest_isochrone_point` is now only the `iso_dist` diagnostic. `assign_photometry_from_catalog`'s `nearest_kwargs` renamed `teff_kwargs`. on fake N-body stars matched to the gd1 catalog, the new Teffs are a median 0.026 dex cooler (16-84%: 0.011-0.041): the nearest-point lookup had partly absorbed the reddening offset by snapping stars back onto the bluer isochrone. decided to live with the reddening for now.
  - with the nearest-point Teffs, catalog stars off the isochrone in G (`iso_dist > 3`) got extrapolated Teffs; the colour fit doesn't use G.
- started regenerating for the colour-Teff fit: `final_datasets.py --which grid` + `--which supplemental`, then the noisy `gmm.py` cases (see the readme TODO). grid pickles were rewritten this day. `via10hr` left out (no table on disk).
- `gd1_2xN_dense` added to `simspect.supplemental_sims` and `gd1_2xN` removed, so `--which supplemental` / `--supplemental True` now do `gd1_2xbinaries` + `gd1_2xN_dense`.
- `gmm.py`: optional fixed `sigma_binary` added to every likelihood function (`binary_excess_array` helper), and the repeat-RV experiment added as `gmm.py`'s last cell: 160 brightest DESI-like GD-1 stars, 34 random repeats, `desi_noise_brightest160_repeat34.fits`. `estimate_sigma_binary` got `'rms_dv'` (default) and `'var_diff'`; `var_diff` hit 0 in up to 27/200 draws at 0.75 pc, through its 2 cov(v_CoM + noise, dv) term. checked: none of the 160 fail the pm/RV error cuts.
- working preference added: panel aspect ratios either very square or very not square.

## 2026-10-29
- `desi_noise_*_supplemental.fits` tables made (gd1_2xN_dense in the supplemental set).

## 2026-10-21
- `particle_spray.py`: **no progenitor potential** (`progenitor_potential=None`). the first version used a fixed Plummer at the *initial* mass (as `animations/particle_spray_vs_nbody.py` still does), but with the release radius following the bound mass, that trapped late releases deep inside a too-massive Plummer. Pal 5 was the extreme case (bound mass 9940 -> 195 Msun, rtid down to 6 pc, peri 7.7 kpc, 13 pericenters): it ran >3x slower than the rest, because gala advances every released particle with one shared adaptive DOP853 step (rtol = atol = 1e-10). the other streams keep 20-33% of their mass, and their Plummer sprays still had 14-72 particles within 20 pc of the progenitor at the present day. with the Plummer each stream took 3-4 min (Pal 5 > 9 min).

## 2026-10-19
- working preference added: make code changes with Edit/Write, not shell commands, so they show up in the editor diff.

## 2026-10-15
- **the real Jarvis+26 GD-1 data put through the pipeline**: `final_datasets.py --which jarvis` -> `jarvis26.pickle`, `gmm.py --jarvis-data True` -> `jarvis26.fits`. reproduces their published decomposition (numbers in the readme). also settled that `DELTA_PHI2` is in degrees.
- `gmm.py` command-line parsing fixed: between commit 3affe0c and this day, `--constrain-widths` and `--supplemental` were parsed as `s=='False'`, i.e. inverted. any table made from the command line with those flags in that window has the opposite setting from what was typed.
- `via_noise_{CoM,binaries}` re-run with `x_err` on the seed-0 pickles, so `desi_noise_*` and `via_noise_*` were then all current (x_err, seed 0, unconstrained). `via10hr_noise_CoM` (30 Sep) was not re-run.
- checked: GD-1 0.75 `via` CoM `x_err[:, v_gsr]` == `rverr_via` exactly; for Jarvis, `x_err[:, v_gsr]` == `V_ERR` and `x_data[:, v_gsr]` == `DELTA_VGSR`.

## 2026-10-14
- **`x_err` added to the GMM likelihood** (per-star errors added in quadrature, so fitted sigmas are deconvolved). GD-1 0.75 `via` CoM, without -> with `x_err` (seed-42 draw): f_cocoon 0.136 -> 0.098, thin sigma_v_gsr 0.81 -> 0.50 (noiseless 0.46), cocoon sigma_v_gsr 1.97 -> 2.07 (noiseless 2.32); v_phi sigmas -> ~0 (unresolved by Gaia). the old noisy fits' 7-14 km/s in v_phi were the error bars.
- **noise seed in `final_datasets.py` 42 -> 0** for the grid (the supplemental loop stays at 42). the seed-42 DESI thin sigma_v_gsr = 0 (GD-1, 0.75, CoM) was a fluke of the draw: with the Jarvis errors (median 3.1 km/s, but 44 stars < 1 km/s) a 0.46 km/s width is resolvable (toy realizations: median estimate 0.45, 5% at 0). in the seed-42 draw the 44 precise stars had noise/error std 0.89, and 500 fresh draws on the same stars put only 3.4% at 0; its profile gave a 95% upper limit of ~0.41 against a true 0.47 (~1.9 sigma low). diagnosed with `claude_profile_desi.py`. after seed 0 the DESI thin dispersion was back near the N-body value.

## 2026-10-13
- DESI RV errors switched to `noise.desi_RVerr_from_jarvis` (the 679 Jarvis+26 `V_ERR`s by brightness rank), from the analytic `desi_RVerr(zmag)`.
- `ontrack_velocity` (|v_gsr| < 30 km/s, to match Jarvis+26) added to the catalog-photometry pool `cut`.

## 2026-10-11
- **supplemental simulations** `gd1_2xbinaries` and `gd1_2xN` added (loaders in `inspect_new_sims.py`, `--which supplemental` in `final_datasets.py`, `--supplemental` in `gmm.py`, `sim` column added to the gmm tables -- older tables don't have it). noiseless fits and the comparison with the grid are in the readme. the grid tables were checked against a refit with the then-current `gmm.py`: exact match.
- `gd1_2xN` checked on all five copies: degenerate CoM split every time (table in the readme). copy 1-4 pickles made.
- `gd1_2xN_dense` (rvir0 = 3 pc, same population as `gd1_2xN`) started on the cluster.

## 2026-10-05
- `final_datasets.py` re-run for all orbits with the `_primaries` catalog keys (pickles 17:29-17:33). the `_primaries` catalog pass was added after the CoM one so the existing CoM draws didn't change (the generator carries over between sims, though, so every sim after the first got a different draw).
- **trim decision:** the trim operates in the space the GMM sees: `CoM` cases use `trim_new`, `binaries` cases `trim_new_primaries`, noiseless and noisy. before this, the noisy `*_binaries` tables were fit with the CoM `trim_new` and CoM-trim catalog photometry. `trim_new_primaries` dropped the std of dv in the matched sample from 10.9 to 7.7 km/s, but didn't change the binary fraction.
- binaries selection effect analysed (`claude_binary_selection_experiment.py`); write-up in the readme's GMM section.

## 2026-10-03
- checked: `gd1` hm rvir_index=0 has all five copies finished.

## 2026-09-30
- **catalog photometry live**: `final_datasets.py` calls `noise.assign_photometry_from_catalog` against the BPW25 Gaia catalog for each stream and saves `catalog_photometry` / `noise_catalog_photometry`; every noisy `gmm.py` case uses those. the `via`, `via10hr` and `desi` tables were refit from these pickles. checked: every matched star passed every downstream noise cut in `gmm.py` (with the nearest-point Teffs).
- **BIC added to `gmm.py`**: single-gaussian fit + `BIC_one` / `BIC_two` columns (tables written before this lack them).
- noted: label switching / degenerate small-N fits in the catalog-noise tables (see the readme bug list).
- pickles written (this set aged out of scratch ~end of December; superseded by the 4 Nov set).

## 2026-09-29
- `build_isochrone_table` / `isochrone_cmd_track`: post-AGB cut on EEP (`max_eep=POST_AGB_EEP = 1409`), because artpop's packaged grid labels all of EEP 808-1710 phase 5, so `max_phase=5` didn't remove post-AGB. before the fix, N-body stars with m0_zams in 0.793574-0.79435 got post-AGB/WD photometry (absolute G ~ +10) and counted as `alive`: 2 stars per pickle (the same two m0 values in every sim, 40 of 193153 `alive` stars across the 20 pickles, 0-2 per sim after `nonrem & unbound`). photometry of every other star was bit-identical with and without the cut. (artpop table went from upper m0 edge 0.79435 to 0.793574, max logTeff 5.05 -> 3.885.)

## 2026-09-25 to 27
- `noiseless_*` tables first fit.

## 2026-09-24
- `final_datasets.py`: binaries get the **primary's** RV/phi2/distance but the **CoM** proper motions (`PM_treatment='CoM'`).

## earlier (undated, recorded in the readme before the log existed)
- **init_displacement units bug fixed**: `init_displacements.txt` is in pc but every consumer wanted kpc, which put the progenitor reference at 16.6 Mpc, where the potential is negligible, so the "orbit" free-streamed in a straight line. `extended_grid_info.__init__` now converts once. anything cached from before is wrong.
- `FINAL_ics_nolmc.csv`: the `name` column for the six run orbits renamed to the short keys (`gd1`, `pa5`, ...). this is what left `get_init_displacements.py`'s `names_to_run` stale.
- `intrinsic_stream_data_v3`: `in_rtid` moved to the core frame (where rtid is measured), as a full 3D radius, applied before trimming; `m0_zams` / `m0_effective` columns added (ID match to `data.0`).
- `straighten_stream_orbit_interp`: the multi-`yval` closures are bound with `y=y` (they used to all interpolate `pm_phi2`); the interpolant now sorts `phi1` first.
- `prog_orbit_track`: `[:-1]` added to drop the duplicated t = 0 point that stalled the chop walk.
- intrinsic-frame `inMW`/`trim` and `outlier_clip` replaced by the observed-frame percentile trim (`trim_new`).
- `gmm.py` switched from fitting `d_phi2` (kpc) to `phi2` (deg); tables from before may be in `d_phi2`.
- a three-component GMM (two narrow components for the thin stream) was tried and dropped in favour of the width-ratio constraint.
- an orbital-phase-dependent 1x/2x rtid boundary was tried and abandoned for 2x everywhere.
- `via10hr` good_RV was rverr < 10 km/s; now < 5 like `via`.
- the isochrone moved from a file downloaded from the MIST web interpolator to `artpop.fetch_mist_iso_cmd`, and from blackbody photometry to isochrone interpolation.
- the figure 1 cell moved from `results.py` to `motivation.py`; `motivation.py` stopped generating sprays itself (`particle_spray.py` does).
- `noise.py`'s `# SCRATCH::::` cells commented out, so `import noise` stopped loading a simulation.
