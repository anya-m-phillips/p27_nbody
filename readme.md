Exploring an N-body grid of mock streams to probe effects of progenitor dynamics on stream dispersion. Questions to answer might include:
- cocoon fraction and dispersion/extent as a function of progenitor density, orbit: Jarvis+26, Carlberg+26 discuss cocoon structures around MW streams, created by 
- structure of the stream velocity dispersion profile: created by the (evolving) cluster dispersion+tidal filling factor and the progenitor orbit/overlap of successive energy ``feathers." see Bovy(2014).
- **detectability**: under realistic survey noise + down-sampling (Via, DESI), when is a cocoon actually recoverable? a one-component model being preferred = a cocoon non-detection (see the `gmm.py` TODO).


**documentation:** this readme describes the project **as it is now**. what changed and when, superseded results, and "used to be" notes are in [`log.md`](log.md) (dated, newest first).


# working preferences
- **every time the code is touched, update both `readme.md` and `log.md`**: the readme so it describes the current state (no dates, no history), and the log with a dated entry for what changed and why (plus any results / checks). this includes changes the user made by hand that turn up along the way.
- **make code changes with the Edit/Write tools, not by running python/sed/shell commands that rewrite files**, so the changes show up in the editor diff view for the user to inspect and approve.
- **panel aspect ratios: either very square or very not square, never in between.** square is the default (set in `scripts/vedant.mplstyle`). go clearly non-square only when the content calls for it, e.g. a panel showing a long stream could be ~8x3. size `figsize` so each panel lands on one or the other.


# workflow

```
get_init_displacements.py   -> data/init_displacements.txt   (ICs for petar; run once, before the sims)
        |
   [petar sims on the cluster]
        |
inspect_new_sims.py  (imported as simspect)  -- load N-body data, observed frames, straightening
noise.py             (imported as noise)     -- isochrone photometry + survey error models
        |                 |
        +---> final_datasets.py  -> <scratch>/p27_data_dicts/<orbit>_<rvir>.pickle   (one per sim)
                        |
                     gmm.py      -> data/gmm_tables/<case>.fits   (+ sanity plots on scratch)
                        |
                        +---> single_component_comparisons_mockobs.py -> data/gmm_tables/<survey>_noise_nominal_multiepoch_scan_<tag>[_dt<mean>yr_sd<sd>yr].fits,
                        |         plots/gd1_<survey>_nominal_multiepoch_scan_<tag>[_dt...].pdf, plots/jarvis_comparison_single_component.pdf
                        |         (pickles + <survey>_noise_* / noiseless_CoM tables; two-epoch RVs, sigma_binary fit.
                        |          needs the petar binaries: from the pickle, or the present-day data.<k>.binary on holystore;
                        |          for_jake/rvir0<X>_<desi|via>.fits, three-epoch mock catalogs for Nibauer;
                        |          last cells: NUTS on those catalogs -> data/gmm_tables/jake_nuts_<tag>.fits, plots/gd1_desi_jake_nuts_<tag>.pdf)
                        |
                     results.py  -> plots/summary_*.pdf etc.
                        |
                   motivation.py -> plots/fig1_motiv_*.pdf   (pickles + gmm table + sprays)
                        ^
particle_spray.py (.slurm) -> data/particle_sprays/*.h5   (gitignored; from the petar ICs + data.tidal, run before motivation.py)
```

`inspect_new_sims.py` and `noise.py` are function libraries (their `__main__` blocks are exploratory leftovers); `final_datasets.py` and `gmm.py` are the programs that produce the data products. `gmm.py` and `results.py` only read the pickles, so nothing downstream of `final_datasets.py` needs the raw simulation data (except `single_component_comparisons_mockobs.py` with current pickles, which don't carry the petar binaries yet; see the data_dict section).

the **supplemental simulations** in `simspect.supplemental_sims` (currently `gd1_2xbinaries` and `gd1_2xN_dense`; see their own section) go through the same pipeline: `python final_datasets.py --which supplemental`, then `python gmm.py --noise None --include-binaries {False,True} --supplemental True`.

the **real Jarvis+26 GD-1 data** goes through the end of the pipeline too, so the sims can be compared to a fit made with the same model (see its own section): `python final_datasets.py --which jarvis`, then `python gmm.py --jarvis-data True`.


# what's in a data_dict (the pickles)
`final_datasets.py` writes one per sim to `/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/<orbit>_<rvir:.2f>.pickle` (`<sim>[_copy<k>]_<rvir>` for the supplemental sims). load with `pickle.load` inside `petar_env` (the `SkyCoord` needs astropy/gala, and `singles`/`binaries` need petar). the real-data `jarvis26.pickle` is different, only three keys: see "the real Jarvis+26 GD-1 data".

it's the dict `paf.intrinsic_stream_data_v3` returns, plus what `final_datasets.build_data_dict` adds. listed from `gd1_0.75.pickle` (nsingles = 11299, nbinaries = 1768). **lengths:**
- **N** = nsingles + nbinaries (13067 here): one row per *system*, ordered `[singles, binaries]`. every key marked N is in this same untrimmed order, so they combine with plain `&`, and row `nsingles + j` is binary `j`.
- **N2** = nsingles + 2 nbinaries (14835): `[singles, binary p1, binary p2]`.
- floats are bare (no astropy units) unless a unit is given in brackets as a Quantity.

**top level, from `intrinsic_stream_data_v3`:**

| key | what | length |
|---|---|---|
| `init_displacement` | t = 0 phase-space position of the progenitor (kpc, km/s) | 6 |
| `nsingles`, `nbinaries` | ints | |
| `IDs` | petar ids, `[singles, p1, p2]`. only lines up with `companions`; for `CoM`/`luminous` binary j's p1 id is `IDs[nsingles + j]` | N2 |
| `pot` | external potential, `[singles, p1, p2]` [pc^2/Myr^2] | N2 |
| `singles` | the raw petar `Particle` object from `data.<k>.single` (core frame, petar units: pc, pc/Myr, Msun; `.mass`, `.pos`, `.vel`, `.id`, `.star.*`, ...) | nsingles |
| `binaries` | the raw petar `Binary` object from `data.<k>.binary` (`.semi` [pc], `.ecc`, `.mass` (total), `.pos`/`.vel` (CoM, core frame), `.p1`/`.p2` (each a `Particle`), ...), rows in binary order | nbinaries |
| `CoM`, `luminous`, `companions` | one subdict per binary treatment, below | |

⚠️ **no pickle on disk has `singles` / `binaries` yet** (they were all written before these keys existed): they appear once `final_datasets.py` is rerun, and nothing else in the pickle changes. (code that needs them, like `single_component_comparisons_mockobs.load_present_day_binaries`, falls back to reading the snapshot.)

**the binary-treatment subdicts** (`CoM`: binaries as their centre of mass; `luminous`: as the brighter component; `companions`: both components, so N2):

| key | what | `CoM` / `luminous` | `companions` |
|---|---|---|---|
| `coords` | intrinsic `StreamFrame` coords (see that section): `phi1`, `phi2` [deg], `pm_phi1`, `pm_phi2` [mas/yr], `r` [kpc, galactocentric], `vr` [km/s, relative to the progenitor] | N | N2 |
| `pos`, `vel` | galactocentric, Quantity [pc], [km/s], shape (·, 3) | N | N2 |
| `mass` | present-day mass, Quantity [Msun]. for `CoM` a binary's row is p1 + p2 | N | N2 |
| `in_rtid` | within `N_rtid_boundary` x rtid of the core (2 rtid in `final_datasets.py`) | N | ⚠️ N (wrong length, see the bugs list) |
| `inMW`, `trim` | intrinsic-frame masks; `trim` has length `inMW.sum()`, use as `x[inMW][trim]`. not used by the pipeline | N, `inMW.sum()` | N2, `inMW.sum()` |
| `phi2_straight`, `r_straight`, `vr_straight`, `pm_phi1_straight`, `pm_phi2_straight` | intrinsic-frame residuals about the progenitor orbit | N | N2 |
| `type`, `L` [Lsun], `R` [Rsun], `m0_zams`, `m0_effective` | BSE type, luminosity, radius; birth mass and BSE's effective initial mass (see the `star.mass0` section). **`luminous` and `companions` only** | N | N2 |

**added by `final_datasets.py`** (all N unless noted; `[_primaries]` = the key exists twice, without the suffix for the CoM trim/coords and with it for the primaries ones, used by `binaries` cases):

| key | what |
|---|---|
| `coords_obs` | CoM `SkyCoord` in the observed great-circle frame (e.g. `GD1Koposov10`), reflex-corrected: `.phi1`, `.phi2`, `.pm_phi1_cosphi2`, `.pm_phi2`, `.distance`, `.radial_velocity` |
| `sc_straighter`, `sc_straighter_primaries` | straightened observed-frame coord dicts: `phi1`, `phi2` [deg], `pm_phi1`, `pm_phi2` [mas/yr], `v_gsr` [km/s], `distance` [kpc], `d_phi2` [kpc], `v_phi1`, `v_phi2` [km/s]; all but `phi1` residuals about the track. `_primaries` = the luminous component's phi2 / v_gsr / distance (so the RVs include binary orbital motion) with CoM PMs |
| `unbound` | `~CoM['in_rtid']` (outside 2 rtid) |
| `nonrem` | `luminous['type'] < 10` |
| `alive` | `m0_zams` on the isochrone's initial-mass range |
| `acceptable_G` | 5 < isochrone mG < 30 (viamock's range) |
| `inMW_na` | all-True placeholder |
| `trim_new`, `trim_new_primaries` | 1-99 percentile trim in every observed dimension, on `sc_straighter` / `sc_straighter_primaries` |
| `phot` | isochrone photometry: `G`, `BP_RP`, `z` (absolute), `mG`, `mz` (apparent), `log_Teff` |
| `noise` | isochrone-based noise: sampled `phi1`, `phi2` [deg], `pm_phi1`, `pm_phi2` [mas/yr], `v_gsr_via`, `v_gsr_via10hr`, `v_gsr_desi` [km/s]; 1 sd errors `rverr_via`, `rverr_via10hr`, `rverr_desi` [km/s], `pm_err_gaia` [**µas/yr**] |
| `ontrack_velocity[_primaries]` | \|`sc_straighter[_primaries]['v_gsr']`\| < 30 km/s (true velocities, Jarvis+26-like) |
| `cut_for_catalog_photometry[_primaries]` | `unbound & trim_new[_primaries] & nonrem & ontrack_velocity[_primaries]`: the pool eligible for catalog photometry |
| `catalog_photometry[_primaries]` | same keys as `phot` plus `iso_dist` and `catalog_index`; NaN outside the matched stars |
| `matched_to_catalog_photometry[_primaries]` | True where a catalog row was assigned. ⚠️ length **`cut.sum()`**, not N (6273 / 6251 here) |
| `noise_catalog_photometry[_primaries]` | same keys and units as `noise`, recomputed from `catalog_photometry[_primaries]`; NaN where unmatched |

- `matched_to_catalog_photometry` is indexed within `cut_for_catalog_photometry`. `gmm.py` expands it with `matched_full = zeros(N); matched_full[cut] = matched`, and `np.isfinite(catalog_photometry['mG'])` gives the same full-length mask.
- `catalog_photometry['catalog_index']` comes out **float64** (the placeholder dict is `np.full(..., nan)`): -1 for unmatched stars inside `cut`, NaN outside `cut`. cast it before using it as an index.
- what each step does is under `final_datasets.py` and its catalog photometry subsection.


# supplemental simulations
one-off GD-1 runs outside the main grid, to separate the effect of binaries and of N from the effect of density. 5 copies each; **only copy 0 goes through the pipeline** (all five `gd1_2xN` copies were checked once, see below).

| sim | rvir0 | N stars | binaries | what changed vs the grid |
|---|---|---|---|---|
| `gd1_2xbinaries` | 1.5 pc | 14999 | 4491 (fbin 0.43) | binary fraction ~doubled (2665 binaries in the grid ICs) |
| `gd1_2xN` | 6.0 pc | 30000 | 5295 (fbin 0.21) | twice the stars, same binary demographics |
| `gd1_2xN_dense` | 3.0 pc | 30000 | 5295 (fbin 0.21) | the same population as `gd1_2xN` (identical `ICs/clust_30k.npz`), just denser |

**`gd1_2xN_dense`:** a 2xN run at rvir0 = 3 pc, at `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/gd1_2xN_dense/<copy>/`, 5 copies. the question: a denser progenitor should make a more substantial cocoon. does the GMM then find a real cocoon, rather than the degenerate split it gives for `gd1_2xN` (see below)? it uses the same stellar population as `gd1_2xN`, so the comparison is density only. it is in `simspect.supplemental_sims` (`'gd1_2xN_dense': {'orbit':'gd1', 'rvir':3.0}`) and `gd1_2xN` is not, so `final_datasets.py --which supplemental` and `gmm.py --supplemental True` do `gd1_2xbinaries` + `gd1_2xN_dense`. the `gd1_2xN` pickles are still on scratch but aren't rebuilt or refit (add it back to the dict for that). noiseless and desi `_supplemental` tables exist; its results aren't written up here yet. watch whether its progenitor is also still bound at the present day: that's what seems to break the `gd1_2xN` fit.

- **the stellar populations are neither `lm` nor `hm`**: each was re-drawn from the IMF. initial cluster masses 11184 / 17494 Msun (from `ICs/*.txt`). don't compare them to the grid as if only one thing changed without keeping that in mind.
- same orbit, age (2700 Myr -> present day = `data.270`), `init_displacement`, potential and petar options as the grid's `gd1` (see `run_gd1.slurm`; the `-t 1000` there was extended by the resume scripts).
- where: `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/gd1_2xbinaries/<copy>/` and `.../gd1_2xN/<copy>/`, with the ICs and slurm scripts alongside. **netscratch, so the 90-day purge applies** -- copy them to holystore. completeness: `gd1_2xbinaries` copies 0, 1, 4 reach the present day (2, 3 stopped at 2320 / 2210 Myr); all five `gd1_2xN` copies do (run to 5000 Myr).
- loading: `retrieve_sim_info` can't find them, so `inspect_new_sims.py` has `supplemental_sims` (name -> orbit, rvir), `supplemental_base_path`, `retrieve_supplemental_sim_info(sim, copy)` and `prepare_nbody_data_supplemental(sim, copies=(0,), ...)`, the analogue of `prepare_nbody_data_anycopy`. to add another one-off run, add a row to `supplemental_sims`.
- pickles: `p27_data_dicts/gd1_2xbinaries_1.50.pickle`, `gd1_2xN_6.00.pickle` (filename is `<sim>_<rvir>`). each supplemental sim gets its own `default_rng(seed=42)`, so its noise draws don't depend on what ran before it (no noisy case has been fit on them, though).
- GMM tables: `data/gmm_tables/noiseless_CoM_supplemental.fits` and `noiseless_binaries_supplemental.fits`, with the extra `sim` column. sanity plots go in the usual `p27_sanity_plots/<case>/` directory as `<sim>_<rvir>.pdf`.
- plotted in the `summary_version = 'gd1'` summary figure (`plots/summary_cocoon_ts_properties[_density]_gd1only.pdf`) in GD-1's colour, star = 2x binaries, diamond = 2xN, filled = CMs, open = with binary orbits.

what the noiseless fits give, against the grid's `gd1` at the same rvir0 (the grid tables were checked against a refit with the same `gmm.py`; they match exactly, so the comparison is like for like):

| | f_cocoon CoM | sigma_phi2,c CoM [deg] | sigma_vgsr,c CoM | f_cocoon binaries | sigma_vgsr,c binaries | sigma_vgsr,thin binaries |
|---|---|---|---|---|---|---|
| grid gd1, 1.5 pc | 0.083 | 0.82 | 1.48 | 0.114 | 6.03 | 0.51 |
| `gd1_2xbinaries` | 0.119 | 1.19 | 2.68 | 0.173 | 11.68 | 0.69 |
| grid gd1, 6 pc | 0.042 | 0.49 | 1.01 | 0.092 | 7.05 | 0.42 |
| `gd1_2xN` | (0.293, see below) | (0.23) | (0.70) | 0.077 | 5.06 | 0.64 |

(velocities in km/s.)
- `gd1_2xbinaries`: the **CoM** cocoon is bigger and hotter too, not just the binaries one -- i.e. doubling the binaries changes the dynamics (binary heating / ejections), not only the RVs. binary fraction among the fitted (unbound, trimmed) systems is 0.258 vs 0.128 in the grid run.
- `gd1_2xN`: **the progenitor survives to the present day** (6087 of 25409 systems unbound), whereas the grid's gd1 6 pc run has fully dissolved (every system unbound). so it's a younger stream with a bound progenitor and strong epicyclic feathers near it -- a different dynamical state, not "the same stream with more stars."
  - ⚠️ its **noiseless CoM fit is degenerate**: the "cocoon" is only 1.2x wider than the thin stream in phi2 and is offset by 0.77 sigma (fails `results.py`'s `okay_mean` < 0.5), because the GMM splits the thin stream into its +phi2 and -phi2 halves (see the sanity plot). so its f_cocoon = 0.29 is not a cocoon fraction. it's still drawn on the summary figure, as the grid points are (no `okay_mean` filter there). the binaries fit (cocoon = binary orbital motion, as usual) is fine.
  - **checked on all five copies, and it's not a fluke of copy 0.** the same gmm.py recipe (noiseless, Powell, same init and selection) gives a degenerate CoM split every time. the binaries fits are stable:

    | copy | unbound | CoM f_c | CoM sigma_c / sigma_thin (phi2, v_phi1, v_phi2, v_gsr) | CoM max \|M_c\|/S_c | binaries f_c | binaries sigma_vgsr,c |
    |---|---|---|---|---|---|---|
    | 0 | 6087 | 0.293 | 1.24, 1.30, 1.15, 1.32 | 0.77 | 0.077 | 5.06 |
    | 1 | 6271 | 0.544 | 1.17, 1.11, 1.05, 1.19 | 0.69 | 0.073 | 6.05 |
    | 2 | 6397 | 0.245 | 1.20, 1.27, 1.20, 1.22 | 0.86 | 0.073 | 6.14 |
    | 3 | 5512 | 0.500 | 1.05, 1.11, 1.06, 1.01 | 0.82 | 0.073 | 5.19 |
    | 4 | 6789 | 0.374 | 1.07, 1.16, 1.01, 1.15 | 0.86 | 0.073 | 5.51 |

    the CoM f_c jumps around between 0.25 and 0.54, while the two components stay within ~1.3x of each other in width. that's two near-equal gaussians splitting one population, not a thin stream plus a cocoon. ΔBIC still favours two components (690-910), which is just the stream not being a single gaussian. by eye, the "cocoon" component is the +phi2 side, concentrated in the tail just leading the bound progenitor (phi1 ≈ -35 to -20), where phi2 steps by ~0.5 deg across the progenitor and v_gsr has strong epicyclic feathers. a degree-5 `poly_straightening` over the whole stream can't follow that step. **so most likely a straightening problem for streams with a surviving progenitor**, not a dynamical cocoon. not tested: masking a window around the progenitor, or straightening each tail separately.
  - the copy 1-4 pickles are `p27_data_dicts/gd1_2xN_copy<k>_6.00.pickle` (`python final_datasets.py --which supplemental --sims gd1_2xN --copies 1 2 3 4`). they aren't in any `gmm.py` table: `--supplemental True` only fits copy 0.


# the real Jarvis+26 GD-1 data
the same two-component mixture, with `x_err`, run on the Jarvis+26 DESI GD-1 members themselves using their own errors. the sims can then be compared to a fit made with exactly the same model, instead of to the published numbers.

```
python final_datasets.py --which jarvis   -> <scratch>/p27_data_dicts/jarvis26.pickle
python gmm.py --jarvis-data True          -> data/gmm_tables/jarvis26.fits  (+ p27_sanity_plots/jarvis26/jarvis26.pdf)
```

**the pickle** (the `args.which in ['all','jarvis']` block at the end of `final_datasets.py`) is built from `data/jarvis26_Table7.fits`, all 679 rows in table order. it is much smaller than a sim pickle and has three keys:

| key | what |
|---|---|
| `sc_straighter` | same key and coordinate names as the sims' straightened CoM dict. `phi1`; `phi2` = `DELTA_PHI2` [deg]; `pm_phi1`/`pm_phi2` = `DELTA_PM_PHI1`/`DELTA_PM_PHI2` [mas/yr]; `v_gsr` = `DELTA_VGSR` [km/s]; `v_phi1`/`v_phi2` [km/s], converted from the PM residuals. all are residuals about **Jarvis' own track**, not our straightening |
| `errors` | the 1 sd uncertainties, under the same keys as `sc_straighter`. there's no noise draw: the data are real, so the noise is already in them. `phi1`/`phi2` = 0 (the table has no position errors), `pm_phi*` = `PM_PHI*_ERR`, `v_gsr` = `V_ERR`, `v_phi*` = the converted PM errors |
| `p_contam` | `1 - P_THIN - P_COCOON`, an astropy `Column` |

- **PM -> v_phi uses a distance track, not per-star distances**: a `CubicSpline` through Jarvis' track distances (phi1 = -90 ... 40 deg, 9.12 ... 19.95 kpc), evaluated at each star's phi1. there's no distance uncertainty in the v_phi errors and `DISTMOD` isn't used. this is a known limitation, accepted for now.

**the fit** (`gmm.py --jarvis-data True`; also `--jarvis_data`, or `jarvis_data=True` in the file):
- forces `noise=None`, `include_binaries=False`, and prints a line if it had to override them. it raises if combined with `--supplemental`.
- `use = p_contam < 0.5`: 608 of the 679 stars.
- `x_data = sc_straighter[keys][use]`, used as it is (it's real data, so no noise is added). `x_err = errors[keys][use]`, so the fitted sigmas are intrinsic, deconvolved widths.
- the initial guesses, Powell, BIC and the optional `--constrain-widths` stage are the same as for the sims.
- the case name is `jarvis26[_constrained]`, not `noiseless_CoM`, so the grid table isn't overwritten. the table has one row: `sim = 'jarvis26'`, `orbit = 'gd1'`, `Rvir0 = NaN`.

**result** (unconstrained; phi2 in deg, velocities in km/s), against the published Jarvis+26 values hardcoded in `results.py`:

| | this model | Jarvis+26 |
|---|---|---|
| f_cocoon | 0.328 | 0.30 ± 0.02 (0.318 = sum P_COCOON / sum(P_THIN + P_COCOON) on the same 608 stars) |
| sigma_phi2, thin | 0.237 | 0.23 ± 0.01 |
| sigma_v_gsr, thin | 2.52 | 2.49 ± 0.28 |
| sigma_phi2, cocoon | 2.15 | 2.18 ± 0.17 |
| sigma_v_gsr, cocoon | 6.13 | 6.13 ± 0.75 |

- the other two dimensions: thin sigma_v_phi1, sigma_v_phi2 = 2.87, 0.87; cocoon 33.7, 3.84. the cocoon's 34 km/s in v_phi1 looks large; it may come from the distance track, but that hasn't been checked.
- ΔBIC (one − two) = 1266, so a cocoon is strongly preferred. every cocoon mean is within 0.25 sigma of zero (passes `okay_mean`).
- **so the model reproduces the published decomposition** (and `DELTA_PHI2` is in degrees).
- **the sims are colder and less cocoon-rich than real GD-1 in every panel except the binaries cocoon v_gsr.** noiseless grid maxima: f_cocoon 0.25 (CoM), cocoon sigma_phi2 1.17 deg, thin sigma_phi2 0.18 deg, thin sigma_v_gsr ≤ 1.06 km/s. the cocoon sigma_v_gsr is 0.7-3.7 km/s for CoM but 4.7-8.3 km/s with binary orbits, which brackets the real 6.1.
- `results.py` draws these values as black dashed lines on the summary figure.

still open: per-star distances (`DISTMOD`) and distance errors in the v_phi conversion, or fitting in proper motions directly / dropping the PM dimensions; whether Jarvis' track matches our straightening closely enough that the residuals are comparable.

# TODOs
the live ones, collected here. details are in each script's section below.
- **catalog photometry, open questions:** which N-body stars are eligible to be matched, dereddening, `p_mem`, and the fact that the noisy sample size *is* the catalog size (46 stars for C-19). see the "catalog photometry" subsection under `final_datasets.py`.
- **refit the noisy tables on the colour-Teff pickles.** the grid pickles were regenerated with `Teff_from_gaia_color_temp_relation` (see the colour -> Teff subsection of photometry), but going by mtime only `via_noise_binaries.fits` has been refit since; check, and run what's missing:
  ```
  python final_datasets.py --which supplemental        # copy 0 of whatever is in simspect.supplemental_sims
  for n in via desi; do for b in False True; do python gmm.py --noise $n --include-binaries $b; done; done
  python gmm.py --noise desi --include-binaries {False,True} --supplemental True
  ```
  only the noisy tables need refitting: noiseless fits never use photometry or noise, and `jarvis26` has no photometry step. the Via tables are the ones that change (Via errors depend on `log_Teff`). the DESI ones should come out the same (Jarvis `V_ERR`s by rank, no Teff dependence); refit them only so every table comes from the same pickles. there's no `via10hr` table on disk.
- **binaries selection effect (for the paper):** the noisy `binaries` cocoons are much hotter in v_gsr than the noiseless one, mostly because the bright catalog-matched stars are ~2x as binary-rich, not because of noise. keeping it in the observed cases; still to do: a "noiseless, observed stars" `gmm.py` case, and checking that the sim's binary fraction vs mass is realistic. see "binaries in the observed samples" in the GMM section.
- **sigma_binary is mis-specified for a heavy-tailed dv:** one gaussian sigma_binary on every star, fit from two-epoch RVs (the current multi-epoch likelihood), drives the stream dispersion to 0 at every rvir0 and f_repeat: a few fast binaries set sigma_binary and it then over-explains everyone else. a per-star mixture (binary with probability f_b, free, else no binary term) recovers sigma_stream to within ~0.1-0.35 km/s in a scratch test; put it in `single_component_comparisons_mockobs.py` (`multiepoch_nll_per_star` is the building block: logaddexp of it with sigma_binary and with 0). alternatives: a heavier-tailed (e.g. student-t) binary term, or a dv distribution from the sim. also: rerun `final_datasets.py` so the pickles carry the petar `binaries`. see `single_component_comparisons_mockobs.py`.
- **NUTS on Jake's catalogs:** the 10-draw runs (DESI, 5 yr, shared and indep) and the 50-draw rvir0 = 3 run with the stars varied too (indep) are done; still to do: look at the divergences / r_hat in the draws runs (indep 3 pc especially), and Via / 1 yr over draws if wanted. see "NUTS on Jake's catalogs" under `single_component_comparisons_mockobs.py`.
- **cocoon populations:** the cocoon's systems are more often binaries, more massive and brighter than the thin stream's (consistent with Weatherford+2025's mass-segregation/ejection selection effect). only shown for GD-1 so far: run the f_bin / <G> / <M_sys> cells on the other orbits, and decide whether the `'fitted'` <M_sys> needs asymmetric/bootstrap errors. see "thin stream vs cocoon populations" under `results.py`.
- **`gmm.py` / model choice:** `gmm.py` fits a single gaussian too and writes `BIC_one` / `BIC_two` columns (see the model choice section). `BIC_one < BIC_two` = one component preferred = a **non-detection of a cocoon**. still to do: make sure every table has the columns, use them in `results.py`, and maybe save the one-component params, `N = len(x_data)` and `result.fun` too. this matters more with catalog photometry, because the noisy fits have N = 46 (c19), 129 (pa5), 144 (jet), 334 (aau), 1578 / 679 (gd1 via / desi), and several of the small-N fits have already come out degenerate (see the GMM results bullet in the bugs list).
- `gmm.py`: AAU at rvir0 = 6 pc gives a weird result; traced to the 2 rtid tidal boundary, and the decision is to live with it rather than relax the cut.
- `gmm.py`: run all noise cases (see the case table below); rewrite the results notebook (`results.py`); add stuff to overleaf.
- supplemental sims (`gd1_2xbinaries`, `gd1_2xN`, `gd1_2xN_dense`): copy them off netscratch; use the other copies; Via noisy cases; decide what to do with the degenerate `gd1_2xN` CoM fit (it happens in all 5 copies; probably straightening around the surviving progenitor); write up the `gd1_2xN_dense` results. see their section.
- `get_init_displacements.py`: rewrite so orbit rows are looked up by the short keys (see its section -- right now it silently writes an empty file).
- observed-frame straightening: not 100% convinced it's successful for M3 or Pal 5, which have notable diverging tails that end up in the cocoon. (M3 is currently dropped from the pipeline anyway.)
- `m3` frame: may settle for "stream on an M3-like orbit" (see the gala gotchas section).
- `animations/`: make the movie scripts array jobs (one sub-job per frame) instead of looping. (`particle_spray_vs_nbody.py` already is.)
- **one noise seed (0) everywhere for the final version:** the grid pickles use seed 0, but the supplemental loop in `final_datasets.py` uses `default_rng(seed=42)`, and the `desi_noise_*_supplemental` tables are fit on those draws. switch it and refit.
- GMM: emcee for posteriors after maximum likelihood (see the GMM "next" section).
- Jarvis+26 real-data fit, still open: `DISTMOD` and distance errors, or fitting in PMs / dropping the PM dimensions; whether their track matches our straightening. see "the real Jarvis+26 GD-1 data".


# map of code
## `/data`
- `FINAL_ics_nolmc.csv`: credit Vedant Chandra, present-day positions, velocities of some promising streams. Age and progenitor mass estimates also given. **the `name` column for the six run orbits has been renamed to the short keys i index by** (`gd1`, `pa5`, `aau`, `m3`, `c19`, `jet`) so that `prog_tab[prog_tab['name']==orbit]` just works. everything else in the table keeps its original name.
- `init_displacements.txt`: generated in `get_init_displacements.py`, Vedant's progenitor locations back-integrated by their stream ages, plus 100 Myr to account for initial expansion due to massive star evolution, rounded to a multiple of 10 so that I can safely output sim snapshots every 10 Myr and get the present day in the final snapshot. **positions in this file are in pc** (that's what `petar.init` wants), velocities in km/s -- see the units bullet in the `paf` notes, this caused a real bug.
- `bpw25_catalogs/<orbit>.fits`: credit Bonaca & Price-Whelan (2025), one per stream, named by the short keys (`gd1`, `pa5`, `aau`, `m3`, `c19`, `jet`); loaded by `noise.load_gaia_catalog(orbit)`. Gaia columns (`phot_g_mean_mag`, `bp_rp`, astrometry) plus `dist`, `Vr`, `FeH`, `p_mem`. **`dist` is mostly masked** -- only 43 of 1578 `gd1` rows have one -- and there's no extinction column, so `bp_rp` is presumably not dereddened. (`stream_catalogs.py` still reads the old long filenames, e.g. `GD-1.fits`.)
- (not in `data/`: `nibauer_gd1.pdf` in the repo root is the Nibauer+26 GD-1 velocity dispersion paper, ApJ 1004:62. gitignored, so it only exists in this checkout.)
- `jarvis26_Table7.fits`: Jarvis+26 table 7, DESI GD-1 members (679 stars). also at `/n/home02/amphillips/data/`. source of `jarvis26.pickle` (see "the real Jarvis+26 GD-1 data") and of the DESI RV errors (`noise.desi_RVerr_from_jarvis`).
- `gmm_tables/<case>.fits`: output of `gmm.py`, one table per case (see `gmm.py` below). what's there:
  - `noiseless_{CoM,binaries}[_supplemental]`: noiseless fits (don't depend on photometry or noise draws).
  - `desi_noise_{CoM,binaries}[_supplemental]`, `via_noise_{CoM,binaries}`: catalog photometry, `x_err`, unconstrained. grid ones on seed-0 draws, supplemental on seed 42. only `via_noise_binaries` has been refit on the current (colour-Teff) pickles; see the TODO.
  - `jarvis26`: the fit to the real Jarvis+26 data.
  - `desi_noise_nominal_multiepoch_scan_{fit_shared,fit_indep,dv_deconvolved,dv_raw,avg,flag_avg,flag_ep0}_dt1yr_sd0.17yr` and `..._{dvcom_deconvolved,dvcom_raw,dvcom_noiseless}`: from `single_component_comparisons_mockobs.py` with `survey = 'desi'`; the two-epoch ones carry the epoch baseline (1 yr ± 2 months) in the name, the one-epoch `dvcom_*` don't. `via_noise_nominal_multiepoch_scan_{avg,flag_avg}_dt1yr_sd0.17yr` are the same with `survey = 'via'`. one per binary treatment (`fit_indep`: only its total width is meaningful; `dvcom_*`: one epoch, dv vs the true CoM v; `dvcom_noiseless`: dv without noise; `avg`: two-epoch means; `flag_*`: RV-variable repeats removed, `flag_ep0` biased low).
  - `desi_noise_nominal_multiepoch_scan[_indep]`: the same numbers under older column names; superseded, nothing reads them.
  - `desi_noise_nominal_repeatfracscan`: the previous version of that experiment; no longer read or written.
  - `jake_nuts_N160_rep35_{shared,indep}_capstd[_ndraw<n>[_varystars][_pooled]]`: NUTS posteriors for sigma_binary and sigma_stream on the `for_jake/` catalogs: one row per survey x rvir0 x baseline (one draw of the repeats), per draw (`_ndraw<n>`), or pooled over the draws (`_pooled`). `_ndraw10`: the 160 brightest, 10 draws of the repeats, every rvir0. `indep_capstd_ndraw50_varystars`: rvir0 = 3 pc only, 50 draws of the 160 stars (from the full catalog) and of the repeats. the Jarvis comparison figure reads the `indep` pooled ones. see "NUTS on Jake's catalogs" under `single_component_comparisons_mockobs.py`.
  - superseded, not read or written by anything: `desi_noise_brightest160_repeat34`, `desi_noise_brightest160_repeatscan`.
- `data_dicts/`: empty. the pickles actually live on scratch (see logistics).
- `for_jake_truth/`: per-star truth for the mock catalogs in `/for_jake` (see "mock catalogs for Jake Nibauer" under `single_component_comparisons_mockobs.py`). not for sending.

## `/for_jake`
`rvir0<X>_<desi|via>.fits`: the three-epoch mock GD-1 RV catalogs sent to Jake Nibauer, from the last cell of `single_component_comparisons_mockobs.py` (see its "mock catalogs for Jake Nibauer" subsection).
- `particle_sprays/`: the gala particle sprays from `particle_spray.py`, `spray_<orbit>_hm_0.75_<copy>.h5` (~30-45 MB each). **in `.gitignore`** because they're large, so they only exist where they were generated -- rerun `sbatch particle_spray.slurm` to rebuild them. `logs/` holds the slurm logs (`out.<array index>`) and has to exist before submitting.

## `/scripts`
- `PETAR_ANALYSIS_FUNCTIONS.py`: a bunch of functions, mostly migrated from other projects; imported as `paf`. see the notes section below for the parts that are actually load-bearing.
- `streamframe.py`: credit Jake Nibauer, transformation to stream frame coordinates as seen from the Galactic center
- `read_mist_models.py`: the MIST team's reader (`ISOCMD`). not on the live path -- the isochrone comes from `artpop` -- but `noise.py` imports it.
- `vedant.mplstyle`: credit Vedant Chandra, plotting style stuff. note it sets `savefig.format : pdf`, `savefig.dpi : 300` and `savefig.bbox : tight`, so movie frames saved as png should pass a lower `dpi` (a 14x14 inch figure at 300 dpi is ~4200 px a side).

## `/memos`
short LaTeX write-ups of the maths behind the code, one directory each (`.tex` + compiled `.pdf`; LaTeX build files are gitignored). build with `pdflatex <name>.tex` (run twice for the cross-references).
- `multiepoch_likelihood/`: the derivation of `multiepoch_nll_per_star` in `single_component_comparisons_mockobs.py`. why the shared stream offset puts sigma_stream^2 on the off-diagonals of each star's covariance, the closed-form inverse (Sherman-Morrison) and log determinant (matrix determinant lemma) the code evaluates, the explicit two-epoch form, the independent / one-epoch / missing-epoch cases, and the same result by integrating out the star's offset.

## `/old`
stuff that is superseded or migrated from older repositories (esp. `~/stream_velocity_structures`). kept for reference, not imported by anything live.
- `DESI_comparison.py`: noise old simulation data like Gaia+DESI and compare to Jarvis+2026 cocoon detection. superseded by `noise.py` + `gmm.py`'s `desi` case.
- `prog_properties_summary.py`: cocoon fractions/dispersions for GD-1 portion of old sim grid
- `cocoons.py`: the hand-cut (`get_cocoon_selection`) cocoon fractions. superseded by the GMM.
- `develop_GMM.py`: the development version of the mixture model. `gmm.py` is the cleaned-up version; things that exist **only** here: `AIC` (`BIC` has been copied into `gmm.py` as `bic`), `cocoon_fraction_constraint`, the commented-out profile likelihood (`nll_fixed_cocoon` / `profile_cocoon_fraction`), and the `max |mu|/sd` mean-bound warning. see the GMM section.

## `/animations`
`velocity_movie.py`, `grid_movie.py`, `grid_movie_long.py`, `grid_rotation_movie.py`: scripts to run with a slurm wrapper for animations. TODO: make these parallel so that the wrapper is submitted as an array job where each sub-job generates one frame. way faster than doing this in a loop.
- `grid_movie.py`, `grid_movie_long.py` and `grid_rotation_movie.py` all build `extended_grid_info(scratch=True)`, i.e. the purged netscratch copy -- switch to `scratch=False` before re-running them (see `extended_grid_info`).
- `grid_movie_long.py` line 25 does `sys.path.append(rotated_pos)` before `rotated_pos` exists, so it dies with `NameError` at import. delete the line.
- frames go to `/n/netscratch/conroy_lab/Lab/amphillips/movies/<movie name>/frame_<index:05d>.png`.

### `particle_spray_vs_nbody.py`
2x2 movie comparing direct N-body (left column) with a gala particle spray (right column) for the same progenitor: `gd1`, `hm`, rvir0 = 0.75 pc, copy 0 (constants at the top). top row is galactocentric x-z (±30 kpc); bottom row is a zoom on the progenitor (`zoom_half_width`, in pc) -- the N-body in the **core frame**, the spray relative to its **integrated** progenitor, because the petar core drifts off the non-interacting orbit (~300 pc at 1.7 Gyr, ~480 pc at present day for this sim). already an array job:

```
python particle_spray_vs_nbody.py --make-spray   # once (~1.5 min): writes spray_gd1_hm_0.75_0.h5 next to the frames
python particle_spray_vs_nbody.py -i 170         # one frame (~12 s); index = petar file index, t = 10*index Myr
```

- frame index == petar file index, 0..`age/10` (0..270 for gd1, 270 = present day; 170 = 1 Gyr before present).
- N-body positions are singles + binary CoMs straight from `data.<k>.single`/`.binary` (already core frame), plus `core.pos[k]` for galactocentric. no `prepare_nbody_data`, which is slow and does streamframe/straightening work a movie doesn't need.
- the spray is integrated at dt = 1 Myr and snapshotted every 10 Myr into gala's HDF5 file, so output index `j` is t = 10 j Myr and lines up with `data.<j>` (the frame asserts this). it starts from the same t = 0 `init_displacement` in the same `BovyMWPotential2014`.
- **exactly 5000 particles** at present day: an `n_particles` array with one release per tail at 2500 evenly spaced steps, rather than `release_every`.
- the Chen DF's progenitor mass follows the N-body **bound mass** from `data.tidal`; the progenitor *potential* is a fixed Plummer at the initial mass with b = rvir0 x 3pi/16 (gala can't evolve it). debatable; the alternative is a constant mass.
- what the frames show: the spray has essentially no particles inside a few tens of pc of the progenitor (Chen releases near the Jacobi radius, ~40 pc here), and the N-body stream is much longer and more diffuse at large scale, with escapers out to ±25 kpc by 1.7 Gyr (probably early escapers from the first ~100 Myr of massive-star mass loss; not checked).

### gala mockstream notes (gala 1.9.1)
- `MockStreamGenerator.run(..., n_particles=...)` takes an **array the length of the time grid** (particles per tail per step), which is how to hit an exact total; an int is combined with `release_every`.
- `prog_mass` may be an **array** over the time grid -- it only sets the DF's release scale, **not** the progenitor potential's mass.
- with `output_every` you must also pass `output_filename`; the file has `stream/pos`, `stream/vel` shaped `(3, n_outputs, n_particles)` (kpc, kpc/Myr), `nbody/pos`/`vel` with the progenitor at index 0, and `stream/time`. **not-yet-released particles are NaN** (h5 `fillvalue`), so mask with `isfinite`. n_outputs = `(ntimes-1)//output_every + 1`, plus one if the last step doesn't land on the cadence.
- `paf.gen_stream` is the older wrapper (int `release_every`, scalar mass, Plummer b = 0.75 pc default); the movie script doesn't use it.

## `/plots`
`summary_<case>.pdf` and friends, made by `results.py` (the `savefig` lines there are mostly commented out, so these are only regenerated by hand). `fig1_motiv_*.pdf` is from `motivation.py`. `jarvis_comparison_single_component*.pdf` and `gd1_desi_nominal_multiepoch_scan.pdf` are from `single_component_comparisons_mockobs.py` (the Jarvis figure's `savefig` is live, the scan plot's commented out; `gd1_desi_*repeat*scan.pdf` are from its previous version), and so are `gd1_<survey>_jake_nuts_<tag>[_ndraw10_5yr].pdf` (the NUTS summaries; live `savefig`s). `orbit_summary.pdf` is from `get_init_displacements.py`. `plots/old/` is from the hand-cut and early-GMM era.

## `get_init_displacements.py`
integrates every orbit in `FINAL_ics_nolmc.csv` for peri/apo/e (`plots/orbit_summary.pdf`), then back-integrates the run orbits by age + 100 Myr and writes `data/init_displacements.txt` (pc, km/s) for `petar.init`.

**stale, and silent about it:** `names_to_run` (and the `nl` label list) still hold the long names (`'ATLAS-Aliqa Uma'`, `'GD-1'`, ...) which no longer match the renamed `name` column, so re-running it right now overwrites `init_displacements.txt` with an *empty* file and no error. the header TODO says so; fix the list before touching it.

## `inspect_new_sims.py` (imported as `simspect`)
functions for loading N-body data, putting it in observed stream frames, and straightening the tracks. contents:
- `prepare_nbody_data(path, include_photometry, i, apo, init_displacement, **kwargs)`: wraps `paf.load_core` + `paf.intrinsic_stream_data_v3`; assumes snapshots every 10 Myr. returns `core, data_dict, CMdict, lumdict, inMW, trim` (plus `G, BP, RP, z` if `include_photometry=True`, which is the slow per-star blackbody path, not the isochrone). extra kwargs go to `intrinsic_stream_data_v3` -- this is how `N_rtid_boundary` gets in.
- `supplemental_sims`, `retrieve_supplemental_sim_info(sim, copy)`, `prepare_nbody_data_supplemental(sim, copies=(0,), verbose, **kwargs)`: the loaders for the one-off runs outside the grid (see the supplemental simulations section). same return signature as `prepare_nbody_data_anycopy`.
- `prepare_nbody_data_anycopy(orbit, stellar_pop, rvir_index, copies, verbose, **kwargs)`: same thing but loops `copy` and returns the first realization that actually finished. an unfinished copy still has a (short) `data.core`, so it gets all the way into `intrinsic_stream_data_v3` before raising `FileNotFoundError` on the missing `data.<file_index>` -- that exception is the signal. dedupes on path because `retrieve_sim_info` ignores `copy` for `circ`. returns `(prepare_nbody_data output), path, apo, age, init_displacement, copy`. use this instead of hand-rolling fallbacks.
- `rotation_matrix(a,b,c)`: rotate by a,b,c about x,y,z. composed `Rz @ Ry @ Rx`, i.e. x rotation applied first -- chosen for movie purposes.
- `streamframe_coords_observed(orbit, data_dict, prog_tab)`: the observed-frame transform. returns `(SkyCoord, frame)`. see the stream frame section.
- `prog_orbit_track` / `observed_orbit_track` / `chop_orbit_track` / `straightened_obscoords_orbit_interp`: the observed-frame orbit straightening. see its own section.
- `poly_straightening(coords, tc=None)`: loops the keys of a coord dict and subtracts a degree-5 polynomial fit via `paf.straighten_stream_polynomial`. applied *after* the orbit-interp straightening to mop up residual tilt. `tc = [mask1, mask2]` selects which stars the polynomial is **fit** to, but it is **subtracted from all of them** -- so outliers don't drag the fit but stay in the dataset. `tc=None` means all-True. copies rather than subtracting in place; note `coords_straighter['phi1']` is the same array object as `coords['phi1']`.
- `clip_coords(coords, [mask1, mask2])`: apply `x[mask1][mask2]` to every key of a coord dict at once.
- `outlier_clip(vr, pmphi1, pmphi2)`: **falling out of use** -- replaced by the percentile trim in `final_datasets.py`. note the code cuts at |pm| < 1 mas/yr and |v| < 75 km/s, not the 1.5 / 100 its docstring says.
- `get_cocoon_selection(coords, cuts)`: the old hand-cut cocoon selection (OR of `|x| > cut` over phi2, pm_phi1, pm_phi2, v_gsr). legacy; only the `__main__` block uses it.
- everything under `if __name__=='__main__':` is the old hand-cut cocoon-separation panels. **the guard is active**, so importing the module does not run it (it does still run `extended_grid_info()` and read `FINAL_ics_nolmc.csv` at import time).

## `noise.py` (imported as `noise`)
survey error models and isochrone photometry. the `__main__` block is entirely commented out (exploratory; it also references `gmm.trim_obstream_percentile`, which now lives in `final_datasets.py` and has a copy in `noise.py`).

the `# SCRATCH::::` cells between `trim_obstream_percentile` and `assign_photometry_from_catalog` (the catalog-photometry prototyping) are commented out, so `import noise` doesn't load a simulation or draw plots.
- `gaia_g_to_lsst_z` / `lsst_z_to_gaia_g`: the RTN-099 sec 1.3.4 polynomial and its exact inverse (0.014 mag rms).
- `desi_RVerr_from_jarvis(mag)`: **the live DESI model.** the 679 Jarvis+26 table 7 `V_ERR`s (0.46-9.98 km/s), assigned rank for rank: brightest finite `mag` gets the smallest error. stars past the 679th, and NaN mags, get NaN (so NaN DESI noise too). output is `len(mag)`, original order. `final_datasets.survey_noise` calls it with apparent G.
- `desi_RVerr(zmag, feh=-2.0)`: an analytic model, `log10(err) = -0.47 + 0.27(z-16) - 0.23[Fe/H]`, from **apparent** z mag. not called.
- `via_RVerr(G, log_Teff, feh=-2.0, exptime_s=3600, nexp=1, **kwargs)`: thin wrapper on `viamock.get_viaspec_errors`, apparent G. viamock's tables cover G in (5, 30); outside that `final_datasets.py` leaves NaN.
- `get_Teff(R, L)`: Stefan-Boltzmann.
- `build_isochrone_table` / `isochrone_photometry` / `gaia_from_isochrone`: ZAMS mass -> Gaia mags + `log_Teff` by interpolating a MIST isochrone. **read the precision subsection before touching them.** see the photometry section.
- `isochrone_cmd_track` / `nearest_isochrone_point` / `Teff_from_gaia_isochrone`: the reverse direction -- a catalog star's (absolute G, BP-RP) -> nearest point on the isochrone -> `log_Teff`. **doesn't set the catalog Teffs**; `nearest_isochrone_point` is only used for the `iso_dist` diagnostic. see the CMD-space lookup subsection of photometry.
- `Teff_from_gaia_color_temp_relation(track, BP_RP, deg=2, G_max=5000.0, return_poly=False)`: **the live catalog Teff.** log Teff from BP-RP alone, via a degree-2 polynomial fit to the isochrone track; NaN outside the fitted colour range. see the colour -> Teff subsection of photometry.
- `load_gaia_catalog(orbit)`: reads `data/bpw25_catalogs/<orbit>.fits`.
- `assign_photometry_from_catalog(gaia_iso, distances, catalog, track, N=None, eligible=None, **teff_kwargs)` -> `(phot, matched)`. `gaia_iso` is the `phot` dict from `isochrone_photometry` (only its G is used, to rank); `distances` a Quantity; `track` from `isochrone_cmd_track`. ranks the N-body stars (finite iso mG, and `eligible` if given) and the catalog (finite G and BP-RP) by apparent G and hands the brightest N N-body stars the catalog's **apparent G and BP-RP exactly**; converts to absolute with the N-body distances; `log_Teff` from `Teff_from_gaia_color_temp_relation` (colour only). `teff_kwargs` (`deg`, `G_max`) go to it (no caller passes any). `phot` has the same keys as `final_datasets`' photdict (`G`, `BP_RP`, `z`, `mG`, `mz`, `log_Teff`) plus `iso_dist` (diagnostic only now) and `catalog_index` (-1 unmatched), all full length in the original order, NaN where unmatched. `log_Teff` is also NaN for matched stars outside the fit's colour range. `N=None` = every usable catalog row; capped at the number of rankable N-body stars (prints if so).
  - **pass `eligible`** (e.g. `unbound & nonrem & trim`), or bound progenitor stars compete for the catalog's rows.
  - measured on `gd1` hm rvir_index=0 with `eligible = unbound & nonrem & trim_new` (6273 stars): all 1578 catalog rows matched, but the matched N-body stars have **isochrone** mG 15.0-23.6 against the catalog's 12.2-20.1 -- the sim has fewer bright stream stars than GD-1, so rank matching brightens each star by ~3 mag, and the catalog-based `log_Teff` runs a median 0.15 dex hotter than the isochrone value for the same star (measured with the old nearest-point Teffs; the colour-fit ones are ~0.026 dex cooler). that's the scheme working as designed, but it means `N` (and the eligible pool) sets how deep into the N-body MS the catalog reaches. `iso_dist` median 1.7, 267 of 1578 > 3.
- `get_gaia_photometry`, `g_phot`: the blackbody + top-hat fallback, duplicated here from `paf` "to illustrate the difference in stellar populations." not used by the pipeline.

`assign_photometry_from_catalog` is on the live path: `final_datasets.py` calls it once per sim, and its output is what the noisy `gmm.py` cases fit. see the catalog photometry subsection under `final_datasets.py`.

## `final_datasets.py`
assembles one data dictionary per simulation and pickles it. imports `paf`, `simspect`, `noise`. no `__main__` guard (it's a script). the per-sim work is `build_data_dict(orbit, prepared, catalog, rng)` and the pickling `dump_data_dict(data_dict, name, rvir)`; `--which {all,grid,supplemental,jarvis}` picks which loop runs (default `all`; read with `parse_known_args`, so it still runs as a notebook). `jarvis` only builds `jarvis26.pickle` from the real Jarvis+26 table, no simulation involved (see "the real Jarvis+26 GD-1 data"). the supplemental loop does `--sims` (default: every `simspect.supplemental_sims` entry) x `--copies` (default `0`), each with its own seed-42 generator. copy 0 is pickled as `<sim>_<rvir>`, and copy k as `<sim>_copy<k>_<rvir>`. what it does, per orbit x rvir:
- **which sims:** `orbits = ['gd1','aau','pa5','jet','c19']` -- **M3 is dropped**, circ is not included. `stellar_pop='hm'` only (`mass_index = 1`). copies tried in order 0..4, except `aau`, which tries 4..0 to dodge a misbehaving copy. (the `if orbit != 'gd1': continue` used to re-run GD-1 alone is still there, commented out.)
- **progenitor removal:** `N_rtid_boundary = 2.0` for every orbit -- "unbound" = outside 2x the petar tidal radius, to make sure the progenitor is fully removed and stars that get recaptured as rtid re-expands aren't counted. (an orbital-phase-dependent 1x/2x choice was tried and abandoned. the orbital phase is still computed, though no longer printed, and it's computed at the wrong time -- see the bugs list.)
- **coordinates:** `coords_obs` (CoM `SkyCoord`), then two straightened versions: `sc` from CoM, and `sc_primaries` from the luminous component with `PM_treatment='CoM'` hardcoded -- i.e. binaries get the **primary's** instantaneous RV/phi2/distance but the **CoM** proper motions. both are then `poly_straightening`-ed. (`streamframe_coords_observed` is called twice in a row with the same arguments, which is harmless but redundant.)
- **trim:** `trim_new = trim_obstream_percentile(sc, p=[1,99], trim_keys=[phi1, d_phi2, v_phi1, v_phi2, v_gsr, distance])` -- a 1-99 percentile clip in *every* observed-frame dimension, AND-ed, computed on the **CoM** coords. `trim_new_primaries` is the same clip computed on `sc_primaries`, i.e. **including binary orbital motion** in phi2 / v_gsr / distance. the trim operates in the space the GMM will see, so `CoM` cases use `trim_new` and `binaries` cases use `trim_new_primaries`, both noiseless and noisy -- which is why the catalog photometry is done once per trim (see that subsection). `sc_primaries` is still poly-straightened with `tc=trim_new`, not `trim_new_primaries`. this replaces both the intrinsic-frame `inMW`/`trim` and the hard `outlier_clip` cuts; `inMW_na` is an all-True placeholder so the two-mask functions still work. the trim only picks which stars `poly_straightening` fits to -- everything is kept in the saved arrays.
- **isochrone photometry** (`phot`): [Fe/H] = -2 MIST v1.2 isochrone from `artpop.fetch_mist_iso_cmd(log_age=log10(12e9), ..., phot_system='UBVRIplus')` -- which actually returns the **12.6 Gyr (log age 10.10), v/vcrit = 0.4** grid isochrone, see the photometry section -- interpolated at `lumdict['m0_zams']` (the luminous component's birth mass for binaries -- the companion's light is not added). absolute G, BP-RP, z (via RTN-099), apparent `mG`/`mz` using the CoM distances, and `log_Teff`. still computed and saved, and it sets the **ranking** for the catalog matching below, but no `gmm.py` case fits it any more.
- **isochrone noise** (`noise`): Gaia DR3 end-of-mission position and PM errors from `pygaia` (total / sqrt(2) per component, in **µas and µas/yr**); DESI RV errors from `mz`; Via RV errors for 1 hr/1 exposure and 10 hr/10 exposures. saved for comparison; unused by `gmm.py`. all three noise sets come from `survey_noise(mG, mz, log_Teff, rng) -> (noise_dict, acceptable_G_viamock)`, defined at the top of the file.
- **catalog photometry** (`catalog_photometry[_primaries]`) + **catalog noise** (`noise_catalog_photometry[_primaries]`): the live path for every noisy case, done once per trim. see the subsection below.
- noise realizations are drawn here, once, from one `default_rng(seed=0)` for the grid (the supplemental loop uses 42, one fresh generator per sim; see the seed TODO), and saved, so every `gmm.py` case sees the same draw. order on the generator, per sim: isochrone noise, then CoM-trim catalog noise, then `_primaries` catalog noise (and `survey_noise` itself draws via, via10hr, desi, pm_phi1, pm_phi2, phi1, phi2 in that order). so changing anything upstream changes every realization after it, and the generator carries over between grid sims.
- **output:** `/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/<orbit>_<rvir:.2f>.pickle` (the grid set was last written 4 Nov 2026). note the filename has no stellar_pop in it, so an `lm` run would overwrite the `hm` pickles. and it's netscratch, so the 90-day purge applies.

every key in the pickle (with lengths, units and the ordering gotchas) is in the "what's in a data_dict" section near the top.

### catalog photometry (the live noisy path)
per sim, after the isochrone photometry, and **twice**: once with `trim_new` (keys as named) and once with `trim_new_primaries` (keys suffixed `_primaries`). it's a loop over `[(trim_new, ''), (trim_new_primaries, '_primaries')]`. the two pools differ slightly, so the catalog rows go to slightly different N-body stars, but the sample size is `len(catalog)` either way.

1. `cut = unbound & trim & nonrem & ontrack_velocity` -- the pool, so the catalog isn't spent on progenitor stars, remnants, trimmed outliers, or stars more than 30 km/s off the v_gsr track (to match Jarvis+26; `ontrack_velocity{sfx}` uses the CoM v_gsr for `''` and the primary's for `_primaries`, so it mostly removes high-amplitude binaries from the binaries cases). every noisy case inherits it through `cut`. noiseless cases (`trim_new & unbound`) do not.
2. `noise.assign_photometry_from_catalog(iso_phot[cut], distances[cut], catalog=noise.load_gaia_catalog(orbit), track=noise.isochrone_cmd_track(isocmd), N=len(catalog))`. subsetting to `cut` beforehand does the job of the `eligible` kwarg. the brightest `N` pool stars in **isochrone** apparent G get the catalog's apparent G and BP-RP exactly, rank for rank; absolute mags come from the N-body CoM distances, `log_Teff` from the catalog BP-RP via the colour-Teff fit (`Teff_from_gaia_color_temp_relation`), and `z`/`mz` from RTN-099.
3. the results are scattered back to full length (NaN elsewhere), and DESI / Via 1 hr / Via 10 hr RV errors and Gaia DR3 position/PM errors are recomputed from the new `mG`/`mz`/`log_Teff` with the same `survey_noise` call as the isochrone noise.

**every catalog row gets used, so the noisy sample *is* the catalog.** the pool is always much bigger than the catalog (4500-11700 stars), and the downstream noise cuts in `gmm.py` (Gaia pm_err < 0.5 mas/yr, Via rverr < 5, DESI rverr < 10) remove essentially nothing -- the brightest catalog G is ~12 and the faintest ~20. the exception is the few matched stars outside the colour-Teff fit's colour range, which get NaN `log_Teff` and so NaN Via errors and fail the Via `good_RV` (gd1 1, aau 17, jet 1; see the colour -> Teff subsection). so the fitted sample size is set by `len(catalog)` and doesn't change with rvir or with via vs via10hr (those two fit the **same stars** with different RV noise):

| orbit | catalog rows = N fit (`via`, `via10hr`) | `desi` N fit |
|---|---|---|
| gd1 | 1578 | 679 |
| aau | 334 | (not run; 308) |
| pa5 | 129 | (not run; 114) |
| jet | 144 | (not run; 144) |
| c19 | 46 | (not run; 43) |

with 17 free parameters, the c19 / pa5 / jet fits are badly underconstrained, which is where the BIC TODO comes in.

things to know / decide:
- **the ranking uses the isochrone mG, so stars off the isochrone can't be matched.** `rankable` in `assign_photometry_from_catalog` needs a finite isochrone mG, i.e. `alive`. that throws out the pool stars with `m0_zams` above the isochrone's 0.793574 Msun edge -- 200-860 per sim, still `nonrem` (at a 2-4 Gyr dynamical age they're mostly living MS/turnoff stars), and **the ones that should be the brightest in the stream**. so the catalog's brightest rows go to the next-brightest stars down. only matters if the catalog stars' kinematics depend on stellar mass (e.g. mass segregation delaying the escape of the heavier stars). ranking the pool by `m0_zams` plus distance modulus, or ranking off-isochrone stars first, would fix it. the ~1000-1400 pool stars below the isochrone's 0.1 Msun low-mass edge are also unrankable, which doesn't matter, since they're faint.
- `alive` and `acceptable_G` in the pickle are still the **isochrone** versions; the catalog block computes its own viamock range locally and doesn't save it (it's redundant anyway: every catalog mG is inside 5-30).
- the catalog is used in full: no `p_mem` cut, and `bp_rp` isn't dereddened. the catalog MS sits ~0.1 mag redder than the isochrone, and Teff comes from colour alone, so that goes straight in: ~0.02-0.03 dex too cool (see the colour -> Teff subsection). accepted for now.
- the rank matching is global in apparent G, with no dependence on phi1, so the catalog's magnitude distribution along the stream isn't reproduced -- just the overall luminosity function.
- `iso_dist > 3` (catalog star off the isochrone in the CMD): ~267/1578 for gd1, ~60-80/334 for aau, ~105/129 for pa5, ~83/144 for jet, ~15/46 for c19. that's most of the Pal 5 catalog. the colour fit doesn't use G, so being off the isochrone in G doesn't matter for Teff, only being outside the fit's colour range does (NaN, see the colour -> Teff subsection).

## `gmm.py`
functions + a main program to fit a two-component (thin stream + cocoon) gaussian mixture to the pickled data, for one **case** at a time. the functions are documented in the GMM section below.

a case is set by three flags at the top of `__main__`, and they build `case_name = <noise>_<binaries>[_constrained]`:

| flag | options | meaning |
|---|---|---|
| `noise` | `None`, `'via'`, `'via10hr'`, `'desi'` | `None` -> `noiseless`. otherwise adds the saved noise draw and applies that survey's selection (below). the inline comment only lists None/via/desi, but `via10hr` works too since it's just a key suffix |
| `include_binaries` | `True` / `False` | `binaries` uses `sc_straighter_primaries` (binary orbital motion in the RVs); `CoM` uses `sc_straighter`. picks the suffix `sfx = '_primaries' if include_binaries else ''`, which selects the coords, `trim_new{sfx}`, and all four catalog products, for noiseless and noisy cases alike |
| `constrain_widths` | `True` / `False` | adds `_constrained`: stage-2 SLSQP with the sigma-ratio constraint. `False` = Powell only, unconstrained |
| `supplemental` | `True` / `False` | fit the supplemental sims instead of the grid; the table is `<case_name>_supplemental.fits` |
| `jarvis_data` | `True` / `False` | fit the real Jarvis+26 GD-1 data (`jarvis26.pickle`) instead of any sim. forces `noise=None`, `include_binaries=False`; the case is `jarvis26[_constrained]`. see "the real Jarvis+26 GD-1 data" |

the flags can also be set from the command line (`--noise None`, `--include-binaries True`, `--constrain-widths`, `--supplemental`, `--make-plots`, `--jarvis-data` (or `--jarvis_data`), each taking `True`/`False`); the values in the file are the defaults. (`--constrain-widths` and `--supplemental` used to be parsed inverted; see the log, 15 Oct, if an old command-line table looks wrong.) the loop is over a flat list of `(sim, orbit, rvir)` runs, where `sim` names the pickle (== `orbit` for the grid).

- **fit dimensions:** `keys = ['phi2', 'v_phi1', 'v_phi2', 'v_gsr']` -- phi2 in **degrees** (not `d_phi2`), and the proper motions converted to **transverse velocities** (km/s, using the CoM distance) before straightening. with noise, the PM noise is converted to v_phi1/v_phi2 with the same (assumed-perfect) distances.
- **noisy cases read the catalog-based products:** `phot = data_dict['catalog_photometry'+sfx]`, `noise_dict = data_dict['noise_catalog_photometry'+sfx]` (likewise `cut_for_*` and `matched_to_*`). the isochrone-based `phot`/`noise` lines are commented out just above them. switching back needs the old `use` lines too (also commented out), since the isochrone path depends on `alive`/`acceptable_G` and has no `matched` mask.
- **selection:**
  - noiseless: `trim_new{sfx} & unbound` (remnants included).
  - noisy: `use = cut & matched_full & good_pm & good_RV & nonrem`, with `good_pm = pm_err_gaia < 0.5 mas/yr` and
    - `via` **and** `via10hr`: `good_RV = rverr < 5 km/s`.
    - `desi`: `good_RV = (brightest 679 stars in catalog mG) & (rverr < 10 km/s)`. 679 = the length of the Jarvis+26 table. since the DESI errors are Jarvis' own, assigned to those same 679 stars (NaN elsewhere, all < 10), this is identical to `np.isfinite(rverr)`.
    - `matched_full` already implies `cut`, which already contains `nonrem`, so it reduces to `matched_full & good_pm & good_RV` (hence the `# do I ...need nonrem????` -- no).
    - `alive` and `acceptable_G` are loaded but unused.
    - in practice the cuts barely bite: N fit = catalog size, less the few NaN-Teff stars in the Via cases (see the catalog photometry subsection).
  - `desi` fits **GD-1 only**.
  - `jarvis_data`: `p_contam < 0.5` (608 of 679), with `x_err` from the pickle's `errors` dict (its real 1 sd uncertainties).
- **orbits:** same five as `final_datasets.py`; no M3.
- **output:** `data/gmm_tables/<case_name>.fits` with columns `sim`, `orbit`, `Rvir0`, `M_ts`, `S_ts` (thin-stream means/sigmas, 4-vectors in `keys` order), `M_c`, `S_c` (cocoon), `f_cocoon`, `M_1comp`, `S_1comp` (the single-gaussian fit), `BIC_one`, `BIC_two`. (old tables can lack `sim` or the BIC columns; see the log.) sanity plots (points coloured by p_cocoon + histograms) to `/n/netscratch/conroy_lab/Lab/amphillips/p27_sanity_plots/<case_name>/<orbit>_<rvir>.pdf` (`jarvis26/jarvis26.pdf` for the real data) when `make_plots=True` -- those directories have to exist already. the `# %%` cells after the table write are quick-look plots.
- the table doesn't record `keys`, so nothing tells you from the file itself whether column 0 is phi2 or d_phi2 -- the "translate back to angle from distance" TODOs in `results.py` suggest some earlier tables were fit in `d_phi2`. regenerate rather than trust an old table if it matters.

**model choice:** each sim also gets a single-component gaussian fit, and both BICs go in the table. the expectation is that for noisy/down-sampled data some sims will prefer one component, which gets called a **non-detection of a cocoon**. see the model choice section below.

## `results.py`
notebook-style cells making the (nice) plots from the tables and pickles (figure 1 is in `motivation.py`, the single-component Jarvis comparison in `single_component_comparisons_mockobs.py`): the main summary figure, `f_cocoon` / cocoon + thin / one-component sigma_phi2, sigma_vgsr vs rvir0 (or rho/rho_0, `x_is_density`), in two versions picked by `summary_version`: `'orbits'` = every grid orbit (ordered by pericenter), CMs only, no supplemental sims, orbits legend only (`plots/summary_cocoon_ts_properties[_density].pdf`); `'gd1'` = GD-1 only, CMs + binary orbits, plus the supplemental sims (reads the `noiseless_*_supplemental.fits` tables), with the CM/binaries and supplemental legends but no orbits legend (`..._gd1only.pdf`). in it the main-grid CM curves are solid with filled markers and the with-binary-orbits curves dotted with **open** markers (white face, edge width 2: the `open_marker` dict), matching the supplemental points (filled = CMs, open = binary orbits); the thin stream vs cocoon population panels for GD-1 (binary fraction, mean G, mean system mass; see the subsection below), and the demo cocoon-separation panels. imports `gmm` for the membership functions. `okay_mean = |M_c| < 0.5 S_c` (in every dimension) is computed to drop fits whose "cocoon" is displaced off-track, but in the summary figure the filter is commented out, so every row is drawn. TODO: rewrite.

the summary figure has two Jarvis+26 references per panel: a grey band = the published values ± errors (hardcoded `jarvis_*` numbers at the top of the cell), and a **black dashed line** = this model fit to their data (read from `data/gmm_tables/jarvis26.fits`, so it follows a refit).

⚠️ `results.py` still reads the **unsuffixed** catalog keys (`catalog_photometry`, `cut_for_catalog_photometry`, ... around lines 163 and 501, and in `demo_data`), so its noisy `binaries` panels use the CoM-trim photometry, unlike `gmm.py`. give it the same `sfx` switch. (the three population cells below do have the `sfx` switch.)

### thin stream vs cocoon populations: f_bin, <G>, <M_sys>
four cells after the Jarvis comparison figure (`plots/jarvis_comparison.pdf`), GD-1 only, all vs rho_prog/rho_0 with the summary figure's x-axis (`x_axis_values`, log scale, inverted, `2^-n` ticks; follows `x_is_density`). so they need the summary-figure cells run first (`gd1_color`, `x_label`, `x_axis_values`, `density_fracs`).
- **membership:** the `noiseless_CoM` table (`case`), refit-free: `gmm.component_membership_probability` on `unbound & trim_new` (with `sfx = '_primaries'` if `case` is a binaries table), p_cocoon > 0.5 = cocoon.
- **style shared by all of them:** thin stream = `gd1_color`, cocoon = `cocoon_fbin_color = "#EE9998"` (defined in the f_bin cell) with alpha-0.7 marker faces/edges, solid lines; "cocoon"/"thin stream" labels in the curve's colour with solid `->` arrows pointing at the midpoint of the segment between the middle two points (indices 1 and 2; geometric mean in x on the log axis). the label offsets are hand-tuned per panel.

| cell | y | stars | error bars | saves |
|---|---|---|---|---|
| f_bin | binary fraction (binaries = rows `[nsingles:]`) | every fitted star | `paf.Prob_of_frac` + `paf.percentile` | `plots/gd1_binary_fractions.pdf`; dotted line at `fbin0 = 0.216` (from IC generation), ylim 0-0.29 |
| <G> | mean apparent catalog G (`catalog_photometry{sfx}['mG']`) | catalog-matched only (`isfinite(mG)`) | standard error of the mean | `plots/gd1_mean_G.pdf`; y inverted (bright up) |
| <M_sys> | system mass, `data_dict['CoM']['mass']`; toggle `mass_stat`: `'mean_log'` (default) = <log10 M_sys>, the same treatment as <G>; `'log_mean'` = log10 <M_sys> | toggle `mass_sample`: `'matched'` (as <G>) or `'fitted'` (every fitted star, remnants included) | `'mean_log'`: SEM of log10 M_sys. `'log_mean'`: SEM in M propagated to log10, sigma_log = SEM / (<M> ln 10) | `plots/gd1_mean_system_mass_<mass_sample>[_meanlog].pdf` (`_meanlog` for `'mean_log'`) |
| combined | the three above as 1x3, <M_sys> / f_bin / <G> (figsize 24x7, wspace 0.25) | reuses the arrays from the three cells (run them first) | as above | `plots/gd1_combined_msys_fbin_meanG_<mass_sample>_<mass_stat>.pdf` (the older `..._<mass_sample>.pdf` files are `log_mean`) |

- **system mass:** `data_dict['CoM']['mass']` is one row per system, and a binary's row is the combined present-day mass: checked on gd1 0.75 that it equals p1 + p2 from the `companions` dict exactly. present-day, not ZAMS.
- `'mean_log'` is the log of the geometric mean: a typical system, not pulled up much by a few heavy ones. `'log_mean'` is dominated by the heavy tail (the mass budget). the <M_sys> numbers in the table below are `'log_mean'`.
- (`'log_mean'` only) the log10 error bar is first-order, fine for `'matched'` (fractional SEM 1-6%), less so for the `'fitted'` cocoon (up to ~28%, so the true interval is asymmetric, e.g. +0.11/-0.14 vs ±0.12 dex), whose masses are also skewed (probably a few heavy remnants; not checked). asymmetric bars or a bootstrap would fix it if that version is shown.
- the catalog photometry assigns the catalog's G values by rank, so the mean G over all matched stars is the same in every sim (all 1578 rows used); the <G> panel only shows how they split between the components.

the numbers (GD-1, `noiseless_CoM`; f_bin is nbin/n, the plotted Prob_of_frac medians are nearly identical; N in brackets):

| rho/rho_0 | f_bin thin | f_bin cocoon | <G> thin (matched) | <G> cocoon (matched) | <M_sys> thin, matched / fitted [Msun] | <M_sys> cocoon, matched / fitted [Msun] |
|---|---|---|---|---|---|---|
| 1 | 0.117 (5699) | 0.165 (788) | 18.76 ± 0.03 (1415) | 18.55 ± 0.08 (163) | 0.543 / 0.347 | 0.673 / 0.56 ± 0.06 |
| 2^-3 | 0.131 (5156) | 0.187 (417) | 18.75 ± 0.03 (1467) | 18.58 ± 0.11 (111) | 0.520 / 0.343 | 0.617 / 0.80 ± 0.17 |
| 2^-6 | 0.153 (7989) | 0.222 (370) | 18.74 ± 0.03 (1503) | 18.64 ± 0.11 (75) | 0.636 / 0.363 | 0.777 / 0.80 ± 0.23 |
| 2^-9 | 0.182 (11699) | 0.241 (398) | 18.73 ± 0.03 (1528) | 18.87 ± 0.14 (50) | 0.724 / 0.412 | 0.801 / 0.82 ± 0.15 |

**the result:** the systems that end up in the cocoon are more often binaries, more massive, and (at all but the lowest density) brighter, i.e. more observable, than the thin stream's. this agrees with the selection effect pointed out by Weatherford+2025: in mass-segregated clusters, ejections come preferentially from the dense core, so the ejected systems tend to be the more massive ones. both f_bin curves rise toward low density (the initial `fbin0` = 0.216, which the cocoon reaches at the two least dense points). the <G> difference shrinks with decreasing density and flips at 2^-9, where the cocoon has only 50 matched stars. caveats: GD-1 only; the matched sample is selected on the luminous component's brightness, so its masses depend on that selection (probably why the matched thin-stream <M_sys> rises toward low density while its <G> doesn't move); and it's the same story as the binaries selection effect in the GMM section, seen from the dynamics side.

## `motivation.py`
the motivation figure (figure 1): every rvir0 = 0.75 pc, `hm` stream at the present day, galactocentric z vs y, coloured by orbit (pericenter order, same colours as the summary plots). three panels:
1. **particle spray** -- the gala spray on the same orbit, read from `data/particle_sprays/` (see `particle_spray.py` below). no progenitor to remove.
2. **thin stream only** -- N-body stars in the `gmm.py` fit with p_cocoon < 0.5, using the `noiseless_CoM` table (`fig1_noise`, `fig1_include_binaries`).
3. **all stars** -- every unbound star (outside 2 rtid). `fig1_include_trimmed = True` includes the stars `trim_new` dropped before the fit, so panels 2 vs 3 = cocoon + trim; `False` = only the fitted stars, so the difference is exactly the cocoon.

saves `plots/fig1_motiv_notrim_with_spray.pdf`. imports `case_label`, `gmm_selection`, `observed`, `observed_errors`, `p_cocoon_from_table` from `results.py`, which runs all of `results.py` on import (TODO in the code). it doesn't generate sprays itself: `make_spray` is still defined there but the call is commented out, and a missing spray file prints a warning and `continue`s -- which also skips that orbit's **N-body** stars in panels 2 and 3, not just the spray. run `particle_spray.slurm` first.

## `particle_spray.py` / `particle_spray.slurm`
make the sprays for `motivation.py`, one orbit per run, since they're too slow to make in the plotting loop.

```
python particle_spray.py --orbit gd1 [--overwrite] [--progress]   # one stream
sbatch particle_spray.slurm                                        # all five as a 0-4 array (gd1 aau pa5 jet c19)
```

- matched to the rvir0 = 0.75 pc, `hm` N-body run: the first finished copy (via `data.tidal` reaching the present day; same copy order as `final_datasets.py`, i.e. 4..0 for `aau`), its t = 0 `init_displacement`, the same `BovyMWPotential2014`, integrated for the same `age` (so it includes the +100 Myr).
- Chen DF, 5000 particles at the present day (one release per tail at 2500 evenly spaced steps), dt = 1 Myr, snapshot every 10 Myr so output index `j` lines up with `data.<j>`. the DF's progenitor mass follows the N-body **bound mass** from `data.tidal`.
- **no progenitor potential** (`progenitor_potential=None`): the particles feel only the MW. a fixed Plummer at the *initial* mass (what `animations/particle_spray_vs_nbody.py` uses) traps late releases deep inside a too-massive potential once the release radius follows the bound mass, and gala advances every released particle with one shared adaptive DOP853 step (rtol = atol = 1e-10), so a few trapped particles slow all of them down -- worst for Pal 5 (bound mass 9940 -> 195 Msun). see the log, 21 Oct.
- output: `data/particle_sprays/spray_<orbit>_hm_0.75_<copy>.h5` (gitignored), the filename `motivation.py` looks for. existing files are skipped unless `--overwrite` (the slurm script currently passes `--overwrite --progress`).
- the slurm script is the same template as the movie jobs (`petar_env`, conroy/itc_cluster/shared/sapphire, intel), 4G / 3 hr per task, logs to `data/particle_sprays/logs/out.<index>`. runtime without the Plummer hasn't been timed (with it: 3-4 min per stream, Pal 5 > 9 min).

## `single_component_comparisons_mockobs.py`
the multi-epoch / sigma_binary experiment: single-component fits to a DESI-like (or Via-like, `survey`) GD-1 sample in which a fraction f_repeat of the stars have a second RV epoch, with sigma_binary either **fit** as a free parameter (Nibauer+26-like) or **estimated** from the repeats' epoch-to-epoch dv. plus the single-component Jarvis comparison figure (`plots/jarvis_comparison_single_component[_with_errs].pdf`) that plots it. imports `measurement_errors`, `binary_excess_array`, `pack_params`, `unpack_params`, `nll_flat`, `N_jarvis` from `gmm`; the multi-epoch likelihood lives in this script, not in `gmm.py`. the cells after the experiment only read tables, so they can be re-run on their own once the table exists (run the imports/settings cells first). the scan plot's `savefig` is commented out; the Jarvis figure's and the NUTS plots' are live. the scan and plot cells run in ~1-2 min (the NUTS cells take longer, see below).

**which measurement: settings cell.**

| setting | options | what |
|---|---|---|
| `survey` | `'desi'` (in the file), `'via'`, `'via10hr'` | whose RV errors and RV selection the mock survey uses, exactly as `gmm.py`'s noise case of that name: `rverr_<survey>` from `noise_catalog_photometry_primaries` (for the fresh per-epoch draws, `x_err` and the flag test) and the saved `v_gsr_<survey>` draw (one-epoch mode). selection `'desi'`: the `N_jarvis` = 679 brightest in catalog G with rverr < 10; `'via'` / `'via10hr'`: rverr < 5 (the whole catalog, so the same 1577 GD-1 stars for both; the one NaN-Teff star drops out). good Gaia PMs and catalog-matched either way. sets the table/plot names (`<survey>_noise_nominal_multiepoch_scan_<tag>`), the `SURVEY` header, the labels, and which `gmm.py` tables (`<survey>_noise_{binaries,CoM}`) the Jarvis figure reads (there's no `via10hr` one, so that figure fails for `'via10hr'`) |
| `sb_mode` | `'fit'`, `'dv'`, `'avg'`, `'flag'` (in the file) | `'fit'`: sigma_binary a free parameter of the two-epoch likelihood (`fit_multiepoch`), jointly with sigma_stream. `'dv'`: sigma_binary estimated from the repeats' dv (`estimate_sigma_binary_dv`), then one gaussian on epoch 0 with it in quadrature on every star (the "old way"). `'avg'`: no sigma_binary; repeat stars' velocity = the mean of their two epochs (error rverr/sqrt 2, `average_epochs`), non-repeats epoch 0; one gaussian. `'flag'`: no sigma_binary; repeat stars that `paf.get_detections` calls RV variable (chi^2 of the epochs about their mean, p < 0.1, hardcoded there) are removed (`flag_rv_variable`), one gaussian on the rest. `'avg'` and `'flag'` need `multi_epoch = True` |
| `flag_average` | `True` (in the file), `False` | `'flag'` only: unflagged repeat stars use their two-epoch mean (`True`) or epoch 0 (`False`); non-repeats always epoch 0 (they can't be flagged) |
| `run_scan` | `True`, `False` (in the file) | top of the scan cell: `False` skips the scan and the plot cells just read the current settings' table |
| `shared_stream` | `True`, `False` (in the file) | `'fit'` only: the shared stream term, or every epoch independent (the Nibauer+26 form, degenerate) |
| `sb_dv_method` | `'deconvolved'` (in the file), `'raw'`, `'noiseless'` | `'dv'` only: which dv estimator. `'noiseless'` only with `multi_epoch = False` (raises otherwise): std of the true binary motion `v_pri - v_com` over the repeats, no noise anywhere |
| `multi_epoch` | `False`, `True` (in the file) | `True`: two mock epochs on Keplerian RV curves; `'dv'` uses dv = epoch 1 - epoch 0 of the repeats. `False`: **one** observed epoch = the sim's instantaneous primary v_gsr + the saved `survey` noise draw (`v_pri + v_noise`); `'dv'` uses dv = observed epoch - **true, noiseless CoM v_gsr** of the repeats (`estimate_sigma_binary_com`), as if the repeats had pinned down each CoM velocity exactly. `'fit'` raises with `multi_epoch=False` (nothing to fit sigma_binary from) |

the file currently has `survey = 'via'`, `multi_epoch = True`, `dt_mean_yr = 5.0`, `dt_sd_yr = 2/12`, `sb_mode = 'flag'`, `shared_stream = True`, `sb_dv_method = 'deconvolved'`, `flag_average = True`, `run_scan = False`. (every table on disk is at the old 1 yr baseline, so with these settings the plot cells have no table to read until the scan is rerun.)

**the epoch baseline is in the table name.** with `multi_epoch = True` the name ends in `dt_tag` = `_dt<dt_mean_yr>yr_sd<dt_sd_yr>yr`, both rounded to 2 decimals and printed with `%g` (1 yr, 2 months -> `_dt1yr_sd0.17yr`), so scans at different baselines don't overwrite each other and the plot cells / Jarvis figure read the one matching the current `dt_*`. the exact values are in the `DTMEANYR` / `DTSDYR` headers. one-epoch (`dvcom_*`) tables have no dt tag, since dt does nothing there. the commented-out scan-plot `savefig` has the tag too.

each combination writes its own table, `<survey>_noise_nominal_multiepoch_scan_<tag><dt_tag>.fits`, with `<tag>` = `fit_shared`, `fit_indep`, `dv_deconvolved`, `dv_raw` (two epochs), `dvcom_deconvolved`, `dvcom_raw`, `dvcom_noiseless` (one epoch, dv vs the CoM), `avg`, `flag_avg`, `flag_ep0` (`scan_tag`). for `'avg'` / `'flag'` the `sigma_binary*` columns are NaN, `S_all` is the two-epoch mean of every star at its CoM v, and there's `N_fit_mean` (stars in the fit); `'flag'` also has `n_flagged_mean`, `n_flagged_bin_mean`, `n_flagged_single_mean` (flagged repeats; true binaries; singles, i.e. false positives) and headers `FLAGTEST`, `FLAGAVG`. the plot cells and the Jarvis figure read the table of the current settings (`scan_name`), and the measurement is in the same columns in every mode, so they work unchanged. all ten DESI tables exist at dt = 1 yr ± 2 months (Via: `avg`, `flag_avg` only).
- **the sample:** GD-1, `survey` noise, binaries case, **the nominal sample of that survey**, i.e. exactly `gmm.py`'s `use` mask for `noise = survey` with `include_binaries=True` (matched, good PMs, nonrem, and `'desi'`: the `N_jarvis` = 679 brightest in catalog mG with rverr < 10, 679 stars at every rvir0; `'via'` / `'via10hr'`: rverr < 5, 1577 stars at 0.75 pc, median rverr 0.49 / 0.18 km/s against DESI's 3.1). the selection is still the primaries one, made with the sim's instantaneous primary v_gsr (`ontrack_velocity_primaries`, `trim_new_primaries`), even though the mock RVs below replace that binary motion with a synthetic one.
- **the mock RVs** (`load_mock_sample(rvir, survey, multi_epoch=True, rng=None, dt_mean_yr, dt_sd_yr)`; `multi_epoch=False` returns only the single-snapshot keys): every binary gets a Keplerian RV curve about its **noiseless CoM v_gsr** (`sc_straighter['v_gsr']`), from `paf.get_orbital_params_from_binaries` (petar `semi`, `ecc`, masses; isotropic inclination, uniform periapsis argument and phase) and `paf.get_rvs`. singles sit at their CoM v_gsr. two epochs, t = 0 and Delta t ~ N(`dt_mean_yr`, `dt_sd_yr`) (every result below is at 1 yr, 2 months). each epoch gets a **fresh** noise draw, N(0, `rverr_<survey>`), not the saved one. returns `t_obs`, `v_epochs_true`, `v_epochs_noise`, `v_epochs_obs` (all (N_obs, 2), days / km/s), `rverr`, `orb_params` (the sample's binaries, `[v0, K, w, phi0, e, P/day]`), plus the old single-snapshot keys.
  - the petar `Binary` object comes from `load_present_day_binaries(dd, rvir)`: `dd['binaries']` if the pickle has it (none do yet), else the first grid copy whose present-day `data.<age/10>.binary` has the same binaries as the pickle (same count, same p1 ids in order; that's copy 0 for all four GD-1 rvir0). row j of the pickle's binaries is `binaries[j - nsingles]`.
  - `calc_amplitude` (and so K) is the orbit of the **heavier** star, `a1 = (M_light/M_tot) a`, whatever the luminous component is.
- **the likelihood** (`multiepoch_nll_per_star`, `nll_multiepoch`, `fit_multiepoch`): star i's epochs are one draw from the stream, shared by its epochs, plus independent binary motion and noise per epoch:

      v_i ~ MVN(mu 1, C_i),   C_i = diag(err_ik^2 + sigma_binary^2) + sigma_stream^2 1 1^T

  evaluated with Sherman-Morrison, any number of epochs. one epoch reduces to N(mu, err^2 + sigma_binary^2 + sigma_stream^2). a star without a second epoch gets `missing_epoch_err` = 1e4 km/s on it (`epoch_arrays(v_epochs, rverr, repeat)` builds the (N, 2) arrays), so it drops out. the within-star scatter (epoch differences) constrains sigma_binary; the between-star scatter of the epoch means is sigma_stream^2 + (sigma_binary^2 + err^2)/2. sigma_binary is applied to **every** star. the full derivation (why sigma_stream^2 is on the off-diagonals, Sherman-Morrison + determinant lemma, the two-epoch form, the special cases) is in `memos/multiepoch_likelihood/`.
  - **`shared_stream`** (with `sb_mode = 'fit'`) picks the likelihood. `True` = the shared stream term above. `False` = every epoch independent, prod_k N(v_ik | mu, err_ik^2 + sigma_binary^2 + sigma_stream^2), the form of Nibauer+26's Eq. 1. (Nibauer+26 put half-normal priors with scale 10 km/s on both sigma_binary and sigma_stream and sample with NUTS; with this likelihood the prior is the only thing that can split the total width. see "what this means for Nibauer+26" below.)
  - **in the independent form sigma_binary and sigma_stream are exactly degenerate**: the likelihood only depends on sigma_binary^2 + sigma_stream^2. checked on the 3 pc sample (f_repeat 0.5): the nll is identical to 1e-6 at every split of the best-fit total (5.67 km/s), from (sigma_stream, sigma_binary) = (0, 5.67) to (5.67, 0). in the scan the split is wherever the optimizer stops, so the `fit_indep` table's `S_vgsr` / `sigma_binary` jump around between f_repeat and draws, while their quadrature sum equals `S_vgsr_nosb` (e.g. 1.5 pc, f_repeat 1: hypot(2.30, 4.05) = 4.66 = `S_vgsr_nosb`). so the `fit_indep` table measures only the total width. the same nll evaluated with the shared term does prefer a split (best of those five at 4.9 / 2.8).
  - **what this means for Nibauer+26** (scripts in `claude_diagnostics/`, no simulation data involved). in the (sigma_stream, sigma_binary) plane the likelihood contours are circles about the origin, so the allowed region is a thin quarter-circle of radius s (the total width), and the posterior along it is the prior's. two half-normals with the **same** scale tau multiply to exp(-s^2/2tau^2), which is rotationally symmetric, so with the polar Jacobian the angle along the arc is **exactly uniform, whatever tau is**; and tau = 10 or 30 km/s is >> s anyway, so the prior is ~flat there. result: identical marginals for the two sigmas, median s/sqrt 2, 16-84% ~[0.25 s, 0.97 s]. their 10 -> 30 km/s test only changes the radial part, which the data fix, so it can't show the dependence; a prior whose shape matters on the scale of s does (log-uniform piles up at the edges; a 1 km/s half-normal on sigma_stream returns ~that prior). (`prior_degeneracy.py`.)
    - **their setup** (Sec. 2 and App. A of the paper, `nibauer_gd1.pdf` in the repo root, gitignored): 160 GD-1 stars (Starkman+25 membership > 0.75) with 195 RVs from DESI, SDSS DR9, LAMOST DR8 and MMT Hectochelle; median errors 0.8 (MMT) to 6.4 km/s (SDSS). the "repeats" are **cross-survey**, not repeat epochs: 34 stars observed by two instruments, one by three, so the baselines span up to ~15-20 yr and sigma_binary is really star-level survey-to-survey scatter (binarity, but also per-star systematics and underestimated errors; they take the reported errors as exact). Δv is relative to the **unperturbed model** track, kept within |Δv| < 30 km/s. three phi1 bins, [-60,-40], [-40,-20] (has the spur and most of the MMT stars), [-20,0] deg, with 35 / 74 / 51 stars. free parameters: one sigma_Δvr per bin, a mean μ_s per survey per bin (zero points, N(0, 10 km/s) priors), and **one global sigma_binary**. their Eq. 1 is exactly the independent form, prod_s N(Δv | μ_s, σ_s,i^2 + sigma_binary^2 + sigma_Δvr^2), with σ_s,i -> 1e4 km/s for an unobserved survey. half-normal priors, scale 10 km/s (30 km/s "consistent"), NUTS in numpyro.
    - so the repeats can't separate the two sigmas: under Eq. 1 the 35 repeat stars count the same as 35 more stars. the paper says the parameter is degenerate "though we still produce a constraint on σ_Δvr thanks to the repeat observations of 35 stars", which Eq. 1 doesn't support. their result: sigma_binary = 2.6 +1.2/-1.4, sigma_Δvr = 1.7 ± 1.2, 4.8 +0.7/-1.2, 3.7 +1.2/-1.8 km/s (16-84%). no corner plot; Fig. 10 shows the 1D sigma_Δvr posteriors.
    - the global sigma_binary still leaves one degenerate direction (sigma_binary up, every sigma_stream,k down in quadrature), now capped at sigma_binary = min_k s_k; the coldest region's arc runs parallel to the sigma_binary axis near the cap, so the posterior piles up toward it. the expected likelihood with **no binary information**, their bin sizes (x 195/160 measurements per star) and the region totals their medians imply (3.1, 5.5, 4.5 km/s) gives, for 4 km/s errors (a guess; 2-3 km/s is about as close), sigma_binary 2.76 +1.15/-1.60 and sigma_stream 1.82 +1.38/-1.24, 4.75 +0.91/-1.23, 3.49 +1.16/-1.78 (`fig10_compare.py`). the medians match by construction, the error bars don't have to. the shapes match their Fig. 10 too: the coldest bin flat from 0 to ~2.5 km/s then falling, the right bin with a tail down to 0. (`three_regions.py` is the first version, N = 50 per region and 2 km/s errors: 2.62 +0.81/-1.48; 1.80, 4.93, 3.75; the same for tau = 30.)
    - so what their data do constrain: each region's total width, the differences sigma_stream,j^2 - sigma_stream,k^2, and bounds sigma_binary <= s_min, sigma_stream,k in [sqrt(s_k^2 - s_min^2), s_k] (the middle bin is ~4.5-5.5 km/s; "middle hotter than left" is robust). not constrained: sigma_binary and the coldest region's dispersion. the tell: going from 50 to 150 stars per region tightens the hot regions but not sigma_binary or region 1.
    - the lower bound on the hot bins rests on the binary scatter being the same in every bin. binary dv is heavy-tailed (see the `fit_shared` result below), so in bins of 35-74 stars the realised binary variance can differ a lot between bins, and a couple of fast binaries inside the 30 km/s cut would read as stream heating. not yet quantified; the sims could give the bin-to-bin spread of the realised binary variance for random 35-74 star subsamples.
    - the hot bins come from the velocity tails: with |Δv| < 15 km/s all three bins drop to ~2-3 km/s (their Fig. 10, bottom), consistent with Gialluca+21's 2.3 ± 0.3 (7 km/s cut). they keep the 30 km/s cut because the tails constrain the dark matter models most. their unperturbed model dispersions are ≲ 1 km/s (Fig. 1).
  - `fit_multiepoch(v, v_err, sigma_binary=None, shared_stream=True)` -> `mu, sigma_stream, sigma_binary, result`. `sigma_binary=None` fits it, a number holds it fixed. theta = [mu, ln sigma_stream(, ln sigma_binary)], ln sigma bounded to [-10, 6]. runs **L-BFGS-B and Nelder-Mead** and keeps the lower nll. **not Powell**: its first line search throws ln sigma_stream onto the lower bound, where the likelihood is flat, and it stays there (nll 8 worse than the truth on synthetic data).
  - checked: matches `scipy.stats.multivariate_normal` per star; on synthetic data with gaussian binary motion (N = 3000, sigma_stream 0.8, sigma_binary 2.0) it recovers 0.81 / 1.96 with every star repeated and 0.89 / 1.94 at f_repeat = 0.2.
- `fit_one_component(x, xe, sb=None)`: single-epoch one gaussian via `gmm.nll_flat`, used for phi2, the noiseless truth, and the "old way" below (`sb` = a fixed sigma_binary array from `gmm.binary_excess_array`).
- **the "old way", sigma_binary estimated, not fit:** `estimate_sigma_binary_dv(v_epochs, rverr, repeat, method)` takes the repeat stars' observed dv = v_1 - v_0 (binaries: two noised RV-curve points; singles: two noise draws). `'deconvolved'` (default) = sqrt((var(dv) - mean(2 rverr^2))/2): noise subtracted, and /2 for a difference of two epochs, so it's the per-epoch binary term alone. `'raw'` = std(dv), which still has sqrt(2) x the noise in it and double counts the errors on top of `x_err`. then one gaussian in (phi2, v_gsr) on **epoch 0 of every star**, with `x_err` and that sigma_binary in quadrature on **every** star (nobody's binary motion is removed). this is `sb_mode = 'dv'`, with `sb_dv_method` picking the estimator. with `multi_epoch = False` it's `estimate_sigma_binary_com(v_obs, v_com, rverr, repeat, method)` instead: dv = the one observed epoch (`v_pri + v_noise`) - the true CoM v; `'deconvolved'` = sqrt(var(dv) - mean(rverr^2)) (one noise term: no factor 2, no /2), `'raw'` = std(dv), `'noiseless'` = std(`v_true - v_com`) with `v_true` = `v_pri`, the epoch without noise (the function's `v_true` kwarg); the fit is then on that observed epoch (with noise) for every star, sigma_binary on every star. also an **oracle**, computed in every mode: the same fit with sigma_binary = the true rms binary motion at epoch 0 over the sample, the best a gaussian-in-quadrature sigma_binary can do.
- **the scan:** repeat fraction `f_repeat` = 0.2, 0.3, ..., 1.0 (`f_repeat_nominal = 0.2` ~ Nibauer+26's 35/160), N_repeat = round(f_repeat x 679). the RV curves and noise are drawn once per rvir0 (`rng_exp`, seed 1); then for each of `n_draws` (currently 10) draws, one random permutation is shared by every f_repeat, so the repeat samples are **nested** (each step adds 10% of the stars). the draws don't depend on the mode, so the four tables are on identical mock data. each draw makes the `sb_mode` measurement (for `'fit'` also the same fit with sigma_binary = 0).
- output: `data/gmm_tables/<survey>_noise_nominal_multiepoch_scan_<tag><dt_tag>.fits`, one row per (rvir, f_repeat): `N_obs`, `f_repeat`, `N_repeat`, `n_draws`, `f_bin`, `n_bin_repeat_mean`, `S_true` / `S_none` / `S_all` (2-vectors in (phi2, v_gsr) order: noiseless CoM; DESI epoch 0 only with binary motion, sigma_binary = 0; DESI at CoM v with sigma_binary = 0, the noise-limited target: both epochs for `'fit'`, epoch 0 for `'dv'`), `S_phi2`, **`S_vgsr[_p16,_p84]`** (the stream sigma_v_gsr of the `sb_mode` measurement) and **`sigma_binary[_p16,_p84]`** (fitted, or estimated from dv), `S_vgsr_nosb[_p16,_p84]` (`'fit'` only, sigma_binary = 0), `sigma_binary_rms` (rms binary motion over every star and epoch), `sigma_binary_rms_diff` (rms(dv_1 - dv_0)/sqrt 2, the part a second epoch sees), `sigma_binary_oracle` / `S_vgsr_sboracle`. header keys `KEYS`, `NESTED`, `SAMPLE`, `DTMEANYR`, `DTSDYR`, `MISSERR`, `SBMODE`, and `LIKELIHD` (`'fit'`) or `SBDVMETH` (`'dv'`). plots: sigma_v_gsr vs f_repeat, one panel per rvir0 (`plots/gd1_<survey>_nominal_multiepoch_scan_<tag><dt_tag>.pdf`, savefig commented out); sigma_binary vs f_repeat.
- superseded tables, same numbers under the older column names, nothing reads them: `desi_noise_nominal_multiepoch_scan.fits` (fit_shared + both dv estimators in one table) and `desi_noise_nominal_multiepoch_scan_indep.fits` (fit_indep).
- the Jarvis comparison figure (`plots/jarvis_comparison_single_component[_with_errs].pdf`, savefig live) draws `S_vgsr` from the current settings' table at the `f_repeat_plot` rows (currently **1.0**), labelled "observational binary treatment". `gd1_color` is recomputed with `results.py`'s orbital-phase ordering (`#E84855`). the annotation offsets are hand-tuned. `fit_shared` and both `dv` tables are ~0 (below), so with those settings that series sits on the x axis.
  - **the Nibauer+26 series** (`show_nibauer = True`): filled diamonds (black edge), always with error bars, labelled "Nibauer+26 method (NUTS, indep. epochs)". y = the pooled NUTS sigma_stream (`nuts_shared_stream = False`), row `nibauer_survey` / `nibauer_baseline` (`'desi'` / `'5yr'`, the only ones in those tables). two versions:
    - `nibauer_focus = True` (in the file): **one point, rvir0 = `nibauer_focus_rvir` = 3 pc**, from `nibauer_focus_table` (`jake_nuts_N160_rep35_indep_capstd_ndraw50_varystars_pooled`: 50 draws of both the 160 stars, from the full catalog, and the 35 repeats). no line. raises if the table has no such row.
    - `nibauer_focus = False`: every rvir0 from `nibauer_table` (`jake_nuts_N160_rep35_indep_capstd_ndraw10_pooled`: 160 brightest, 10 draws of the repeats), joined by a dash-dot line.
    - error bar per `nibauer_err`: `'pooled'` (in the file) = 16-84% of the posterior pooled over the draws (every draw's 8000 samples concatenated, equal weight), so it has both each draw's posterior width and the scatter between draws; `'draws'` = 16-84% of the 50 per-draw medians, the between-draw scatter only. the point is the pooled median either way (the median of the per-draw medians, `sigma_stream_draw_med`, isn't plotted).
  - **the cell after it: per-draw NUTS medians**, sigma_stream vs sigma_binary, one point per draw of `nibauer_focus_table`'s per-draw table (same survey / baseline / rvir0 settings, so run the Jarvis cell first), with that draw's marginal 16-84% as x and y error bars, coloured by the number of binaries among the draw's repeats. the pooled median / 16-84% as a big square; truths (median over draws of the 160-star noiseless CoM dispersion, horizontal; of the rms true binary motion, vertical) as dashed lines with their 16-84% over draws as bands; 1:1 dotted. its `savefig` (`plots/gd1_<survey>_<per-draw table>_rvir<X>_<baseline>_medians.pdf`) is commented out. **not the joint posterior** (that's the next cell): it uses only the tables' per-draw summaries, so the correlation within a draw isn't in it (with independent epochs each draw's posterior is a thin arc, sigma_stream^2 + sigma_binary^2 ~ const). what it shows: in the rvir0 = 3, 50-draw run the two medians are equal in every draw (sigma_stream / sigma_binary = 1.01, 16-84% 0.98-1.08; correlation 0.88), so the points lie along 1:1 at a distance that varies between draws (sqrt(sigma_stream^2 + sigma_binary^2) of the medians 1.65-5.55 km/s). i.e. the likelihood fixes the total width, the prior splits it evenly, and the between-draw spread of sigma_stream is mostly the spread of the total width (which binaries got sampled).
  - **the cell after that: joint and marginal posteriors** from the saved samples (`data/nuts_samples/<per-draw table>_<survey>_rvir<X>_<baseline>.npz`, same `nibauer_*` settings). corner layout (9x9, 4:1 panels): the pooled joint posterior as a grey 2D histogram (80 bins, `n_bins_joint`), with its 68% / 95% highest-density regions shaded and outlined in the GD-1 colour (computed on the histogram smoothed by `contour_smooth` = 2.5 bins); each draw's own 68% contour faint blue (`show_draw_contours`); truths dashed, 1:1 dotted. marginals on top / right: every draw's (faint blue steps) and the pooled one (filled), all unit-area densities, axis clipped at `marg_ymax` = 2.5x the pooled peak; pooled median solid, 16-84% dotted. `savefig` (`..._joint.pdf`) commented out, as in the medians cell. what it shows (rvir0 = 3, indep): **every draw's posterior is a quarter-circle arc**, sigma_stream^2 + sigma_binary^2 ≈ its total width, with radii from ~1 to ~6 km/s between draws (pooled radius 3.67 [1.68, 5.54]; the truth point (3.96, 0.50) sits at radius 4.0). the pooled joint is a superposition of those arcs: its 68% region is everything inside radius ~4.2 km/s (minus a small hole at the origin) plus a patch at 5-5.6 km/s, and it contains the truth only because the arcs pass through every split. the marginals are broad and nearly flat out to ~4 km/s, with the sigma_binary ones cut off at each draw's `sb_cap`. the dark band at radius ~1 km/s in the grey image is mostly geometry: every draw puts 8000 samples on an arc of length ~(pi/2) r with a similar radial width (0.5-0.8 km/s), so samples per unit area go as p(r)/r, and the 8 draws with totals 0.8-1.45 km/s (14% of the samples) look ~5x darker than the same number at r ~ 5-6. those 8 picked 160 stars with slow binaries (true rms 1.3-3.4 vs 4-6 km/s for the largest totals; the fitted total tracks the true rms binary motion, correlation 0.77, not the number of binaries, 0.24), and they also undershoot their own true total sqrt(S_true^2 + rms^2) (fitted/true 0.33-0.84 vs a median 0.93 over all draws; not checked why).
  - **the cell after that: the same plot for one draw**, i.e. one Nibauer+26-like measurement. reads draw `single_draw` from the same `.npz` (no new sampling: each draw there is its own fit), so run the cell above first. `single_draw = 'median'` (in the file) = the draw whose median total width is the median over draws (draw 8; totals range ~1-6 km/s), or an int. adds this draw's truths (dashed), its true total width as a grey dashed arc, and its `sb_cap` (dash-dot, top marginal). `savefig` (`..._draw<k>_joint.pdf`) commented out. draw 8 (59 binaries, 10 among the repeats): sigma_stream 2.59 +0.88/-1.65, sigma_binary 2.44 +0.99/-1.75 km/s; truth 0.51 / 4.92. the posterior is one arc of radius ~3.5 km/s, inside the true total of 4.95, and the cap (4.51) doesn't bind. both marginals peak near the arc's ends (sigma_stream ~3.5, sigma_binary ~3.3): that's the projection of a roughly uniform-in-angle arc onto each axis (density ∝ 1/sqrt(r^2 - x^2), peaked at x ≈ r), not information about the split. the measurement has no stream width, so **x = the binary-CoM single-component sigma_phi2** (`<survey>_noise_CoM`'s `S_1comp`) at the same rvir0. with `survey = 'via'` that's the Via sample's width against a DESI-sample sigma_stream (the two widths differ by ≤ 0.02 deg).

**result (`fit_shared`): with one gaussian sigma_binary on every star, sigma_stream goes to 0** at every rvir0 and every f_repeat (median over draws; at f_repeat ≥ 0.5 the 16-84% range is 0 too, except 0.75 pc). sigma_binary absorbs everything. (km/s; 10 draws.)

| rvir0 | truth | 2 epochs, CoM ("all") | 1 epoch, binaries in | f_bin | rms dv / epoch-to-epoch | sigma_stream at 0.2 / 0.5 / 1.0 | sigma_binary at 0.2 / 0.5 / 1.0 | sigma_stream with sb = 0, at 1.0 |
|---|---|---|---|---|---|---|---|---|
| 0.75 | 0.84 | 0.64 | 2.16 | 0.23 | 3.15 / 3.22 | 0.00 / 0.12 / 0.00 | 2.23 / 2.61 / 2.48 | 1.62 |
| 1.5 | 0.67 | 0.65 | 4.98 | 0.25 | 4.21 / 4.57 | 0.00 / 0.00 / 0.00 | 4.68 / 4.87 / 4.66 | 3.19 |
| 3 | 0.52 | 0.35 | 4.87 | 0.30 | 4.26 / 5.36 | 0.00 / 0.00 / 0.00 | 4.60 / 4.87 / 4.95 | 2.05 |
| 6 | 0.39 | 0.48 | 2.90 | 0.32 | 4.01 / 4.32 | 0.00 / 0.00 / 0.00 | 2.99 / 3.35 / 3.34 | 1.80 |

- why: the binary motion is heavy-tailed. 68-77% of the sample are singles with dv = 0, and a few binaries move by tens of km/s between epochs. the likelihood is dominated by those few large epoch differences, so it fits a sigma_binary of 2-5 km/s, and applying that to **every** star already over-explains the between-star scatter, leaving nothing for the stream. the code is fine: the same fit recovers sigma_stream when the binary motion really is gaussian. it's the same mis-specification as the old in-quadrature approach (see the TODO).
- holding sigma_binary at 0 instead leaves the binaries in the stream width (1.6-3.2 km/s at f_repeat = 1).
- scratch test, **not in the script**: a per-star mixture (a star is a binary with probability f_b and gets C = diag(err^2 + sigma_b^2) + sigma_stream^2 11^T, else diag(err^2) + sigma_stream^2 11^T; f_b free) gives sigma_stream 0.67-0.90 / 0.87-0.95 / 0.60-0.65 / 0.72-0.77 km/s at 0.75 / 1.5 / 3 / 6 pc for f_repeat 0.2-1.0 (truths 0.84 / 0.67 / 0.52 / 0.39), with f_b ~ 0.02-0.06 and sigma_b ~ 15-21 km/s: it finds the few fast binaries and leaves the rest alone. 5 draws, unconverged scatter not quantified.
- "all" differs from the noiseless truth (0.35 vs 0.52 at 3 pc, 0.48 vs 0.39 at 6 pc): the noise draw on 679 stars with median rverr 3.1 km/s.

**the old way (`dv_deconvolved`, `dv_raw`: sigma_binary estimated from the repeats' dv, in quadrature on every star, epoch 0) also gives ~0**, and so does the oracle with the exact sigma_binary (km/s; median over 10 draws; f_repeat 0.2 / 0.5 / 1.0):

| rvir0 | truth | oracle sigma_binary -> sigma_v_gsr | `dv_deconvolved` sigma_binary | `dv_deconvolved` S_vgsr | `dv_raw` sigma_binary | `dv_raw` S_vgsr |
|---|---|---|---|---|---|---|
| 0.75 | 0.84 | 3.42 -> 0.00 | 3.75 / 3.28 / 3.22 | 0.00 / 0.00 / 0.00 | 8.43 / 7.89 / 7.93 | 0 / 0 / 0 |
| 1.5 | 0.67 | 4.24 -> 2.62 | 5.55 / 5.53 / 4.80 | 0.00 / 0.00 / 1.33 | 10.4 / 10.1 / 9.4 | 0 / 0 / 0 |
| 3 | 0.52 | 4.23 -> 2.41 | 3.54 / 5.01 / 5.34 | 3.34 / 0.24 / 0.00 | 8.0 / 9.6 / 10.0 | 0 / 0 / 0 |
| 6 | 0.39 | 3.87 -> 0.00 | 3.72 / 4.50 / 4.03 | 0.00 / 0.00 / 0.00 | 8.0 / 9.0 / 8.6 | 0 / 0 / 0 |

- the deconvolved estimator gets sigma_binary about right (within ~1 km/s of the oracle at f_repeat ≥ 0.5); the raw one is ~2x too big and always gives 0.
- **the oracle failing means it isn't the sigma_binary estimate, it's subtracting binaries in quadrature at all.** the stream variance (0.15-0.7 km^2/s^2) is a few % of the binary term (~10-30 km^2/s^2), and the random cross terms in the realised sample are bigger than it: the noise x binary term 2 <e b> alone scatters by ~2 x 3.1 x 4/sqrt(679) ~ 1 km^2/s^2. the gaussian fit also weights stars by 1/(err^2 + ...), so it measures a weighted width, not the variance sigma_binary describes, and the brighter (smaller-error) stars are the more binary-rich ones. which way it misses (0, or 2.4-2.6 km/s) depends on the draw.
- **one epoch, dv vs the true CoM** (`dvcom_*`, the sim's instantaneous binary motion + the saved noise; f_repeat 0.2 / 0.5 / 1.0): `dvcom_deconvolved` sigma_binary 2.5-3.6 km/s against a true rms of 2.37-3.06, S_vgsr 0 at 0.75 / 3 / 6 pc and 2.03 / 1.44 / 0.65 at 1.5 pc (truth 0.67); the oracle gives 0.00 / 1.64 / 1.53 / 0.00. `dvcom_raw` sigma_binary 5.1-6.0, S_vgsr 0 everywhere. `dvcom_noiseless` (sigma_binary from the true binary motion of the repeats): sigma_binary 2.0-3.2, and at f_repeat = 1 it equals the true rms and S_vgsr equals the oracle exactly (0.00 / 1.64 / 1.54 / 0.00 at 0.75 / 1.5 / 3 / 6 pc, truths 0.84 / 0.67 / 0.52 / 0.39); at f_repeat 0.2 / 0.5, S_vgsr 0.00 / 0.00, 2.17 / 2.01, 1.94 / 1.63, 0.43 / 0.00. so even a perfectly known sigma_binary from the repeats doesn't recover the stream (the oracle point). the deconvolved estimate runs 0.3-0.7 km/s above the true rms at f_repeat = 1 at every rvir0: checked, dv = b + noise exactly, and it's the realised noise variance of the saved draw sitting 0.8-4.4 km^2/s^2 above <rverr^2> (a few rverr near 10 km/s make that sample variance noisy), plus a 2 cov(b, noise) of -0.3 to +1.1. not a construction error. the binary motion here (rms 2.4-3.1) is smaller than on the Keplerian curves (3.2-4.3), because the sample was selected with |v_pri| < 30 km/s.
- so with DESI errors and this binary population, a sub-km/s stream dispersion only comes out of a method that identifies the individual binaries: the per-star mixture (scratch test above), or flagging stars whose two epochs disagree and fitting the rest (`'flag'`, below).

**averaging the epochs (`avg`) and removing RV-variable stars (`flag_avg`, `flag_ep0`)** (stream sigma_v_gsr, km/s, median over 10 draws, f_repeat 0.2 / 0.5 / 1.0):

| rvir0 | truth | `avg` | `flag_avg` | `flag_ep0` | flagged at f_repeat = 1: binaries / all binaries, singles / all singles |
|---|---|---|---|---|---|
| 0.75 | 0.84 | 1.83 / 1.97 / 1.62 | 1.64 / 1.37 / 1.02 | 1.68 / 1.42 / 0.74 | 28 / 156, 44 / 523 |
| 1.5 | 0.67 | 4.01 / 3.66 / 3.19 | 3.63 / 3.01 / 0.88 | 3.70 / 3.15 / 0.58 | 36 / 167, 55 / 512 |
| 3 | 0.52 | 4.56 / 3.61 / 2.05 | 4.59 / 3.48 / 1.32 | 4.67 / 3.69 / 1.23 | 49 / 202, 44 / 477 |
| 6 | 0.39 | 2.62 / 1.97 / 1.80 | 2.54 / 1.65 / 1.09 | 2.58 / 1.66 / 0.61 | 38 / 217, 40 / 462 |

- `avg` is exactly `fit_shared`'s `S_vgsr_nosb` (two-epoch fit, sigma_binary = 0): with equal errors per star, the shared-stream likelihood with no binary term is the same as fitting the epoch mean. averaging over a ~1 yr baseline barely touches the binary motion.
- `flag` is the only method so far that doesn't go to 0, and at f_repeat = 1 it gets within ~0.2-0.8 km/s of the truth. but it's still well above the truth at f_repeat ≤ 0.5 (every unflagged non-repeat binary is still in the sample), and the test only catches 18-24% of the binaries (DESI errors ~3 km/s, p < 0.1), while removing ~8-11% of the singles (the expected 10% false-positive rate).
- **`flag_ep0` is biased low; use `flag_avg`.** the flag keeps stars whose two epochs agree, and for a single star the epoch-0 noise is correlated with the epoch difference (n_0 = mean - diff/2), so the kept stars' epoch-0 noise is smaller than rverr says (by ~20% in variance for a cut at p = 0.1) and the fit over-subtracts it. the two-epoch mean is independent of the difference (equal errors), so `flag_avg` has no such bias. `flag_ep0`'s lower numbers at f_repeat = 1 are that bias, not a better measurement.
- (the previous version of the experiment, perfect CoM correction for the repeats + an estimated sigma_binary in quadrature, is in the log, 15 Nov; its table `desi_noise_nominal_repeatfracscan.fits` is still on disk but nothing reads or writes it.)

### mock catalogs for Jake Nibauer (last cell)
`jake_catalog` + the loop under `### assemble catalogs for jake ###`: three-epoch mock GD-1 RV catalogs to send to Jake Nibauer, so he can run his own measurement on the sims. needs only the imports/settings and function cells (not the scan or plot cells); ~40 s.
- output: `for_jake/rvir0<X>_<desi|via>.fits` (`%g`, so `rvir00.75_desi.fits`, `rvir01.5_via.fits`, ...), one per rvir0 x survey. the truth for each (`is_binary`, `phi2_true`, `v_com`, noiseless `v_true_k`, `orb_params`, `mG`) is `data/for_jake_truth/<same>_truth.fits`, **not for sending**.
- stars: the binaries case's catalog-matched stars (`cut_for_catalog_photometry_primaries` & matched): progenitor (2 rtid), remnants, the 1-99 percentile trim and the |v_gsr,pri| < 30 km/s cut are applied; **no error-based cuts**. `via` = all of them, 1578 = the BPW25 catalog; `desi` = the 679 brightest in catalog G (the Jarvis V_ERRs), a subset of the via stars. one via star per rvir0 has a NaN Via error (NaN `log_Teff`), so NaN `vr_k` and `vr_k_err`.
- epochs: t_0 = 0, t_1 ~ N(1 yr, 2 months), t_2 ~ N(5 yr, 6 months) (`jake_epochs_yr`). singles at their noiseless CoM v_gsr, binaries on Keplerian curves about it (as `load_mock_sample`; orbits drawn once per rvir0 and shared by both surveys); each epoch gets a fresh N(0, rverr_<survey>) draw, same error every epoch. seed `default_rng([jake_seed = 2026, rvir_index])`.
- columns: `phi1`, `phi2` [deg] with the saved Gaia position noise draw (`_err` = pygaia total position uncertainty / sqrt 2), `vr_0..2` + `_err` [km/s], `t_0..2` [yr], `phi2_track` [deg], `vr_track` [km/s]. **`phi2` and `vr_k` are residuals about our track** (as gmm.py fits); add the `_track` columns for absolute GD1Koposov10 phi2 / reflex-corrected v_gsr. phi1 is absolute. header: `RVIR0`, `SURVEY`, `NEPOCH`, `SEED`, `T<k>MEAN`/`T<k>SD`, `COMMENTS`. rows ordered by phi1.

### NUTS on Jake's catalogs (last four cells)
posteriors for sigma_binary and the stream's v_gsr dispersion from the `for_jake/` catalogs, sampled with numpyro's NUTS and Nibauer+26's priors. no stream width. reads only the catalogs and their truth files, so it needs the imports/settings and function cells but not the scan or Jake's catalog cell. ~2.5 min for all 16 rows of the single-draw cell. the file currently has `nuts_shared_stream = False`; tables exist for both settings (`_shared_` and `_indep_` in the name).
- **sample:** per rvir0 x survey, the `nuts_n_obs` = 160 brightest stars in catalog G (`mG` from the truth file) with finite errors, of which `nuts_n_repeat` = 35, chosen at random, have a second epoch (Nibauer+26: 160 stars, 35 repeats). the same 35 stars for both baselines and both surveys (`default_rng([nuts_seed = 7, rvir_index])`). binaries in the 160: 32 / 38 / 43 / 53 at 0.75 / 1.5 / 3 / 6 pc (DESI; Via within 1), and in the 35 repeats 5 / 7 / 12 / 14.
- **baselines** (`nuts_baselines`): every star has `vr_0`; the repeats' second epoch is `vr_1` (`'1yr'`, ~1 yr later) or `vr_2` (`'5yr'`, ~5 yr). one table row each.
- **likelihood:** `multiepoch_nll_per_star` called with `xp=jax.numpy` (the function takes an array-module argument; numpy is the default), on the (N, 2) arrays from `epoch_arrays`, so it's exactly the ML fits' likelihood. `nuts_shared_stream = True` = the shared stream term; `False` = Nibauer+26's Eq. 1 (degenerate, see "what this means for Nibauer+26"). needs `numpyro.enable_x64()`: missing epochs have err^2 = 1e8.
- **priors:** sigma_stream ~ HalfNormal(10 km/s); sigma_binary ~ HalfNormal(10 km/s) **truncated at `sb_cap`**, the repeats' epoch-to-epoch scatter (`sigma_binary_cap`): `nuts_sb_cap = 'std'` (in the file) caps it at std(dv)/sqrt 2, `'max'` at max |dv|/sqrt 2, with dv = second epoch - `vr_0` (observed, noise included). std(dv)/sqrt 2 is the per-epoch scatter a difference of two epochs implies, so under the model it's an upper bound on sigma_binary. mu ~ N(0, 10 km/s).
- sampler: 4 chains x (1000 warmup + 2000 samples), sequential, `jax.random.PRNGKey(7)`.
- output: `data/gmm_tables/jake_nuts_<tag>.fits`, `<tag>` = `N<n_obs>_rep<n_repeat>_<shared|indep>_cap<std|max>`. columns `survey`, `Rvir0`, `baseline`, `epoch2`, `dt_mean_yr`, `dt_sd_yr`, `N_obs`, `N_repeat`, `n_bin`, `n_bin_repeat`, `sb_cap`, `S_true` (std of the noiseless CoM v_gsr over the 160), `sigma_binary_rms` (rms true binary motion over the 160 at the two epochs used), and for each of `sigma_stream`, `sigma_binary`, `mu`: the posterior median, `_p16`, `_p84`, `_err` = (p84 - p16)/2, `_mean`, `_std`; `rhat_max`, `n_divergent`. headers `NOBS`, `NREPEAT`, `PRIOR`, `MUPRIOR`, `SBCAP`, `LIKELIHD`, `SEED`, `NWARMUP`, `NSAMPLE`, `NCHAIN`. plot: `plots/gd1_<nuts_plot_survey>_jake_nuts_<tag>.pdf`, sigma_binary and sigma_stream vs rvir0, median and 16-84%, 1 yr (red circles) vs 5 yr (blue squares), with the truths as crosses.
- every run converged: r_hat ≤ 1.005, one divergent transition in all 16 runs.

result, DESI, `shared`, one draw of the repeats (km/s; posterior median [16, 84]):

| rvir0 | truth sigma_stream | sigma_stream 1 yr | sigma_stream 5 yr | sigma_binary 1 yr (cap; true rms) | sigma_binary 5 yr (cap; true rms) |
|---|---|---|---|---|---|
| 0.75 | 0.81 | 2.39 [2.15, 2.62] | 2.00 [1.67, 2.31] | 0.68 [0.30, 1.06] (1.47; 2.97) | 1.34 [1.00, 1.61] (1.74; 2.43) |
| 1.5 | 0.65 | 4.09 [3.79, 4.41] | 4.05 [3.76, 4.36] | 1.76 [1.45, 2.06] (2.26; 4.74) | 2.03 [1.82, 2.18] (2.25; 4.17) |
| 3 | 0.43 | 0.52 [0.15, 1.09] | 2.64 [2.25, 3.00] | 6.22 [5.91, 6.56] (11.88; 5.05) | 3.46 [3.15, 3.79] (5.29; 4.74) |
| 6 | 0.44 | 0.64 [0.24, 0.99] | 0.42 [0.13, 0.78] | 1.42 [1.22, 1.60] (2.04; 1.72) | 1.89 [1.72, 2.06] (3.23; 2.00) |

- **the posteriors are tight, and wrong wherever the 35 repeats aren't representative of the 160.** sigma_binary is set by the binaries that happen to be among the repeats (5-14 of them); the binary motion of the 125 non-repeats that sigma_binary doesn't cover goes into sigma_stream. at 0.75 / 1.5 pc the repeats under-sample the fast binaries, so sigma_binary is low and sigma_stream is 2-4 km/s (truth 0.65-0.81); at 3 pc with the 1 yr baseline they over-sample them and sigma_stream comes out near the truth; at 6 pc, where the binary motion is small (rms ~2 km/s), both baselines recover it. the 16-84% ranges don't cover this, because the choice of repeats is not in the posterior.
- **the cap barely matters for DESI**: rerun with sb_cap = 1e4 (effectively no cap), every sigma_stream and sigma_binary moves by ≤ 0.15 km/s except 0.75 pc 5 yr (sigma_stream 1.69 [0.75, 2.19], sigma_binary 1.66). for Via it binds: the repeats' dv scatter is small (Via errors ~0.1 km/s for the brightest stars), e.g. 0.75 pc 1 yr cap 0.18 km/s, and the sigma_binary posterior piles up against it.
- Via (same stars, so the same repeats): sigma_stream 3.17 / 3.12, 4.47 / 4.48, 2.97 / 2.81, 0.77 / 0.42 (1 yr / 5 yr) at 0.75 / 1.5 / 3 / 6 pc. no better than DESI, since the error isn't the limitation.
- 1 yr vs 5 yr changes little except at 3 pc. the baseline changes which binaries look variable, but with 35 repeats the draw of repeats matters more.

**over random draws** (the last two cells): the same fit (always 160 stars, 35 of them repeats) for `nuts_n_draws` random draws, for `nuts_draw_rvirs` x `nuts_draw_surveys` x `nuts_draw_baselines`. in the file: **50 draws, rvir0 = 3 pc only** (its one-component width best matches the GD-1 thin stream), `['desi']`, `{'5yr': 2}`, `nuts_draw_vary_stars = True`. draw d uses jax key `nuts_seed + d`. needs the single-draw cell's definitions (run it with its loop commented out); ~1.5 min per 10 fits.
- `nuts_draw_vary_stars = False`: the same 160 brightest stars every draw, only the 35 repeats change. draws come from the single-draw cell's per-rvir0 rng, so draw 0 is that cell's choice of repeats.
- `nuts_draw_vary_stars = True`: each draw picks its 160 stars uniformly at random from the **full catalog** of the survey (`jake_nuts_pool`: every star of Jake's catalog with finite errors, i.e. every catalog-matched star with an RV error for that survey: 679 for DESI, the Jarvis-error stars; 1577 of 1578 for Via), then 35 of those at random as the repeats. so the sample is no longer the brightest 160.
- output: `data/gmm_tables/jake_nuts_<nuts_draw_tag>.fits` (one row per survey x rvir0 x baseline x draw, same columns as above plus `draw`) and `..._pooled.fits` (one row per survey x rvir0 x baseline), with `nuts_draw_tag = <nuts_tag>_ndraw<n>[_varystars]`. pooled = every draw's samples concatenated with equal weight (the posterior marginalized over the draws; `<p>`, `_p16`, `_p84`, `_err`), plus `<p>_draw_med`, `_draw_p16`, `_draw_p84` (spread of the per-draw medians), `n_bin_min` / `_max`, `n_bin_repeat_min` / `_max`, `S_true` and `sigma_binary_rms` (medians over draws: they vary when the stars do), `S_true_pool` (noiseless CoM dispersion over the whole pool), `vary_stars`, `N_pool`, `rhat_max`, `n_divergent` (summed). older pooled tables (`_ndraw10`) lack the `vary_stars` / `N_pool` / `n_bin_min` / `S_true_pool` columns and have `n_bin` instead. **raw samples**: `data/nuts_samples/jake_nuts_<nuts_draw_tag>_<survey>_rvir<X>_<baseline>.npz` (`nuts_samples_path`; gitignored, ~4 MB each), with `sigma_stream`, `sigma_binary`, `mu` as (n_draws, 8000) float32 arrays and per-draw `S_true`, `sigma_binary_rms`, `sb_cap`, `n_bin`, `n_bin_repeat`. only the rvir0 = 3, 50-draw, varied-stars run has one so far (rerun with fixed seeds: it reproduced the table exactly). plot: `plots/gd1_desi_jake_nuts_<nuts_draw_tag>_5yr.pdf` (each draw faint, pooled big square, truth cross).
- pooled, DESI, 5 yr, 10 draws of the repeats with the 160 brightest (`_ndraw10`) (km/s; median [16, 84]; binaries among the repeats 5-9 / 5-13 / 4-12 / 8-15 at 0.75 / 1.5 / 3 / 6 pc):

| rvir0 | truth | sigma_stream `shared` | sigma_stream `indep` (Nibauer+26 Eq. 1) | sigma_binary `shared` | sigma_binary `indep` |
|---|---|---|---|---|---|
| 0.75 | 0.81 | 2.05 [1.46, 2.46] | 2.11 [1.62, 2.39] | 1.30 [0.50, 1.87] | 0.96 [0.32, 1.69] |
| 1.5 | 0.65 | 4.13 [2.66, 4.58] | 4.16 [3.69, 4.58] | 1.85 [0.70, 3.51] | 1.12 [0.36, 2.79] |
| 3 | 0.43 | 2.30 [0.54, 2.68] | 2.44 [1.67, 4.07] | 2.70 [0.70, 4.51] | 1.45 [0.42, 3.94] |
| 6 | 0.44 | 0.87 [0.26, 1.35] | 1.12 [0.49, 1.58] | 1.33 [0.48, 1.76] | 0.97 [0.34, 1.55] |

- marginalized over the repeats, both likelihoods over-estimate the stream dispersion at every rvir0 (by 1.3-3.5 km/s at 0.75-3 pc, 0.4-0.7 at 6 pc), and `indep` sits at or above `shared`, with a lower sigma_binary: in Eq. 1 the repeats carry no information on the split, so less of the width goes to sigma_binary.
- **rvir0 = 3 pc, 50 draws of the stars and the repeats** (`indep_capstd_ndraw50_varystars`, DESI, 5 yr; pool = the 679 DESI stars): pooled sigma_stream **2.24 [0.66, 4.33]** km/s (per-draw medians 2.56 [1.17, 4.03]) against a truth of 0.50 (median over draws of the 160-star noiseless CoM dispersion, 0.42-0.63; 0.52 over the whole pool); sigma_binary 2.04 [0.62, 3.93] against a true rms binary motion of 3.96. binaries per draw 35-59 of 160, 4-17 of them among the repeats. so varying the stars as well widens the interval a lot (it was 2.44 [1.67, 4.07] with the brightest 160 and 10 draws of the repeats) and it now just reaches down to the truth, but the median is still ~4x the truth. convergence: one draw with r_hat 1.05, 6 of 50 draws with divergent transitions (438 in total, at most 235 in one draw).
- ⚠️ **convergence isn't clean in the draws runs**, unlike the single-draw runs: `indep` 3 pc has r_hat up to 1.10 and 1231 divergent transitions over the 10 draws (4 draws with 219-404 each), 1.5 pc one draw with 182; `shared` 6 pc r_hat 1.15, 74 divergent; `shared` 3 pc 28. probably the degenerate arc in `indep`; not looked into. treat those pooled intervals as approximate.
- not done: `nuts_shared_stream = False` on the single-draw cell's DESI/Via x 1 yr/5 yr grid is on disk (`jake_nuts_N160_rep35_indep_capstd.fits`) but not written up here; the draws cells for Via or the 1 yr baseline.

## `/claude_diagnostics`
one-off diagnostic scripts, not part of the pipeline: nothing imports them, and they don't write anything. run them from the repo root (`python claude_diagnostics/<script>.py`). the two that use the repo put it on `sys.path` by absolute path, so they work from anywhere.

### `claude_binary_selection_experiment.py`
only reads a pickle. for one sim it prints (1) the binary fraction and binary orbital-motion spread in v_gsr for the full unbound sample, the via-matched and desi-top-679 samples, and in bins down the brightness ranking of the eligible pool; and (2) noiseless GMM fits on each of those selections, for both binary treatments. written to explain the binaries selection effect -- see "binaries in the observed samples" in the GMM section. run across orbits/rvir to check the effect holds beyond GD-1 rvir0 = 0.75.

### `claude_profile_desi.py`
only reads a pickle. `python claude_diagnostics/claude_profile_desi.py <rvir>` (GD-1, `desi`, CoM; a few minutes). rebuilds `gmm.py`'s DESI selection and prints a noiseless fit on the same stars (the truth), the noisy fit with `x_err`, 1D `mle_sigma` on the p_thin > 0.5 stars, and a **profile likelihood** in the thin-stream sigma_v_gsr: that one parameter fixed on a grid, everything else refit (Powell, warm-started). delta nll = 1.35 is the one-sided 95% upper limit. written to diagnose a DESI thin sigma_v_gsr = 0 that turned out to be a fluke of the noise draw (log, 14 Oct). the profile loop is the template for an upper-limit helper in `gmm.py`. it hardcodes the thin component's v_gsr index in theta for n = 2, K = 4.

### `prior_degeneracy.py`
synthetic data only, a few seconds. grid posterior in (sigma_stream, sigma_binary) for the independent-epoch (Nibauer+26) likelihood under several priors (equal and unequal half-normal scales, flat, log-uniform, a 1 km/s half-normal on sigma_stream), printing the 16/50/84% of each marginal. shows the prior alone sets the split, and that changing the half-normal scale can't reveal it. see "what this means for Nibauer+26" in the `single_component_comparisons_mockobs.py` section.

### `three_regions.py`
no data at all, ~10 s. the Nibauer+26 setup (one global sigma_binary, one sigma_stream per GD-1 region) with the expected (Asimov) likelihood for N stars per region and equal errors, i.e. no binary information. given the region totals, prints the posteriors for tau = 10 / 30 and N = 50 / 150, plus one region alone. reproduces the Nibauer+26 values; same section as above.

### `fig10_compare.py`
no data at all, a few seconds. `three_regions.py` with the actual Nibauer+26 bin sizes (35 / 74 / 51 stars x 195/160 measurements per star) for errors of 2, 3, 4 km/s. prints sigma_binary, each bin's sigma_stream, and each bin's posterior shape on 0-8 km/s in 0.5 km/s steps, for comparison by eye with their Fig. 10. same section as above.


# notes on `paf` (PETAR_ANALYSIS_FUNCTIONS.py)
~2500 lines, a lot of it leftover from the old grid (the `#<-- LEFTOVER FROM OLD GRID` markers are honest -- `define_paths*`, `define_apocenters`, `unpack_escaper_dict` etc. are superseded by `extended_grid_info`). the notes below are for the parts the current pipeline actually touches.

## conventions to keep straight
- **`i` vs `file_index`.** the new grid writes snapshots every 10 Myr, so `file_index = int(i/10)` where `i` is simulation time in Myr. functions are not consistent about which one they want, and nothing validates it:
  - sim time in Myr: `intrinsic_stream_data_v3(path, i, ...)` (computes `file_index` itself), `prog_position(init_displacement, i)`, `xform_to_core_frame` (it converts internally via `file_naming_convention`).
  - file index: `CM_to_galcen_frame`, `core_to_galcen_frame`, `is_dissolved` (indexes `tidal.n[i]`), and the `file_index=` kwargs.
  - both: `load_coords_v2(path, i, ..., file_index=...)` takes time as `i` *and* index as `file_index`; `straighten_stream_orbit_interp(coords, yval, core, i, ...)` wants `i` as a time when `use_core=False` and a file index when `use_core=True`.
  - `file_naming_convention="every integer"` is the *default* on `load_particle`/`xform_to_core_frame`/`clip_outside_rtid`, which is the OLD grid's convention. for the new grid pass `"every 10"` (or pass the file index directly).
- **`data.core` is written at the same 10 Myr cadence as the snapshots**, so `core.pos[file_index]` lines up. (checked against `m3/lm/0.75/0`: 501 snapshots, core time column steps by 10.)
- **the sims keep running past the present day.** the present day is `data.<age/10>` (`age` from `retrieve_sim_info`, which already includes the +100 Myr), but e.g. `gd1/hm/0.75/0` has snapshots up to `data.1500` (15 Gyr) and 1501 rows in `data.core`/`data.tidal`. so "last file in the directory" is **not** the present day.
- **units.** petar outputs pc, pc/Myr, Msun. `StreamFrame` wants kpc, kpc/Myr and returns deg, mas/yr, kpc, km/s. most `paf` functions hand back astropy quantities, but the streamframe coord dicts are bare floats.
- **`init_displacement` is kpc, km/s** everywhere it is *consumed* (`prog_position`, `integrate_prog_orbit`, `straighten_stream_orbit_interp` all do `init_displacement[:3] * u.kpc`), but `data/init_displacements.txt` *writes* it in pc because that's what petar wants. `extended_grid_info.__init__` does the pc->kpc conversion once, at the bottom, for the six real orbits (`circ` was always already in kpc). **this was a silent 1000x bug** -- it put the progenitor reference at 16.6 Mpc instead of 16.6 kpc, where the potential is negligible, so the "orbit" free-streamed in a straight line. sanity check if you ever touch it: `paf.prog_position(init_displacement, age)` must reproduce that orbit's row in `FINAL_ics_nolmc.csv` (it does, to 0.0000 kpc, for all six).
  - **`init_displacement` is the *initial* (t = -(age+100) Myr) phase space position, not the present day.** the present day is the orbit's row in `FINAL_ics_nolmc.csv`. easy to mix up -- see the orbital phase bug.
- **potential is `gp.BovyMWPotential2014(units=galactic)` everywhere** -- `prog_position`, `integrate_prog_orbit`, `straighten_stream_orbit_interp`, `simspect.prog_orbit_track`. petar was run with `external_mode='galpy'` to match. don't change one without the others.
- **`interrupt_mode='bse'`** is assumed by every loader (stellar evolution on).
- **`star.mass0` is NOT the ZAMS mass except at snapshot 0.** see its own section below -- this one is easy to get wrong silently because the field is literally labelled "initial stellar mass."

## `star.mass0` is an *effective* initial mass, and BSE rewrites it
if you want the IMF, or want to paint isochrone photometry onto stars by their birth mass, **use `m0_zams`** (off the `luminous`/`companions` subdicts), or read `mass0` out of `data.0` and match to the present day by `id` yourself. the present-day `mass0` is a different quantity.

`mass0` (`M0`/`mass` in the Fortran, `m0` in petar's C++) is the ZAMS mass *whose evolutionary track the star is currently being interpolated along*. SSE is a set of analytic fits parameterized by (M0, Z, age), so whenever the real mass `mt` changes, BSE re-parameterizes the star onto the track of a different initial mass to keep the fitting formulae valid. `mt` (petar's `mass`) is the real current mass; `mass0` is a fitting coordinate.

this is stock Hurley, not a petar quirk:
- Hurley, Pols & Tout (2000), MNRAS 315, 543 -- SSE. the mass-loss section is where M0 gets reset for MS/HG stars so they follow a lower-mass track.
- Hurley, Tout & Pols (2002), MNRAS 329, 897 -- BSE, for the mass-transfer / rejuvenation case.

the actual rewrite sites, in `~/software/PeTar/bse-interface/bse/`:

| site | trigger |
|---|---|
| `evolv1.f:152-157`, `evolv1.f:213-215` | wind mass loss, single star: `if(kw.le.2.or.kw.eq.7) mass = mt` |
| `evolv2.f:712-717` | same, in the binary evolver |
| `evolv2.f:2041-2043` | Roche-lobe mass transfer -- **both** donor and accretor, if MS or naked-He MS |
| `evolv2.f:2047-2063` | HG stars, with a revert if the reduced M0 would leave the core too massive for the new track |
| `hrdiag.f:787` | **remnant formation: at WD birth `mass = mt`, i.e. M0 is overwritten with the remnant mass.** this is the big one |
| `hrdiag.f:446` | core exposed -> naked He star |
| `mix.f`, `merge.f` | collisions/mergers; M0 is passed by reference and rewritten |

petar just passes `m0` by reference into `evolv1_`/`evolv2_`/`merge_`/`mix_` (`bse_interface.h:1423, 1529, 1649`) and writes back whatever comes out. the field comment `///> Initial stellar mass in solar units` (`bse_interface.h:259`) is true only at t=0.

measured on `gd1/hm/rvir=6/copy 0`, snapshot 0 vs. 2700 Myr (14999 -> 14839 particles):

| present-day type | N | N with changed `mass0` | median delta |
|---|---|---|---|
| 0 (low-mass MS) | 12714 | 1661 | 0.0000 |
| 1 (MS) | 1050 | 1050 | -0.0007 |
| 2-6 (giants) | 37 | 37 | -0.060 |
| 10-14 (remnants) | 1038 | 1038 | **-1.098** |

3786 of 14839 survivors changed (3589 down, 197 up). summed `mass0` over survivors falls 9022 -> 7094 Msun and the count above 1 Msun falls 1618 -> 841, so the present-day `mass0` histogram is depleted at the high-mass end and piled up at low mass -- it does not look like an IMF, and no amount of adding merged stars back fixes it.

**ID bookkeeping is fine, so the ID match works.** checked on the same run: IDs are unique in both snapshots, nothing is renumbered, and there are **zero** present-day IDs absent from snapshot 0 -- 160 IDs disappear (mergers) and none appear. so

```python
order = np.argsort(ids0)
zams_index = order[np.searchsorted(ids0[order], ids1)]
assert (ids0[zams_index] == ids1).all()
m0_zams = m00[zams_index]
```

recovers the birth mass of every surviving star, and `np.concatenate([m0_zams, m00[~np.isin(ids0, ids1)]])` reproduces `m00` exactly (verified with `np.sort`). `intrinsic_stream_data_v3` does exactly this internally (see below); the hand-rolled version is still what you want if you're working straight off `load_particle` output, or if you need the stars that have *disappeared*, which the subdicts by construction don't contain.

⚠️ naming: the field is **`star.mass0`** in the petar python reader (`petar/bse.py`) and **`m0`** in the C++ (`bse_interface.h`). there is no `star.m0` and no top-level `particle.mass0` -- both are `AttributeError`.

## `extended_grid_info(scratch=False)`
holds paths, apocenters, stream ages (Vedant's), and init displacements as attributes. grid axes are orbit x stellar_pop (`lm`/`hm`) x rvir (`0.75, 1.5, 3, 6` pc, indexed 0-3) x copy. use `retrieve_sim_info(orbit, stellar_pop, rvir_index, copy) -> path, apo, age, init_displacement` rather than assembling paths by hand.

**`scratch=False` (the default, and what every script uses) is the holystore copy** at `/n/holystore01/LABS/itc_lab/Users/amphillips/extended_grid/`, which is complete and not purged. `scratch=True` is netscratch, which is actively losing files -- the plain `data.<i>` snapshots age out before the `.single`/`.binary` siblings, so a purged directory *looks* populated and `prepare_nbody_data_anycopy` reports every copy "unfinished". the source comment at the top of `__init__` has the details.

the six `*_init_displacement` literals are pasted verbatim from `init_displacements.txt` (i.e. in pc), and a loop at the end of `__init__` converts the position components to kpc. if you paste a new orbit in, paste the raw pc numbers and add its key to that loop -- don't pre-convert.

**not every copy finished.** every `(orbit, stellar_pop, rvir_index)` has at least one that did, but which one varies; use `prepare_nbody_data_anycopy` rather than assuming `copy=0`. availability at `stellar_pop='lm'`, present-day snapshot present (`Y`) or not (`.`) (the pipeline runs `hm`, which hasn't been tabulated here):

| orbit | rvir_index=0 | 1 | 2 | 3 |
|---|---|---|---|---|
| gd1 | Y . . . Y | Y Y . Y . | Y . Y Y Y | Y Y Y Y Y |
| pa5 | Y . Y . . | Y Y Y Y . | . Y Y Y . | Y Y Y Y Y |
| aau | Y Y Y Y Y | . Y Y Y Y | Y Y Y Y Y | Y Y Y Y Y |
| m3  | Y Y Y Y Y | Y Y Y . Y | Y Y Y Y Y | Y Y Y Y Y |
| c19 | Y . . . Y | . . Y Y Y | Y Y Y Y Y | Y Y Y Y Y |
| jet | Y Y . . Y | . . . Y Y | Y . Y Y Y | Y Y Y Y Y |

the one `hm` cell checked so far: `gd1` rvir_index=0 has all five copies finished.

sharp edges:
- `retrieve_sim_info` is a chain of bare `if`s with no `else`, so a typo'd orbit string gives `UnboundLocalError` on `base_paths`.
- `circ` is the odd one out (not used by `final_datasets.py`/`gmm.py`): it comes from the OLD grid under `conroy_lab/Lab/amphillips/finished_grid/`, has one realization per rvir (no `copy` subdir), `circ_lm_paths[0]` is flagged unfinished, and `retrieve_sim_info` returns `age = 10000*10 = 100000` -- which is why `prepare_nbody_data_anycopy` hardcodes `i=30000` for circ instead of using the returned age.

## loading petar data
`load_core(path)`, `load_tidal(path)`, `load_particle(path, i)` are thin `petar.*` wrappers. `is_dissolved(path, i, threshold=100)` = fewer than 100 stars inside the tidal radius at index `i`.

`data.tidal` (via `load_tidal`) is a cheap way to get the cluster's evolution without loading snapshots: columns `time` [Myr], `rtid` [pc], `mass` (**bound** mass, Msun), `n` (bound count), `pot`, at the same 10 Myr cadence. e.g. `gd1/hm/0.75/0`: rtid 48.8 / 37.1 / 25.7 pc, bound mass 9940 / 3384 / 2388 Msun, n 14999 / 9054 / 6205 at t = 0 / 1700 / 2700 Myr (at t = 0 everything is bound, so `tidal.mass[0]` is the initial cluster mass).

**which frame a file is in** (this is the thing that bites):
- `data.<i>` (all particles) is in the **CM frame**; the offset to galactocentric lives in the file header (`petar.PeTarDataHeader(...).pos_offset/vel_offset`). use `CM_to_galcen_frame`.
- `data.<i>.single` / `data.<i>.binary` are in the **core frame**; the offset is `core.pos[file_index]`. use `core_to_galcen_frame`.
- `xform_to_core_frame` goes the other way for the all-particles file (CM -> galcen -> core), and also returns `rrel`; `clip_outside_rtid` pairs it with `tidal.rtid` to get a bound mask (currently unused).

`fix_core_vel(core)` differentiates the core *position* track to get a core velocity, because the raw `core.vel` has spurious jumps. its docstring says it only works at 1 Myr cadence, but it does `np.gradient(core_x, times)` against the actual `core.time` array, so it is fine at 10 Myr too -- the caveat is stale. **caveat:** the corrected velocity is only used for the *progenitor reference coordinate*. putting singles/binaries into the galactocentric frame (`core_to_galcen_frame`, and the equivalent block inside `load_coords_v2`) still adds the raw uncorrected `core.vel`, so those velocities inherit exactly the jumps `fix_core_vel` exists to avoid. probably fine since velocities get differenced against the progenitor downstream, but know that it's there.

`correct_core(core)` is a near-duplicate of `fix_core_vel` that also returns the position, and drops units on the velocity (`np.array` of quantities silently strips them). nothing calls it -- dead code, delete or ignore.

## `load_coords_v2(path, i, ...)`
loads all/single/binary particles, puts each in the galactocentric frame, and runs each through `StreamFrame` against a progenitor reference coordinate. returns `particle_data, streamframe_data` -- two lists ordered `[all, singles, binaries]` per the `load_*` flags. (the docstring listing six return values is stale.)

the branching is all about **where the progenitor reference comes from**:
- `check_dissolved=True`: use the live core if the cluster survives, else fall back to integrating from the core position at `tdis_estimate`.
- `check_dissolved=False, use_core=True`: always integrate forward from the core at `tdis_estimate`.
- `check_dissolved=False, use_core=False`: ignore the core entirely, integrate the progenitor orbit from `init_displacement`. **this is what `intrinsic_stream_data_v3` uses by default**, and it's the right choice once the cluster has dissolved.

latent bug: the all-particles filename is built as `path+"data."+str(file_index)` (the raw kwarg), not `file_i`, so calling with `file_index=None` and `load_all=True` tries to open `data.None`. always pass `file_index`.

## the `StreamFrame` (intrinsic) coordinate system
credit Jake Nibauer, `scripts/streamframe.py`. origin is the **Galactic center**; x_hat is the progenitor's position direction, z_hat its angular momentum direction. so:
- `phi1`, `phi2` [deg] -- longitude/latitude along the progenitor orbit plane
- `r` [kpc] -- **galactocentric** radius, not a heliocentric distance
- `vr` [km/s], `pm_phi1`, `pm_phi2` [mas/yr] -- velocities are taken **relative to the progenitor** (`DeltaV = v - v_prog`) before projecting, and `pm_phi1` is *not* multiplied by cos(phi2)

this is the "intrinsic"/God's-eye frame. contrast with `streamframe_coords_observed`, which is heliocentric, reflex-corrected, and uses real ICRS great circles -- the key names overlap but mean different things (`r` vs `distance`, `vr` vs `v_gsr`). the current pipeline uses the intrinsic frame only for `in_rtid`; everything fit is in the observed frame.

## `intrinsic_stream_data_v3(path, i, core, apo, init_displacement, use_core=False, binary_treatments=[...], N_rtid_boundary=1.)`
the main entry point. returns one big nested dict.

top level: `init_displacement`, `nsingles`, `nbinaries`, `IDs`, `pot`, the raw petar `singles` / `binaries` objects of that snapshot (core frame, petar units), plus one subdict per binary treatment. the full key list is in the "what's in a data_dict" section.

**ordering gotcha:** top-level `IDs` and `pot` are concatenated as `[singles, binary_p1, binary_p2]`, so they always have length `nsingles + 2*nbinaries`. the subdict arrays are concatenated as `[singles, binaries]`, where "binaries" is one row per binary for `CoM`/`luminous` but two for `companions`. so **`IDs` only lines up with the `companions` subdict**; for `CoM`/`luminous` you have to slice `[:nsingles+nbinaries]`.

the three treatments:
- `CoM` -- each binary as a single center-of-mass point. no `L`/`R`/`type`/`m0_*`.
- `luminous` -- each binary represented by its brighter component (`binaries.p1.star.lum >= binaries.p2.star.lum`).
- `companions` -- both components kept separately; this is the only one that double-counts.

each subdict has: `coords` (the `StreamFrame` dict), `pos`/`vel` (galactocentric, with units), `mass`, `inMW`, `trim`, `in_rtid`, and `phi2_straight`/`r_straight`/`vr_straight`/`pm_phi1_straight`/`pm_phi2_straight`. `luminous`/`companions` also get `L`, `R`, `type`, `m0_zams`, `m0_effective`.

### the two initial-mass columns
`luminous` and `companions` carry **both** flavours of initial mass:
- **`m0_zams`** -- the genuine birth mass. read out of `data.0` (`p0.star.mass0`, which at t=0 is exactly `p0.mass`) and matched to the present-day particle by `id`. **this is the one for the IMF and for isochrone photometry**, and what `final_datasets.py` uses.
- **`m0_effective`** -- the present-day `star.mass0`, i.e. whichever ZAMS track BSE is currently interpolating the star along. use it if you want to reproduce what BSE itself thinks the star is, e.g. when cross-checking `L`/`R`/`type` against an SSE track.

both are bare floats in Msun (no astropy units, unlike `mass`). neither is the *current* mass -- that's `mass`. on `gd1/hm/rvir=6/copy 0` at 2700 Myr, the `companions` subdict has 3786 of 14839 entries where they disagree, 1489 vs 841 above 1 Msun, and max 111.1 vs 61.0 Msun.

implementation notes:
- **the ID match is done per component, before the `luminous` selection** (`sid`, `bp1id`, `bp2id` each get their own `searchsorted`, then `np.where(luminous_mask, ...)`), so the `IDs` ordering gotcha does *not* bite here. don't "simplify" this by slicing top-level `IDs`.
- each lookup `assert`s `ids0[index] == id`, so a present-day particle with no snapshot-0 counterpart fails loudly.
- it costs one extra read of `data.0` per call, via `load_particle(path, 0, file_naming_convention="every integer")`. the lookup is rebuilt inside the `binary_treatment` loop (twice when both `luminous` and `companions` are requested) -- cheap, doesn't matter.

### `in_rtid` (for cutting the progenitor out)
`in_rtid` = distance from the core <= `N_rtid_boundary * tidal.rtid[file_index]`, so `~in_rtid` is the "remove the progenitor" mask. `N_rtid_boundary` defaults to 1; **`final_datasets.py` uses 2.0** so stars that will be recaptured as rtid re-expands are also excluded. it is computed on the **core-frame** `singles.pos`/`binaries.pos` (the frame `rtid` is measured in), as a full 3D radius, and it is applied **before** any trimming.

still wrong: **it is the wrong length for the `companions` subdict.** `b_r` comes from `binaries.pos`, one row per binary, so `in_rtid` is always `nsingles + nbinaries` long. that matches `CoM` and `luminous`, but `companions` is `nsingles + 2*nbinaries`, so indexing a companions array with it raises `IndexError`. either build `[in_rtid_s, in_rtid_b, in_rtid_b]` inside the companions branch or don't use `in_rtid` there.

also minor: `load_tidal(path)` is called inside the `binary_treatment` loop, so the tidal file gets read once per treatment.

## mock binary RV curves (the "Mock observations and detection fractions" block)
- `get_orbital_params_from_binaries(binaries, v_sys_vals, rng=None)` -> `(N_bin, 6)` array `[v0 (km/s), K (km/s), w, phi0, e, P (day)]` from a petar `Binary` object and systemic velocities (km/s, one per binary). isotropic inclination (`draw_inclinations`), `w` uniform in [0, 2 pi), `phi0` uniform in [0, 1). ~1 s per 2000 binaries (`calc_amplitude` loops over binaries). `get_orbital_params(paths, times, n)` is the older version that loads a snapshot itself and puts the binaries in the intrinsic frame.
- `calc_amplitude(binary, inclination)`: K of the **heavier** star (`a1 = (M_light/M_tot) a`), not necessarily the luminous one.
- `radial_velocity(t, params)`: Keplerian RV, `M = 2 pi (t/P - phi0)`, i.e. `phi0` is the phase as a fraction of the orbit. `t` and `P` in the same units (days).
- `get_rvs(params_all, obstimes_all, verbose, add_noise, noise_level, rng)`: loops `radial_velocity` over binaries; `obstimes_all` is (N_bin, n_epochs). `add_noise` only takes one scalar `noise_level`.
- `solve_kepler`: Newton from E = M; converges to < 2e-10 for every binary in the gd1 0.75 pc snapshot (e up to 0.998).
- `get_obstimes(N, rng)`: three epochs (0, +30-90 d, +3-5 yr). not used by the two-epoch experiment, which draws its own.
- `get_detections`: chi^2 binary detection; its `bool_arr='detet'` branch is spelled that way.

## masks: `inMW` and `trim` (intrinsic frame)
from `trim_coords_percentile(coords, low=1, high=99, apo=apo)`:
1. `inMW = coords['r'] <= 1.5*apo` -- drops stars flung well outside the orbit. **this is the only thing `apo` is used for**, so passing the wrong orbit's apocenter silently gives you the wrong outlier mask.
2. `trim` = 1st-99th percentile clip in phi1 and phi2, **computed on the already-`inMW`-selected array**.

so `trim` has length `inMW.sum()`, not `len(inMW)`, which is why the older code is always `x[inMW][trim]` and why the two are not interchangeable or combinable with `&`.

**the current pipeline doesn't use these** -- `final_datasets.py` replaced them with `trim_obstream_percentile` in the observed frame, where `trim_new` is full length and combines with `&`. they're still computed and still in the pickles, and `simspect.__main__` still uses them.

## `straighten_stream_orbit_interp(coords, yval, core, i, ...)`
the *intrinsic*-frame straightening that produces the `*_straight` residuals in `intrinsic_stream_data_v3`. integrates the progenitor orbit +/- `Dt` Myr around time `i`, pushes that orbit through `StreamFrame` with the same reference coordinate, interpolates `yval` against `phi1`, and subtracts it. `yval` can be a string or a list of strings.

`Dt` is found by an adaptive search (`Dt_start=240` from v3): it grows `Dt` until the orbit's phi1 range covers the data, but shrinks it if the orbit wraps around in phi1 (`dphi1 > 100` deg between samples), halving the step whenever it flips between those two cases.

**it returns the iteration count as the last value -- check it.** the loop caps at `max_iters=100` and if it exhausts them it just falls through and straightens anyway, with no warning. an `itr` of 100 means the straightening never converged. the returned residuals are for *all* particles (untrimmed).

the interpolant sorts `orbit_coords['phi1']` first (`np.interp` requires increasing `xp` and doesn't check), and the multi-`yval` closures are bound with `y=y` (without it they'd all interpolate the last key).

`straighten_stream_polynomial(phi1, y, degree=5, trim_criteria=[mask1, mask2], return_poly_fn=True)` is the polynomial alternative, used by `simspect.poly_straightening`. note `trim_criteria` has no working default -- it unpacks two masks unconditionally, so leaving it `None` is a `TypeError` (`poly_straightening` always passes one).

## photometry
### isochrone interpolation: ZAMS mass -> Gaia mags (`noise.py`) -- the live path
three functions, no state, all float64. the point of doing it this way rather than by blackbody is that **the isochrone is a free parameter**: swap the `artpop.fetch_mist_iso_cmd` call in `final_datasets.py` and the whole stellar population changes. it is **[Fe/H] = -2, requested at 12 Gyr (actually 12.6 Gyr, v/vcrit = 0.4 -- see the artpop subsection below)** right now, which does *not* match the 2.4-4 Gyr dynamical ages -- a known caveat, justified by all the massive star evolution being over by a few Gyr anyway.

- `build_isochrone_table(iso, bands=ISO_BANDS, max_phase=5, max_eep=POST_AGB_EEP, mass_col='initial_mass')` -> `(m0_grid, {band: values})`. `ISO_BANDS = ('Gaia_G_EDR3', 'Gaia_BP_EDR3', 'Gaia_RP_EDR3', 'log_Teff')` -- `log_Teff` rides along because viamock needs it. it sorts on `mass_col` and **drops any non-increasing node** rather than trusting the file, because `np.interp` requires increasing `xp` and does not check. it prints if it drops anything.
- `isochrone_photometry(m0_query, m0_grid, table)` -> `(phot, on_iso)`. **off-isochrone stars get `NaN`, not a clamped edge value**, and `on_iso` flags them (two-sided: the sim IMF reaches below the isochrone's low-mass end as well as above its turnoff). deliberate: a forgotten mask then shows up as a hole in the CMD rather than a fake pile-up at the tip of the AGB. `final_datasets.py` saves `on_iso` as `alive`.
- `gaia_from_isochrone(m0_query, iso, ...)` -> `(G, BP, RP, on_iso)`, the one-call wrapper.

mags are **absolute**; `final_datasets.py` applies `paf.m_from_M` with the CoM distances.

query with **`m0_zams`**, and interpolate against the isochrone's `initial_mass` column, NOT `star_mass` -- see the `star.mass0` section.

#### precision: this is the part with no headroom
the evolved sequence is nearly degenerate in initial mass. measured on the 12 Gyr, [Fe/H]=-2 file downloaded from the MIST web interpolator (1460 rows; see the artpop subsections below for the isochrone that's actually used -- same scale, different edges):

| phase | N | m0 range |
|---|---|---|
| 0 MS | 203 | 0.1022 - 0.7884 |
| 2 RGB | 151 | 0.7887 - 0.8025 |
| 3 CHeB | 102 | 0.80249 - 0.80416 |
| 4 EAGB | 101 | 0.80416 - 0.80438 |
| 5 TPAGB | 601 | 0.804376 - 0.804398 |
| 6 post-AGB | 302 | 0.804398 - 0.80521 |

so **everything above the turnoff occupies 0.0168 Msun**, and the tightest node spacing is **1.58e-11 Msun**, on the TP-AGB. that is ~1e5 x float64 eps at 0.8 Msun, so it is resolvable, but:
- **never round, bin, or cast to float32 anywhere on the mass path.**
- **linear, not cubic.** on the RGB/AGB segments dM/dmag is ~1e-9; a spline overshoots enormously between nodes and invents points off the isochrone.
- `initial_mass`, not `star_mass`, is also what makes the map invertible at all -- `star_mass` turns over once winds start.

checked (on the downloaded file): interpolation at the nodes reproduces the table bit-exactly, the tightest node pair comes back cleanly distinguished (0.012 mag in G), shuffle-invariant, NaN/`on_iso` correct on both sides.

#### artpop's isochrone is not the one on the MIST web interpolator
`artpop.fetch_mist_iso_cmd(log_age, feh, phot_system, v_over_vcrit=0.4)` doesn't interpolate anything: it reads artpop's packaged MIST v1.2 grid file (here `~/.artpop/mist/MIST_v1.2_vvcrit0.4_UBVRIplus/MIST_v1.2_feh_m2.00_afe_p0.0_vvcrit0.4_UBVRIplus.iso.cmd`, 107 ages) and picks one age with `age_index`. so relative to the downloaded file:
- **v/vcrit = 0.4 by default** (the downloaded file is 0.0). `final_datasets.py` and the `noise.py` scratch cell both have `v_over_vcrit` commented out. `0.0` is the other packaged option.
- **the age snaps to the grid**: `log_age = log10(12e9) = 10.079` comes back as **log age 10.10 = 12.6 Gyr**. `[Fe/H]` must match a grid file (-2.00 does).
- the header `Zinit` differs too (1.43e-4 vs 1.89e-4 for the downloaded file, both at [Fe/H] = -2.00).

these are presumably why the upper mass differs: the artpop isochrone also has 1460 rows but **tops out at 0.79435 Msun, not 0.805** (not checked which of the three is responsible).

#### EEPs vs `phase`
the MIST primary EEPs are defined in **table II of the MIST README** ([README_tables.pdf](http://waps.cfa.harvard.edu/MIST/README_tables.pdf), "Primary Equivalent Evolutionary Points"). artpop hardcodes the same boundaries in `SSP.select_phase` (`artpop/stars/populations.py`, which cites that table): MS 202-453, RGB 454-605, CHeB 606-706, EAGB 707-807, TP-AGB 808-1408, post-AGB 1409-1710, WD cooling > 1710. the in-between primary EEPs are RGB tip = 605 and ZACHeB = 631.

the `phase` column is a coarser label, and **the two files label post-AGB differently**: in the downloaded file `phase == 6` is exactly `EEP >= 1409` (302 rows), but artpop's packaged grid labels **all of EEP 808-1710 phase 5**, so its `phase` can't separate TP-AGB from post-AGB. its other phase boundaries fall on primary EEPs (0 | 2 at 454, 2 | 3 at 605, 3 | 4 at 707, 4 | 5 at 808). the artpop isochrone split by primary EEP (it starts at EEP 251, i.e. 0.1 Msun, not at the ZAMS EEP):

| EEP range | stage | N | artpop `phase` | m0 range |
|---|---|---|---|---|
| 251-453 | MS (to TAMS) | 203 | 0 | 0.100000 - 0.777822 |
| 454-604 | SGB + RGB (TAMS to RGB tip) | 151 | 2 | 0.778112 - 0.791740 |
| 605-630 | He flash (RGB tip to ZACHeB) | 26 | 3 | 0.791745 - 0.791765 |
| 631-706 | CHeB (ZACHeB to TACHeB) | 76 | 3 | 0.791769 - 0.793341 |
| 707-807 | EAGB | 101 | 4 | 0.793342 - 0.793553 |
| 808-1408 | TP-AGB | 601 | **5** | 0.793553 - 0.793574 |
| 1409-1710 | post-AGB (to WDCS) | 302 | **5** | 0.793574 - 0.794351 |

so `max_phase=5` does **not** remove post-AGB here, which is why `build_isochrone_table` also cuts `EEP < max_eep` (default `noise.POST_AGB_EEP = 1409`; see the features bullet below). `isochrone_cmd_track` uses the same constant.

#### two things that are features, not bugs
- **interpolating in mass gets the relative numbers of giants right for free.** a phase is sampled in proportion to the initial-mass interval it occupies, which is exactly its lifetime x IMF weight.
- **post-AGB (the proto-WD tail) is dropped by the `max_eep=POST_AGB_EEP` default** (`EEP < 1409`), since those stars are already cut from the sim side by `nonrem` (`type < 10`). `max_phase=5` is also still the default but only does anything on the downloaded MIST file, where post-AGB is phase 6 -- on the artpop isochrone it's phase 5 (table above). on the downloaded file the two cuts remove the identical 302 rows. with the EEP cut the artpop table is 1158 nodes, upper m0 edge **0.793574** (0.79435 without it), max logTeff 3.885 (5.05 without it). `max_eep=None, max_phase=None` keeps everything, and moves the upper m0 edge (so the `alive` count depends on it). without the cut, the 2 or so N-body stars per pickle with m0_zams in 0.793574-0.79435 get post-AGB/WD photometry (absolute G ~ +10) and count as `alive` (log, 29 Sep).

caveats that are real selections, not rounding details:
- the isochrone is old, so it stops at 0.793574 Msun (artpop, after the post-AGB cut; 0.804398 for the downloaded file) and every more massive sim star is thrown out (`alive`). at `hm` the discarded stars are the luminous ones, so quote how many got dropped alongside any luminosity-weighted number.
- for binaries only the luminous component's photometry is used; the companion's light isn't added.

### colour -> Teff: polynomial fit to the isochrone (`noise.py`) -- the live path
step 5 of `assign_photometry_from_catalog`: a catalog star's `log_Teff` (for viamock) comes from its **BP-RP alone**, through a least-squares polynomial fit of log Teff vs BP-RP to the isochrone track. colour only, so it doesn't depend on the N-body distance each catalog star was painted onto (the CMD nearest-point lookup in the next subsection does, through absolute G).

```python
track    = noise.isochrone_cmd_track(isocmd)
log_Teff = noise.Teff_from_gaia_color_temp_relation(track, BP_RP)                     # (N,)
log_Teff, poly = noise.Teff_from_gaia_color_temp_relation(track, BP_RP, return_poly=True)
```

- `Teff_from_gaia_color_temp_relation(track, BP_RP, deg=2, G_max=5000.0, return_poly=False)`. refits every call (~1000 points, cheap). `Polynomial.fit` maps the colour onto [-1, 1], so the fit is well conditioned. `G_max` keeps only track points with absolute G < `G_max`; the default 5000 (or `None`) means **no cut**, i.e. the whole track down to the 0.1 Msun end of the MS.
- **why a fit, not interpolation:** log Teff is nearly a single function of BP-RP, but the isochrone's MS, RGB, HB and AGB overlap in colour at slightly different Teff (log g), so the nodes aren't a function of colour and `np.interp` on colour-sorted nodes would zigzag between branches. the fit averages over them.
- **outside the fitted colour range -> NaN**, no extrapolation (BP-RP 0.23-2.54 on the artpop [Fe/H] = -2 isochrone; the blue end is the model HB). those stars then get NaN Via errors and fail `good_RV` in `gmm.py`. catalog stars affected: gd1 1 / 1578, **aau 17 / 334** (BP-RP down to -0.08), jet 1 / 144, pa5 0, c19 0. probably blue stragglers / hot HB stars / contaminants; clamping to the fit's edge is the alternative if they should stay.
- fit quality, degree 2, residuals about the isochrone points: ~0.006 dex rms (max ~0.017) with no cut; ~0.005 dex rms (max ~0.014) with `G_max = 4`. degree 3-5 gains little (0.0037-0.0031 rms at `G_max = 4`). the difference between cuts: at BP-RP > 1.2 the lower MS runs ~0.02-0.03 dex cooler than the giants, and with no cut it pulls the fit toward itself there. with no cut: `log Teff = 3.951 - 0.288 (BP-RP) + 0.0476 (BP-RP)^2`.
- **reddening goes straight into Teff.** the catalog `bp_rp` isn't dereddened, and the catalog MS sits ~0.1 mag redder than the isochrone. the slope is ~-0.2 dex/mag near BP-RP ~ 0.7, so that's ~0.02 dex too cool (a median 0.026 dex cooler than the nearest-point Teffs on the gd1 catalog, which partly absorbed the offset by snapping stars onto the bluer isochrone). **accepted for now.**
- the scratch cell just above the function in `noise.py` (commented out) is where the fit was developed and plotted (fit + residual panel, coloured by G).

### CMD-space lookup: catalog (G, BP-RP) -> nearest isochrone point (`noise.py`) -- only the `iso_dist` diagnostic
the reverse of the mass interpolation: real catalog stars never sit exactly on the isochrone, so each one is snapped to the nearest point on the track. `assign_photometry_from_catalog` only calls it for `iso_dist` (default metric); the colour-Teff fit above sets `log_Teff`. (pickles that haven't been regenerated since the colour-Teff switch, which includes the supplemental ones, still have nearest-point Teffs; see the refit TODO.)

```python
track   = noise.isochrone_cmd_track(isocmd)                       # EEP-ordered polyline, post-AGB dropped
nearest = noise.nearest_isochrone_point(G_abs, BP_RP, track)      # dict of (N,) arrays
logTeff = nearest['log_Teff']      # == noise.Teff_from_gaia_isochrone(track, G_abs, BP_RP)
```

- `isochrone_cmd_track(iso, max_eep=POST_AGB_EEP, G_col=..., BP_col=..., RP_col=..., carry=('log_Teff',))` -> dict `EEP`, `G`, `BP_RP` (absolute) + each of `carry`, ordered by EEP. takes the **raw** artpop/MIST table, not `build_isochrone_table`'s output, because:
  - **it's ordered by EEP, not `initial_mass`.** consecutive EEPs are consecutive points along the sequence. on the TP-AGB many nodes share an identical `initial_mass`; `build_isochrone_table` drops those as non-increasing, which is right for mass interpolation but would leave holes in the CMD track.
  - **the post-AGB cut is on EEP.** `POST_AGB_EEP = 1409` is the MIST primary EEP where post-AGB starts (see the EEPs vs `phase` subsection), since artpop's `phase` can't distinguish it. without the cut, faint blue catalog stars snap onto a WD cooling track at logTeff ~ 5. `max_eep=None` keeps everything. with the default the track is 1158 nodes, max logTeff 3.885.
  - MIST phase boundaries are continuous in the CMD (the He flash is smoothed over; every phase-change step is < 0.05 mag), so the whole thing is one polyline with no breaks.
- `nearest_isochrone_point(G, BP_RP, track, sigma_G=0.2, sigma_color=0.05, chunk=1024)` -> dict with every track key evaluated at the matched point, plus `dist`.
  - **nearest point on the piecewise-linear track** (projection onto each segment, clipped to [0,1]), not the nearest node -- node spacing is very uneven. every track column is linearly interpolated at the projection, so the interpolated `EEP` tells you which stage a star landed on.
  - the metric is `sqrt((dcolor/sigma_color)^2 + (dG/sigma_G)^2)`, i.e. `dist` is in "sigmas". the defaults weight colour 4x more than G: a catalog star painted onto an N-body star carries that star's distance (~0.2 mag spread along the stream) in its absolute G, while BP-RP is what actually sets Teff. only the ratio matters for which point is nearest; the absolute scale only matters for reading `dist`.
  - `dist` > ~3 = not on this isochrone (blue stragglers, BHB stars bluer than the model HB, contaminants) -- the natural "matched" flag for step 6.
  - NaN inputs come back NaN. vectorized in chunks of `chunk` query stars (memory ~ chunk x nodes); 1578 stars take ~0.03 s.
- `Teff_from_gaia_isochrone(track, G, BP_RP, **kwargs)`: one-liner returning just `log_Teff`.

checked: querying at the track's own nodes and at segment midpoints returns `dist` = 0 (to 1e-14) and the exact `log_Teff`. on the `gd1` catalog (1578 stars, the 1535 missing `dist` filled with the median as a stand-in for the N-body distances): median `dist` 1.2, 239 stars > 3, logTeff 3.66-3.88. MS stars project mostly horizontally (Teff set by colour), SGB/RGB stars land on the RGB, the two BHB stars land on the HB. the catalog MS sits ~0.1 mag redder than the isochrone -- probably reddening (`bp_rp` isn't dereddened) and/or [Fe/H] -- and with colour weighted heavily that shifts the assigned Teffs cooler. decide whether to deredden before matching.

### the blackbody fallback (in `paf` and duplicated in `noise.py`, not the live path)
synthetic photometry from a **blackbody**, with **top-hat filters** -- no real response curves. `define_photometric_bands()` (Gaia G/BP/RP from the DR2 paper + a 100 nm z at 900 nm), `integrated_mag(nu_min, nu_max, T, R, d=10)` (AB, **cgs bare numbers**), `get_gaia_photometry`, `get_z_photometry` (scalar only, hence slow loops), `m_from_M(M, dist)` (dist with astropy units -- this one *is* still used by the pipeline). the 1/nu normalization at the `np.trapezoid` lines is marked "suggestion that idk why works" in the source and has never been validated against a real zero point -- don't trust these absolute mags.


# stream frame coordinates ("observed" frames)
`streamframe_coords_observed(orbit, data_dict, prog_tab)` in `inspect_new_sims.py`. `data_dict` is one of the binary-treatment subdicts (`CoM`, `luminous`, `companions`) so it has `pos`/`vel` with astropy units. pipeline is:

`galcen pos/vel` -> `paf.galcen_to_ICRS` -> `gala.coordinates.reflex_correct` -> `.transform_to(<great circle frame>)`

**returns `sc, selected_streamframe` where `sc` is a `SkyCoord`.** so pull attributes off it: `sc.phi1`, `sc.phi2`, `sc.pm_phi1_cosphi2`, `sc.pm_phi2`, `sc.distance`, `sc.radial_velocity` -- all with units, unlike the coord dicts. returning the frame object as well is what lets the progenitor orbit be pushed through the *same* frame downstream. these are for *all* particles, untrimmed.

when it gets flattened to a dict (in `straightened_obscoords_orbit_interp`) the keys are, as bare floats:

| key | unit | |
|---|---|---|
| `phi1`, `phi2` | deg | |
| `pm_phi1`, `pm_phi2` | mas/yr | `pm_phi1` = `pm_phi1_cosphi2` |
| `v_gsr` | km/s | reflex-corrected line-of-sight velocity |
| `distance` | kpc | heliocentric |
| `d_phi2` | kpc | phi2 [rad] x distance |
| `v_phi1`, `v_phi2` | km/s | pm [rad/s] x distance [km], i.e. transverse velocities |

all but `phi1` are residuals about the orbit track. contrast with the intrinsic dict: `distance` is heliocentric where `r` is galactocentric, and `v_gsr` is a real reflex-corrected line-of-sight velocity where `vr` is relative to the progenitor.

`paf.galcen_to_ICRS(pos, vel)` wants shape (N,3) with units (it transposes internally). a single (3,) vector works too and gives a scalar coord, which is how the progenitor origins get built.

## how each frame is defined
| `orbit` | source | how |
|---|---|---|
| `gd1` | Koposov+2010 | `gc.GD1Koposov10`, built into gala |
| `pa5` | Price-Whelan+2018 | `gc.Pal5PriceWhelan18`, built into gala |
| `c19` | Ibata+2024 (see also Mohammed+2026) | `from_pole_ra0`, pole + alpha_0=354.356 deg |
| `jet` | Do+2026 | `from_pole_ra0`, pole + origin |
| `aau` | Shipp+2018 table 1 (the ATLAS half) | `from_endpoints` + `ra0` from the progenitor |
| `m3` | Yang+2023 sec 4.4 | `from_endpoints` + progenitor origin, `priority='origin'` |
| `circ` | -- | no frame; not a real stream |

phi1 zero points are set by the present-day progenitor sky position read out of `data/FINAL_ics_nolmc.csv` (except `jet`, where Do+26 give an origin directly), looked up as `prog_tab[prog_tab['name']==orbit]`. there is no `else` branch -- an unrecognized `orbit` string leaves `selected_streamframe` undefined and you get an `UnboundLocalError` at the `transform_to` line rather than a useful message. (the commented-out `else: raise` was firing for `gd1` because it only paired with the last `if` -- it'd need to be an `if/elif` chain.)

## gala `GreatCircleICRSFrame` gotchas
version in `petar_env` is gala 1.9.1. the modern API is pole+origin only; passing `ra0=` or `rotation=` straight to the constructor raises. use the classmethods.

- **`from_pole_ra0` and `from_endpoints` are not actually different animals.** both just compute a pole and an origin and hand them to the same constructor. if you build a frame one way and then rebuild it from its own `.pole` and `.origin.ra` the other way, you get bit-identical phi1/phi2. so if one "works" and the other doesn't, suspect a typo before suspecting gala.
- **the `priority` kwarg matters a lot.** the pole and origin have to be orthogonal. if your origin is off the great circle, gala emits only a `RuntimeWarning` (easy to miss in a notebook cell) and silently fixes it. default is `priority='origin'`, which *moves the pole* -- i.e. throws away the great circle you defined by endpoints. `priority='pole'` keeps the great circle and projects the origin onto it. passing `ra0=` instead of `origin=` also keeps the pole.
- for `m3` this is not academic: the `M3` row of `FINAL_ics_nolmc.csv` sits ~17 deg off the great circle through the Yang+23 endpoints, so `priority` changes where the endpoints land by >12 deg in phi2. also, `data/bpw25_catalogs/m3.fits` spans ra 189-269 deg with ~7.5 deg of phi2 scatter in every endpoint-derived frame, so one great circle may just not describe the whole M3 stream. TODO in the code: may end up settling for "stream on an M3-like orbit," since that run is really a high-e / low-pericenter test. (M3 is currently dropped from `final_datasets.py`/`gmm.py`.)

## straightening in the observed frame
TODO: I am not 100% convinced that this is successful for M3 or Pal 5, which have notable diverging tails that get added to cocoons.

four functions in `inspect_new_sims.py`, chained by `straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab, Dt=500, lumdict=None, PM_treatment='CoM')`. **returns a dict** (keys in the table above).

**binary handling:** with `lumdict=None` everything is CoM. with `lumdict` given, phi2 / v_gsr / distance come from the luminous component (so the RVs include binary orbital motion), and the proper motions (and hence `d_phi2`, `v_phi1`, `v_phi2`) come from the CoM if `PM_treatment='CoM'` or from the luminous component if `'primary'`. the progenitor orbit track is always the same. `final_datasets.py` uses `PM_treatment='CoM'`.

1. `prog_orbit_track(w0, Dt)` -- integrate the progenitor +/- `Dt` Myr in mwp2014 and concatenate `[backward reversed, forward]`. **note the `[:-1]`**: both integrations contain t=0, so the naive concatenation duplicates the progenitor and puts an exact `dphi1 = 0` step at the midpoint, which stalls any sign-based walk. with it dropped the progenitor sits at index exactly `Dt`, which is also `n//2`. `w0` is the present-day row of `FINAL_ics_nolmc.csv`.
2. `observed_orbit_track(pos, vel, obs_streamframe)` -- same galcen -> ICRS -> reflex_correct -> great circle pipeline as the data. **the reflex correction on the orbit is not optional**: it doesn't change phi1/phi2 but it absolutely changes `v_gsr` and the proper motions.
3. `chop_orbit_track(phi1, jump_threshold=45)` -- walk outward from the progenitor and cut at the first bad step on each side, so what comes back is contiguous, single-valued and increasing in phi1 (ready for `np.interp`). two kinds of bad step:
   - `|dphi1| > jump_threshold`: the +/-180 seam of the great circle frame. these are ~356 deg while real steps at 1 Myr sampling are < 7 deg/Myr for gd1/pa5/aau/c19/jet, so the threshold is not delicate at all.
   - a change of direction: the orbit genuinely doubling back in phi1. **m3 needs this** -- near its low pericenter it hits 31 deg/Myr, so no jump threshold can catch its turnaround.
4. subtract `np.interp(data_phi1, orbit_phi1[idx], orbit_y[idx])` from each key.

the nice property is that the chopped chunk **stops depending on `Dt`** once `Dt` is big enough to reach a seam or a turnaround on both sides. measured on `lm`, `rvir_index=0`, `copy=0` with `jump_threshold=45` -- chunk length at Dt = 200 / 400 / 800:

| orbit | data phi1 range | 200 | 400 | 800 | phi2 residual std |
|---|---|---|---|---|---|
| gd1 | -113.7 to 30.9 | 401 | 608 | 652 | 0.334 |
| pa5 | -24.3 to 26.4 | 332 | 332 | 332 | 0.234 |
| aau | -38.2 to 23.7 | 401 | 703 | 850 | 0.214 |
| m3 | -50.7 to 77.2 | 398 | 479 | 479 | 0.786 |
| c19 | -29.7 to 23.7 | 401 | 508 | 508 | 0.188 |
| jet | -28.2 to 23.8 | 401 | 720 | 720 | 0.120 |

coverage of the trimmed data was 1.000 for all six (measured at Dt=400 with the old intrinsic `inMW`/`trim`; the default is 500). the phi2 residual std column is against the orbit track only, before `poly_straightening`.

`straightened_obscoords_orbit_interp` calls plain `np.interp` with no `left`/`right`, so anything outside the track's phi1 range gets **silently clamped to the edge value** rather than flagged. fine while coverage is complete, but if you change `Dt` or the frame, pass `left=np.nan, right=np.nan` and check for NaNs.


# cocoon separation by gaussian mixture (`gmm.py`)
a two-component gaussian mixture fit to the straightened observed-frame residuals, in place of the old hand-cut `get_cocoon_selection`, so cocoon membership is a fitted probability per star instead of a by-hand threshold per orbit -- the manual cuts didn't standardize across orbits. the functions are n-component general; the main program runs n=2. the convention throughout is **narrowest first, cocoon last**: components `0 .. n-2` are the thin stream and the final (widest) one is the cocoon.

the data vector is `keys = ['phi2','v_phi1','v_phi2','v_gsr']`, stacked `(N, 4)` in that order, so **every `mu` and `sigma` row is a 4-vector in that same order** (and so are the `M_*`/`S_*` table columns). mixing up the column order silently gives a nonsense fit. ⚠️ several docstrings and comments in `gmm.py` (`sigma_ratio_constraint`, `sort_components`, `constraint_dims`) still say `pm_phi1`/`pm_phi2` -- the index positions are the same, the quantities are now transverse velocities.

a persistent issue has been that the thin stream itself is non-gaussian, so when the cocoon doesn't contain enough stars the optimizer just fits two gaussians to the thin stream. this is a particular issue for M3 and Pal 5, with close-in orbits and significant width variations along the stream. the fix is the width-ratio constraint (`constrain_widths=True`); an earlier three-component version that gave the thin stream two narrow gaussians was dropped in favour of it.

**the n-1 convention.** everything that takes or returns `component_fractions` uses length `n_components - 1` -- the last weight is implicit, pinned by `sum(Q) = 1`. so `n_components` is always inferred as `len(component_fractions) + 1`, and **the cocoon fraction is `1 - fracs.sum()`**. `means` and `sigmas` are the full `(n_components, K)`. (the main program appends the implicit weight after sorting, so there `fracs_fit` is full length and `fracs_fit[:-1]` goes back into the functions.)

## the likelihood
each star independently is drawn from one of the n components:

    L = prod_i [ sum_j Q_j p_j(x_i) ],   sum_j Q_j = 1

**the mixture weights go inside the product over stars.** `sum_j Q_j prod_i p_j` is a different and wrong model (the *whole stream* drawn from one component).

everything stays in logs: the inner product over the 4 dims is `norm.logpdf(...).sum(axis=1)`, the outer over stars is the final `np.sum`, and the per-star mixing over components is `scipy.special.logsumexp(axis=0)`, which is exact even when every term is ~ -3000.

- `component_likelihood(x_data, mu, sigma) -> (N,)` -- ln p per star under ONE component. `-inf` for sigma <= 0.
- `gmm_negative_loglikelihood(component_fractions, means, sigmas, data) -> scalar` -- **`data` is last** so `minimize`'s `args` can supply it. out-of-bounds guard returns `+np.inf` (we're minimizing). raises if `len(means)`/`len(sigmas)` disagree with `len(component_fractions) + 1`.
- `component_responsibilities(x_data, fracs, means, sigmas) -> (n_components, N)` -- `R[j,i] = Q_j p_j(x_i) / sum_k Q_k p_k(x_i)`. **order-agnostic** -- `j` only means something physical if you sorted first.
- `component_membership_probability(x_data, fracs, means, sigmas, component=0, sort_dim=None) -> (N,)` -- the indexable one. `component` is an int (negatives work), a sequence, or a slice; sequences/slices *sum* the responsibilities. `sort_dim=-1` makes a mislabelled (unsorted) fit raise.
- `membership_probability(x_data, fracs, means, sigmas, sort_dim=-1) -> (N,)` -- `p_thin`, i.e. `component=slice(0,-1)`. `p_cocoon = 1 - p_thin`. enforces the sorted check.

## measurement errors in the likelihood (`x_err`)
every likelihood function (`component_likelihood`, `gmm_negative_loglikelihood`, `nll_flat`, `component_responsibilities`, `component_membership_probability`, `membership_probability`) takes an optional `x_err`: an `(N, K)` array of per-star 1 sd uncertainties, in the same column order as `keys`. each star's variance in component j becomes `sigma_jk^2 + e_ik^2` (a gaussian convolved with a gaussian error is a gaussian with the variances added), so **the fitted sigmas are intrinsic, error-deconvolved widths**. it's the same idea as `mle_sigma` (Vedant's 1D version), and the 1D one-component case reproduces `mle_sigma` exactly. it's also extreme deconvolution (Bovy, Hogg & Roweis 2011) with diagonal errors.
- `nll_flat(theta, x_data, min_components=2, x_err=None)`. `x_err` comes after `min_components` so `args=(x_data,)` / `(x_data, 1)` work without errors. with errors: `args=(x_data, 2, x_err)` / `(x_data, 1, x_err)`.
- `x_err=None` is the error-free likelihood, and noiseless cases use it.
- **`sigma_binary`**: every likelihood function above also takes an optional, *fixed* (not fit) binary-orbital-motion dispersion, added in quadrature like `x_err`: `sigma_jk^2 + e_ik^2 + sb_ik^2`. it's anything that broadcasts to (N, K) in `keys` order; build it with `gmm.binary_excess_array(sigma_binary, keys, N, apply_to=mask)` (that number in the v_gsr column for the `mask` stars, 0 elsewhere). `nll_flat(theta, x_data, min_components, x_err, sigma_binary)`, so `args=(x_data, 1, x_err, sb)`. `None` (default) leaves the likelihood unchanged. pass the same `sigma_binary` to the membership functions as to the fit. **nothing uses it at the moment**: the repeat-RV experiment now fits sigma_binary as a free parameter with its own multi-epoch likelihood in `single_component_comparisons_mockobs.py` (see its section).
- `gmm.measurement_errors(noise_dict, noise, distances, keys)` builds the full-length array: phi2 = 0 (Gaia position errors are ~1e-5 deg, and they aren't saved); v_phi1 = v_phi2 = `pm_err_gaia` converted to km/s with the CoM distance (same as the noise draw, no distance error); v_gsr = `rverr_<noise>`. `gmm.py` main and `results.py` (`observed_errors` -> `p_cocoon_from_table(x_data, row, x_err)`) both call it, so p_cocoon in the demo panels and Figure 1 uses the same errors as the fit. checked: noise draw / error has std 1.00 +- 0.02 in each velocity dimension.
- **`x_err` is always the 1 sd uncertainty, never the sampled noise.** for the sims the saved noise draw only goes into `x_data` (and the sanity-plot scatter), and `x_err` comes from the error columns (`rverr_<noise>`, `pm_err_gaia`). for the real Jarvis+26 data nothing is added to `x_data`, and `x_err` is the pickle's `errors` dict. (checked both; see the log, 15 Oct.)
- **memberships need `x_err` too.** with the table's intrinsic sigmas and no `x_err`, every noisy star looks too far off-track and gets pushed into the cocoon.
- **with `x_err`, the v_phi1 / v_phi2 sigmas go to ~0**, because the Gaia transverse-velocity errors (median ~7 km/s) are ~10x the true scatter of these stars (~0.6-0.9 km/s). so a ~0 there means "not resolved by Gaia" (an upper limit), not a measurement of 0, and the `sigma_ratio_constraint` on v_phi2 can't mean anything in that regime. without `x_err` the noisy fits' 7-14 km/s in v_phi are the error bars. (before/after numbers in the log, 14 Oct.)
- every noisy table in `data/gmm_tables/` is fit with `x_err`, unconstrained; `jarvis26` with the real errors. `results.py` broadens with `x_err` too.

### a width at ~0 = unresolved, not "approaches the errors"
- **ignoring the errors**, the fitted width tends to the *observed* spread `sqrt(sigma^2 + e^2)`, i.e. ~e when sigma << e. that measures the instrument, not the stream.
- **with `x_err`**, the fit asks how much width is left after the errors are accounted for. for an unresolved width the likelihood is flat from sigma = 0 up to some upper limit, then rises steeply. the best fit often sits exactly at 0, since a variance can't go negative: with equal errors the estimate is `max(0, s^2 - e^2)`, and s^2 scatters by ~`e^2 sqrt(2/N)`. so sigma = e is *excluded* (it predicts sqrt(2) e of observed spread), and a ~0 should be quoted as an upper limit. rough resolution limit for equal errors: sigma >~ e (2/N)^(1/4).
- in the code, "0.000" means ln sigma drifted down that flat floor; the exact value means nothing.
- **the v_phi zeros are genuinely unresolved** (~7 km/s Gaia errors, N ~ 1600 -> resolvable down to ~1.4 km/s, true ~0.35).
- **the DESI thin sigma_v_gsr is resolvable** (GD-1, 0.75, CoM): with the Jarvis errors (median 3.1 km/s, but 44 stars < 1 km/s, which carry most of the information), a 0.46 km/s width gives a median estimate of 0.45 in toy realizations, with 5% at 0. so a 0 there is a noise-draw fluke, not a resolution limit (it happened with seed 42; see the log, 14 Oct). on the current seed-0 draw it's near the N-body value.
- real GD-1 is much hotter than the sim: Jarvis+26 find ~2.5 km/s for the thin stream (`mle_sigma` on their P_THIN > 0.5 stars gives 2.37, against a naive std of 5.85 and a median V_ERR of 3.3; this model fit to their data gives 2.52, see "the real Jarvis+26 GD-1 data"), vs 0.46 km/s for the simulated GD-1 at 0.75 pc. their value being close to the DESI errors is a coincidence; it's a deconvolved width.
- to do: a profile-likelihood helper in `gmm.py` for upper limits (the loop in `claude_diagnostics/claude_profile_desi.py`), and refitting the small-N cases on several noise seeds to see how much of each table entry is the draw.

diagonal sigma only: the model gives up correlations *within* a component (e.g. phi2 with v_gsr along the track). those should be small post-straightening, but that's an assumption worth checking on the residuals.

a positive `ln L` is not a bug: these are log *densities*. it also means `ln L` is only comparable at fixed units. **relatedly: never threshold `component_likelihood`'s output** -- threshold a responsibility.

## fitting it: the packing layer
`scipy.optimize.minimize` wants **one flat 1-D array**, and `args` has to be a tuple (`args=(x_data,)` -- no trailing comma and scipy splats the array). so `pack_params` / `unpack_params` flatten to a **length `(n-1) + 2nK`** vector (17 for n=2, K=4), and `nll_flat(theta, x_data)` is what gets handed to `minimize`. layout:

    [alpha_1 ... alpha_{n-1},  mu_1, ln sigma_1,  mu_2, ln sigma_2,  ...]

the packing **reparameterizes to make the fit unconstrained**: `sigma -> log(sigma)`, and the weights go through a softmax with the last component's alpha pinned at 0 (`alpha_j = ln(Q_j / Q_last)`; for n=2 this is exactly `logit`). this keeps the optimizer off the `inf` walls and fixes the conditioning between dimensions with very different scales.

guards: `unpack_params` **infers `n_components` from `len(theta)`** and raises if the length isn't of the form `(n-1) + 2nK`; `nll_flat` reads `K` off `x_data.shape[1]`.

**`nll_flat(theta, x_data, min_components=1)` is the single-gaussian null model** -- the thing the BIC TODO needs -- with no special-casing: `pack_params(np.array([]), means (1,K), sigmas (1,K))` gives a length-2K theta, the validity guards are vacuously true on the empty fractions array, `Qs = [1.0]`, and the single-row `logsumexp` is the identity. `min_components` defaults to 2 so an *accidentally* emptied fractions array still fails loudly. call it with `args=(x_data, 1)`. (`unpack_params`' docstring calls this `nll_flat_anyn`, which no longer exists.)

`sort_components(fracs, means, sigmas, sort_dim=-1)` handles **label switching** -- returns components ordered **narrowest -> widest** in `sort_dim` (default -1 = v_gsr), with the returned fractions being the n-1 narrowest.

## constraints and bounds
**both only do anything for `SLSQP`, `trust-constr` or `COBYLA`.** Powell and Nelder-Mead accept `constraints=`, emit only a `RuntimeWarning`, then ignore it and return `success=True`. (Powell *does* honour `bounds`.) a constraint function is called as `g(theta)` only -- `args` go to the objective. `constraints=` takes a list.

### `sigma_ratio_constraint(min_ratio, n_components, K, dims, thin_component=0)`
requires the cocoon to be at least `min_ratio` times **wider** than the thin component in each of `dims`. a much better handle than bounding the cocoon fraction (which just clips the answer at the wall): constraining the **separation** says "only call it a cocoon if it's much wider." exactly **linear in theta** (`(ln sigma_cocoon)_k - (ln sigma_thin)_k >= ln R`), one `LinearConstraint` row per dimension, so a **per-dimension ratio** is free: `min_ratio` is a scalar or a list aligned elementwise with `dims`. from Jarvis+26 the GD-1 cocoon is ~10x wider in phi2 but only ~3x in radial velocity, so one scalar is either too weak in phi2 or too strong in v_gsr.

`gmm.py` currently uses `constraint_dims = [0, 2, 3]` (phi2, v_phi2, v_gsr), `min_ratio = [10.0, 5.0, 5.0]`, v_phi1 unconstrained.

⚠️ `thin_component=0` only, and "cocoon = last component" is hardcoded -- during the fit that's a definition. with n=3 you'd need a second constraint for the middle component.

### `pack_bounds` (available, not currently used)
`pack_bounds(fracs_bounds, means_bounds, sigmas_bounds, n_components, K)` pushes natural-space bounds through the same reparameterization as `pack_params`. means pass straight through; sigma bounds become `ln sigma` (so positivity is free); frac bounds are `logit` and only honest for n=2 (it raises otherwise; its docstring points at `cocoon_fraction_constraint`, which is only in `old/develop_GMM.py`). `()` is not the unbounded spelling -- `None` is.

**`gmm.py` runs with `bounds=None`.** the ±1 sd mean box that `develop_GMM.py` used isn't applied, and there's no active-bound warning. instead `results.py` drops fits with `|M_c| >= 0.5 S_c` in any dimension -- i.e. off-track "cocoons" (diverging tails, feather blobs) are filtered after the fact rather than prevented. if they're prevented with bounds again, remember a mean sitting exactly on the wall is the optimizer reporting the bound, not a fit, and scipy says nothing about it; also pass `bounds` to *both* stages or SLSQP silently clips Powell's `x0` back into the box.

## initial guesses
**do not start any two components identical** -- that's an exact saddle point (zero gradient wrt the weights; confirmed, L-BFGS-B returns f_thin = 0.5 with bit-identical sigmas). and scale off `x_data.std(axis=0)` since the dimensions have mixed units.

`gmm.py` uses `sigma_0 = [0.1, 10.0] * sd`, both `mu = 0`, and `f_1 = 0.9` -- the symmetry is broken twice over. the 0.9 is a soft prior that the stream is mostly thin; starting at 0.5 splits the thin stream itself in two and calls the wider half a cocoon.

⚠️ **multimodality.** the likelihood surface is multimodal, `minimize` only finds a local optimum, and **no `success` flag warns about it**, so the reported cocoon fraction is conditional on the starting basin. the independent check is physical: cocoon fraction should *decrease* with increasing rvir (decreasing initial density) -- which is exactly what the `results.py` summary panels show, so a non-monotonic line there is the signal. still worth doing: scan the inits over a grid and keep the best `result.fun`. the profile likelihood in `f_cocoon` (commented out in `old/develop_GMM.py`) is the honest version of this check -- flat toward `f_cocoon -> 0` means report an upper limit, not a value; delta-nll ~0.5 is 1-sigma.

## optimizer: the two-stage recipe
**SLSQP from a cold start collapses this likelihood onto a single component.** so:

    # stage 1: Powell, to land in the right basin (ignores constraints, honours bounds)
    result_free = minimize(nll_flat, theta0, args=(x_data,), method='Powell',
                           bounds=bounds, options={'maxiter': 100000, 'maxfev': 100000})
    # stage 2 (constrain_widths=True only): re-fit from there WITH the constraint
    result = minimize(nll_flat, result_free.x, args=(x_data,), method='SLSQP',
                      bounds=bounds, constraints=constraints, options={'maxiter': 5000})

the constrained nll is necessarily `>=` the free one, and the gap is how hard the data resist being told the cocoon must be much wider -- worth printing/saving, `gmm.py` currently doesn't.

benchmarked on synthetic data (n=2, K=4, N=30000, unconstrained):

| method | wall time | final nll | fitted f_thin (true 0.7495) |
|---|---|---|---|
| Nelder-Mead | 312 s | 34121.9 | 0.6915 |
| L-BFGS-B | 64 s | 33168.1 | 0.7492 |
| Powell | 5 s | 33168.1 | 0.7491 |

**all three returned `success=True`.** so: **check `result.fun`, never `result.success`** -- which is also `True` for an ignored constraint, a clipped `x0`, and a fit pinned on a bound. with the good optimum, recovered sigmas match truth to <1% and per-star label recovery is 99.7%.

## model choice: BIC
`bic(lnL0, k, N) = -2 lnL0 + k ln N`, copied from `old/develop_GMM.py`. lower is preferred. in the main loop, right after the two-component fit:
- the single gaussian is fit with **Powell** on `nll_flat` with `args=(x_data, 1)` (the `min_components=1` route; see the packing layer section), starting from `mu = 0`, `sigma = sd`, i.e. the same data scaling as the two-component guesses. its maximum likelihood is analytic (sample mean and std), and Powell reproduces it to 1e-3 in nll on synthetic data, so the fit is just a consistency check with the same machinery.
- `BIC_one = bic(-result_one.fun, len(result_one.x), N)`, `BIC_two = bic(-result.fun, len(result.x), N)`, with `N = len(x_data)`. k = 8 and 17.
- in the `constrain_widths=True` cases, `result` is the SLSQP (constrained) fit, so `BIC_two` uses the constrained likelihood; k is still counted as 17 (the inequality constraints don't remove parameters).
- on synthetic data (N = 500) it picks one component for a single gaussian (BIC 7356 vs 7412) and two for a 70/30 thin + cocoon mix (10068 vs 7756).
- `ΔBIC = BIC_one - BIC_two`: > 0 favours a cocoon; the usual scale treats |ΔBIC| > ~10 as strong. the two-component optimum is only a local one (see multimodality above), so a two-component fit stuck in a bad basin can make `BIC_two` look worse than it really is.

note `old/develop_GMM.py`'s `AIC` is actually AICc (includes `2k(k+1)/(N-k-1)`; numerically irrelevant at N ~ 1e4). for the noisy cases N gets small enough (hundreds, for DESI) that the choice of criterion can matter.

## binaries in the observed samples: a selection effect -- for the paper
**the noisy `binaries` cocoons have a much larger sigma_v_gsr than the noiseless `binaries` cocoon, and it's mostly not noise -- it's which stars get observed.** the bright, catalog-matched stars are about twice as likely to be binaries as the stream as a whole. this is a real effect (a bright, magnitude-limited stream sample *would* be binary-rich), so the decision is to **keep it in the observed cases**. but understand it, quote it, and don't present noisy-vs-raw differences in the `binaries` cases as the effect of measurement noise.

GD-1, hm, rvir0 = 0.75 pc, cocoon sigma_v_gsr [km/s] from the tables:

| | noiseless | via | desi |
|---|---|---|---|
| `CoM` | 2.32 | 1.97 | 4.62 |
| `binaries` | 6.40 | 10.28 | 12.19 |

i.e. binaries add ~4 km/s raw but ~8 km/s observed.

**why.** `dv` = binary orbital motion in v_gsr (`sc_straighter_primaries['v_gsr'] - sc_straighter['v_gsr']`), same pickle, `binaries`-case trim:

| sample | N | binary fraction | std of dv over binaries [km/s] |
|---|---|---|---|
| `trim_new_primaries & unbound` (noiseless) | 6469 | 0.118 | 6.8 |
| catalog-matched (`via`) | 1578 | 0.216 | 7.7 |
| catalog-matched, top 679 in mG (`desi`) | 679 | 0.236 | 9.6 |
| eligible pool, isochrone-brightness ranks 679-1578 | 899 | 0.201 | 5.4 |
| eligible pool, ranks 1578-3313 | 1735 | 0.111 | 5.9 |
| eligible pool, ranks 3313-5049 (faintest rankable) | 1736 | 0.064 | 6.5 |

the binary fraction falls steadily with brightness. the dv spread is noisier, since a std is dominated by a few very tight binaries: ranks 4000-5049 alone give 1.7 km/s.

- the catalog matching gives photometry to the brightest pool stars, ranked by the isochrone mG of the **luminous component's** `m0_zams`. a binary's luminous component is the heavier of two stars, so binaries move up the ranking: median luminous m0_zams 0.47 Msun in the matched sample vs 0.37 in the full unbound sample. the full unbound sample is dominated by faint, mostly single M dwarfs.
- heavier systems also have larger orbital velocities at fixed separation, so the dv spread grows with brightness too.
- dispersions add in **quadrature**, not linearly, and the GMM puts nearly all the binary variance into the cocoon. so binaries don't add a fixed amount; roughly `sigma_c,bin^2 ~ sigma_c,CoM^2 + f_bin <dv^2> / f_cocoon`. that predicts a binary term of 6.0 / 10.3 / 11.8 km/s (noiseless / via / desi), against 6.0 / 10.1 / 11.3 km/s measured as `sqrt(sigma_c,bin^2 - sigma_c,CoM^2)`.
- `trim_new_primaries` (vs the CoM trim) helped only a little: it drops the std of dv in the matched sample from 10.9 to 7.7 km/s by clipping the extreme tails, but it doesn't touch the binary fraction.

**the check that isolates it:** a **noiseless** GMM fit (Powell, same initial guess as `gmm.py`) on exactly the stars each survey case keeps:

| noiseless fit on | f_cocoon `CoM` | sigma_v_gsr `CoM` | f_cocoon `binaries` | sigma_v_gsr `binaries` |
|---|---|---|---|---|
| all unbound | 0.127 | 2.32 | 0.152 | 6.40 |
| via-matched 1578 | 0.107 | 2.09 | 0.154 | **9.17** |
| desi top 679 | 0.134 | 1.86 | 0.175 | **11.08** |

so with zero noise, selection alone takes the `binaries` cocoon to 9.2 / 11.1 km/s, against 10.3 / 12.2 with noise. noise adds only the last ~1 km/s. the `CoM` cocoon hardly moves with the selection. all of this comes from `claude_diagnostics/claude_binary_selection_experiment.py` (`--orbit`, `--rvir`; defaults gd1 / 0.75; ~1 min). it builds the matched mask the same way `gmm.py` does, but leaves out the pm/RV error cuts, which don't remove anything. then it fits `sc_straighter{sfx}` without noise.

**to come back to:**
- add a case to `gmm.py`: `noise=None`, but with the catalog-photometry `use` mask, i.e. "noiseless, observed stars". then noisy vs this separates noise from selection, and this vs all unbound is the selection effect alone.
- is the sim's binary fraction vs primary mass realistic? it's set by how the petar binary ICs were generated / paired; check before leaning on the size of the effect in the paper.
- the companion's light isn't added to the luminous component (see the photometry caveats), so binaries are ranked slightly too faint. adding it would push the matched binary fraction up a little more.
- not the same "cocoon" in every case: in `via_noise_CoM` the cocoon is partly the faint stars with large Gaia PM errors (cocoon sigma_v_phi1 13.5 vs thin 7.8 km/s), while in `via_noise_binaries` the split moves to v_gsr (cocoon sigma_v_phi1 7.1 is *smaller* than the thin 8.4). so differencing the two cocoons' sigma_v_gsr is a bit apples-to-oranges.
- only measured for GD-1 rvir0 = 0.75. check it holds across orbits/rvir; the small catalogs (c19 46, pa5 129, jet 144) will be noisier.

## next
maximum likelihood first, then emcee for posteriors. `nll_flat` negated is the right sign for `log_prob`, the softmax/log-sigma parameterization is already unconstrained, and the `inf` guards are what emcee expects for a rejected step -- add a prior and it's ready. a sampler needs `sort_components` applied **per sample** or the marginals come out as mush.


# bugs / stale things that remain
in rough order of how much they'd hurt:
- `results.py` reads the unsuffixed catalog keys, so its noisy `binaries` panels don't match `gmm.py`'s selection (see the `results.py` section).
- **GMM label switching / degenerate fits in the catalog-noise tables** (examples below are from the first catalog-noise fits, before `x_err`; not rechecked on the current tables). `sort_components` orders by v_gsr width only, so when the two components have nearly equal v_gsr sigma the label is a coin flip: `pa5` rvir 6 in `via_noise_CoM` / `via10hr_noise_CoM` has the "thin" component at sigma_v_phi1 ~ 34 km/s and the "cocoon" at ~5.7, with v_gsr sigmas 1.08 vs 1.14, so the reported `f_cocoon` = 0.61 is really the thin-stream fraction. and some small-N fits have the thin component collapse onto a handful of stars: `via10hr_noise_CoM` `jet` rvir 1.5 (`f_cocoon` = 1.000) and `c19` rvir 0.75 (`f_cocoon` = 0.935, thin sigma_phi2 = 0.01 deg). **all three pass `results.py`'s `okay_mean` filter**, so they'll show up in the summary plots. the BIC TODO is the principled fix; sorting on a combined width (e.g. the product of sigma/sd over dims), or requiring agreement across dims, would catch the pa5 case.
- **in the `binaries` cases the "cocoon" is mostly the binary-orbital-motion population**: `via_noise_binaries` / `desi_noise_binaries` cocoons have sigma_v_gsr ~ 6-20 km/s against ~1 km/s for `CoM`, with thin-stream-like sigma_phi2. expected physics, but it means `f_cocoon` in those tables isn't a spatial/kinematic cocoon fraction, so don't compare it directly to the CoM tables. and the noisy `binaries` cocoons are hotter than the noiseless one mostly through selection (bright stars are more often binaries), not noise -- see "binaries in the observed samples" in the GMM section.
- in `gmm.py`'s plotting block, `cut = 3*sigmas_fit[-1][jj]` reuses the name of the catalog-photometry `cut` mask. harmless now because `cut` is reloaded from the pickle at the top of each rvir iteration, but it'll bite if that load ever moves.
- **pickle filenames don't include stellar_pop** -- an `lm` run of `final_datasets.py` would silently overwrite the `hm` pickles. (`gmm.py`'s `mass_index = 1` line says "LOW mass" but 1 is `hm`, and it isn't used there anyway.)
- **pickles live on netscratch**, which is purged on a 90-day clock (the grid set was written 4 Nov 2026, so it ages out around early February 2027; the supplemental ones are older, 11-29 Oct, and go around mid-January). copy them somewhere permanent before they age out.
- **`star.mass0` at the present day is not the ZAMS mass** -- use `m0_zams`. fails silently. see its own section.
- **the GMM likelihood is multimodal** and `minimize` only finds a local optimum -- `success=True` either way. check `result.fun` and the `f_cocoon` vs rvir trend.
- **GMM tables don't record which `keys` were fit**, and `results.py` has "translate back to angle from distance" TODOs that only make sense for a `d_phi2` fit -- the current `gmm.py` fits `phi2` in degrees. tables made before the switch may not match the current code.
- **`in_rtid` is the wrong length for the `companions` subdict** (fine for `CoM` and `luminous`, which is all the pipeline uses).
- **`get_init_displacements.py` silently writes an empty file** -- `names_to_run` still uses the long stream names.
- **anything cached from before the init_displacement units fix is wrong**, including the intrinsic `coords`. regenerate.
- stale comments/docstrings in `gmm.py`: `pm_phi1`/`pm_phi2` where it's now `v_phi1`/`v_phi2` (see the GMM section), `nll_flat_anyn` in `unpack_params`, `cocoon_fraction_constraint` in `pack_bounds` (only in `old/`), the `noise` options comment missing `via10hr`. `outlier_clip`'s docstring thresholds don't match its code.
- the main program reads the cocoon as index `[1]` / `[-1]` -- fine at n=2, wrong if n=3 comes back (the cocoon could be two of three components).
- `straightened_obscoords_orbit_interp` clamps instead of flagging outside the orbit track's phi1 range.
- `straighten_stream_polynomial`'s `trim_criteria=None` default is a `TypeError`, not a default.
- no `else` branch in either `retrieve_sim_info` (`UnboundLocalError`) or `streamframe_coords_observed` for an unrecognized orbit string.
- `load_coords_v2` builds the all-particles filename from the raw `file_index` kwarg, so `file_index=None` + `load_all=True` opens `data.None`.
- `core_to_galcen_frame` still adds the raw `core.vel` rather than the `fix_core_vel` version.
- `correct_core` is dead code that also strips units.
- `animations/grid_movie_long.py` `NameError`s at import (`sys.path.append(rotated_pos)`), and the three grid movie scripts use the purged `scratch=True` grid.
- the `m3` great circle doesn't describe the whole stream (phi2 residual std 0.786 vs ~0.2 for everything else).
- the isochrone age (requested 12 Gyr, actually 12.6 Gyr, and v/vcrit = 0.4 by artpop's default) doesn't match the dynamical ages; the isochrone cut (`alive`) throws out every star above 0.793574 Msun (artpop isochrone).
- `straighten_stream_orbit_interp_arbitrary_frame` in `paf` is a commented-out stub.


# logistics:
## simulation data:
- storage (default, `extended_grid_info(scratch=False)`): `/n/holystore01/LABS/itc_lab/Users/amphillips/extended_grid/` -- complete, not purged.
- scratch: `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/` -- 90-day purge, actively losing snapshots. don't use. **exception:** the supplemental sims (`gd1_2xbinaries/`, `gd1_2xN/`) only live here so far -- copy them to holystore.
- the circular orbit simulation is from the [Phillips+26](https://iopscience.iop.org/article/10.3847/1538-4357/ae680b) grid, stored at `/n/holystore01/LABS/conroy_lab/Lab/amphillips/finished_grid/` in the directories that begin with 0-7 (see `extended_grid_info`).

## derived data products:
- pickles (`final_datasets.py`): `/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/<orbit>_<rvir:.2f>.pickle` -- scratch, see the purge warning. the real-data pickle is `jarvis26.pickle` in the same directory (cheap to rebuild with `--which jarvis`).
- GMM sanity plots (`gmm.py`): `/n/netscratch/conroy_lab/Lab/amphillips/p27_sanity_plots/<case_name>/` -- one subdir per case, must exist before running.
- movie frames (`animations/`): `/n/netscratch/conroy_lab/Lab/amphillips/movies/<movie name>/`. `particle_spray_vs_nbody/` also holds the cached spray (`spray_gd1_hm_0.75_0.h5`); scratch, so regenerate with `--make-spray` if it's been purged.
- particle sprays for `motivation.py` (`particle_spray.py`): `data/particle_sprays/` in the repo directory but **gitignored** (large), so not in version control -- regenerate with `sbatch particle_spray.slurm`.
- GMM tables: `data/gmm_tables/` in the repo.
- Jarvis+26 table 7 (DESI GD-1 members, for the mag-distribution comparison and the photometry-matching TODO): `data/jarvis26_Table7.fits` in the repo (identical copy at `/n/home02/amphillips/data/jarvis26_Table7.fits`).
- Bonaca & Price-Whelan (2025) Gaia member catalogs: `data/bpw25_catalogs/<orbit>.fits`.

## conda environment:
`petar_env`, stored at `~/.conda/envs/petar_env`, containing standard packages like numpy, scipy, matplotlib, astropy, but in particular gala (1.9.1), plus `petar`, `sklearn`, `pygaia`, `artpop` (MIST isochrone fetching), and `jax` (0.4.35) + `numpyro` (0.19.0, for the NUTS cell in `single_component_comparisons_mockobs.py`). `viamock` (Via RV errors) is **not installed** -- `noise.py` and `final_datasets.py` put `~/software/viamock` on `sys.path`.

nothing outside this env can import gala or petar, so run scripts with `~/.conda/envs/petar_env/bin/python` (or activate the env) rather than the system python. importing `inspect_new_sims` (which `noise`, `gmm`, `final_datasets`, `results` all do) builds `extended_grid_info()` and reads `FINAL_ics_nolmc.csv` at import time.

## petar documentation
the `README.md` from Long Wang's [PeTar github](https://github.com/lwang-astro/PeTar) is useful.
