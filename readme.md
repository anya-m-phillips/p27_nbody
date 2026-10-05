Exploring an N-body grid of mock streams to probe effects of progenitor dynamics on stream dispersion. Questions to answer might include:
- cocoon fraction and dispersion/extent as a function of progenitor density, orbit: Jarvis+26, Carlberg+26 discuss cocoon structures around MW streams, created by 
- structure of the stream velocity dispersion profile: created by the (evolving) cluster dispersion+tidal filling factor and the progenitor orbit/overlap of successive energy ``feathers." see Bovy(2014).
- **detectability**: under realistic survey noise + down-sampling (Via, DESI), when is a cocoon actually recoverable? a one-component model being preferred = a cocoon non-detection (see the `gmm.py` TODO).


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
                     results.py  -> plots/summary_*.pdf etc.
```

`inspect_new_sims.py` and `noise.py` are function libraries (their `__main__` blocks are exploratory leftovers); `final_datasets.py` and `gmm.py` are the programs that produce the data products. `gmm.py` and `results.py` only read the pickles, so nothing downstream of `final_datasets.py` needs the raw simulation data.

# TODOs
the live ones, collected here. details are in each script's section below.
- ~~**`noise.py` / photometry:** star-by-star catalog photometry~~ **done** (30 Sep): `final_datasets.py` calls `noise.assign_photometry_from_catalog` against the BPW25 Gaia catalog for each stream and saves a second photometry + noise set (`catalog_photometry`, `noise_catalog_photometry`), and every noisy `gmm.py` case now uses those. open questions it left behind are in the "catalog photometry" subsection under `final_datasets.py`: which N-body stars are eligible to be matched, dereddening, `p_mem`, and the fact that the noisy sample size now *is* the catalog size (46 stars for C-19).
- ~~**`final_datasets.py`:** re-run all orbits with the `_primaries` catalog keys~~ **done** (5 Oct, pickles 17:29-17:33).
- **binaries selection effect (for the paper):** the noisy `binaries` cocoons are much hotter in v_gsr than the noiseless one, mostly because the bright catalog-matched stars are ~2x as binary-rich, not because of noise. keeping it in the observed cases; still to do: a "noiseless, observed stars" `gmm.py` case, and checking that the sim's binary fraction vs mass is realistic. see "binaries in the observed samples" in the GMM section.
- **`gmm.py` / model choice:** the BIC is **in** (30 Sep): `gmm.py` fits a single gaussian too and writes `BIC_one` / `BIC_two` columns (see the model choice section). `BIC_one < BIC_two` = one component preferred = a **non-detection of a cocoon**. still to do: re-run every case so the tables have the columns, use them in `results.py`, and maybe save the one-component params, `N = len(x_data)` and `result.fun` too. this matters more with catalog photometry, because the noisy fits have N = 46 (c19), 129 (pa5), 144 (jet), 334 (aau), 1578 / 679 (gd1 via / desi), and several of the small-N fits have already come out degenerate (see the GMM results bullet in the bugs list).
- `gmm.py`: AAU at rvir0 = 6 pc gives a weird result; traced to the 2 rtid tidal boundary, and the decision is to live with it rather than relax the cut.
- `gmm.py`: run all noise cases (see the case table below); rewrite the results notebook (`results.py`); add stuff to overleaf.
- `get_init_displacements.py`: rewrite so orbit rows are looked up by the short keys (see its section -- right now it silently writes an empty file).
- observed-frame straightening: not 100% convinced it's successful for M3 or Pal 5, which have notable diverging tails that end up in the cocoon. (M3 is currently dropped from the pipeline anyway.)
- `m3` frame: may settle for "stream on an M3-like orbit" (see the gala gotchas section).
- `animations/`: make the movie scripts array jobs (one sub-job per frame) instead of looping. (`particle_spray_vs_nbody.py` already is.)
- GMM: emcee for posteriors after maximum likelihood (see the GMM "next" section).


# map of code
## `/data`
- `FINAL_ics_nolmc.csv`: credit Vedant Chandra, present-day positions, velocities of some promising streams. Age and progenitor mass estimates also given. **the `name` column for the six run orbits has been renamed to the short keys i index by** (`gd1`, `pa5`, `aau`, `m3`, `c19`, `jet`) so that `prog_tab[prog_tab['name']==orbit]` just works. everything else in the table keeps its original name.
- `init_displacements.txt`: generated in `get_init_displacements.py`, Vedant's progenitor locations back-integrated by their stream ages, plus 100 Myr to account for initial expansion due to massive star evolution, rounded to a multiple of 10 so that I can safely output sim snapshots every 10 Myr and get the present day in the final snapshot. **positions in this file are in pc** (that's what `petar.init` wants), velocities in km/s -- see the units bullet in the `paf` notes, this caused a real bug.
- `bpw25_catalogs/<orbit>.fits`: credit Bonaca & Price-Whelan (2025), one per stream, named by the short keys (`gd1`, `pa5`, `aau`, `m3`, `c19`, `jet`); loaded by `noise.load_gaia_catalog(orbit)`. Gaia columns (`phot_g_mean_mag`, `bp_rp`, astrometry) plus `dist`, `Vr`, `FeH`, `p_mem`. **`dist` is mostly masked** -- only 43 of 1578 `gd1` rows have one -- and there's no extinction column, so `bp_rp` is presumably not dereddened. (`stream_catalogs.py` still reads the old long filenames, e.g. `GD-1.fits`.)
- `jarvis26_Table7.fits`: Jarvis+26 table 7, DESI GD-1 members (679 stars). also at `/n/home02/amphillips/data/`.
- `gmm_tables/<case>.fits`: output of `gmm.py`, one table per case (see `gmm.py` below). as of 30 Sep, the `via`, `via10hr` and `desi` tables come from the catalog-photometry pickles. the `noiseless_*` tables are older (25-27 Sep) and predate today's pickles. the noisy `*_binaries` tables were fit with the CoM `trim_new` and CoM-trim catalog photometry, which `gmm.py` no longer does (see the trim bullet under `final_datasets.py`). re-run everything before comparing noisy with noiseless.
- `data_dicts/`: empty. the pickles actually live on scratch (see logistics).

## `/scripts`
- `PETAR_ANALYSIS_FUNCTIONS.py`: a bunch of functions, mostly migrated from other projects; imported as `paf`. see the notes section below for the parts that are actually load-bearing.
- `streamframe.py`: credit Jake Nibauer, transformation to stream frame coordinates as seen from the Galactic center
- `read_mist_models.py`: the MIST team's reader (`ISOCMD`). no longer on the live path -- the isochrone now comes from `artpop` -- but `noise.py` still imports it.
- `vedant.mplstyle`: credit Vedant Chandra, plotting style stuff. note it sets `savefig.format : pdf`, `savefig.dpi : 300` and `savefig.bbox : tight`, so movie frames saved as png should pass a lower `dpi` (a 14x14 inch figure at 300 dpi is ~4200 px a side).

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
`summary_<case>.pdf` and friends, made by `results.py` (the `savefig` lines there are mostly commented out, so these are only regenerated by hand). `orbit_summary.pdf` is from `get_init_displacements.py`. `plots/old/` is from the hand-cut and early-GMM era.

## `get_init_displacements.py`
integrates every orbit in `FINAL_ics_nolmc.csv` for peri/apo/e (`plots/orbit_summary.pdf`), then back-integrates the run orbits by age + 100 Myr and writes `data/init_displacements.txt` (pc, km/s) for `petar.init`.

**stale, and silent about it:** `names_to_run` (and the `nl` label list) still hold the long names (`'ATLAS-Aliqa Uma'`, `'GD-1'`, ...) which no longer match the renamed `name` column, so re-running it right now overwrites `init_displacements.txt` with an *empty* file and no error. the header TODO says so; fix the list before touching it.

## `inspect_new_sims.py` (imported as `simspect`)
functions for loading N-body data, putting it in observed stream frames, and straightening the tracks. contents:
- `prepare_nbody_data(path, include_photometry, i, apo, init_displacement, **kwargs)`: wraps `paf.load_core` + `paf.intrinsic_stream_data_v3`; assumes snapshots every 10 Myr. returns `core, data_dict, CMdict, lumdict, inMW, trim` (plus `G, BP, RP, z` if `include_photometry=True`, which is the slow per-star blackbody path, not the isochrone). extra kwargs go to `intrinsic_stream_data_v3` -- this is how `N_rtid_boundary` gets in.
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

the `# SCRATCH::::` cells between `trim_obstream_percentile` and `assign_photometry_from_catalog` (the catalog-photometry prototyping) are commented out now, so `import noise` no longer loads a simulation or draws plots.
- `gaia_g_to_lsst_z` / `lsst_z_to_gaia_g`: the RTN-099 sec 1.3.4 polynomial and its exact inverse (0.014 mag rms).
- `desi_RVerr(zmag, feh=-2.0)`: `log10(err) = -0.47 + 0.27(z-16) - 0.23[Fe/H]`, from **apparent** z mag. no Teff dependence (TODO in the docstring).
- `via_RVerr(G, log_Teff, feh=-2.0, exptime_s=3600, nexp=1, **kwargs)`: thin wrapper on `viamock.get_viaspec_errors`, apparent G. viamock's tables cover G in (5, 30); outside that `final_datasets.py` leaves NaN.
- `get_Teff(R, L)`: Stefan-Boltzmann.
- `build_isochrone_table` / `isochrone_photometry` / `gaia_from_isochrone`: ZAMS mass -> Gaia mags + `log_Teff` by interpolating a MIST isochrone. **read the precision subsection before touching them.** see the photometry section.
- `isochrone_cmd_track` / `nearest_isochrone_point` / `Teff_from_gaia_isochrone`: the reverse direction -- a catalog star's (absolute G, BP-RP) -> nearest point on the isochrone -> `log_Teff` (for viamock). see the CMD-space lookup subsection of photometry.
- `load_gaia_catalog(orbit)`: reads `data/bpw25_catalogs/<orbit>.fits`.
- `assign_photometry_from_catalog(gaia_iso, distances, catalog, track, N=None, eligible=None, **nearest_kwargs)` -> `(phot, matched)`. `gaia_iso` is the `phot` dict from `isochrone_photometry` (only its G is used, to rank); `distances` a Quantity; `track` from `isochrone_cmd_track`. ranks the N-body stars (finite iso mG, and `eligible` if given) and the catalog (finite G and BP-RP) by apparent G and hands the brightest N N-body stars the catalog's **apparent G and BP-RP exactly**; converts to absolute with the N-body distances; `log_Teff` from `nearest_isochrone_point`. `phot` has the same keys as `final_datasets`' photdict (`G`, `BP_RP`, `z`, `mG`, `mz`, `log_Teff`) plus `iso_dist` and `catalog_index` (-1 unmatched), all full length in the original order, NaN where unmatched. `N=None` = every usable catalog row; capped at the number of rankable N-body stars (prints if so).
  - **pass `eligible`** (e.g. `unbound & nonrem & trim`), or bound progenitor stars compete for the catalog's rows.
  - measured on `gd1` hm rvir_index=0 with `eligible = unbound & nonrem & trim_new` (6273 stars): all 1578 catalog rows matched, but the matched N-body stars have **isochrone** mG 15.0-23.6 against the catalog's 12.2-20.1 -- the sim has fewer bright stream stars than GD-1, so rank matching brightens each star by ~3 mag, and the catalog-based `log_Teff` runs a median 0.15 dex hotter than the isochrone value for the same star. that's the scheme working as designed, but it means `N` (and the eligible pool) sets how deep into the N-body MS the catalog reaches. `iso_dist` median 1.7, 267 of 1578 > 3.
- `get_gaia_photometry`, `g_phot`: the blackbody + top-hat fallback, duplicated here from `paf` "to illustrate the difference in stellar populations." not used by the pipeline.

`assign_photometry_from_catalog` is on the live path now: `final_datasets.py` calls it once per sim, and its output is what the noisy `gmm.py` cases fit. see the catalog photometry subsection under `final_datasets.py`.

## `final_datasets.py`
assembles one data dictionary per simulation and pickles it. imports `paf`, `simspect`, `noise`. no `__main__` guard (it's a script). what it does, per orbit x rvir:
- **which sims:** `orbits = ['gd1','aau','pa5','jet','c19']` -- **M3 is dropped**, circ is not included. `stellar_pop='hm'` only (`mass_index = 1`). copies tried in order 0..4, except `aau`, which tries 4..0 to dodge a misbehaving copy. (the `if orbit != 'gd1': continue` used to re-run GD-1 alone is still there, commented out.)
- **progenitor removal:** `N_rtid_boundary = 2.0` for every orbit -- "unbound" = outside 2x the petar tidal radius, to make sure the progenitor is fully removed and stars that get recaptured as rtid re-expands aren't counted. (an orbital-phase-dependent 1x/2x choice was tried and abandoned. the orbital phase is still computed, though no longer printed, and it's computed at the wrong time -- see the bugs list.)
- **coordinates:** `coords_obs` (CoM `SkyCoord`), then two straightened versions: `sc` from CoM, and `sc_primaries` from the luminous component with `PM_treatment='CoM'` hardcoded -- i.e. binaries get the **primary's** instantaneous RV/phi2/distance but the **CoM** proper motions (decided 24 Sep). both are then `poly_straightening`-ed. (`streamframe_coords_observed` is called twice in a row with the same arguments, which is harmless but redundant.)
- **trim:** `trim_new = trim_obstream_percentile(sc, p=[1,99], trim_keys=[phi1, d_phi2, v_phi1, v_phi2, v_gsr, distance])` -- a 1-99 percentile clip in *every* observed-frame dimension, AND-ed, computed on the **CoM** coords. `trim_new_primaries` is the same clip computed on `sc_primaries`, i.e. **including binary orbital motion** in phi2 / v_gsr / distance. decided 5 Oct: the trim operates in the space the GMM will see, so `CoM` cases use `trim_new` and `binaries` cases use `trim_new_primaries`, both noiseless and noisy -- which is why the catalog photometry is done once per trim (see that subsection). `sc_primaries` is still poly-straightened with `tc=trim_new`, not `trim_new_primaries`. this replaces both the intrinsic-frame `inMW`/`trim` and the hard `outlier_clip` cuts; `inMW_na` is an all-True placeholder so the two-mask functions still work. the trim only picks which stars `poly_straightening` fits to -- everything is kept in the saved arrays.
- **isochrone photometry** (`phot`): [Fe/H] = -2 MIST v1.2 isochrone from `artpop.fetch_mist_iso_cmd(log_age=log10(12e9), ..., phot_system='UBVRIplus')` -- which actually returns the **12.6 Gyr (log age 10.10), v/vcrit = 0.4** grid isochrone, see the photometry section -- interpolated at `lumdict['m0_zams']` (the luminous component's birth mass for binaries -- the companion's light is not added). absolute G, BP-RP, z (via RTN-099), apparent `mG`/`mz` using the CoM distances, and `log_Teff`. still computed and saved, and it sets the **ranking** for the catalog matching below, but no `gmm.py` case fits it any more.
- **isochrone noise** (`noise`): Gaia DR3 end-of-mission position and PM errors from `pygaia` (total / sqrt(2) per component, in **µas and µas/yr**); DESI RV errors from `mz`; Via RV errors for 1 hr/1 exposure and 10 hr/10 exposures. saved for comparison; unused by `gmm.py`. all three noise sets come from `survey_noise(mG, mz, log_Teff, rng) -> (noise_dict, acceptable_G_viamock)`, defined at the top of the file.
- **catalog photometry** (`catalog_photometry[_primaries]`) + **catalog noise** (`noise_catalog_photometry[_primaries]`): the live path for every noisy case, done once per trim. see the subsection below.
- noise realizations are drawn here, once, from one `default_rng(seed=42)`, and saved, so every `gmm.py` case sees the same draw. order on the generator, per sim: isochrone noise, then CoM-trim catalog noise, then `_primaries` catalog noise (and `survey_noise` itself draws via, via10hr, desi, pm_phi1, pm_phi2, phi1, phi2 in that order). so changing anything upstream changes every realization after it. the `_primaries` pass was added last (5 Oct) so the existing CoM draws didn't change -- but the generator carries over between sims, so every sim after the first one does get a different draw.
- **output:** `/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/<orbit>_<rvir:.2f>.pickle`. note the filename has no stellar_pop in it, so an `lm` run would overwrite the `hm` pickles. and it's netscratch, so the 90-day purge applies.

what's in each pickle -- the `intrinsic_stream_data_v3` dict (so `CoM`/`luminous`/`companions` subdicts, `nsingles`, `nbinaries`, `IDs`, ...) plus:

| key | what | length |
|---|---|---|
| `coords_obs` | CoM `SkyCoord` in the observed frame | N = nsingles + nbinaries |
| `sc_straighter` | straightened coord dict, CoM | N |
| `sc_straighter_primaries` | straightened coord dict, primary RV/phi2/distance + CoM PMs | N |
| `unbound` | `~CMdict['in_rtid']` at 2 rtid | N |
| `nonrem` | `lumdict['type'] < 10` | N |
| `alive` | on the isochrone's initial-mass range | N |
| `acceptable_G` | 5 < mG < 30 (viamock's range) | N |
| `inMW_na`, `trim_new` | all-True placeholder; the percentile trim on CoM coords (used by `CoM` cases) | N |
| `trim_new_primaries` | the percentile trim on `sc_primaries`, i.e. with binary orbital motion (used by `binaries` cases) | N |
| `phot` | isochrone photometry: `G`, `BP_RP`, `z`, `mG`, `mz`, `log_Teff` (absolute unless prefixed `m`) | N |
| `noise` | isochrone-based noise: sampled `phi1`, `phi2` [deg], `pm_phi1`, `pm_phi2` [mas/yr], `v_gsr_via`, `v_gsr_via10hr`, `v_gsr_desi` [km/s]; errors `rverr_via`, `rverr_via10hr`, `rverr_desi` [km/s], `pm_err_gaia` [**µas/yr**] | N |
| `cut_for_catalog_photometry` | `unbound & trim_new & nonrem`, the pool eligible for catalog photometry | N |
| `catalog_photometry` | same keys as `phot` plus `iso_dist` and `catalog_index`; NaN outside the matched stars | N |
| `matched_to_catalog_photometry` | True where a catalog row was assigned | ⚠️ **`cut.sum()`**, not N |
| `noise_catalog_photometry` | same keys and units as `noise`, recomputed from `catalog_photometry`; NaN where unmatched | N |
| `cut_for_catalog_photometry_primaries`, `catalog_photometry_primaries`, `matched_to_catalog_photometry_primaries`, `noise_catalog_photometry_primaries` | the same four, with `trim_new_primaries` in place of `trim_new` in the cut | as above (`matched_*` is `cut_primaries.sum()`) |

all of these except the `matched_to_catalog_photometry*` keys are in the *same* (untrimmed, pre-`[inMW][trim]`) CoM/luminous ordering, so they combine with plain `&`. `matched_to_catalog_photometry` is indexed within `cut_for_catalog_photometry`. `gmm.py` expands it with `matched_full = zeros(N); matched_full[cut] = matched`, and `np.isfinite(catalog_photometry['mG'])` gives the same full-length mask. `catalog_photometry['catalog_index']` comes out **float64** (the placeholder dict is `np.full(..., nan)`): -1 for unmatched stars inside `cut`, NaN outside `cut`. cast it before using it as an index.

### catalog photometry (the live noisy path)
per sim, after the isochrone photometry, and **twice**: once with `trim_new` (keys as named) and once with `trim_new_primaries` (keys suffixed `_primaries`). it's a loop over `[(trim_new, ''), (trim_new_primaries, '_primaries')]`. the two pools differ slightly, so the catalog rows go to slightly different N-body stars, but the sample size is `len(catalog)` either way.

1. `cut = unbound & trim & nonrem` -- the pool, so the catalog isn't spent on progenitor stars, remnants or trimmed outliers.
2. `noise.assign_photometry_from_catalog(iso_phot[cut], distances[cut], catalog=noise.load_gaia_catalog(orbit), track=noise.isochrone_cmd_track(isocmd), N=len(catalog))`. subsetting to `cut` beforehand does the job of the `eligible` kwarg. the brightest `N` pool stars in **isochrone** apparent G get the catalog's apparent G and BP-RP exactly, rank for rank; absolute mags come from the N-body CoM distances, `log_Teff` from the nearest isochrone point, and `z`/`mz` from RTN-099.
3. the results are scattered back to full length (NaN elsewhere), and DESI / Via 1 hr / Via 10 hr RV errors and Gaia DR3 position/PM errors are recomputed from the new `mG`/`mz`/`log_Teff` with the same `survey_noise` call as the isochrone noise.

**every catalog row gets used, so the noisy sample *is* the catalog.** the pool is always much bigger than the catalog (4500-11700 stars), and on 30 Sep every matched star passed every downstream noise cut in `gmm.py` (Gaia pm_err < 0.5 mas/yr, Via rverr < 5, DESI rverr < 10) -- the brightest catalog G is ~12 and the faintest ~20. so the fitted sample size is set by `len(catalog)` alone and doesn't change with rvir or with via vs via10hr (those two now fit the **same stars** with different RV noise):

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
- the catalog is used in full: no `p_mem` cut, and `bp_rp` isn't dereddened (see the CMD-space lookup subsection -- the catalog MS sits ~0.1 mag redder than the isochrone, which pushes the assigned `log_Teff` cooler).
- the rank matching is global in apparent G, with no dependence on phi1, so the catalog's magnitude distribution along the stream isn't reproduced -- just the overall luminosity function.
- `iso_dist > 3` (catalog star off the isochrone, so its Teff is extrapolated): ~267/1578 for gd1, ~60-80/334 for aau, ~105/129 for pa5, ~83/144 for jet, ~15/46 for c19. that's most of the Pal 5 catalog, so the Pal 5 Via errors rest on extrapolated Teffs.

## `gmm.py`
functions + a main program to fit a two-component (thin stream + cocoon) gaussian mixture to the pickled data, for one **case** at a time. the functions are documented in the GMM section below.

a case is set by three flags at the top of `__main__`, and they build `case_name = <noise>_<binaries>[_constrained]`:

| flag | options | meaning |
|---|---|---|
| `noise` | `None`, `'via'`, `'via10hr'`, `'desi'` | `None` -> `noiseless`. otherwise adds the saved noise draw and applies that survey's selection (below). the inline comment only lists None/via/desi, but `via10hr` works too since it's just a key suffix |
| `include_binaries` | `True` / `False` | `binaries` uses `sc_straighter_primaries` (binary orbital motion in the RVs); `CoM` uses `sc_straighter`. picks the suffix `sfx = '_primaries' if include_binaries else ''`, which selects the coords, `trim_new{sfx}`, and all four catalog products, for noiseless and noisy cases alike |
| `constrain_widths` | `True` / `False` | adds `_constrained`: stage-2 SLSQP with the sigma-ratio constraint. `False` = Powell only, unconstrained |

- **fit dimensions:** `keys = ['phi2', 'v_phi1', 'v_phi2', 'v_gsr']` -- phi2 in **degrees** (not `d_phi2`), and the proper motions converted to **transverse velocities** (km/s, using the CoM distance) before straightening. with noise, the PM noise is converted to v_phi1/v_phi2 with the same (assumed-perfect) distances.
- **noisy cases read the catalog-based products:** `phot = data_dict['catalog_photometry'+sfx]`, `noise_dict = data_dict['noise_catalog_photometry'+sfx]` (likewise `cut_for_*` and `matched_to_*`). the isochrone-based `phot`/`noise` lines are commented out just above them. switching back needs the old `use` lines too (also commented out), since the isochrone path depends on `alive`/`acceptable_G` and has no `matched` mask.
- **selection:**
  - noiseless: `trim_new{sfx} & unbound` (remnants included, as before).
  - noisy: `use = cut & matched_full & good_pm & good_RV & nonrem`, with `good_pm = pm_err_gaia < 0.5 mas/yr` and
    - `via` **and** `via10hr`: `good_RV = rverr < 5 km/s` (via10hr used to be < 10).
    - `desi`: `good_RV = (brightest 679 stars in catalog mG) & (rverr < 10 km/s)`. 679 = the length of the Jarvis+26 table. note the top 679 are picked *before* the pm/RV cuts, so this can come out under 679 (it doesn't for gd1: all 679 pass).
    - `matched_full` already implies `cut`, which already contains `nonrem`, so it reduces to `matched_full & good_pm & good_RV` (hence the `# do I ...need nonrem????` -- no).
    - `alive` and `acceptable_G` are loaded but no longer used.
    - in practice none of the cuts bite: N fit = catalog size (see the table in the catalog photometry subsection).
  - `desi` fits **GD-1 only**.
- **orbits:** same five as `final_datasets.py`; no M3.
- **output:** `data/gmm_tables/<case_name>.fits` with columns `orbit`, `Rvir0`, `M_ts`, `S_ts` (thin-stream means/sigmas, 4-vectors in `keys` order), `M_c`, `S_c` (cocoon), `f_cocoon`, `BIC_one`, `BIC_two` (tables written before 30 Sep don't have the BIC columns). sanity plots (points coloured by p_cocoon + histograms) to `/n/netscratch/conroy_lab/Lab/amphillips/p27_sanity_plots/<case_name>/<orbit>_<rvir>.pdf` when `make_plots=True` -- those directories have to exist already. the `# %%` cells after the table write are quick-look plots.
- the table doesn't record `keys`, so nothing tells you from the file itself whether column 0 is phi2 or d_phi2 -- the "translate back to angle from distance" TODOs in `results.py` suggest some earlier tables were fit in `d_phi2`. regenerate rather than trust an old table if it matters.

**model choice:** each sim also gets a single-component gaussian fit, and both BICs go in the table. the expectation is that for noisy/down-sampled data some sims will prefer one component, which gets called a **non-detection of a cocoon**. see the model choice section below.

## `results.py`
notebook-style cells making the (nice) plots from the tables and pickles: `f_cocoon` / cocoon + thin sigma_phi2, sigma_vgsr vs rvir0 per orbit (ordered by pericenter), CoM vs binaries overlay, binary fractions in the cocoon vs thin stream (`plots/binary_fractions.pdf`), and the demo cocoon-separation panels. imports `gmm` for the membership functions. rows are filtered with `okay_mean = |M_c| < 0.5 S_c` in every dimension, i.e. a fit whose "cocoon" is displaced off-track is dropped from the plot rather than shown. TODO: rewrite.

⚠️ `results.py` still reads the **unsuffixed** catalog keys (`catalog_photometry`, `cut_for_catalog_photometry`, ... around lines 163 and 501, and in `demo_data`), so its noisy `binaries` panels use the CoM-trim photometry, unlike `gmm.py`. give it the same `sfx` switch.

## `claude_binary_selection_experiment.py`
a one-off diagnostic, not part of the pipeline; only reads a pickle. for one sim it prints (1) the binary fraction and binary orbital-motion spread in v_gsr for the full unbound sample, the via-matched and desi-top-679 samples, and in bins down the brightness ranking of the eligible pool; and (2) noiseless GMM fits on each of those selections, for both binary treatments. written to explain the binaries selection effect -- see "binaries in the observed samples" in the GMM section. run across orbits/rvir to check the effect holds beyond GD-1 rvir0 = 0.75.


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

the one `hm` cell checked so far: `gd1` rvir_index=0 has all five copies finished (3 Oct 2026).

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

top level: `init_displacement`, `nsingles`, `nbinaries`, `IDs`, `pot`, plus one subdict per binary treatment.

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

the interpolant sorts `orbit_coords['phi1']` first (`np.interp` requires increasing `xp` and doesn't check), and the multi-`yval` closures are bound with `y=y` (they used to all interpolate `pm_phi2`).

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
- **post-AGB (the proto-WD tail) is dropped by the `max_eep=POST_AGB_EEP` default** (`EEP < 1409`), since those stars are already cut from the sim side by `nonrem` (`type < 10`). `max_phase=5` is also still the default but only does anything on the downloaded MIST file, where post-AGB is phase 6 -- on the artpop isochrone it's phase 5 (table above). on the downloaded file the two cuts remove the identical 302 rows. with the EEP cut the artpop table is 1158 nodes, upper m0 edge **0.793574** (was 0.79435), max logTeff 3.885 (was 5.05). `max_eep=None, max_phase=None` keeps everything, and moves the upper m0 edge (so the `alive` count depends on it).
  - before this cut went in (fixed 29 Sep 2026), N-body stars with m0_zams in 0.793574-0.79435 got post-AGB/WD photometry (absolute G ~ +10) and counted as `alive`. that was **2 stars per pickle** (the same two m0 values in every sim, 40 of 193153 `alive` stars across the 20 pickles, 0-2 per sim after `nonrem & unbound`). photometry of every other star is bit-identical with and without the cut. the 30 Sep pickles were made after the fix, so they don't have these stars.

caveats that are real selections, not rounding details:
- the isochrone is old, so it stops at 0.793574 Msun (artpop, after the post-AGB cut; 0.804398 for the downloaded file) and every more massive sim star is thrown out (`alive`). at `hm` the discarded stars are the luminous ones, so quote how many got dropped alongside any luminosity-weighted number.
- for binaries only the luminous component's photometry is used; the companion's light isn't added.

### CMD-space lookup: catalog (G, BP-RP) -> nearest isochrone point (`noise.py`)
the reverse of the mass interpolation, for the catalog-photometry TODO: real catalog stars never sit exactly on the isochrone, so to get a `log_Teff` for viamock each one is snapped to the nearest point on the track. step 5 of `assign_photometry_from_catalog`.

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

### the blackbody fallback (in `paf` and duplicated in `noise.py`, no longer the live path)
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

## binaries in the observed samples: a selection effect (5 Oct) -- for the paper
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

so with zero noise, selection alone takes the `binaries` cocoon to 9.2 / 11.1 km/s, against 10.3 / 12.2 with noise. noise adds only the last ~1 km/s. the `CoM` cocoon hardly moves with the selection. all of this comes from `claude_binary_selection_experiment.py` (`--orbit`, `--rvir`; defaults gd1 / 0.75; ~1 min). it builds the matched mask the same way `gmm.py` does, but leaves out the pm/RV error cuts, which don't remove anything. then it fits `sc_straighter{sfx}` without noise.

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
- **GMM label switching / degenerate fits in the catalog-noise tables** (30 Sep). `sort_components` orders by v_gsr width only, so when the two components have nearly equal v_gsr sigma the label is a coin flip: `pa5` rvir 6 in `via_noise_CoM` / `via10hr_noise_CoM` has the "thin" component at sigma_v_phi1 ~ 34 km/s and the "cocoon" at ~5.7, with v_gsr sigmas 1.08 vs 1.14, so the reported `f_cocoon` = 0.61 is really the thin-stream fraction. and some small-N fits have the thin component collapse onto a handful of stars: `via10hr_noise_CoM` `jet` rvir 1.5 (`f_cocoon` = 1.000) and `c19` rvir 0.75 (`f_cocoon` = 0.935, thin sigma_phi2 = 0.01 deg). **all three pass `results.py`'s `okay_mean` filter**, so they'll show up in the summary plots. the BIC TODO is the principled fix; sorting on a combined width (e.g. the product of sigma/sd over dims), or requiring agreement across dims, would catch the pa5 case.
- **in the `binaries` cases the "cocoon" is mostly the binary-orbital-motion population**: `via_noise_binaries` / `desi_noise_binaries` cocoons have sigma_v_gsr ~ 6-20 km/s against ~1 km/s for `CoM`, with thin-stream-like sigma_phi2. expected physics, but it means `f_cocoon` in those tables isn't a spatial/kinematic cocoon fraction, so don't compare it directly to the CoM tables. and the noisy `binaries` cocoons are hotter than the noiseless one mostly through selection (bright stars are more often binaries), not noise -- see "binaries in the observed samples" in the GMM section.
- in `gmm.py`'s plotting block, `cut = 3*sigmas_fit[-1][jj]` reuses the name of the catalog-photometry `cut` mask. harmless now because `cut` is reloaded from the pickle at the top of each rvir iteration, but it'll bite if that load ever moves.
- **pickle filenames don't include stellar_pop** -- an `lm` run of `final_datasets.py` would silently overwrite the `hm` pickles. (`gmm.py`'s `mass_index = 1` line says "LOW mass" but 1 is `hm`, and it isn't used there anyway.)
- **pickles live on netscratch**, which is purged on a 90-day clock (the current set was written 30 Sep 2026, so it ages out around the end of December). copy them somewhere permanent before they age out.
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
- scratch: `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/` -- 90-day purge, actively losing snapshots. don't use.
- the circular orbit simulation is from the [Phillips+26](https://iopscience.iop.org/article/10.3847/1538-4357/ae680b) grid, stored at `/n/holystore01/LABS/conroy_lab/Lab/amphillips/finished_grid/` in the directories that begin with 0-7 (see `extended_grid_info`).

## derived data products:
- pickles (`final_datasets.py`): `/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/<orbit>_<rvir:.2f>.pickle` -- scratch, see the purge warning.
- GMM sanity plots (`gmm.py`): `/n/netscratch/conroy_lab/Lab/amphillips/p27_sanity_plots/<case_name>/` -- one subdir per case, must exist before running.
- movie frames (`animations/`): `/n/netscratch/conroy_lab/Lab/amphillips/movies/<movie name>/`. `particle_spray_vs_nbody/` also holds the cached spray (`spray_gd1_hm_0.75_0.h5`); scratch, so regenerate with `--make-spray` if it's been purged.
- GMM tables: `data/gmm_tables/` in the repo.
- Jarvis+26 table 7 (DESI GD-1 members, for the mag-distribution comparison and the photometry-matching TODO): `data/jarvis26_Table7.fits` in the repo (identical copy at `/n/home02/amphillips/data/jarvis26_Table7.fits`).
- Bonaca & Price-Whelan (2025) Gaia member catalogs: `data/bpw25_catalogs/<orbit>.fits`.

## conda environment:
`petar_env`, stored at `~/.conda/envs/petar_env`, containing standard packages like numpy, scipy, matplotlib, astropy, but in particular gala (1.9.1), plus `petar`, `sklearn`, `pygaia`, `artpop` (MIST isochrone fetching). `viamock` (Via RV errors) is **not installed** -- `noise.py` and `final_datasets.py` put `~/software/viamock` on `sys.path`.

nothing outside this env can import gala or petar, so run scripts with `~/.conda/envs/petar_env/bin/python` (or activate the env) rather than the system python. importing `inspect_new_sims` (which `noise`, `gmm`, `final_datasets`, `results` all do) builds `extended_grid_info()` and reads `FINAL_ics_nolmc.csv` at import time.

## petar documentation
the `README.md` from Long Wang's [PeTar github](https://github.com/lwang-astro/PeTar) is useful.
