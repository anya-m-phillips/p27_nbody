Exploring an N-body grid of mock streams to probe effects of progenitor dynamics on stream dispersion. Questions to answer might include:
- cocoon fraction and dispersion/extent as a function of progenitor density, orbit: Jarvis+26, Carlberg+26 discuss cocoon structures around MW streams, created by 
- structure of the stream velocity dispersion profile: created by the (evolving) cluster dispersion+tidal filling factor and the progenitor orbit/overlap of successive energy ``feathers." see Bovy(2014).


# map of code, (not so) briefely:
`/data`:
- `FINAL_ics_nolmc.csv`: credit Vedant Chandra, present-day positions, velocities of some promising streams. Age and progenitor mass estimates also given. **the `name` column for the six run orbits has been renamed to the short keys i index by** (`gd1`, `pa5`, `aau`, `m3`, `c19`, `jet`) so that `prog_tab[prog_tab['name']==orbit]` just works. everything else in the table keeps its original name.
- `init_displacements.txt`: generated in `get_init_displacements.py`, Vedant's progenitor locations back-integrated by their stream ages, plus like 100 Myr to account for initial expansion due to massive star evolution, rounded to a multiple of 10 so that I can safely output sim snapshots every 10 Myr and get the present day in the final snapshot. **positions in this file are in pc** (that's what `petar.init` wants), velocities in km/s -- see the units bullet below, this caused a real bug.
- `ATLAS-Aliqa Uma.fits`, `C-19.fits`, `Jet.fits`, `M3.fits`: credit Bonaca & Price-Whelan (2025); TODO: pull GD-1 and Pal 5 as well. 

`/scripts`:
- `PETAR_ANALYSIS_FUNCTIONS.py`: a bunch of functions, mostly migrated from other projects; imported as `paf`. see the notes section below for the parts that are actually load-bearing.
- `streamframe.py`: credit Jake Nibauer, transformation to stream frame coordinates as seen from the Galactic center
- `vedant.mplstyle`: credit Vedant Chandra, plotting style stuff

`/old`: stuff that is migrated from older repositories (esp. `~/stream_velocity_structures`)
- `DESI_comparison.py`: noise old simulation data like Gaia+DESI and compare to Jarvis+2026 cocoon detection. Will eventually systematize and do for all new sims
- `prog_properties_summary.py`: cocoon fractions/dispersions for GD-1 portion of old sim grid

scripts in top directory for now:
- `get_init_displacements.py`: used to generate init_displacements.txt, fed as inputs to `petar.init ...` for displacing progenitors from the galactic center at the initial condition. **stale:** its `names_to_run` list still holds the long names (`'ATLAS-Aliqa Uma'`, `'GD-1'`, ...) which no longer match the renamed `name` column, so re-running it right now writes an *empty* file with no error. update that list before touching it again.
- `inspect_new_sims.py`: writing a bunch of functions to process sim data alongside `paf`; in particular transforming to ``observed" stream frame coordinates (based on ICRS coordinates; correcting for solar reflex motion, etc). see the stream frame section below. contents:
  - `prepare_nbody_data(path, include_photometry, i, apo, init_displacement)`: wraps `paf.load_core` + `paf.intrinsic_stream_data_v3`; assumes snapshots every 10 Myr. returns `core, data_dict, CMdict, lumdict, inMW, trim` (plus `G, BP, RP, z` if `include_photometry=True`, which is slow -- it loops `paf.get_gaia_photometry` per star). `data_dict` has `CoM`/`luminous`/`companions` subdicts for the three binary treatments; `inMW`/`trim` masks live on `CMdict`.
  - `prepare_nbody_data_anycopy(orbit, stellar_pop, rvir_index, copies, **kwargs)`: same thing but loops `copy` and returns the first realization that actually finished. an unfinished copy still has a (short) `data.core`, so it gets all the way into `intrinsic_stream_data_v3` before raising `FileNotFoundError` on the missing `data.<file_index>` -- that exception is the signal. dedupes on path because `retrieve_sim_info` ignores `copy` for `circ`. returns `(prepare_nbody_data output), path, apo, age, init_displacement, copy`. use this instead of hand-rolling fallbacks: at `rvir_index=1` alone, jet is missing copies 0,1,2 and c19 is missing 0,1.
  - `rotation_matrix(a,b,c)`: rotate by a,b,c about x,y,z. composed `Rz @ Ry @ Rx`, i.e. x rotation applied first -- chosen for movie purposes.
  - `streamframe_coords_observed(orbit, data_dict, prog_tab)`: the main one. see below. **returns a `SkyCoord`, not a dict, as of the current version.**
  - `prog_orbit_track` / `observed_orbit_track` / `chop_orbit_track` / `straightened_obscoords_orbit_interp`: the observed-frame orbit straightening. see its own section below.
  - `poly_straightening(coords, tc)`: loops the keys of a coord dict and subtracts a degree-5 polynomial fit via `paf.straighten_stream_polynomial`. applied *after* the orbit-interp straightening to mop up whatever residual tilt is left. it copies rather than subtracting in place, so it's safe to call on a dict you still want. note `coords_straighter['phi1']` is the same array object as `coords['phi1']`, not a copy.
  - `clip_coords(coords, [inMW, trim])`: apply the masks to every key of a coord dict at once.
  - `outlier_clip(vr, pmphi1, pmphi2)`: hard cuts at 100 km/s and 1.5 mas/yr, to stop one or two wack stars from dominating the cocoon dispersion.
  - `get_cocoon_selection(coords, cuts)`: OR of `|x| > cut` over phi2, pm_phi1, pm_phi2, v_gsr -- order matters, it zips against the `*_cuts` lists at the bottom of the script.
  - `desi_RVerr(zmag, feh)`, `add_noise(...)`: survey error models. `add_noise` is still a stub.
  - everything under `if __name__=='__main__':` loops all orbits and makes the cocoon-separation panels, so the function defs can be imported elsewhere without running it. comment that line out and un-indent if working interactively in the notebook cells. **right now that guard is commented out**, so importing this module runs the whole loop.

- `cocoons.py` summarizing cocoon fractions and properties as a function of progenitor properties. functions are called from `inspect_new_sims.py` as `simspect.[...]`
- `velocity_movie.py`, `grid_movie.py`, `grid_movie_long.py` : scripts to run with a slurm wrapper for animations. TODO: make all of these parallel so that the wrapper is submitted as an array job where each sub-job generates one frame. indescribably faster than doing this in a loop. 
- `nfc_plots.py` was used for plotting in preparation for a conference; will likely abandon soon

# notes on `paf` (PETAR_ANALYSIS_FUNCTIONS.py)
~2400 lines, a lot of it leftover from the old grid (the `#<-- LEFTOVER FROM OLD GRID` markers are honest -- `define_paths*`, `define_apocenters`, `unpack_escaper_dict` etc. are superseded by `extended_grid_info`). the notes below are for the parts the current pipeline actually touches.

## conventions to keep straight
- **`i` vs `file_index`.** the new grid writes snapshots every 10 Myr, so `file_index = int(i/10)` where `i` is simulation time in Myr. functions are not consistent about which one they want, and nothing validates it:
  - sim time in Myr: `intrinsic_stream_data_v3(path, i, ...)` (computes `file_index` itself), `prog_position(init_displacement, i)`, `xform_to_core_frame` (it converts internally via `file_naming_convention`).
  - file index: `CM_to_galcen_frame`, `core_to_galcen_frame`, `is_dissolved` (indexes `tidal.n[i]`), and the `file_index=` kwargs.
  - both: `load_coords_v2(path, i, ..., file_index=...)` takes time as `i` *and* index as `file_index`; `straighten_stream_orbit_interp(coords, yval, core, i, ...)` wants `i` as a time when `use_core=False` and a file index when `use_core=True`.
  - `file_naming_convention="every integer"` is the *default* on `load_particle`/`xform_to_core_frame`/`clip_outside_rtid`, which is the OLD grid's convention. for the new grid pass `"every 10"` (or pass the file index directly).
- **`data.core` is written at the same 10 Myr cadence as the snapshots**, so `core.pos[file_index]` lines up. (checked against `m3/lm/0.75/0`: 501 snapshots, core time column steps by 10.)
- **units.** petar outputs pc, pc/Myr, Msun. `StreamFrame` wants kpc, kpc/Myr and returns deg, mas/yr, kpc, km/s. most `paf` functions hand back astropy quantities, but the streamframe coord dicts are bare floats.
- **`init_displacement` is kpc, km/s** everywhere it is *consumed* (`prog_position`, `integrate_prog_orbit`, `straighten_stream_orbit_interp` all do `init_displacement[:3] * u.kpc`), but `data/init_displacements.txt` *writes* it in pc because that's what petar wants. `extended_grid_info.__init__` now does the pc->kpc conversion once, at the bottom, for the six real orbits (`circ` was always already in kpc). **this was a silent 1000x bug** -- it put the progenitor reference at 16.6 Mpc instead of 16.6 kpc, where the potential is negligible, so the "orbit" free-streamed in a straight line. sanity check if you ever touch it: `paf.prog_position(init_displacement, age)` must reproduce that orbit's row in `FINAL_ics_nolmc.csv` (it does, to 0.0000 kpc, for all six).
- **potential is `gp.BovyMWPotential2014(units=galactic)` everywhere** -- `prog_position`, `integrate_prog_orbit`, `straighten_stream_orbit_interp`. petar was run with `external_mode='galpy'` to match. don't change one without the others.
- **`interrupt_mode='bse'`** is assumed by every loader (stellar evolution on).

## `extended_grid_info(scratch=True)`
holds paths, apocenters, stream ages (Vedant's), and init displacements as attributes. grid axes are orbit x stellar_pop (`lm`/`hm`) x rvir (`0.75, 1.5, 3, 6` pc, indexed 0-3) x copy. use `retrieve_sim_info(orbit, stellar_pop, rvir_index, copy) -> path, apo, age, init_displacement` rather than assembling paths by hand.

the six `*_init_displacement` literals are still pasted verbatim from `init_displacements.txt` (i.e. in pc), and a loop at the end of `__init__` converts the position components to kpc. if you paste a new orbit in, paste the raw pc numbers and add its key to that loop -- don't pre-convert.

**not every copy finished.** every `(orbit, stellar_pop, rvir_index)` has at least one that did, but which one varies; use `prepare_nbody_data_anycopy` rather than assuming `copy=0`. availability at `stellar_pop='lm'`, present-day snapshot present (`Y`) or not (`.`):

| orbit | rvir_index=0 | 1 | 2 | 3 |
|---|---|---|---|---|
| gd1 | Y . . . Y | Y Y . Y . | Y . Y Y Y | Y Y Y Y Y |
| pa5 | Y . Y . . | Y Y Y Y . | . Y Y Y . | Y Y Y Y Y |
| aau | Y Y Y Y Y | . Y Y Y Y | Y Y Y Y Y | Y Y Y Y Y |
| m3  | Y Y Y Y Y | Y Y Y . Y | Y Y Y Y Y | Y Y Y Y Y |
| c19 | Y . . . Y | . . Y Y Y | Y Y Y Y Y | Y Y Y Y Y |
| jet | Y Y . . Y | . . . Y Y | Y . Y Y Y | Y Y Y Y Y |

sharp edges:
- `scratch=False` prints a warning and `return`s from `__init__` *before assigning any attributes*, so you get an object that AttributeErrors on first use rather than a useful failure. also the storage path hardcoded there is `conroy_lab/Lab/...`, which does not match the `itc_lab/Users/...` path in the logistics section below -- worth reconciling when doing the storage TODO.
- `retrieve_sim_info` is a chain of bare `if`s with no `else`, so a typo'd orbit string gives `UnboundLocalError` on `base_paths`.
- `circ` is the odd one out: it comes from the OLD grid under `finished_grid/`, has one realization per rvir (no `copy` subdir), `circ_lm_paths[0]` is flagged unfinished, and `retrieve_sim_info` returns `age = 10000*10 = 100000` -- which is why `inspect_new_sims.py` hardcodes `i=30000` for circ instead of using the returned age.

## loading petar data
`load_core(path)`, `load_tidal(path)`, `load_particle(path, i)` are thin `petar.*` wrappers. `is_dissolved(path, i, threshold=100)` = fewer than 100 stars inside the tidal radius at index `i`.

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

this is the "intrinsic"/God's-eye frame. contrast with `streamframe_coords_observed` below, which is heliocentric, reflex-corrected, and uses real ICRS great circles -- the key names overlap but mean different things (`r` vs `distance`, `vr` vs `v_gsr`).

## `intrinsic_stream_data_v3(path, i, core, apo, init_displacement, use_core=False, binary_treatments=[...])`
the main entry point. returns one big nested dict.

top level: `init_displacement`, `nsingles`, `nbinaries`, `IDs`, `pot`, plus one subdict per binary treatment.

**ordering gotcha:** top-level `IDs` and `pot` are concatenated as `[singles, binary_p1, binary_p2]`, so they always have length `nsingles + 2*nbinaries`. the subdict arrays are concatenated as `[singles, binaries]`, where "binaries" is one row per binary for `CoM`/`luminous` but two for `companions`. so **`IDs` only lines up with the `companions` subdict**; for `CoM`/`luminous` you have to slice `[:nsingles+nbinaries]` (this is what the docstring note is getting at).

the three treatments:
- `CoM` -- each binary as a single center-of-mass point. no `L`/`R`/`type`.
- `luminous` -- each binary represented by its brighter component (`binaries.p1.star.lum >= binaries.p2.star.lum`).
- `companions` -- both components kept separately; this is the only one that double-counts.

each subdict has: `coords` (the `StreamFrame` dict), `pos`/`vel` (galactocentric, with units), `mass`, `inMW`, `trim`, `in_rtid`, and `phi2_straight`/`r_straight`/`vr_straight`/`pm_phi1_straight`/`pm_phi2_straight`. `luminous`/`companions` also get `L`, `R`, `type`.

### `in_rtid` (for cutting the progenitor out)
`in_rtid` = distance from the core <= `tidal.rtid[file_index]`, i.e. still bound-ish, so `~in_rtid` is the "remove the progenitor" mask. it is computed on the **core-frame** `singles.pos`/`binaries.pos` (which is the frame those files are already in, and the frame `rtid` is measured in), and it is applied **before** `inMW`/`trim` -- so the usage is `in_rtid[inMW][trim]`, same as everything else.

it used to sum over `range(2)`, i.e. a projected cylindrical radius with z dropped, which over-counted bound stars by ~4% (6935 vs 6672 singles for `gd1/lm/1.5/0`). that's fixed -- it's `range(3)` now.

still wrong: **it is the wrong length for the `companions` subdict.** `b_r` comes from `binaries.pos`, one row per binary, so `in_rtid` is always `nsingles + nbinaries` long. that matches `CoM` and `luminous`, but `companions` is `nsingles + 2*nbinaries`, so `in_rtid[inMW]` raises `IndexError: boolean index did not match`. it fails loudly rather than silently, but it does fail -- either build `[in_rtid_s, in_rtid_b, in_rtid_b]` inside the companions branch or just don't ask for `in_rtid` there.

also minor: `load_tidal(path)` is called inside the `binary_treatment` loop, so the tidal file gets read once per treatment.

## masks: `inMW` and `trim` (and why it's always `[inMW][trim]`)
from `trim_coords_percentile(coords, low=1, high=99, apo=apo)`:
1. `inMW = coords['r'] <= 1.5*apo` -- drops stars flung well outside the orbit, i.e. **this is the only thing `apo` is used for**, so passing the wrong orbit's apocenter silently gives you the wrong outlier mask.
2. `trim` = 1st-99th percentile clip in phi1 and phi2, **computed on the already-`inMW`-selected array**.

so `trim` has length `inMW.sum()`, not `len(inMW)`, which is why every plotting line is `x[inMW][trim]` and why the two are not interchangeable or combinable with `&`.

## `straighten_stream_orbit_interp(coords, yval, core, i, ...)`
integrates the progenitor orbit +/- `Dt` Myr around time `i`, pushes that orbit through `StreamFrame` with the same reference coordinate, interpolates `yval` against `phi1`, and subtracts it -- so the output is a residual about the orbit track. `yval` can be a string or a list of strings (`intrinsic_stream_data_v3` passes `['phi2','r','vr','pm_phi1','pm_phi2']`).

`Dt` is found by an adaptive search, not fixed: it grows `Dt` until the orbit's phi1 range covers the data, but shrinks it if the orbit wraps around in phi1 (`dphi1 > 100` deg between samples), halving the step whenever it flips between those two cases. `Dt_start=240` from v3.

**it returns the iteration count as the last value -- check it.** the loop caps at `max_iters=100` and if it exhausts them it just falls through and straightens anyway, with no warning. an `itr` of 100 means the straightening never converged and the residuals are suspect. the returned residuals are for *all* particles (untrimmed), so apply `inMW`/`trim` afterwards. (the `itr=100` you used to get for every orbit was the init_displacement units bug, not this function.)

two fixes went in here:
- the interpolant now **sorts `orbit_coords['phi1']`** first. `np.interp` requires increasing `xp` and does not check -- the orbit chunk is ordered in *time*, so phi1 comes out increasing or decreasing depending on which way the progenitor runs, and the decreasing case silently returned garbage.
- the multi-`yval` branch had a **late-binding closure**: all five returned `interp_orbit` functions ended up interpolating `pm_phi2`. bound with `y=y` now. the *residuals* were always fine (they're evaluated inside the loop), only the returned callables were wrong.

`straighten_stream_polynomial(phi1, y, degree=5, trim_criteria=[inMW, trim], return_poly_fn=True)` is the polynomial alternative. note `trim_criteria` has no working default -- it unpacks `inMW, trim = trim_criteria` unconditionally, so leaving it `None` is a `TypeError`.

## photometry stuff [MAJOR WIP ⚠️]
synthetic photometry from a **blackbody**, with **top-hat filters** -- no real Gaia/SDSS response curves, so treat colors as approximate.
- `define_photometric_bands()` -- Gaia G (330-1050 nm), BP (330-680), RP (630-1050) from the DR2 paper, plus a 100 nm-wide z centered at 900 nm.
- `integrated_mag(nu_min, nu_max, T, R, d=10)` -- AB magnitude (the -48.60 zero point), `d` in pc, so `d=10` gives an absolute mag. **T and R must be cgs bare numbers** (K and cm) since `B(nu,T)` uses cgs constants -- this is why `inspect_new_sims.py` does `.cgs.value` first. note the 1/nu weighting normalization at the `np.trapezoid` lines is marked "suggestion that idk why works. " in the source and has not been validated against a real photometric zero point; worth checking against a known star before trusting absolute mags (colors are probably safer).
- `get_gaia_photometry(Teff, Radius, distance) -> [G, BP, RP]`, `get_z_photometry(...) -> z`. both are scalar-only, hence the slow per-star loop in `prepare_nbody_data(include_photometry=True)`. vectorizing these is the obvious win if photometry ever becomes a bottleneck.
- `m_from_M(M, dist)` for the distance modulus (dist needs astropy units here, unlike the above -- inconsistent, but that's how it is).

# stream frame coordinates ("observed" frames)
`streamframe_coords_observed(orbit, data_dict, prog_tab)` in `inspect_new_sims.py`. `data_dict` is one of the binary-treatment subdicts (`CoM`, `luminous`, `companions`) so it has `pos`/`vel` with astropy units. pipeline is:

`galcen pos/vel` -> `paf.galcen_to_ICRS` -> `gala.coordinates.reflex_correct` -> `.transform_to(<great circle frame>)`

**returns `sc, selected_streamframe` where `sc` is a `SkyCoord`, not a dict** (it used to be a dict; the dict-building block is commented out in the source). so pull attributes off it: `sc.phi1`, `sc.phi2`, `sc.pm_phi1_cosphi2`, `sc.pm_phi2`, `sc.distance`, `sc.radial_velocity` -- all with units, unlike the intrinsic coord dicts. returning the frame object as well is what lets the progenitor orbit be pushed through the *same* frame downstream. these are for *all* particles -- apply `inMW`/`trim` yourself.

when it does get flattened to a dict (in `straightened_obscoords_orbit_interp`) the keys are `phi1`, `phi2` [deg], `pm_phi1` (= `pm_phi1_cosphi2`), `pm_phi2` [mas/yr], `v_gsr` [km/s], `distance` [kpc], bare floats. contrast with the intrinsic dict: `distance` is heliocentric where `r` is galactocentric, and `v_gsr` is a real reflex-corrected line-of-sight velocity where `vr` is relative to the progenitor.

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
| `circ` | -- | no frame; not a real stream, so it's skipped in the `__main__` loop |

phi1 zero points are set by the present-day progenitor sky position read out of `data/FINAL_ics_nolmc.csv` (except `jet`, where Do+26 give an origin directly), looked up as `prog_tab[prog_tab['name']==orbit]` -- which works because the table's `name` column was renamed to the short keys. there is still no `else` branch -- an unrecognized `orbit` string leaves `selected_streamframe` undefined and you get a `NameError` at the `transform_to` line rather than a useful message.

## straightening in the observed frame
four functions in `inspect_new_sims.py`, chained by `straightened_obscoords_orbit_interp(orbit, CMdict, prog_tab, Dt=500)`:

1. `prog_orbit_track(w0, Dt)` -- integrate the progenitor +/- `Dt` Myr in mwp2014 and concatenate `[backward reversed, forward]`. **note the `[:-1]`**: both integrations contain t=0, so the naive concatenation duplicates the progenitor and puts an exact `dphi1 = 0` step at the midpoint, which stalls any sign-based walk. with it dropped the progenitor sits at index exactly `Dt`, which is also `n//2`.
2. `observed_orbit_track(pos, vel, obs_streamframe)` -- same galcen -> ICRS -> reflex_correct -> great circle pipeline as the data. **the reflex correction on the orbit is not optional**: it doesn't change phi1/phi2 (positions are untouched) but it absolutely changes `v_gsr` and the proper motions, and if you skip it the residuals are offset by the solar motion.
3. `chop_orbit_track(phi1, jump_threshold=45)` -- walk outward from the progenitor and cut at the first bad step on each side, so what comes back is contiguous, single-valued and increasing in phi1 (ready for `np.interp`). two kinds of bad step:
   - `|dphi1| > jump_threshold`: the +/-180 seam of the great circle frame. these are ~356 deg while real steps at 1 Myr sampling are < 7 deg/Myr for gd1/pa5/aau/c19/jet, so the threshold is not delicate at all.
   - a change of direction: the orbit genuinely doubling back in phi1. **m3 needs this** -- near its low pericenter it hits 31 deg/Myr, so no jump threshold can catch its turnaround, and without the direction test the segment comes back non-monotonic and `np.interp` quietly returns nonsense.
4. subtract `np.interp(data_phi1, orbit_phi1[idx], orbit_y[idx])` from each of phi2, pm_phi1, pm_phi2, v_gsr, distance.

the nice property is that the chopped chunk **stops depending on `Dt`** once `Dt` is big enough to reach a seam or a turnaround on both sides. measured on `lm`, `rvir_index=0`, `copy=0` with `jump_threshold=45` -- chunk length at Dt = 200 / 400 / 800:

| orbit | data phi1 range | 200 | 400 | 800 | phi2 residual std |
|---|---|---|---|---|---|
| gd1 | -113.7 to 30.9 | 401 | 608 | 652 | 0.334 |
| pa5 | -24.3 to 26.4 | 332 | 332 | 332 | 0.234 |
| aau | -38.2 to 23.7 | 401 | 703 | 850 | 0.214 |
| m3 | -50.7 to 77.2 | 398 | 479 | 479 | 0.786 |
| c19 | -29.7 to 23.7 | 401 | 508 | 508 | 0.188 |
| jet | -28.2 to 23.8 | 401 | 720 | 720 | 0.120 |

coverage of the trimmed data is 1.000 for all six (measured at Dt=400; the script's default is 500). gd1 and aau are the ones still growing at 800, and they're bounded by the seam anyway. the phi2 residual std column is against the orbit track only, before `poly_straightening`.

`straightened_obscoords_orbit_interp` calls plain `np.interp` with no `left`/`right`, so anything outside the track's phi1 range gets **silently clamped to the edge value** rather than flagged. fine at the moment since coverage is complete, but if you change `Dt` or the frame, pass `left=np.nan, right=np.nan` and check for NaNs instead of trusting it.

## gala `GreatCircleICRSFrame` gotchas
version in `petar_env` is gala 1.9.1. the modern API is pole+origin only; passing `ra0=` or `rotation=` straight to the constructor raises. use the classmethods.

- **`from_pole_ra0` and `from_endpoints` are not actually different animals.** both just compute a pole and an origin and hand them to the same constructor. if you build a frame one way and then rebuild it from its own `.pole` and `.origin.ra` the other way, you get bit-identical phi1/phi2. so if one "works" and the other doesn't, suspect a typo before suspecting gala (ask me how I know).
- **the `priority` kwarg matters a lot.** the pole and origin have to be orthogonal. if your origin is off the great circle, gala emits only a `RuntimeWarning` (easy to miss in a notebook cell) and silently fixes it. default is `priority='origin'`, which *moves the pole* -- i.e. throws away the great circle you defined by endpoints. `priority='pole'` keeps the great circle and projects the origin onto it, which is what you want when the endpoints are the literature definition and the progenitor is only setting the phi1 zero point. passing `ra0=` instead of `origin=` also keeps the pole.
- for `m3` this is not academic: the `M3` row of `FINAL_ics_nolmc.csv` sits ~17 deg off the great circle through the Yang+23 endpoints, so `priority` changes where the endpoints land by >12 deg in phi2. also, `data/M3.fits` spans ra 189-269 deg with ~7.5 deg of phi2 scatter in every endpoint-derived frame, so one great circle may just not describe the whole M3 stream. see the TODO in the code -- may end up settling for "stream on an M3-like orbit," since that run is really a high-e / low-pericenter test.

# bugs / stale things that remain
in rough order of how much they'd hurt:
- **`in_rtid` is the wrong length for the `companions` subdict**, so `in_rtid[inMW]` raises `IndexError` there. fine for `CoM` and `luminous`. details above.
- **`get_init_displacements.py` silently produces nothing** -- its `names_to_run` list still uses the long stream names that were renamed out of `FINAL_ics_nolmc.csv`.
- **anything cached from before the init_displacement units fix is wrong** -- not just the `*_straight` residuals but the intrinsic `coords` themselves, since the progenitor reference position and velocity were both bogus. regenerate.
- `straightened_obscoords_orbit_interp` clamps instead of flagging outside the orbit track's phi1 range.
- `straighten_stream_polynomial`'s `trim_criteria=None` default is a `TypeError`, not a default.
- no `else` branch in either `retrieve_sim_info` (`UnboundLocalError`) or `streamframe_coords_observed` (`NameError`) for an unrecognized orbit string.
- `extended_grid_info(scratch=False)` returns from `__init__` before assigning anything, and its hardcoded storage path disagrees with the one in logistics below.
- `load_coords_v2` builds the all-particles filename from the raw `file_index` kwarg, so `file_index=None` + `load_all=True` opens `data.None`.
- `core_to_galcen_frame` still adds the raw `core.vel` rather than the `fix_core_vel` version (see above).
- `correct_core` is dead code that also strips units.
- the `m3` great circle still doesn't describe the whole stream (see the gala gotchas section); the phi2 residual std of 0.786 vs ~0.2 for everything else is that showing up.
- photometry is still blackbody + top-hat filters with an unvalidated normalization.
- `straighten_stream_orbit_interp_arbitrary_frame` in `paf` is a commented-out stub.
- `add_noise` is a stub, and `viamock` isn't in the env yet.

# logistics:
## simulation data:
- scratch: `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/` (stopped keeping track of when my 90 days is up though; will be removed during monthly maintanence eventually)
- storage: `/n/holystore01/LABS/itc_lab/Users/amphillips/extended_grid/` TODO: add a storage path to `paf.extended_grid_info()` 

note that the circular orbit simulation is from the [Phillips+26](https://iopscience.iop.org/article/10.3847/1538-4357/ae680b) grid, so it is stored at `/n/holystore01/LABS/conroy_lab/Lab/amphillips/finished_grid/` in the directories that begin with 0-7 (see `PETAR_ANALYSIS_FUNCTIONS.py`'s `get_extended_grid_info`)

## conda environment:
`petar_env`, stored at `~/.conda/envs/petar_env`, containing standard packages like numpy, scipy, matplotlib, astropy, but in particular gala (1.9.1), plus `petar`, `sklearn`, `pygaia`. TODO: add `viamock` for the Via velocity errors.

nothing outside this env can import gala or petar, so run scripts with `~/.conda/envs/petar_env/bin/python` (or activate the env) rather than the system python. 

## petar documentation
the `README.md` from Long Wang's [PeTar github](https://github.com/lwang-astro/PeTar) is useful. 