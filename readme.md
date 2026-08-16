Exploring an N-body grid of mock streams to probe effects of progenitor dynamics on stream dispersion. Questions to answer might include:
- cocoon fraction and dispersion/extent as a function of progenitor density, orbit: Jarvis+26, Carlberg+26 discuss cocoon structures around MW streams, created by 
- structure of the stream velocity dispersion profile: created by the (evolving) cluster dispersion+tidal filling factor and the progenitor orbit/overlap of successive energy ``feathers." see Bovy(2014).


# map of code, briefely:
`/data`:
- `FINAL_ics_nolmc.csv`: credit Vedant Chandra, present-day positions, velocities of some promising streams. Age and progenitor mass estimates also given.
- `init_displacements.txt`: generated in `get_init_displacements.py`, Vedant's progenitor locations back-integrated by their stream ages, plus like 100 Myr to account for initial expansion due to massive star evolution, rounded to a multiple of 10 so that I can safely output sim snapshots every 10 Myr and get the present day in the final snapshot.
- `ATLAS-Aliqa Uma.fits`, `C-19.fits`, `Jet.fits`, `M3.fits`: credit Bonaca & Price-Whelan (2025); TODO: pull GD-1 and Pal 5 as well. 

`/scripts`:
- `PETAR_ANALYSIS_FUNCTIONS.py`: a bunch of functions, mostly migrated from other projects; imported as `paf`. see the notes section below for the parts that are actually load-bearing.
- `streamframe.py`: credit Jake Nibauer, transformation to stream frame coordinates as seen from the Galactic center
- `vedant.mplstyle`: credit Vedant Chandra, plotting style stuff

`/old`: stuff that is migrated from older repositories (esp. `~/stream_velocity_structures`)
- `DESI_comparison.py`: noise old simulation data like Gaia+DESI and compare to Jarvis+2026 cocoon detection. Will eventually systematize and do for all new sims
- `prog_properties_summary.py`: cocoon fractions/dispersions for GD-1 portion of old sim grid

scripts in top directory for now:
- `get_init_displacements.py`: used to generate init_displacements.txt, fed as inputs to `petar.init ...` for displacing progenitors from the galactic center at the initial condition. 
- `inspect_new_sims.py`: writing a bunch of functions to process sim data alongside `paf`; in particular transforming to ``observed" stream frame coordinates (based on ICRS coordinates; correcting for solar reflex motion, etc). see the stream frame section below. contents:
  - `prepare_nbody_data(path, include_photometry, i, apo, init_displacement)`: wraps `paf.load_core` + `paf.intrinsic_stream_data_v3`; assumes snapshots every 10 Myr. returns `core, data_dict, CMdict, lumdict, inMW, trim` (plus `G, BP, RP, z` if `include_photometry=True`, which is slow -- it loops `paf.get_gaia_photometry` per star). `data_dict` has `CoM`/`luminous`/`companions` subdicts for the three binary treatments; `inMW`/`trim` masks live on `CMdict`.
  - `rotation_matrix(a,b,c)`: rotate by a,b,c about x,y,z. composed `Rz @ Ry @ Rx`, i.e. x rotation applied first -- chosen for movie purposes.
  - `streamframe_coords_observed(orbit, data_dict)`: the main one. see below.
  - `desi_RVerr(zmag, feh)`, `add_noise(...)`: survey error models. `add_noise` is still a stub.
  - everything under `if __name__=='__main__':` loops all orbits and plots phi1/phi2, so the function defs can be imported elsewhere without running it. comment that line out and un-indent if working interactively in the notebook cells.
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
- **potential is `gp.BovyMWPotential2014(units=galactic)` everywhere** -- `prog_position`, `integrate_prog_orbit`, `straighten_stream_orbit_interp`. petar was run with `external_mode='galpy'` to match. don't change one without the others.
- **`interrupt_mode='bse'`** is assumed by every loader (stellar evolution on).

## `extended_grid_info(scratch=True)`
holds paths, apocenters, stream ages (Vedant's), and init displacements as attributes. grid axes are orbit x stellar_pop (`lm`/`hm`) x rvir (`0.75, 1.5, 3, 6` pc, indexed 0-3) x copy. use `retrieve_sim_info(orbit, stellar_pop, rvir_index, copy) -> path, apo, age, init_displacement` rather than assembling paths by hand.

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

each subdict has: `coords` (the `StreamFrame` dict), `pos`/`vel` (galactocentric, with units), `mass`, `inMW`, `trim`, and `phi2_straight`/`r_straight`/`vr_straight`/`pm_phi1_straight`/`pm_phi2_straight`. `luminous`/`companions` also get `L`, `R`, `type`.

## masks: `inMW` and `trim` (and why it's always `[inMW][trim]`)
from `trim_coords_percentile(coords, low=1, high=99, apo=apo)`:
1. `inMW = coords['r'] <= 1.5*apo` -- drops stars flung well outside the orbit, i.e. **this is the only thing `apo` is used for**, so passing the wrong orbit's apocenter silently gives you the wrong outlier mask.
2. `trim` = 1st-99th percentile clip in phi1 and phi2, **computed on the already-`inMW`-selected array**.

so `trim` has length `inMW.sum()`, not `len(inMW)`, which is why every plotting line is `x[inMW][trim]` and why the two are not interchangeable or combinable with `&`.

## `straighten_stream_orbit_interp(coords, yval, core, i, ...)`
integrates the progenitor orbit +/- `Dt` Myr around time `i`, pushes that orbit through `StreamFrame` with the same reference coordinate, interpolates `yval` against `phi1`, and subtracts it -- so the output is a residual about the orbit track. `yval` can be a string or a list of strings (`intrinsic_stream_data_v3` passes `['phi2','r','vr','pm_phi1','pm_phi2']`).

`Dt` is found by an adaptive search, not fixed: it grows `Dt` until the orbit's phi1 range covers the data, but shrinks it if the orbit wraps around in phi1 (`dphi1 > 100` deg between samples), halving the step whenever it flips between those two cases. `Dt_start=240` from v3.

**it returns the iteration count as the last value -- check it.** the loop caps at `max_iters=100` and if it exhausts them it just falls through and straightens anyway, with no warning. an `itr` of 100 means the straightening never converged and the residuals are suspect. the returned residuals are for *all* particles (untrimmed), so apply `inMW`/`trim` afterwards.

TODO already noted in `inspect_new_sims.py`: `/old/prog_properties_summary.py` has a `poly_straightening()` that reportedly worked better than the `paf` versions, and `/old/DESI_comparison.py: get_GD1_coords_nbody()` has another orbit-interp implementation.

## photometry stuff [MAJOR WIP ⚠️]
synthetic photometry from a **blackbody**, with **top-hat filters** -- no real Gaia/SDSS response curves, so treat colors as approximate.
- `define_photometric_bands()` -- Gaia G (330-1050 nm), BP (330-680), RP (630-1050) from the DR2 paper, plus a 100 nm-wide z centered at 900 nm.
- `integrated_mag(nu_min, nu_max, T, R, d=10)` -- AB magnitude (the -48.60 zero point), `d` in pc, so `d=10` gives an absolute mag. **T and R must be cgs bare numbers** (K and cm) since `B(nu,T)` uses cgs constants -- this is why `inspect_new_sims.py` does `.cgs.value` first. note the 1/nu weighting normalization at the `np.trapezoid` lines is marked "suggestion that idk why works. " in the source and has not been validated against a real photometric zero point; worth checking against a known star before trusting absolute mags (colors are probably safer).
- `get_gaia_photometry(Teff, Radius, distance) -> [G, BP, RP]`, `get_z_photometry(...) -> z`. both are scalar-only, hence the slow per-star loop in `prepare_nbody_data(include_photometry=True)`. vectorizing these is the obvious win if photometry ever becomes a bottleneck.
- `m_from_M(M, dist)` for the distance modulus (dist needs astropy units here, unlike the above -- inconsistent, but that's how it is).

# stream frame coordinates ("observed" frames)
`streamframe_coords_observed(orbit, data_dict)` in `inspect_new_sims.py`. `data_dict` is one of the binary-treatment subdicts (`CoM`, `luminous`, `companions`) so it has `pos`/`vel` with astropy units. pipeline is:

`galcen pos/vel` -> `paf.galcen_to_ICRS` -> `gala.coordinates.reflex_correct` -> `.transform_to(<great circle frame>)`

returns a plain dict (like Jake's, but different keys): `phi1`, `phi2` [deg], `pm_phi1` (this is `pm_phi1_cosphi2`), `pm_phi2` [mas/yr], `distance` [kpc], `v_gsr` [km/s]. all bare floats, units stripped. note these are for *all* particles -- apply `inMW`/`trim` yourself when plotting.

`paf.galcen_to_ICRS(pos, vel)` wants shape (N,3) with units (it transposes internally). a single (3,) vector works too and gives a scalar coord, which is how the progenitor origins get built.

## how each frame is defined
| `orbit` | source | how |
|---|---|---|
| `gd1` | Koposov+2010 | `gc.GD1Koposov10`, built into gala |
| `pa5` | Price-Whelan+2018 | `gc.Pal5PriceWhelan18`, built into gala |
| `c19` | Ibata+2024 (see also Mohammed+2026) | `from_pole_ra0`, pole + alpha_0=354.356 deg |
| `jet` | Do+2026 | `from_pole_ra0`, pole + origin |
| `aau` | Shipp+2018 table 1 (the ATLAS half) | `from_endpoints` + `ra0` from the progenitor |
| `m3` | Yang+2023 sec 4.4 | `from_endpoints` + progenitor origin, `priority='pole'` |
| `circ` | -- | no frame; not a real stream, so it's skipped in the `__main__` loop |

phi1 zero points are set by the present-day progenitor sky position read out of `data/FINAL_ics_nolmc.csv` (except `jet`, where Do+26 give an origin directly). there is no `else` branch -- an unrecognized `orbit` string leaves `sc` undefined and you get a `NameError` at the `coords_stream['phi1']` line rather than a useful message.

## gala `GreatCircleICRSFrame` gotchas
version in `petar_env` is gala 1.9.1. the modern API is pole+origin only; passing `ra0=` or `rotation=` straight to the constructor raises. use the classmethods.

- **`from_pole_ra0` and `from_endpoints` are not actually different animals.** both just compute a pole and an origin and hand them to the same constructor. if you build a frame one way and then rebuild it from its own `.pole` and `.origin.ra` the other way, you get bit-identical phi1/phi2. so if one "works" and the other doesn't, suspect a typo before suspecting gala (ask me how I know).
- **the `priority` kwarg matters a lot.** the pole and origin have to be orthogonal. if your origin is off the great circle, gala emits only a `RuntimeWarning` (easy to miss in a notebook cell) and silently fixes it. default is `priority='origin'`, which *moves the pole* -- i.e. throws away the great circle you defined by endpoints. `priority='pole'` keeps the great circle and projects the origin onto it, which is what you want when the endpoints are the literature definition and the progenitor is only setting the phi1 zero point. passing `ra0=` instead of `origin=` also keeps the pole.
- for `m3` this is not academic: the `M3` row of `FINAL_ics_nolmc.csv` sits ~17 deg off the great circle through the Yang+23 endpoints, so `priority` changes where the endpoints land by >12 deg in phi2. also, `data/M3.fits` spans ra 189-269 deg with ~7.5 deg of phi2 scatter in every endpoint-derived frame, so one great circle may just not describe the whole M3 stream. see the TODO in the code -- may end up settling for "stream on an M3-like orbit," since that run is really a high-e / low-pericenter test.

# logistics:
## simulation data:
- scratch: `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/` (stopped keeping track of when my 90 days is up though; will be removed during monthly maintanence eventually)
- storage: `/n/holystore01/LABS/itc_lab/Users/amphillips/extended_grid/` TODO: add a storage path to `paf.extended_grid_info()` 

note that the circular orbit simulation is from the [Phillips+26](https://iopscience.iop.org/article/10.3847/1538-4357/ae680b) grid, so it is stored at `/n/holystore01/LABS/conroy_lab/Lab/amphillips/finished_grid/` in the directories that begin with 0-7 (see `PETAR_ANALYSIS_FUNCTIONS.py`'s `get_extended_grid_info`)

## conda environment:
`petar_env`, stored at `~/.conda/envs/petar_env`, containing standard packages like numpy, scipy, matplotlib, astropy, but in particular gala (1.9.1), plus `petar`, `sklearn`, `pygaia`. TODO: add `viamock` for the Via velocity errors.

nothing outside this env can import gala or petar, so run scripts with `~/.conda/envs/petar_env/bin/python` (or activate the env) rather than the system python. 

## petar documentation
the `README.md` from Long Wang's (PeTar github)[https://github.com/lwang-astro/PeTar] is useful. 