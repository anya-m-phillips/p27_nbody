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

`/animations`:
- `velocity_movie.py`, `grid_movie.py`, `grid_movie_long.py` : scripts to run with a slurm wrapper for animations. TODO: make all of these parallel so that the wrapper is submitted as an array job where each sub-job generates one frame. way faster than doing this in a loop. 

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

- `cocoons.py` summarizing cocoon fractions and properties as a function of progenitor properties. functions are called from `inspect_new_sims.py` as `simspect.[...]`. eyeballing hard boundaries to decide what gets counted as the cocoon. 
- `develop_GMM.py`: replacing the hard cuts in `cocoons.py` with an n-component (narrow ... narrow + cocoon) gaussian mixture, so cocoon membership is a fitted probability per star instead of a by-hand threshold per orbit. this is the point of the whole thing -- the manual cuts didn't standardize across orbits. the fitting layer is n-component general; the script currently runs 3. see the GMM section below. For messing around in the development stage. 

- `gmm.py`: where we will actually run the GMM. option for whether to constrain the cocoon widths in phi2 and v_gsr or not. Also an option to make plots. TODO: run again and experiment with outlier clipping, e.g., clipping outside 3-5 sigma (of the whole dataset) before fitting rather than picking hard lines of what to cut out. my cuts seem a little arbitrary and are based on looking at GD-1 orbit data. TODO: do the "trimming" step in the observed stream frame. 

- `results.py` making (nice) plots.

- `nfc_plots.py` was used for plotting in preparation for a conference; will likely abandon soon

- `noise.py`: adding survey noise and re-fitting the GMM. currently the live place where **MIST isochrone photometry** is being painted onto stars by initial mass, replacing the blackbody + top-hat scheme in `paf` (see the photometry section). reads the isochrone with `read_mist_models.ISOCMD`; `gaia_g_to_lsst_z` / `lsst_z_to_gaia_g` are the RTN-099 sec 1.3.4 polynomial and its exact inverse. the interpolation itself is `build_isochrone_table` / `isochrone_photometry` / `gaia_from_isochrone`, defined here rather than in `paf` -- **read the precision subsection before touching them**, the post-turnoff sequence spans 0.017 Msun and the tightest node spacing is 1.6e-11. `USE_ISOCHRONE = False` in the main block reverts to the blackbody path for comparison. the block at the bottom reconstructs the IMF from two snapshots -- **that block only works because it takes `mass0` from snapshot 0**, see the `star.mass0` section. for the per-star photometry use `lumdict['m0_zams']`, which `intrinsic_stream_data_v3` now provides directly.

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
- **`star.mass0` is NOT the ZAMS mass except at snapshot 0.** see its own section below -- this one is easy to get wrong silently because the field is literally labelled "initial stellar mass."

## `star.mass0` is an *effective* initial mass, and BSE rewrites it
if you want the IMF, or want to paint isochrone photometry onto stars by their birth mass, **read `mass0` out of `data.0` and match to the present day by `id`.** the present-day `mass0` is a different quantity.

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

recovers the birth mass of every surviving star, and `np.concatenate([m0_zams, m00[~np.isin(ids0, ids1)]])` reproduces `m00` exactly (verified with `np.sort`). this is what `noise.py` does when it reconstructs the IMF from two snapshots.

**as of 21 Sep 2026 you mostly don't have to do this by hand.** `intrinsic_stream_data_v3` now does the lookup internally and hands back *both* columns on the `luminous` and `companions` subdicts -- `m0_zams` (genuine birth mass) and `m0_effective` (the BSE-modified one). see the "two initial-mass columns" section under `intrinsic_stream_data_v3`. the hand-rolled version above is still what you want if you're working straight off `load_particle` output, or if you need the stars that have *disappeared*, which the subdicts by construction don't contain.

⚠️ naming: the field is **`star.mass0`** in the petar python reader (`petar/bse.py`) and **`m0`** in the C++ (`bse_interface.h`). there is no `star.m0` and no top-level `particle.mass0` -- both are `AttributeError`.

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
- `CoM` -- each binary as a single center-of-mass point. no `L`/`R`/`type`/`m0_*`.
- `luminous` -- each binary represented by its brighter component (`binaries.p1.star.lum >= binaries.p2.star.lum`).
- `companions` -- both components kept separately; this is the only one that double-counts.

each subdict has: `coords` (the `StreamFrame` dict), `pos`/`vel` (galactocentric, with units), `mass`, `inMW`, `trim`, `in_rtid`, and `phi2_straight`/`r_straight`/`vr_straight`/`pm_phi1_straight`/`pm_phi2_straight`. `luminous`/`companions` also get `L`, `R`, `type`, `m0_zams`, `m0_effective`.

### the two initial-mass columns (added 21 Sep 2026)
`luminous` and `companions` carry **both** flavours of initial mass, so you can pick:

- **`m0_zams`** -- the genuine birth mass. read out of `data.0` (`p0.star.mass0`, which at t=0 is exactly `p0.mass`) and matched to the present-day particle by `id`. **this is the one you want for the IMF and for isochrone photometry.**
- **`m0_effective`** -- the present-day `star.mass0`, i.e. whichever ZAMS track BSE is currently interpolating the star along. use it if you want to reproduce what BSE itself thinks the star is, e.g. when cross-checking `L`/`R`/`type` against an SSE track.

both are bare floats in Msun (no astropy units, unlike `mass`). neither is the *current* mass -- that's `mass`.

the two are genuinely different: on `gd1/hm/rvir=6/copy 0` at 2700 Myr, the `companions` subdict has 3786 of 14839 entries where they disagree, 1489 vs 841 above 1 Msun, and max 111.1 vs 61.0 Msun. see the `star.mass0` section for why.

implementation notes:
- **the ID match is done per component, before the `luminous` selection** (`sid`, `bp1id`, `bp2id` each get their own `searchsorted`, then `np.where(luminous_mask, ...)`), so the `IDs` ordering gotcha above does *not* bite here -- binaries get the right ZAMS mass regardless of which component is brighter. don't "simplify" this by slicing top-level `IDs`.
- each lookup `assert`s `ids0[index] == id`, so a present-day particle with no snapshot-0 counterpart fails loudly rather than silently grabbing a neighbour's mass. (there are none in the runs checked: IDs are unique, nothing is renumbered, and no new IDs appear -- only mergers remove them.)
- it costs one extra read of `data.0` per call, via `load_particle(path, 0, file_naming_convention="every integer")` -- deliberately not `load_coords_v2`, since only the `id` and `star.mass0` columns are wanted and the streamframe/core machinery is irrelevant at t=0.
- the lookup is rebuilt inside the `binary_treatment` loop, so it runs twice when both `luminous` and `companions` are requested. it's an argsort of ~15k ints, so it doesn't matter, but the source comment noting it is correct.

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
TODO: i think it would be way easier to just pull a MIST isochrone with gaia + the z-band to grab photometry given masses. would not be totally honest since a caveat of this project is that the stellar+dynamical ages of the streams are the same, however all of the massive stellar evolution stuff should be over by a few Gyr. the tophat photometric bandbass integration is easy to add bugs to. **this is what `noise.py` is doing now** -- and the thing it wants is the ZAMS mass, so interpolate against the isochrone's `initial_mass` column (not `star_mass`) using `lumdict['m0_zams']`, NOT `lumdict['m0_effective']` and NOT `lumdict['mass']` (current mass). see the `star.mass0` section above.

**as of 21 Sep 2026 that TODO is done** -- the isochrone interpolation lives in `noise.py` and is what the CMDs are made from now. the blackbody functions below are still in `paf` and still reachable (`USE_ISOCHRONE = False`), but nothing current uses them. see the next section; the `paf` bullets after it are the fallback path.

### isochrone interpolation: ZAMS mass -> Gaia mags (`noise.py`)
three functions, no state, all float64. the point of doing it this way rather than by blackbody is that **the isochrone is a free parameter**: swap the file + `age_index` at the top of `noise.py` and the whole stellar population changes, which is the thing the blackbody scheme couldn't do.

- `build_isochrone_table(iso, bands=ISO_BANDS, max_phase=5, mass_col='initial_mass')` -> `(m0_grid, {band: mags})`. `iso` is one age's rows, i.e. `isocmd.isocmds[age_ind]`. it sorts on `mass_col` and **drops any non-increasing node** rather than trusting the file, because `np.interp` requires increasing `xp` and does not check -- same trap as `straighten_stream_orbit_interp`, and it fails silently in exactly the same way. it prints if it drops anything (it drops nothing on the 12 Gyr file).
- `isochrone_photometry(m0_query, m0_grid, table)` -> `(phot, on_iso)`. **off-isochrone stars get `NaN`, not a clamped edge value**, and `on_iso` flags them. deliberate: a forgotten mask then shows up as a hole in the CMD rather than a fake pile-up at the tip of the AGB (contrast `straightened_obscoords_orbit_interp`, which *does* silently clamp).
- `gaia_from_isochrone(m0_query, iso, ...)` -> `(G, BP, RP, on_iso)`, the one-call wrapper.

mags are **absolute**, same convention as `g_phot(..., dpc=10)`, so `gaia_g_to_lsst_z` / `m_from_M` downstream need no change.

**`alive` is now two-sided.** it used to be `m0s <= max(m0_iso)`; it's `on_iso` now, because the sim's IMF reaches below the isochrone's low-mass end (0.102 Msun) as well as above its turnoff.

#### precision: this is the part with no headroom
the whole reason to be careful is that the evolved sequence is nearly degenerate in initial mass. on the 12 Gyr, [Fe/H]=-2 file (1460 rows, `initial_mass` strictly increasing over all of them, float64 in the file already):

| phase | N | m0 range |
|---|---|---|
| 0 MS | 203 | 0.1022 - 0.7884 |
| 2 RGB | 151 | 0.7887 - 0.8025 |
| 3 CHeB | 102 | 0.80249 - 0.80416 |
| 4 EAGB | 101 | 0.80416 - 0.80438 |
| 5 TPAGB | 601 | 0.804376 - 0.804398 |
| 6 post-AGB | 302 | 0.804398 - 0.80521 |

so **everything above the turnoff occupies 0.0168 Msun**, and the tightest node spacing is **1.58e-11 Msun**, on the TP-AGB. that is ~1e5 x float64 eps at 0.8 Msun, so it is resolvable, but:
- **never round, bin, or cast to float32 anywhere on the mass path.** interpolate against the raw `initial_mass` column and query with the raw `m0_zams` floats.
- **linear, not cubic.** on the RGB/AGB segments dM/dmag is ~1e-9; a spline overshoots enormously between nodes and invents points off the isochrone. linear is monotonic and can't.
- `initial_mass`, not `star_mass`, is also what makes the map invertible at all -- `star_mass` turns over once winds start.

checked: interpolation at the nodes reproduces the table **bit-exactly** (max |err| = 0 in all three bands), and the tightest node pair (1.58e-11 Msun apart) differs by 0.012 mag in G and 0.013 in BP-RP and comes back cleanly distinguished, with the midpoint strictly between. also verified shuffle-invariant, and NaN/`on_iso` correct on both sides of the range.

#### two things that are features, not bugs
- **interpolating in mass gets the relative numbers of giants right for free.** a phase is sampled in proportion to the initial-mass interval it occupies, which is exactly its lifetime x IMF weight. no extra weighting needed -- and it's why the giants come out rare despite the RGB having 151 of the 1460 rows.
- **`max_phase=5` is the default**, dropping post-AGB (phase 6, 302 rows). those are the proto-WD tail, and they're already cut from the sim side by `nonrem` (`type < 10`). set `ISO_MAX_PHASE = None` to keep them. note this moves `max_m0` from 0.80521 to 0.80440, so the `alive` count depends on it.

⚠️ the old `iso_cutoff = -700` magic number in the plotting cells is a *crude* stand-in for this -- index 760 of 1460 corresponds to m0 = 0.804388, i.e. partway through the TP-AGB, so it throws away most of phase 5 and all of 6 but keeps some. use `max_phase` instead; the `iso_cutoff` lines are only still there for overplotting the raw isochrone.

### the blackbody fallback (in `paf`, no longer the live path)
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

TODO: I am not 100% convinced that this is successful for M3 of Pal 5, which have notable diverging tails that get added to cocoons. 

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

# cocoon separation by gaussian mixture (`develop_GMM.py`)

adding __constraints__ seems to help quite a bit. a persistent issue has been that the thin stream component itself is non-gaussian, meaning that when the highly dispersed cocoon doesn't contain enough stars the optimizer is clearly just fitting two gaussians to the thin stream (really just by eye). This is a particular issue for M3 and Pal 5, with close-in orbits and significant width variations along the stream in multiple dimensions. The constraints have been to require that the cocoon component be at least $R$ times wider than the thin stream.

**$R$ is now per-dimension, which resolves the "possibly it should only be in one dimension" question** -- it doesn't have to be a single number, and it shouldn't be. from Jarvis+26 the GD-1 cocoon is ~10x wider in phi2 but only ~3x wider in radial velocity, so one scalar $R$ is either too weak in phi2 or actively excludes the right answer in v_gsr. the script runs `min_ratio = [10.0, 5.0]` on `dims = [0, 3]` (phi2, v_gsr) and leaves the proper motions unconstrained. see the constraints/bounds section below.

**the model is back to two components** (the `# , f_2` / `# , mu_3` bits are commented out in both the interactive cell and the grid loop). the earlier three-component version gave the thin stream two narrow gaussians to soak up its non-gaussian shape; the conclusion was that the constraints do that job more directly and honestly, so there's no reason to carry the extra 10 parameters. the fitting layer below is still n-component general and the n=3 path still works -- it's the *settings* that are n=2, so anything that says "the thin stream is components 0..n-2" degenerates to "component 0" right now.

an **n-component** gaussian mixture fit to the straightened observed-frame residuals, in place of the hard `get_cocoon_selection` cuts. the convention throughout is **narrowest first, cocoon last**: components `0 .. n-2` are the thin stream and the final (widest) one is the cocoon.

the data vector is `keys = ['phi2','pm_phi1','pm_phi2','v_gsr']`, stacked `(N, 4)` in that order, so **every `mu` and `sigma` row is a 4-vector in that same order** and mixing up the column order silently gives a nonsense fit. built on `sc_straighter` (orbit-interp straightening, then `poly_straightening`), masked by `unbound & ol_clip`. eventually the proper motions should become transverse velocities *before* straightening; not done yet.

**the n-1 convention.** everything that takes or returns `component_fractions` uses length `n_components - 1` -- the last weight is implicit, pinned by `sum(Q) = 1`, because only n-1 of them are actually free. so `n_components` is always inferred as `len(component_fractions) + 1`, and **the cocoon fraction is `1 - fracs.sum()`**, never an element of the array. `means` and `sigmas` are the full `(n_components, K)`.

## the likelihood
the model is: each star independently is drawn from one of the n components.

    L = prod_i [ sum_j Q_j p_j(x_i) ],   sum_j Q_j = 1

**the mixture weights go inside the product over stars.** the other ordering, `sum_j Q_j prod_i p_j`, is a different and wrong model -- it says the *whole stream* is drawn from a single component, and it collapses onto whichever one wins globally. (the first draft had it that way.)

everything stays in logs. there are two products and they become two sums:
- **inner**, over the 4 phase space dims: `norm.logpdf(x, mu, sigma).sum(axis=1)`. broadcasting does the loop.
- **outer**, over stars: the final `np.sum`.

the *one* place a sum-of-probabilities survives is the per-star mixing over components, and that's what `logsumexp` is for -- it factors out the largest term, so it's exact even when every term is ~ -3000. no intermediate ever leaves log space. (this was `np.logaddexp` when there were only two components; `scipy.special.logsumexp` with `axis=0` over an `(n_components, N)` array is the n-component version.)

- `component_likelihood(x_data, mu, sigma) -> (N,)` -- ln p per star under ONE component. returns `-inf` for sigma <= 0 so an optimizer stepping out of bounds gets rejected instead of NaN. unchanged from the 2-component version; it never knew about the mixture.
- `gmm_negative_loglikelihood(component_fractions, means, sigmas, data) -> scalar` -- the minimization objective. **`data` is last** so `minimize`'s `args` can supply it. its out-of-bounds guard returns `+np.inf`, not `-inf` -- getting that sign wrong rewards the optimizer for leaving the valid region. it raises if `len(means)` / `len(sigmas)` disagree with `len(component_fractions) + 1`, which is the one wiring mistake that would otherwise be silent.
- `component_responsibilities(x_data, fracs, means, sigmas) -> (n_components, N)` -- the full responsibility matrix, `R[j,i] = Q_j p_j(x_i) / sum_k Q_k p_k(x_i)`. columns sum to 1 by construction. **deliberately order-agnostic** -- it labels components in whatever order you hand them over, so `j` only means something physical if you sorted first.
- `component_membership_probability(x_data, fracs, means, sigmas, component=0, sort_dim=None) -> (N,)` -- **the indexable one.** `component` takes an int (negatives work, so `-1` is the cocoon once sorted), a sequence (`[0,1]` = the two narrow components together), or a slice (`slice(0,-1)` = everything but the cocoon, whatever n is). a sequence or slice *sums* the responsibilities, which is the right operation -- those are mutually exclusive posteriors, so "belongs to any of this set" is their sum. `sort_dim=None` takes the components as given; pass `sort_dim=-1` when the index is supposed to carry the narrow->cocoon meaning and you want a mislabelled fit to raise.
- `membership_probability(x_data, fracs, means, sigmas, sort_dim=-1) -> (N,)` -- thin wrapper on the above with `component=slice(0,-1)`: the probability of belonging to the thin stream, i.e. any of the n-1 narrower components. `p_cocoon = 1 - p_thin`. **this is the actual cocoon separation**, a soft per-star weight rather than a boolean cut. it enforces the sorted check, so run `sort_components` first.

the diagonal-sigma model *is* the multivariate one with `cov = np.diag(sigma**2)` -- verified identical to 1e-10 back when both existed (`gmm_likelihood_multivariate` is the commented-out block at the top of the script). the only thing given up is correlations *within* a component (e.g. phi2 with v_gsr along the track). those should be small post-straightening, but that's an assumption worth checking on the residuals rather than assuming.

a positive `ln L` is not a bug: these are log *densities*, and with sigma_phi2 ~ 0.1 deg the density exceeds 1. it also means `ln L` is not comparable across different unit choices -- only differences at fixed units mean anything. **relatedly: `component_likelihood` returns log densities, not probabilities, so a `p > 0.5` cut on its output is meaningless** -- if you want a per-star threshold, threshold a responsibility from `component_responsibilities` / `component_membership_probability`.

## fitting it: the packing layer
`scipy.optimize.minimize` wants **one flat 1-D array** as the objective's first argument. it cannot take `(fracs, means, sigmas)` -- numpy makes that a ragged object array and it dies with `ValueError: setting an array element with a sequence`. and `args=x_data` (no trailing comma) is not a tuple, so scipy iterates the array and splats it as separate arguments. it has to be `args=(x_data,)`.

so `pack_params` / `unpack_params` flatten to a **length `(n-1) + 2nK`** vector (17 for the n=2, K=4 fit running now; 27 for n=3), and `nll_flat(theta, x_data)` is what actually gets handed to `minimize`. layout:

    [alpha_1 ... alpha_{n-1},  mu_1, ln sigma_1,  mu_2, ln sigma_2,  ...]

the mu/ln-sigma pairs stay blocked **by component** rather than all-mus-then-all-sigmas, which keeps the n=2 layout bit-identical to the old one, so an old `result.x` still unpacks.

the packing also **reparameterizes to make the fit unconstrained**:
- `sigma -> log(sigma)`, so any real maps back to sigma > 0.
- the weights go through a **multinomial logit (softmax)**: `alpha_j = ln(Q_j / Q_last)` with the last component's alpha pinned at 0 as the reference. n-1 free reals then map onto the interior of the simplex -- all `Q_j > 0` and `sum(Q_j) = 1` automatically. for n=2 this reduces to exactly `logit`/`expit`, which is what it used to be. the inverse subtracts `logsumexp(alpha)` rather than dividing by a sum of exponentials, so it doesn't overflow for large |alpha|.

the optimizer then can't step into invalid territory and hit the `inf` walls, which is what wrecks a finite-difference gradient, and it fixes the conditioning: phi2 is ~0.1 deg while v_gsr is ~10 km/s, and in logs those are comparable steps.

two guards worth knowing about:
- `unpack_params` **infers `n_components` from `len(theta)`** via `divmod(len(theta) + 1, 2K + 1)`, and raises if the length isn't of that form. so a K/n mismatch fails loudly instead of misreshaping.
- `nll_flat` reads `K` off `x_data.shape[1]` rather than defaulting it, so `theta` can never be unpacked against the wrong number of phase space dimensions.

`nll_flat(theta, x_data, min_components=1)` is also **the single-gaussian null model** for the AIC/BIC comparison, and it needs no special-casing anywhere: `pack_params` turns a `(np.array([]), (1,K), (1,K))` triple into a length-2K theta with no weights in it, the validity guards in `gmm_negative_loglikelihood` are vacuously true on an empty fractions array, so `Qs = [1.0]`, `ln Q = 0`, and the single-row `logsumexp` is the identity. it comes back as the plain single-gaussian log likelihood. `min_components` stays at 2 by default so an *accidentally* emptied fractions array still fails loudly. note the call is `args=(x_data, 1)` -- the trailing-comma rule again.

`sort_components(fracs, means, sigmas, sort_dim=-1)` handles **label switching** -- the likelihood is exactly invariant under relabelling components, so the fit has no idea which one you meant to call "thin." it returns them ordered **narrowest -> widest**, and `sort_dim` picks which phase space dimension does the ordering (a component can be widest in one coordinate and not another). **the default is `sort_dim=-1`, i.e. v_gsr**, not phi2 as in the earliest 2-component version. the implicit weight is re-derived after sorting, so the returned fractions are the n-1 narrowest and the cocoon fraction is `1 - sum(returned)`. note the **sigma-ratio constraint does not sort for you** -- it hardcodes "cocoon = last component," so during the fit that's a definition rather than an observation, and `sort_components` is still what makes the index meaningful afterwards.

## constraints and bounds: the degeneracy fix
this is the part that made the fits usable, and it is two separate mechanisms that scipy treats differently. **both only do anything for `SLSQP`, `trust-constr` or `COBYLA`.** Powell and Nelder-Mead accept `constraints=` , emit only `RuntimeWarning: Method Powell cannot handle constraints`, and then ignore it and hand back a normal-looking result with `success=True` -- same genre of trap as the gala `priority` warning, and you end up thinking you constrained a fit you didn't. (Powell *does* honour `bounds`, though. only constraints get ignored.)

three general things:
1. a constraint function is called as `g(theta)` **only**. `args=(x_data,)` goes to the OBJECTIVE, not to the constraints -- if a constraint needs the data, close over it.
2. `constraints=` takes a **list** of constraint objects, one entry per constraint. `[]` or `None` is unconstrained.
3. for the legacy dict form the convention is `g(theta) >= 0`, i.e. "feasible when non-negative". `LinearConstraint`/`NonlinearConstraint(g, lb, ub)` state their bounds explicitly, so the dict form is the one that's easy to get backwards.

### `sigma_ratio_constraint(min_ratio, n_components, K, dims, thin_component=0)`
requires the cocoon to be at least `min_ratio` times **wider** than the thin component. this is a much better handle than bounding the cocoon fraction: bounding `f_cocoon` just clips the answer at whatever wall you put up, so the optimizer reports the bound rather than a fit. what actually goes wrong for a puffy progenitor is that the two components stop being *distinguishable* and the second one soaks up the non-gaussian wings of the thin stream. constraining the **separation** says what you actually mean -- "only call it a cocoon if it's much wider."

it is exactly **linear in theta**, since `ln sigma` is stored there directly: `sigma_cocoon/sigma_thin >= R` is `(ln sigma_cocoon)_k - (ln sigma_thin)_k >= ln R`, one `LinearConstraint` row per dimension, no reparameterization needed.

that row structure is also why the **per-dimension ratio is free**: `LinearConstraint` takes a *vector* lower bound, one entry per row, so `ln R` just stops being a scalar that broadcasts. `min_ratio` is either one number for every requested dimension or a list **aligned elementwise with `dims`**:

    constraint_dims = [0, 3]     # phi2, v_gsr  (keys order is phi2, pm_phi1, pm_phi2, v_gsr)
    min_ratio       = [10.0, 5.0]

a scalar is expanded to `len(dims)` inside the function rather than left to broadcast, so `len(ratios) == len(rows)` always holds and a wrong-length list raises instead of silently constraining a subset. `ratios <= 0` raises too (it's a width ratio, and `ln R` is the bound).

⚠️ **`thin_component=0` only, which matters if you go back to n=3.** the constraint compares the cocoon against exactly one thin component. with two narrow components the middle one is unconstrained relative to the cocoon, so you'd want a second constraint object with `thin_component=1`, or to apply the ratio against the widest thin component.

### bounding the means near zero (`pack_bounds`)
the other half of the fix. the sigma-ratio constraint stops the second component from being a slightly-wider copy of the thin stream; bounding `mu` stops it from **wandering off-track instead** and soaking up some *offset* blob of stars -- a diverging tail, the far end of an epicyclic feather -- as if it were a cocoon. the residuals are centred on the orbit track by construction, so a real cocoon shares the thin stream's mean: it is **wider, not displaced**.

`pack_bounds(fracs_bounds, means_bounds, sigmas_bounds, n_components, K)` takes bounds in **natural** parameter space and pushes them through the same monotonic reparameterization `pack_params` applies, so `(lo, hi) -> (T(lo), T(hi))` with no reordering. the three groups are very unequal in how useful they are:
- **means: the only group where a bound is exactly what it looks like.** `mu` sits in theta unchanged, so the pair passes straight through. nothing to get backwards.
- **sigmas: free, so this group costs nothing and buys nothing.** `sigma > 0` becomes `ln sigma > -inf`, i.e. unbounded -- positivity is already in the parameterization.
- **fracs: `logit`, and only honest for n=2.** `alpha_j = ln(Q_j/Q_last)` couples *all* the weights, so a box in Q space is not a box in alpha space for n>2 and `pack_bounds` raises rather than pretending. that n=2 case does cover bounding the cocoon fraction, since `f_cocoon = 1 - f_1` there; for n>2 use `cocoon_fraction_constraint(lo, hi, n_components)`, a `NonlinearConstraint` on `Q_last = 1/(1 + sum_j exp(alpha_j))` which is exact for any n.

the mean box is **+/- 1 sd of the data, per dimension**:

    mu_halfwidth = 1.0 * sd      # sd = x_data.std(axis=0)
    means_bound  = np.broadcast_to(np.column_stack([-mu_halfwidth, mu_halfwidth]),
                                   (ncomponents, len(sd), 2))

it has to be per-dimension, not one pair broadcast over the group: `sd` is ~0.1 deg in phi2 but ~10 km/s in v_gsr, so a single `(-1, 1)` would be ~10 sd in phi2 and ~0.1 sd in v_gsr. the `(n_components, K, 2)` form is what says "same rule, different number in each coordinate," and `_bound_pairs` reshapes it to `(nK, 2)`.

two sharp edges here:
- **`()` is not the unbounded spelling.** `_bound_pairs` reads an empty tuple as *zero* pairs and raises `got 0 (lo, hi) pairs, expected 8`. `None` -- either a whole group, or one side of a pair -- is what means unbounded.
- **an active mean bound is silent.** a `mu` sitting exactly on the wall is not a fitted mean, it's the optimizer reporting the bound, the same failure mode as bounding the cocoon fraction. scipy says nothing about it, so the script now prints a `WARNING mu[comp j, key] is ON its +/-x bound` per offending entry (and otherwise prints `max |mu|/sd`). that warning means the data wanted an offset component, which is the straightening leaving structure behind -- see the m3/pa5 diverging-tail TODO -- rather than evidence of a cocoon. verified on synthetic data both ways: a centred cocoon fits interior at `max |mu|/sd = 0.07`, while a deliberately offset wide blob pins `mu` on the phi2 and v_gsr walls and trips it.

## initial guesses: this is the part that actually matters
**do not start any two components identical.** identical components are an exact saddle point: if `p_j == p_k` then every star's responsibility is the prior weight regardless of the weight, the gradient wrt those weights is exactly zero, and the components can never split. confirmed empirically in the 2-component case -- from a `mu=0, sigma=1` both start, L-BFGS-B returns `f_thin = 0.5` with the two sigma vectors bit-identical. `sigma = 1` is also meaningless across mixed units (1 deg in phi2 is the whole stream; 1 km/s in v_gsr is nothing), so scale off `x_data.std(axis=0)`.

what the script does now, in both the interactive cell and the grid loop: `sigma_0 = [0.1, 10.0] * sd`, both `mu = 0`, and `f_1 = 0.9`. so the symmetry is broken **twice over** -- by a 100x gap in the width scales *and* by the starting weight. that 0.9 is doing real work and is not a neutral choice: it's a soft prior that the stream is mostly thin, and starting at `f_1 = 0.5` is the known failure mode below.

⚠️ **the multimodality caveat still stands.** starting at `f_1 = 0.5` splits *the thin stream itself* in two and calls the wider half a cocoon, which is why `f_1` starts at 0.9. the constraints now block the worst version of that (a "cocoon" only 1.5x wider can't satisfy `min_ratio`), but the likelihood surface is still genuinely multimodal, `minimize` still only ever finds a local optimum, and **there is no `success` flag that will warn about it**. so the reported cocoon fraction is still conditional on the starting basin and that should be stated in any writeup.

the independent check is physical, not numerical: cocoon fraction should *decrease* with increasing rvir (decreasing initial cluster density). the grid loop at the bottom of the script plots exactly that (`f_cocoon`, `sigma_phi2,cocoon`, `sigma_vgsr,cocoon` vs rvir, one line per orbit, ordered by pericenter), so a non-monotonic line there is the signal that fits in that row landed in different basins. still worth doing: scan the width-scale / fraction inits over a grid and keep the best `result.fun`.

## optimizer choice: the two-stage recipe
**SLSQP from a cold start collapses this particular likelihood onto a single component** -- both sigmas equal, thin weight -> 0. so it runs in two stages, and this is the actual recipe:

    # stage 1: Powell, to land in the right basin
    result_free = minimize(nll_flat, theta0, args=(x_data,), method='Powell',
                           bounds=bounds, options={'maxiter': 100000, 'maxfev': 100000})
    # stage 2: re-fit from there, WITH the constraint
    result = minimize(nll_flat, result_free.x, args=(x_data,), method='SLSQP',
                      bounds=bounds, constraints=constraints, options={'maxiter': 5000})

Powell finds the basin fast (the benchmark below) but ignores constraints; SLSQP respects them but can't find the basin. bounds and constraints coexist happily in stage 2.

**pass `bounds` to stage 1 as well, even though Powell can't use the constraints.** Powell does honour bounds, and if it isn't bounded it can walk a mean far outside the box -- at which point SLSQP *silently clips* `x0` back into the box and throws away the basin Powell was run to find. that clip is not an error and not a warning.

the constrained nll is necessarily `>=` the free one, and the gap is the useful number: it's how hard the data resist being told the cocoon must be much wider. the script prints both plus the difference.

benchmarked on synthetic data with 17 free parameters (n=2, K=4), N=30000, unconstrained:

| method | wall time | final nll | fitted f_thin (true 0.7495) |
|---|---|---|---|
| Nelder-Mead | 312 s | 34121.9 | 0.6915 |
| L-BFGS-B | 64 s | 33168.1 | 0.7492 |
| Powell | 5 s | 33168.1 | 0.7491 |

**all three returned `success=True`.** Nelder-Mead stopped ~950 nll units short of the optimum with a 6% error in `f_thin` and said it converged -- simplex methods degrade badly above ~10 dimensions, and n=3 would be 27 parameters, so it only gets worse. so: **check `result.fun`, never `result.success`.** with the good optimum, recovered sigmas match truth to <1% and per-star label recovery is 99.7%. (the `fatol`/`xatol` options still in the grid loop's `minimize` call are Nelder-Mead options and are ignored by Powell -- harmless, but they're not doing anything.)

`result.success` is now doubly useless: it's `True` for an ignored constraint, `True` for a clipped `x0`, and `True` for a fit pinned against an active mean bound. the three things worth actually reading are `result.fun`, the free-vs-constrained nll gap, and the mean-bound warning.

## model choice: AIC / BIC
`AIC(lnL0, k, N)` and `BIC(lnL0, k, N)`, with `k = len(theta0)` (i.e. `(n-1) + 2nK`, so 17 for n=2) and `N = len(x_data)`. **note `AIC` is actually AICc** -- it includes the small-sample correction `2k(k+1)/(N-k-1)` on top of `-2lnL + 2k`. with `N ~ 10^4` and `k = 17` that term is ~0.1, so it's numerically irrelevant here, but the name says AIC and the formula says AICc.

they take a **log likelihood**, not the nll, hence the sign flip at the call site (`L0 = -gmm_negative_loglikelihood(...)`). the intended comparison is against the single-gaussian null via `nll_flat(..., min_components=1)` -- see the packing layer section. the "removed stuff testing a single-component model here" comment marks where that used to live, so the null side of the comparison is **not currently wired up**; the printed AIC/BIC is a number for the n=2 fit alone, which on its own says nothing.

## profile likelihood in the cocoon fraction (written, commented out)
`nll_fixed_cocoon(psi, x_data, f_cocoon, n_components)` and `profile_cocoon_fraction(x_data, f_grid, means_0, sigmas_0, ...)` sit commented out above `sort_components`. they're the honest way to ask *"does the likelihood actually want a big cocoon, or is the optimizer just putting it there?"* -- pin `f_cocoon` on a grid, re-fit every other parameter at each grid point, look at the curve. it is independent of where the fit starts, of the reparameterization, and of the choice of method, so it settles the question a single `result.x` cannot.

`psi` has length `(n-2) + 2nK`: the mu/ln-sigma half is exactly `pack_params`' layout, but the thin components share the remaining `1 - f_cocoon` via their own softmax with `beta_1` pinned at 0, so there are n-2 free weight reals. **for n=2 that is zero** -- nothing left to split -- and `psi` is just the mu/ln-sigma block. `warm_start=True` walks the grid from the best-fitting end outward, seeding each point with the previous solution, which keeps the profile inside one basin (set it `False` to hunt for the second minimum at every point independently).

how to read the curve:
- **rising monotonically toward `f_cocoon -> 0`**: the data genuinely prefer a big cocoon. not an optimizer problem -- look at the *model* (are the two fitted sigmas actually distinct, or is component 2 soaking up non-gaussian wings?).
- **flat below some `f_cocoon`**: the cocoon weight is unconstrained down there, any small value fits as well as any other, and the number your fit reports is arbitrary. **report an upper limit, not a value.**
- **a second local minimum**: the multimodality above, made visible. the delta-nll between minima tells you how much it matters.

delta-nll of ~0.5 is the 1-sigma interval on one parameter, ~2 is 2-sigma.

## next
maximum likelihood first, then emcee for posteriors (per the header comment). the pieces are there: `nll_flat` negated is the right sign for `log_prob`, the softmax/log-sigma parameterization is already unconstrained so the sampler doesn't need bounds, and the `-inf` / `+inf` guards are what emcee expects for a rejected step -- add a prior and it's ready. sampling would also expose the multimodality directly, which point estimation hides. note that a sampler needs the label-switching convention applied *per sample* (`sort_components` on each draw) or the marginals come out as mush.

# bugs / stale things that remain
in rough order of how much they'd hurt:
- **`star.mass0` at the present day is not the ZAMS mass** -- BSE rewrites it on mass loss, mass transfer and especially remnant formation. it fails *silently*: you get a plausible-looking mass array that just isn't an IMF. use `m0_zams` off the `luminous`/`companions` subdicts, not `m0_effective`; if you're working off raw `load_particle` output, read `star.mass0` from snapshot 0 and match by `id` yourself. see its own section above.
- **`in_rtid` is the wrong length for the `companions` subdict**, so `in_rtid[inMW]` raises `IndexError` there. fine for `CoM` and `luminous`. details above.
- **`get_init_displacements.py` silently produces nothing** -- its `names_to_run` list still uses the long stream names that were renamed out of `FINAL_ics_nolmc.csv`.
- **anything cached from before the init_displacement units fix is wrong** -- not just the `*_straight` residuals but the intrinsic `coords` themselves, since the progenitor reference position and velocity were both bogus. regenerate.
- **the GMM likelihood is multimodal** and `minimize` only finds a local optimum -- no warning is emitted and `success=True` either way, so the reported cocoon fraction is conditional on the starting basin. check `result.fun`, and check the `f_cocoon` vs rvir trend for monotonicity. see the GMM section.
- **the grid loop at the bottom of `develop_GMM.py` is still bare unconstrained, unbounded Powell** -- no `bounds`, no `constraints`, `n=2` from a cold start. so the `f_cocoon` vs rvir trend plots, which are the one independent physical check on the whole method, are made with the *old* fits and get none of the benefit of the constraint/bound work. the two-stage recipe needs pushing into that loop before those panels mean anything.
- the grid loop reads the cocoon dispersions as `sigmas_fit[-1]`, with its own inline warning that **this is wrong if >2 components are fit** -- the cocoon could be two of three components. fine at n=2, a trap if n=3 comes back.
- `unpack_params`' docstring refers to **`nll_flat_anyn`**, which no longer exists -- it's `nll_flat(..., min_components=1)` now (`develop_GMM.py:421`).
- ~~the exploratory plotting cells cut on `p1 > 0.5` where `p1` is a `component_likelihood` output~~ **fixed** -- those cells now build `p1, p2` from `component_membership_probability` on sorted components, so the 0.5 threshold is a real responsibility. the general warning still stands: `component_likelihood` returns log *densities*, so never threshold its output.
- `straightened_obscoords_orbit_interp` clamps instead of flagging outside the orbit track's phi1 range.
- `straighten_stream_polynomial`'s `trim_criteria=None` default is a `TypeError`, not a default.
- no `else` branch in either `retrieve_sim_info` (`UnboundLocalError`) or `streamframe_coords_observed` (`NameError`) for an unrecognized orbit string.
- `extended_grid_info(scratch=False)` returns from `__init__` before assigning anything, and its hardcoded storage path disagrees with the one in logistics below.
- `load_coords_v2` builds the all-particles filename from the raw `file_index` kwarg, so `file_index=None` + `load_all=True` opens `data.None`.
- `core_to_galcen_frame` still adds the raw `core.vel` rather than the `fix_core_vel` version (see above).
- `correct_core` is dead code that also strips units.
- the `m3` great circle still doesn't describe the whole stream (see the gala gotchas section); the phi2 residual std of 0.786 vs ~0.2 for everything else is that showing up.
- ~~photometry is still blackbody + top-hat filters with an unvalidated normalization~~ **fixed** -- `noise.py` interpolates a MIST isochrone in `initial_mass` now, and MIST mags come with real response curves and a real zero point. the blackbody path in `paf` is still there behind `USE_ISOCHRONE = False` and still has the unvalidated 1/nu normalization, so don't trust *that* one's absolute mags if you go back to it. the live caveat is different now: the isochrone age is a free parameter that does **not** match the dynamical age, which is the whole reason it's swappable.
- the isochrone is old, so it stops at ~0.805 Msun and every more massive sim star is thrown out (`alive` / `on_iso`). that's a real selection, not a rounding detail -- at `hm` the discarded stars are the luminous ones, so quote how many got dropped alongside any luminosity-weighted number.
- `straighten_stream_orbit_interp_arbitrary_frame` in `paf` is a commented-out stub.
- `add_noise` is a stub, and `viamock` isn't in the env yet.

# logistics:
## simulation data:
- scratch: `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/` (90 days is running out though, and some sims have started being deleted.)
- storage: `/n/holystore01/LABS/itc_lab/Users/amphillips/extended_grid/` TODO: add a storage path to `paf.extended_grid_info()` 

note that the circular orbit simulation is from the [Phillips+26](https://iopscience.iop.org/article/10.3847/1538-4357/ae680b) grid, so it is stored at `/n/holystore01/LABS/conroy_lab/Lab/amphillips/finished_grid/` in the directories that begin with 0-7 (see `PETAR_ANALYSIS_FUNCTIONS.py`'s `get_extended_grid_info`)

## conda environment:
`petar_env`, stored at `~/.conda/envs/petar_env`, containing standard packages like numpy, scipy, matplotlib, astropy, but in particular gala (1.9.1), plus `petar`, `sklearn`, `pygaia`. TODO: add `viamock` for the Via velocity errors.

nothing outside this env can import gala or petar, so run scripts with `~/.conda/envs/petar_env/bin/python` (or activate the env) rather than the system python. 

## petar documentation
the `README.md` from Long Wang's [PeTar github](https://github.com/lwang-astro/PeTar) is useful. 