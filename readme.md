an N-body grid of mock streams to probe effects of progenitor dynamics on stream dispersion

# map of code, briefely:
`/data`:
- `FINAL_ics_nolmc.csv`: credit Vedant Chandra, present-day positions, velocities of some promising streams. Age and progenitor mass estimates also given.
- `init_displacements.txt`: generated in `get_init_displacements.py`, Vedant's progenitor locations back-integrated by their stream ages, plus like 100 Myr to account for initial expansion due to massive star evolution, rounded to a multiple of 10 so that I can safely output sim snapshots every 10 Myr and get the present day in the final snapshot.
- `ATLAS-Aliqa Uma.fits`, `C-19.fits`, `Jet.fits`, `M3.fits`: credit Bonaca & Price-Whelan (2025); TODO: pull GD-1 and Pal 5 as well. 

`/scripts`:
- `PETAR_ANALYSIS_FUNCTIONS.py`: a bunch of functions, mostly migrated from other projects; imported as `paf`
- `streamframe.py`: credit Jake Nibauer, transformation to stream frame coordinates as seen from the Galactic center
- `vedant.mplstyle`: credit Vedant Chandra, plotting style stuff

`/old`: stuff that is migrated from older repositories (esp. `~/stream_velocity_structures`)
- `DESI_comparison.py`: noise old simulation data like Gaia+DESI and compare to Jarvis+2026 cocoon detection. Will eventually systematize and do for all new sims
- `prog_properties_summary.py`: cocoon fractions/dispersions for GD-1 portion of old sim grid

scripts in top directory for now:
- `get_init_displacements.py`: used to generate init_displacements.txt, fed as inputs to `petar.init ...` for displacing progenitors from the galactic center at the initial condition. 
- `inspect_new_sims.py`: writing a bunch of functions to process sim data alongside `paf`; in particular transforming to ``observed" stream frame coordinates (based on ICRS coordinates; correcting for solar reflex motion, etc)
- `velocity_movie.py`, `grid_movie.py`, `grid_movie_long.py` : scripts to run with a slurm wrapper for animations. TODO: make all of these parallel so that the wrapper is submitted as an array job where each sub-job generates one frame. indescribably faster than doing this in a loop. 
- `nfc_plots.py` was used for plotting in preparation for a conference; will likely abandon soon

SIMULATION DATA LOCATION:
scratch: `/n/netscratch/conroy_lab/Lab/amphillips/extended_grid/`
storage: `/n/holystore01/LABS/itc_lab/Users/amphillips/extended_grid/` TODO: add a storage path to `paf.extended_grid_info()` 