#--------------------------------------------------------#
#   Chen particle spray stream for one orbit, matched     #
#   to the rvir0 = 0.75 pc, hm N-body stream: same t=0    #
#   phase-space position and potential                    #
#   (BovyMWPotential2014) as petar, with the DF's         #
#   progenitor mass following the N-body bound mass       #
#   track from data.tidal. no progenitor potential, so    #
#   nothing gets trapped around a progenitor whose        #
#   mass doesn't track the DF's. saved as hdf5 (every     #
#   frame_dt snapshot, not just the present day) to       #
#   data/particle_sprays/, which is where motivation.py   #
#   looks for them. moved out of motivation.py so the     #
#   orbits can run in parallel (particle_spray.slurm).    #
#--------------------------------------------------------#
# usage:
#   python particle_spray.py --orbit gd1               # skips if the file already exists
#   python particle_spray.py --orbit gd1 --overwrite
print("loading packages")
import sys, os
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"
import argparse
import numpy as np

import astropy.units as u

import gala.dynamics as gd
import gala.potential as gp
from gala.dynamics import mockstream as ms
from gala.units import galactic

sys.path.append(script_path)
import PETAR_ANALYSIS_FUNCTIONS as paf
print("done")

print("defining stuff")
spray_dir = repo_path+"/data/particle_sprays/"

stellar_pop = 'hm'
rvir_index = 0 #<-- 0.75 pc
rvir0 = 0.75 #pc
N_spray = 5000 #<-- total particles at the present day, both tails
spray_dt = 1 #Myr
frame_dt = 10 #Myr, matches petar's snapshot cadence

grid_info = paf.extended_grid_info(scratch=False)

def find_finished_copy(orbit, stellar_pop, rvir_index, copy_options):
    """
    retrieve_sim_info()'s `copy` just picks a subdirectory; not every copy ran to the present
    day. try them in copy_options order (the same order final_datasets.py uses) and take the
    first whose data.tidal reaches the present day -- cheap, unlike prepare_nbody_data_anycopy.
    """
    for copy in copy_options:
        path, apo, age, init_displacement = grid_info.retrieve_sim_info(
            orbit=orbit, stellar_pop=stellar_pop, rvir_index=rvir_index, copy=copy)
        tidal = paf.load_tidal(path)
        if tidal.time.max() >= age:
            return path, apo, age, init_displacement, copy
    raise FileNotFoundError("no finished copy for %s rvir_index=%i"%(orbit, rvir_index))

def nbody_bound_mass_track(path, t):
    """bound mass from petar's data.tidal (10 Myr cadence), linearly interpolated onto t [Myr]."""
    tidal = paf.load_tidal(path)
    return np.interp(t, tidal.time, tidal.mass)*u.Msun

def make_spray(path, init_displacement, age, spray_file,
               dt=spray_dt, output_every=int(frame_dt/spray_dt), N=N_spray, progress=False):
    """
    Chen-type particle spray from the same t=0 progenitor phase-space position as petar, in
    the same potential (BovyMWPotential2014). writes every frame_dt snapshot to spray_file.
    the progenitor mass only sets the DF's release scale; the stream particles feel the
    MW potential alone (progenitor_potential=None).
    """
    n_steps = int(age/dt)
    t = np.arange(n_steps+1)*dt
    n_per_tail = int(N/2)
    release_steps = np.round(np.linspace(0, n_steps, n_per_tail)).astype(int)
    assert len(np.unique(release_steps))==n_per_tail, "more releases than steps; lower dt"
    n_particles = np.zeros(n_steps+1, dtype=int)
    n_particles[release_steps] = 1

    prog_mass = nbody_bound_mass_track(path, t)

    prog_w0 = gd.PhaseSpacePosition(init_displacement[:3]*u.kpc, init_displacement[3:]*(u.km/u.s))
    df = ms.ChenStreamDF(lead=True, trail=True)
    H = gp.Hamiltonian(gp.BovyMWPotential2014(units=galactic))
    gen = ms.MockStreamGenerator(df, H, progenitor_potential=None)

    stream, prog = gen.run(prog_w0, prog_mass, dt=dt*u.Myr, n_steps=n_steps,
                           n_particles=n_particles,
                           output_every=output_every, output_filename=spray_file, overwrite=True,
                           progress=progress)
    print("wrote %s: %i particles at the present day"%(spray_file, len(stream.x)))

print("done. main program below")
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--orbit', required=True, choices=['gd1','aau','pa5','jet','c19','m3'])
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--progress', action='store_true', help="gala's progress bar (messy in slurm logs)")
    args = parser.parse_args()
    orbit = args.orbit

    os.makedirs(spray_dir, exist_ok=True)

    copy_options = [4,3,2,1,0] if orbit=='aau' else [0,1,2,3,4] #<-- same order final_datasets.py uses
    path, apo, age, init_displacement, copy = find_finished_copy(orbit, stellar_pop, rvir_index, copy_options)
    print("%s: copy %i, age %i Myr, %s"%(orbit, copy, age, path))

    # same filename motivation.py looks for
    spray_file = spray_dir+"spray_%s_%s_%.2f_%i.h5"%(orbit, stellar_pop, rvir0, copy)
    if os.path.exists(spray_file) and not args.overwrite:
        print("%s already exists; pass --overwrite to regenerate"%spray_file)
        sys.exit(0)

    print("particle spraying %s"%orbit)
    make_spray(path, init_displacement, age, spray_file, progress=args.progress)
