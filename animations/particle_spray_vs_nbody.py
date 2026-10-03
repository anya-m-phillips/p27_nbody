#----------------------------------------------------#
#   making a movie w a framerate of 10 Myr that      #
#   shows large-scale stream formation and zooms     #
#   in on the site of the progenitor so that we can  #
#   see that in ps there is no progenitor. for       #
#   the particle spray just use gala to generate     #
#   a stream on a GD-1-like orbit and for the nbody  #
#   use a dense GD-1-like progenitor. integrate      #
#   the gala stream over the same time/locations     #
#   as the N-body stream. ps zoom-in will be based   #
#   on the integerated progenitor orbit and nbody    #
#   zoom in will be centered on the origin of        #
#   the core frame from petar since the n-body       #
#   progenitors drift a little due to idk stellar    #
#   evolution? differential tidal drag? or somethin  #
#----------------------------------------------------#
# usage:
#   python particle_spray_vs_nbody.py --make-spray   # once, before any frames (writes SPRAY_FILE)
#   python particle_spray_vs_nbody.py -i 170          # one frame; index = petar file index, t = 10*index Myr
# %%
################ import packages
print("importing packages...")
import sys, os
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/scripts"

import petar
import numpy as np
import h5py


import astropy.units as u

import gala.dynamics as gd
import gala.potential as gp
from gala.dynamics import mockstream as ms
from gala.units import galactic

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

plt.style.use(script_path+'/vedant.mplstyle')

sys.path.append(script_path)
import PETAR_ANALYSIS_FUNCTIONS as paf
import argparse

print("done")
# %%

### WHICH SIMULATION
orbit = 'gd1'
stellar_pop = 'hm'
rvir_index = 0 # <-- 0.75 pc
copy = 0
rvir0 = 0.75*u.pc

grid_info = paf.extended_grid_info(scratch=False)
path, apo, age, init_displacement = grid_info.retrieve_sim_info(
    orbit=orbit, stellar_pop=stellar_pop, rvir_index=rvir_index, copy=copy
)
# init_displacement is the t=0 phase space position of the progenitor in kpc, km/s (see readme)

### TIME STEPPING
# petar snapshots every 10 Myr, so frame index == petar file index, t = 10*index Myr.
# the spray is integrated at 1 Myr and snapshotted every 10 Myr so its output index lines up.
frame_dt = 10 # Myr
spray_dt = 1 # Myr
n_steps = int(age/spray_dt) # <-- age already includes the +100 Myr, so the last frame is the present day
i_list = np.arange(0, int(age/frame_dt)+1) # 0 ... 270 for gd1

### PARTICLE SPRAY SETUP
N_spray = 5000 # <-- total at the present day, both tails
output_dir = '/n/netscratch/conroy_lab/Lab/amphillips/movies/particle_spray_vs_nbody/'
spray_file = output_dir+f'spray_{orbit}_{stellar_pop}_{rvir0.value:.2f}_{copy}.h5'

### PLOTTING
wide_lim = 30 # kpc
zoom_half_width = 50 # pc, i.e. 20 pc to an axis

#------------------------------------------------#
#                  FUNCTIONS                     #
#------------------------------------------------#
def nbody_bound_mass_track(path, t):
    """
    bound mass from petar's data.tidal (10 Myr cadence), linearly interpolated onto
    the spray integration times t [Myr]. returned with astropy units.
    """
    tidal = paf.load_tidal(path)
    return np.interp(t, tidal.time, tidal.mass)*u.Msun

def make_spray(path, init_displacement, n_steps, dt=spray_dt, output_every=int(frame_dt/spray_dt), N=N_spray):
    """
    run a Chen-type particle spray from the same t=0 progenitor position/velocity as petar,
    in the same potential (BovyMWPotential2014). writes gala's snapshot file to spray_file.

    release schedule: N/2 particles per tail, released one at a time at N/2 integration
    steps spread as evenly as possible over [0, n_steps], so the present day has exactly N.
    the df's progenitor mass follows the N-body bound mass; the progenitor *potential* is
    a fixed Plummer with the initial mass and b set by rvir0 (gala can't evolve it).
    """
    t = np.arange(n_steps+1)*dt
    n_per_tail = int(N/2)
    release_steps = np.round(np.linspace(0, n_steps, n_per_tail)).astype(int)
    assert len(np.unique(release_steps))==n_per_tail, "more releases than steps; lower dt"
    n_particles = np.zeros(n_steps+1, dtype=int)
    n_particles[release_steps] = 1

    prog_mass = nbody_bound_mass_track(path, t)
    plummer_b = rvir0*3*np.pi/16 # <-- virial radius of a plummer sphere is (16/3pi) b

    prog_w0 = gd.PhaseSpacePosition(init_displacement[:3]*u.kpc, init_displacement[3:]*(u.km/u.s))
    prog_pot = gp.PlummerPotential(m=prog_mass[0], b=plummer_b, units=galactic)
    df = ms.ChenStreamDF(lead=True, trail=True)
    H = gp.Hamiltonian(gp.BovyMWPotential2014(units=galactic))
    gen = ms.MockStreamGenerator(df, H, progenitor_potential=prog_pot)

    os.makedirs(output_dir, exist_ok=True)
    stream, prog = gen.run(prog_w0, prog_mass, dt=dt*u.Myr, n_steps=n_steps,
                           n_particles=n_particles,
                           output_every=output_every, output_filename=spray_file, overwrite=True,
                           progress=True)
    print(f"\nwrote {spray_file}: {len(stream.x)} particles at the present day")

def load_spray_snapshot(it):
    """
    returns (stream_pos (3,n) kpc, prog_pos (3,) kpc, time Myr) at output index it.
    particles not yet released are NaN in gala's file, so they're dropped here.
    """
    with h5py.File(spray_file, 'r') as f:
        stream_pos = f['stream']['pos'][:, it, :]
        prog_pos = f['nbody']['pos'][:, it, 0]
        time = f['stream']['time'][it]
    released = np.all(np.isfinite(stream_pos), axis=0)
    return stream_pos[:, released], prog_pos, time

def load_nbody_snapshot(path, file_index, core):
    """
    singles + binary CoMs (no binary orbital motion) from data.<file_index>.single/.binary,
    which are in the CORE frame. returns (core-frame pos (N,3) pc, galactocentric pos (N,3) kpc).
    """
    singles = petar.Particle(interrupt_mode='bse', external_mode='galpy')
    singles.loadtxt(path+f"data.{file_index}.single")
    binaries = petar.Binary(member_particle_type=petar.Particle, G=petar.G_MSUN_PC_MYR,
                            interrupt_mode='bse', external_mode='galpy')
    binaries.loadtxt(path+f"data.{file_index}.binary")

    pos_core = np.vstack([singles.pos, binaries.pos]) # pc, core frame
    pos_galcen = (pos_core + core.pos[file_index])/1e3 # kpc
    return pos_core, pos_galcen

#------------------------------------------------#
#              MAIN PROGRAM BELOWIDK             #
#------------------------------------------------#
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', '-i', type=int)
    parser.add_argument('--make-spray', action='store_true',
                        help='generate the particle spray snapshot file (run once, before the frames)')
    args = parser.parse_args()

    if args.make_spray:
        make_spray(path, init_displacement, n_steps)
        sys.exit()

    if not os.path.exists(spray_file):
        raise FileNotFoundError(f"{spray_file} doesn't exist; run with --make-spray first")

    it = args.index
    file_index = i_list[it]
    print(it, "t =", file_index*frame_dt, "Myr")

    ### N-body
    core = paf.load_core(path)
    nb_core, nb_galcen = load_nbody_snapshot(path, file_index, core)

    ### particle spray
    ps_galcen, ps_prog, ps_time = load_spray_snapshot(it)
    assert np.isclose(ps_time, file_index*frame_dt), "spray snapshot time doesn't match the petar snapshot"
    ps_rel = (ps_galcen - ps_prog[:, None])*1e3 # pc, relative to the integrated progenitor

    fig, axs = plt.subplots(2, 2, figsize=(14, 14))
    point_kw = dict(s=1, c='k', lw=0, rasterized=True)
    point_kw_zoom = dict(s=50, c='gold',edgecolor='k', lw=0.5, rasterized=True)

    # top row: zoomed out, galactocentric
    axs[0,0].scatter(nb_galcen[:,0], nb_galcen[:,2], **point_kw)
    axs[0,1].scatter(ps_galcen[0], ps_galcen[2], **point_kw)
    for ax in axs[0]:
        ax.set_xlim(-wide_lim, wide_lim)
        ax.set_ylim(-wide_lim, wide_lim)
        ax.set_xlabel(r'$x$ [kpc]')
        ax.set_aspect('equal')

        for width, height in [(2, 2), (3, 2), (0.5, 0.5), (30, 0.5)]:
            ax.add_patch(Ellipse(xy=[0, 0], width=width, height=height, angle=0,
                                facecolor='none', edgecolor='black', linewidth=1))


    axs[0,0].set_ylabel(r'$z$ [kpc]')

    # bottom row: zoomed in on the progenitor. N-body in the core frame, spray around the integrated orbit
    axs[1,0].scatter(nb_core[:,0], nb_core[:,2], **point_kw)
    axs[1,1].scatter(ps_rel[0], ps_rel[2], **point_kw_zoom)
    for ax in axs[1]:
        ax.plot(0, 0, marker='+', ms=20, mew=1.5, color='C3')
        ax.set_xlim(-zoom_half_width, zoom_half_width)
        ax.set_ylim(-zoom_half_width, zoom_half_width)
        ax.set_xlabel(r'$\Delta x$ [pc]')
        ax.set_aspect('equal')
    axs[1,0].set_ylabel(r'$\Delta z$ [pc]')

    axs[0,0].set_title(r'direct $N$-body')
    axs[0,1].set_title('particle spray')
    fig.suptitle(f'$t = {file_index*frame_dt/1e3:.2f}$ Gyr', y=0.93)

    filename = f"frame_{it:05d}.png"
    plt.savefig(output_dir+filename, dpi=150)
    plt.close()

# %%
