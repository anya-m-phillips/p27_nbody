#------------------------------------------------------------#
#   print initial phase space displacement from zero for     #
#   the progenitor cluster when initializing petar sims      #
#------------------------------------------------------------#
# %%
import sys
# sys.path.append("/Users/anyaphillips/Downloads/software/petar_install/include") # for my machine
# script_path = '/Users/anyaphillips/Desktop/harvard/research/stream_velocity_structures/scripts'
script_path = "/n/home02/amphillips/p27_nbody/scripts" # for cannon

import petar
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
import astropy.units as u
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from astropy.table import Table
from scipy.stats import binned_statistic_2d


import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import gala.coordinates as gc
import gala.dynamics as gd
import gala.potential as gp
from gala.dynamics import mockstream as ms

from gala.units import galactic

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import pickle


from tqdm import tqdm


sys.path.append(script_path)
# from analyze_petar import analyze_petar
import PETAR_ANALYSIS_FUNCTIONS as paf
# from get_escapers import rms, make_key, dedupe_true_copies
import astropy.constants as const
from streamframe import StreamFrame
from scipy.stats import binned_statistic

import matplotlib.pyplot as plt
plt.style.use(script_path+'/vedant.mplstyle')
# %config InlineBackend.figure_format='retina'

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit
from scipy.stats import binom

import io

from scipy.ndimage import gaussian_filter1d
time_cmap = paf.define_time_cmap()
# %%
### load present-day prog locations + stream ages
#  based on particle spray streams from vedant
dat = Table.read('/n/home02/amphillips/p27_nbody/data/FINAL_ics_nolmc.csv')
### integrate some orbits, store names, pericenters, eccentricities. 
positions = np.vstack([dat['x'], dat['y'], dat['z']]).T * u.kpc
velocities = np.vstack([dat['vx'], dat['vy'], dat['vz']]).T * (u.km / u.s)
names = np.array(dat['name'].value)

mwp = gp.BovyMWPotential2014(units=galactic)
H = gp.Hamiltonian(mwp)


peris, apos, eccs = [], [], []

for k in tqdm(range(len(dat))):
    pos = positions[k]
    vel = velocities[k]
    w = gd.PhaseSpacePosition(pos, vel)

    orbit = H.integrate_orbit(w, dt=1*u.Myr, n_steps=10000)

    # if names[k] == 'gaia_11':
    #     orbit.plot()

    p, a, e = orbit.pericenter(), orbit.apocenter(), orbit.eccentricity()
    peris.append(p.to(u.kpc).value)
    apos.append(a.to(u.kpc).value)
    eccs.append(e.value)


peris = np.array(peris)
apos = np.array(apos)
eccs = np.array(eccs)
# %%
# summarize orbits:
lmc, hmc, simcolors = paf.define_simcolors()
fig, ax = plt.subplots(figsize=(8,8))
ax.scatter(peris, eccs, c=simcolors[0])


nl = ['GD-1','Palomar 5', 'M3', 'ATLAS-Aliqa Uma', 'Jet', 'C-19']
for x, y, name in zip(peris, eccs, names):
    plt.annotate(name, (x,y), xytext=(5,5),
                 textcoords='offset points', fontsize=12 if name in nl else 8,
                 color='k' if name in nl else '0.7')
    # plt.annotate(name, (x,y), textcoords='offset points', fontsize=12, color='0.7')
    # if name=='GD-1' or name=='Palomar 5':
    #     ax.scatter(x,y,marker="*", s=300, c=hmc[1], edgecolor='k')
    if name=='M3' or name=='ATLAS-Aliqa Uma' or name=='Jet' or name=='C-19' or name=='GD-1' or name=='Palomar 5':
        ax.scatter(x, y, marker="*", s=300, c=hmc[-1], edgecolor='k')
    
ax.set_xlim(0,21)
ax.set_ylim(0,1)
ax.set_xlabel(r'$r_{\rm{peri}}~\rm[kpc]$', fontsize=25)
ax.set_ylabel(r'$e$', fontsize=25)
plt.savefig('plots/orbit_summary.pdf', dpi=300, bbox_inches='tight')
# %%
#### now given ages, back-integrate to get ICs for petar
### integrate some orbits, store names, pericenters, eccentricities. 
positions = np.vstack([dat['x'], dat['y'], dat['z']]).T * u.kpc
velocities = np.vstack([dat['vx'], dat['vy'], dat['vz']]).T * (u.km / u.s)
names = np.array(dat['name'].value)
ages = np.array(dat['age_gyr'].value)
mprogs = np.array(dat['m_prog'])

mwp = gp.BovyMWPotential2014(units=galactic)
H = gp.Hamiltonian(mwp)

names_to_run = ['ATLAS-Aliqa Uma', 'Jet', 'M3', 'C-19', 'GD-1','Palomar 5']

outfile = open('data/init_displacements.txt', 'w') #<-- write these to a file to copy from later. 
def tprint(*a, **k):  # print to console and to init_displacements.txt
    print(*a, **k)
    print(*a, **k, file=outfile)

for k in tqdm(range(len(dat))):
    name = names[k]
    age_Myr = np.round(ages[k]*u.Gyr.to(u.Myr), -1).astype(int) # round to the nearest 10 so that we are guaranteed a "present day" output. 
    mprog = mprogs[k]
    if name in names_to_run:
        tprint(name)
        tprint("apocenter:", apos[k])


        pos = positions[k]
        vel = velocities[k]
        w = gd.PhaseSpacePosition(pos, vel)

        # back-integrate progenitor orbit
        orbit = H.integrate_orbit(w, dt=-1*u.Myr, 
                                  n_steps=age_Myr+100 # add 100 for stellar evolution effects. 
                                  )

        ### print initial condition, want x,y,z,vx,vy,vz in pc, km/s
        x=orbit.pos.x[-1].to(u.pc).value
        y=orbit.pos.y[-1].to(u.pc).value
        z=orbit.pos.z[-1].to(u.pc).value
        vx = orbit.vel.d_x[-1].to(u.km/u.s).value
        vy = orbit.vel.d_y[-1].to(u.km/u.s).value
        vz = orbit.vel.d_z[-1].to(u.km/u.s).value
        tprint("\t","%.20f,%.20f,%.20f,%.20f,%.20f,%.20f"%(x,y,z,vx,vy,vz))
        tprint("\t M_prog: %.2f Msun"%mprog)
        tprint("\t N=%i"%(mprog/0.5))
        tprint("\t age=%i Myr"%(age_Myr+100))
        tprint("\n")

        fig = orbit.plot()
        fig.suptitle(name)

outfile.close()
# %%