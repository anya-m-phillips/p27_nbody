# %%
script_path = "/n/home02/amphillips/stream_velocity_structures/scripts" # for cannon

print("loading packages.....")
import argparse
import petar
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import gala.coordinates as gc
import gala.dynamics as gd
import gala.potential as gp
from gala.dynamics import mockstream as ms
from gala.units import galactic

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec
plt.style.use(script_path+'/vedant.mplstyle')

from tqdm import tqdm

import sys

sys.path.append(script_path)
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
import astropy.constants as const
from scipy.stats import binned_statistic
from scipy.spatial import cKDTree
print("loaded")

time_cmap = paf.define_time_cmap()
init_displacements = paf.define_init_displacements()
# %%

# define helper functions
def dispersion_5_95(data):
    low, high = np.percentile(data, q=[5,95])
    return high-low

def get_initial_dispersion(path, core):
    particle_data0, streamframe_data0 = paf.load_coords_v2(
        path, 0, core, load_all=False
    )
    singles0, binaries0 = particle_data0
    single_coords0, binary_coords0 = streamframe_data0


    core_v0 = core.vel[0] * u.pc/u.Myr
    core_vx, core_vy, core_vz = core_v0
    core_v = np.sqrt(core_vx**2 + core_vy**2 + core_vz**2).to(u.km/u.s)

    bvx, bvy, bvz = binaries0.vel.T * (u.pc/u.Myr)
    bv = np.sqrt(bvx**2 + bvy**2 + bvz**2).to(u.km/u.s)
    svx, svy, svz = singles0.vel.T * (u.pc/u.Myr)
    sv = np.sqrt(svx**2 + svy**2 + svz**2).to(u.km/u.s)
    v = np.concatenate([sv, bv]) - core_v
    sigv = dispersion_5_95(v)

    return sigv

def get_local_densities(x,y,z, n_neighbors): ### <-- want to add this to paf eventually
    points = np.column_stack([x,y,z])
    tree = cKDTree(points)
    distances, _ = tree.query(points, k=n_neighbors+1)
    r_n = distances[:,-1] # distance to 200th nearest neighbor of each point
    volumes = (4/3) * np.pi * r_n**3
    number_densities = n_neighbors / volumes
    return number_densities


def make_frame(n, path, i):
    apo = apocenters[n]
    init_displacement = init_displacements[n]

    particle_data, streamframe_data = paf.load_coords_v2(path, i, use_core=False,
                                                         init_displacement=init_displacement,
                                                            # tdis_estimate=int(i-1),
                                                            load_all=False,
                                                            check_dissolved=False
                                                            )
    singles, binaries = particle_data

    spos_galcen, svel_galcen = paf.core_to_galcen_frame(path, singles, i, core=core)
    bpos_galcen, bvel_galcen = paf.core_to_galcen_frame(path, binaries, i, core=core)
    xs,ys,zs = spos_galcen.T.to(u.pc).value # * u.pc.to(u.kpc)
    Rs = np.sqrt(xs**2 + ys**2)
    xb,yb,zb = bpos_galcen.T.to(u.pc).value # * u.pc.to(u.kpc)
    Rb = np.sqrt(xb**2 + yb**2)

    # z, R = np.concatenate([zs, zb]), np.concatenate([Rs, Rb])

    x_all, y_all, z_all = np.concatenate([xs, xb]), np.concatenate([ys, yb]), np.concatenate([zs, zb])
    number_densities = get_local_densities(x_all, y_all, z_all, n_neighbors=150)

    single_coords, binary_coords = streamframe_data

    # concatenate single + binary com data
    sb_coords = paf.concat_sb_coords(single_coords, binary_coords)
    # trim step
    inMW, trim = paf.trim_coords_percentile(sb_coords, low=1, high=99, apo=apo)

    phi1 = sb_coords['phi1']


    ############# for 3D velocity dispersions -- also just going to fit a polynomial i think! 
    sb_vels = np.concatenate([singles.vel, binaries.vel])
    sb_vels = sb_vels * (u.pc/u.Myr).to(u.km/u.s) + core.vel[i]*(u.pc/u.Myr).to(u.km/u.s)  # km/s in galcen frame
    vx, vy, vz = sb_vels.T 
    v = np.sqrt(vx**2+vy**2+vz**2) 


    ### the straightening:
    bins = np.arange(-100, 101, 1) # degrees



    # _, vr_straight = paf.straighten_stream_binmeds(sb_coords['phi1'], sb_coords['vr'],  bins=bins,
    #                                                 trim_criteria=[inMW, trim],
    #                                                 sigma_kernel=5)
    sc, interp_orbit, orbit_coords, orbit_w, iter = paf.straighten_stream_orbit_interp(
        sb_coords, ['r', 'vr'], core, i, trim_criteria=[inMW, trim], return_orbit_chunk=True,
        correct_core_vel=True, use_core=False, init_displacement=init_displacement
    )
    r_straight, vr_straight = sc # <'straight coords'

    # orbit_phi1 = orbit_coords['phi1']
    # orbit_vx, orbit_vy, orbit_vz = orbit_w[:,3:].T *(u.kpc/u.Myr).to(u.km/u.s)
    # orbit_v = np.sqrt(orbit_vx**2 + orbit_vy**2 + orbit_vz**2) 
    # vx_straight = vx - np.interp(phi1, orbit_phi1, orbit_vx)
    # vy_straight = vy - np.interp(phi1, orbit_phi1, orbit_vy)
    # vz_straight = vz - np.interp(phi1, orbit_phi1, orbit_vz)
    # v_straight = v - np.interp(phi1, orbit_phi1, orbit_v) # why the FUCK ;jkl isn't this working ???? 
   
    # _, v_straight = paf.straighten_stream_binmeds(sb_coords['phi1'], v, bins=bins,
    #                                                 trim_criteria=[inMW, trim],
    #                                                 sigma_kernel=5)
    



    fig = plt.figure(figsize=[21,10])
    # gs = GridSpec(4,3,wspace=0.2, hspace=0.1, width_ratios=[1,0.5, 1])#, figsize=[25,9])
    gs = GridSpec(3,2, hspace=0.1, width_ratios = [2,1])


    ax1 = fig.add_subplot(gs[0,0])
    lo, hi = np.percentile(phi1[inMW][trim], [2,98])
    lim = np.max([np.abs(lo), np.abs(hi)])

    ax1.tick_params(labelbottom=False)

    ax1.set_ylim(-0.3,0.3)
    ax1.set_xlim(-lim, lim)

    ax1.scatter(phi1[inMW][trim], r_straight[inMW][trim], 
                c=vr_straight[inMW][trim], cmap='coolwarm', vmin=-2, vmax=2,
                s=.5, rasterized=True)
    ax1.set_ylabel(r'$\Delta r~[\rm kpc]$')



    ######## panel 2 vr vs phi1 ############
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    # just take 1-99% range for stream x-axis
    ax2.scatter(phi1[inMW][trim], vr_straight[inMW][trim],
                c='k',
                s=.5, rasterized=True)

    ax2.set_ylabel(r'$\Delta v_{r}\ \rm[km~s^{-1}]$')
    ax2.tick_params(labelbottom=False)

    ax2.set_ylim(-2,2)



    ######## dispersions ############
    ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
    ax3.set_ylabel(r'$\sigma_{v_r}\ \rm[km~s^{-1}]$')

    ### generate bin centers so that there are 200 systems per bin:
    N_per_bin = 150
    n_bins = max(1, len(phi1[inMW][trim])//N_per_bin)
    q = np.linspace(0, 1, n_bins+1)
    edges = np.quantile(phi1, q)
    edges = np.unique(edges) # helpful in cases where phi1 values repeat a lot..
    bins=edges.copy()

    bin_centers = (edges[1:]+edges[:-1])/2
    bin_widths = (edges[1:] - edges[:-1])

    sigvr_binned, binedges, binnumber = binned_statistic(x=phi1[inMW][trim],
                                        values=vr_straight[inMW][trim], 
                                        statistic=dispersion_5_95, 
                                        bins=bins)  


    sel = bin_widths<=5
    ax3.plot(bin_centers[sel], sigvr_binned[sel], c='k', label=r'$\sigma_{v_r}$')
    # ax3.set_yscale('log')

    ax3.set_ylim(0, 6)



    ax8 = fig.add_subplot(gs[:,1])
    ax8.set_xlabel(r'$R\ \rm[kpc]$')
    ax8.set_ylabel(r'$z\ \rm[kpc]$')
    # ax8.plot(R, z, c='cornflowerblue', lw=0.1)
    ax8.scatter(np.concatenate([Rs,Rb])[inMW][trim]*u.pc.to(u.kpc), 
                np.concatenate([zs,zb])[inMW][trim]*u.pc.to(u.kpc), 
                # c=vr_straight[inMW][trim], 
                c='k',
                cmap=time_cmap.reversed(), s=.1, rasterized=True)#, vmin=-1, vmax=1, s=.1)
    ax8.set_xlim(8,30)
    ax8.set_ylim(-22,22)


    ax8.text(0.05, 0.95, "%i Myr"%i, fontsize=40, va='top', ha='left', transform=ax8.transAxes)

    return fig, ax8

# %%


#------------------------------------------------------#
#               main program belowwww                  #
#------------------------------------------------------#

n=8 # idk guysss
# path = '/n/netscratch/conroy_lab/Lab/amphillips/finished_grid/8_gd1_rvir0.75_lm/'
path = '/n/holystore01/LABS/conroy_lab/Lab/amphillips/finished_grid/optimized_dense_sims/8_gd1_rvir0.75_lm/'
apocenters = paf.define_apocenters()
core=petar.Core(interrupt_mode='bse', external_mode='galpy')
core.loadtxt(path+"data.core")



# i_list = np.arange(2750, 10805, 5)
i_list = np.arange(3000, 15000, 5)
# i_list = [5300]

# mwp = gp.BovyMWPotential2014(units=galactic)

# H = gp.Hamiltonian(mwp)
# Dt = 1
# orbit = H.integrate_orbit()


it=0
p = '/n/netscratch/conroy_lab/Lab/amphillips/movies/rafa_velocities/'
for i in tqdm(i_list):
    fig, ax8 = make_frame(n, path, i)

    filename = f"frame_{it:05d}.png"
    plt.savefig(p+filename)
    plt.close()
    it+=1 

print(".....idk man")
# %%
