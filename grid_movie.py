# %%
################ import packages
print("importing packages...")
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
script_path = repo_path+"/script"

import petar
import numpy as np

import astropy.units as u
import astropy.constants as const


import matplotlib.pyplot as plt
# %matplotlib inline
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
plt.style.use(script_path+'/vedant.mplstyle')
# %config InlineBackend.figure_format='retina'

from tqdm import tqdm

sys.path.append(script_path)
sys.path.append(repo_path)
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
from inspect_new_sims import prepare_nbody_data, rotation_matrix
import argparse

print("done")
# %%

### DEFINE STUFF ABOUT THE GRID: 
grid_info = paf.extended_grid_info(scratch=True) 
# %%
orbits = ['circ','gd1','aau','pa5','jet','m3','c19']

lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
theta_list = np.arange(0, 361, 1)*u.degree.to(u.radian)
# it = 0

i_list = np.arange(2360, -10, -10) # <-- set by the youngest stream, AAU for now
#%%
# print(len(i_list))
#%%
theta = 0*u.degree.to(u.radian) # <--- not rotating azimuthally for now. 
# step 2: prepare to rotate everyone
phi =  5*u.degree.to(u.radian) 
R = rotation_matrix(a=phi, b=0, c=theta)

#------------------------------------------------#
#              MAIN PROGRAM BELOWIDK             #
#------------------------------------------------#  
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', '-i', type=int)    

    args = parser.parse_args()

    print(args.index)

    it = args.index #<--index
    i = i_list[it]

    fig, ax = plt.subplots()
    # go through each dictionary
    e0,e1,e2 = np.array([]), np.array([]), np.array([])
    n_sim = np.array([])

    for ii, orbit in tqdm(enumerate(orbits)):
        path, apo, age, init_displacement = grid_info.retrieve_sim_info(
            orbit=orbit, stellar_pop='lm', rvir_index=0, copy=0
        )

        loading_time = int(age-i)

        # print(orbit, age)
        # print("loading at", loading_time)


        core, data_dict, CMdict, lumdict, inMW, trim = prepare_nbody_data(
            path = path,
            include_photometry=False,
            i=loading_time if orbit !="circ" else int(30000-i), #age if orbit != "circ" else 30000,
            init_displacement = init_displacement
        )
        # dicts.append(data_dict)


        x, y, z = data_dict['CoM']['pos'].T.to(u.kpc)
        pos = data_dict['CoM']['pos'].to(u.kpc)

        rotated_pos = (R@pos.T).T

        e0_, e1_, e2_ = rotated_pos.T
        e0 = np.append(e0, e0_.to(u.kpc).value)
        e1 = np.append(e1, e1_.to(u.kpc).value)
        e2 = np.append(e2, e2_.to(u.kpc).value)
        n_sim = np.append(n_sim, 
                        np.array([ii]*len(e0_)))

    # reorder all of the points along the axis pointing out of the page.         
    reordered = np.argsort(e1)[::-1]
    # cmap = LinearSegmentedColormap.from_list('cmap', cc)
    n_orbits = len(orbits)
    cmap = mcolors.ListedColormap(cc[:n_orbits])
    # discrete norm: one color band per orbit, boundaries on the half-integers
    bounds = np.arange(n_orbits + 1) - 0.5   # [-0.5, 0.5, ..., n_orbits-0.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    ob = ax.scatter(e0[reordered], e2[reordered],
            c=n_sim[reordered],
            cmap=cmap,
            norm=norm,
            s=1)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right','5%',pad=0.05)
    cbar = fig.colorbar(ob, cax=cax, ticks=np.arange(n_orbits))
    cbar.ax.set_yticklabels(orbits)   # bottom -> top follows the orbits list order

    ax.set_xlim(-40,40)
    ax.set_ylim(-40,40)
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_xticklabels([])


    dir='/n/netscratch/conroy_lab/Lab/amphillips/movies/moving_grid_rotation/'
    filename = f"frame_{it:05d}.png"
    plt.savefig(dir+filename)
    plt.close()
