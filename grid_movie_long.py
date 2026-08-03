# %%
################ import packages
print("importing packages...")
import sys
script_path = "/n/home02/amphillips/p27_nbody/scripts" # for cannon

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
from streamframe import StreamFrame
import PETAR_ANALYSIS_FUNCTIONS as paf
import argparse

print("done")
# %%

### DEFINE STUFF ABOUT THE GRID: 
grid_info = paf.extended_grid_info(scratch=True) 
# %%
def rotation_matrix(a,b,c):
    """
    rotate angles a,b,c about x,y,z axes
    https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
    """
    Rx = np.array([
        [1,0,0],
        [0, np.cos(a), -np.sin(a)],
        [0, np.sin(a), np.cos(a)]
    ])
    Ry = np.array([
        [np.cos(b), 0, np.sin(b)],
        [0,1,0],
        [-np.sin(b), 0, np.cos(b)]
    ])
    Rz = np.array([
        [np.cos(c), -np.sin(c), 0],
        [np.sin(c), np.cos(c), 0],
        [0,0,1]
    ])

    # R = Rx @ Ry @ Rz
    R = Rz @ Ry @ Rx # <-- do x rotation first, _then_ z rotation.. mostly for movie purposes. 
    return R

# main helper function: 
def prepare_nbody_data(path = grid_info.gd1_lm_paths[0]+"0/", 
                       include_photometry=False, 
                       i=grid_info.gd1_age,
                       apo = grid_info.gd1_apo,
                       init_displacement = grid_info.gd1_init_displacement,
                       ):
    """
    this is expecting N-body data that outputs every 10 Myr. 
    """
    core = paf.load_core(path)
    file_index = int(i/10)

    data_dict = paf.intrinsic_stream_data_v3(
        path, i, core, apo, init_displacement,
        use_core=False,
        binary_treatments=['CoM','companions','luminous']
    )

    CMdict = data_dict['CoM']
    lumdict = data_dict['luminous']
    inMW, trim = CMdict['inMW'], CMdict['trim']


    if include_photometry==True:
        # do the calculations -- ORDERED AS ALL_PARTICLES
        primary_Ls = lumdict["L"][inMW][trim] #* u.Lsun
        primary_Rs = lumdict['R'][inMW][trim] #* u.Rsun
        primary_Teffs = (primary_Ls / (4*np.pi*primary_Rs**2 * const.sigma_sb))**(1/4)

        R_cgs = primary_Rs.cgs.value
        Teff_cgs = primary_Teffs.cgs.value
        G, BP, RP, z = [],[],[],[]
        for Tval, Rval in tqdm(zip(Teff_cgs, R_cgs)):
            Gval, BPval, RPval = paf.get_gaia_photometry(Tval, Rval, 10) # 10 pc. 
            G.append(Gval)
            BP.append(BPval)
            RP.append(RPval)

            zval = paf.get_z_photometry(Tval, Rval, 10) # 10pc
            z.append(zval)
        G = np.array(G)
        BP = np.array(BP)
        RP = np.array(RP)
        z = np.array(z)

        return core, data_dict, CMdict, lumdict, inMW, trim, G, BP, RP, z
    
    else:
        return core, data_dict, CMdict, lumdict, inMW, trim
# %%

orbits = ['circ','gd1','aau','pa5','jet','m3','c19']

lm_colors, hm_colors, simcolors = paf.define_simcolors()
reordered_colors = hm_colors + lm_colors[::-1]
cc = reordered_colors[:-1]
theta_list = np.arange(0, 361, 1)*u.degree.to(u.radian)
# it = 0

i_list = np.arange(3950, -10, -10) # <-- set by the youngest stream, AAU for now
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

        if loading_time<0: #<-- don't load sims that haven't started yet!
            continue
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


    dir='/n/netscratch/conroy_lab/Lab/amphillips/movies/evolving_grid_long/'
    filename = f"frame_{it:05d}.png"
    plt.savefig(dir+filename)
    plt.close()

# %%
