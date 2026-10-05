#--------------------------------------------------------------#
#  why are the noisy "binaries" cocoons so much hotter in      #
#  v_gsr than the noiseless one? (5 Oct) -- see "binaries in   #
#  the observed samples" in the readme's GMM section.          #
#--------------------------------------------------------------#
# two checks on one pickle:
#   1. binary fraction + spread of binary orbital motion in v_gsr
#      for the full unbound sample vs the catalog-matched samples
#   2. NOISELESS GMM fits on exactly the stars each survey case
#      keeps, so selection is separated from noise.
#
# usage:
#   python claude_binary_selection_experiment.py                  # gd1, rvir0 = 0.75
#   python claude_binary_selection_experiment.py --orbit pa5 --rvir 3
import sys
repo_path = "/n/home02/amphillips/p27_nbody"
sys.path.append(repo_path)
import argparse
import pickle
import numpy as np
from scipy.optimize import minimize

import gmm #<-- main program is behind __main__, so this only gets the functions


datapath = '/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/'
keys = ['phi2','v_phi1','v_phi2','v_gsr'] #<-- same as gmm.py
N_jarvis = 679 #<-- desi case keeps the brightest N_jarvis catalog-matched stars


def fit_gmm(sc, use):
    """noiseless two-component fit, same initial guess + optimizer as gmm.py (unconstrained, Powell)."""
    x_data = np.column_stack([sc[k][use] for k in keys])
    sd = x_data.std(axis=0)
    theta0 = gmm.pack_params(np.array([0.9]),
                             np.zeros((2, len(keys))),
                             np.array([0.1*sd, 10.0*sd]))
    result = minimize(gmm.nll_flat, theta0, args=(x_data,), method='Powell',
                      options={'maxiter': 100000, 'maxfev': 100000})
    fracs, means, sigmas = gmm.unpack_params(result.x)
    fracs, means, sigmas = gmm.sort_components(fracs, means, sigmas)
    return 1 - fracs.sum(), sigmas #<-- cocoon fraction, (2, K) sigmas narrowest first


def survey_selections(data_dict, sfx):
    """
    the full unbound sample and the stars the via / desi cases keep, for one
    binary treatment. (gmm.py's pm/RV error cuts are left out: on 30 Sep none of
    them removed a catalog-matched star.)
    """
    N = data_dict['nsingles'] + data_dict['nbinaries']
    raw = data_dict['trim_new'+sfx] & data_dict['unbound']

    matched_full = np.zeros(N, dtype=bool)
    matched_full[data_dict['cut_for_catalog_photometry'+sfx]] = data_dict['matched_to_catalog_photometry'+sfx]

    top = np.zeros(N, dtype=bool)
    top[np.argsort(data_dict['catalog_photometry'+sfx]['mG'])[:N_jarvis]] = True #<-- nans sort to the end

    return {
        'all unbound': raw,
        'via-matched': matched_full,
        'desi top %d'%N_jarvis: matched_full & top,
    }


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--orbit', default='gd1')
    parser.add_argument('--rvir', type=float, default=0.75)
    args = parser.parse_args()

    with open(datapath+'%s_%.2f.pickle'%(args.orbit, args.rvir), 'rb') as handle:
        data_dict = pickle.load(handle)

    ns, nb = data_dict['nsingles'], data_dict['nbinaries']
    isbin = np.zeros(ns+nb, dtype=bool)
    isbin[ns:] = True
    # binary orbital motion in v_gsr (+ a tiny difference between the two poly straightenings)
    dv = data_dict['sc_straighter_primaries']['v_gsr'] - data_dict['sc_straighter']['v_gsr']
    m0_lum = np.asarray(data_dict['luminous']['m0_zams'])

    print("%s, rvir0 = %.2f pc"%(args.orbit, args.rvir))

    #### 1. who's in each sample
    print("\n1. binary fraction and binary orbital motion, binaries-case trim (trim_new_primaries)")
    sels_primaries = survey_selections(data_dict, '_primaries')
    for name, use in sels_primaries.items():
        b = use & isbin
        print(f"  {name:14s} N={use.sum():6d}  f_bin={b.sum()/use.sum():.3f}  "
              f"std dv={dv[b].std():6.2f} km/s  frac |dv|>5 km/s={np.mean(np.abs(dv[b])>5):.3f}  "
              f"median m0_zams (luminous) of binaries={np.median(m0_lum[b]):.3f}")

    # binary fraction down the brightness ranking of the eligible pool
    print("  eligible pool, ranked by isochrone mG:")
    pool = data_dict['cut_for_catalog_photometry_primaries']
    mG = data_dict['phot']['mG']
    rankable = np.where(pool & np.isfinite(mG))[0]
    ranked = rankable[np.argsort(mG[rankable])]
    n_catalog = sels_primaries['via-matched'].sum()
    edges = [0, N_jarvis, n_catalog, (n_catalog+len(ranked))//2, len(ranked)]
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = ranked[lo:hi]
        b = isbin[s]
        print(f"    ranks {lo:5d}-{hi:5d}: f_bin={b.mean():.3f}  std dv(binaries)={dv[s][b].std():6.2f} km/s")

    #### 2. noiseless GMM on each selection
    print("\n2. NOISELESS gmm fits on each selection (cocoon fraction, cocoon sigma_v_gsr)")
    for sfx, case in [('', 'CoM'), ('_primaries', 'binaries')]:
        sc = data_dict['sc_straighter'+sfx]
        for name, use in survey_selections(data_dict, sfx).items():
            fc, sigmas = fit_gmm(sc, use)
            print(f"  {case:8s} {name:14s} N={use.sum():5d}  f_cocoon={fc:.3f}  "
                  f"S_thin={np.round(sigmas[0], 2)}  S_cocoon={np.round(sigmas[-1], 2)}")
