"""
a one-off diagnostic, not part of the pipeline: profile likelihood in the
thin-stream sigma_vgsr for the desi_noise_CoM GD-1 fit, with x_err.

    python claude_profile_desi.py 0.75      # rvir0; ~a few minutes

rebuilds gmm.py's desi selection (catalog-matched, top 679 in mG, rverr < 10)
from the pickle, then prints:
  - a noiseless fit on the same stars (the truth to compare against),
  - the noisy fit with errors,
  - 1D mle_sigma on the p_thin > 0.5 stars,
  - the profile: thin sigma_vgsr fixed on a grid, every other parameter
    refit (Powell, warm-started), delta nll vs the best fit. a profile that's
    flat down to 0 = unresolved width; quote the value where delta nll = 1.35
    as the one-sided 95% upper limit.
written 14 Oct to show that a thin sigma_vgsr = 0 from the seed-42 noise draw
was a ~2 sigma low fluctuation, not a resolution limit.
"""
import sys, pickle
sys.path.insert(0, '/n/home02/amphillips/p27_nbody')
import numpy as np
import astropy.units as u
from scipy.optimize import minimize
import gmm

rvir = float(sys.argv[1]) if len(sys.argv) > 1 else 0.75
keys = ['phi2', 'v_phi1', 'v_phi2', 'v_gsr']
d = pickle.load(open('/n/netscratch/conroy_lab/Lab/amphillips/p27_data_dicts/gd1_%.2f.pickle' % rvir, 'rb'))
dist = d['coords_obs'].distance
sc = d['sc_straighter']
cut = d['cut_for_catalog_photometry']
matched = np.zeros(len(cut), bool); matched[cut] = d['matched_to_catalog_photometry']
phot = d['catalog_photometry']
nd = dict(d['noise_catalog_photometry'])
nd['v_gsr'] = nd['v_gsr_desi']
rverr = nd['rverr_desi']
nd['v_phi1'] = dist.to(u.km).value * (nd['pm_phi1']*u.mas/u.yr).to(u.radian/u.s).value
nd['v_phi2'] = dist.to(u.km).value * (nd['pm_phi2']*u.mas/u.yr).to(u.radian/u.s).value
cf = (u.microarcsecond/u.yr).to(u.mas/u.yr)
top = np.zeros(len(cut), bool); top[np.argsort(phot['mG'])[:679]] = True
use = cut & matched & (nd['pm_err_gaia']*cf < 0.5) & top & (rverr < 10) & d['nonrem']

x = np.column_stack([sc[k][use] + nd[k][use] for k in keys])
xt = np.column_stack([sc[k][use] for k in keys])           # the TRUE values of the same stars
xe = gmm.measurement_errors(nd, 'desi', dist, keys)[use]
print("N = %i, median rverr %.2f, N(rverr<1) = %i" % (len(x), np.median(xe[:, 3]), np.sum(xe[:, 3] < 1)))

opts = {'maxiter': 100000, 'maxfev': 100000}
def fit2(x, e):
    sd = x.std(axis=0)
    th0 = gmm.pack_params(np.array([0.9]), np.zeros((2, 4)), np.array([0.1*sd, 10*sd]))
    return minimize(gmm.nll_flat, th0, args=(x, 2, e), method='Powell', options=opts)

# truth for these same stars: noiseless fit
r_true = fit2(xt, None)
f, m, s = gmm.sort_components(*gmm.unpack_params(r_true.x))
print("noiseless fit, same stars: f_c %.3f  thin %s  cocoon %s" % (1-f.sum(), np.round(s[0], 3), np.round(s[1], 2)))

r = fit2(x, xe)
f, m, s = gmm.sort_components(*gmm.unpack_params(r.x))
print("noisy fit w/ errors:       f_c %.3f  thin %s  cocoon %s   nll %.3f" % (1-f.sum(), np.round(s[0], 4), np.round(s[1], 2), r.fun))

# which unsorted component is thin in v_gsr, and where its ln sigma_vgsr lives in theta
_, _, s_un = gmm.unpack_params(r.x)
j = int(np.argmin(s_un[:, 3]))
idx = 1 + j*8 + 4 + 3

# 1D view: just the stars the fit calls thin
p_thin = gmm.membership_probability(x, f, m, s, x_err=xe)
th = p_thin > 0.5
print("1D mle_sigma on p_thin>0.5 stars (N=%i): noisy+errors %.3f   true values (no noise) std %.3f"
      % (th.sum(), gmm.mle_sigma(x[th, 3], xe[th, 3]), xt[th, 3].std()))

print("\nprofile: thin sigma_vgsr fixed, everything else refit")
best = r.fun
theta_free = np.delete(r.x, idx)
for sig in [1e-3, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5]:
    def obj(tf, sig=sig):
        return gmm.nll_flat(np.insert(tf, idx, np.log(sig)), x, 2, xe)
    rp = minimize(obj, theta_free, method='Powell', options=opts)
    theta_free = rp.x  # warm start the next grid point
    print("  sigma_thin,vgsr = %5.3f   delta nll = %7.3f" % (sig, rp.fun - best))
