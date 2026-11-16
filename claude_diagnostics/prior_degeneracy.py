"""
a one-off diagnostic, not part of the pipeline; reads nothing (synthetic data).

    python claude_diagnostics/prior_degeneracy.py      # a few seconds

with the independent-epoch likelihood (Nibauer+26 as written),
    v_ik ~ N(mu, err_ik^2 + sigma_binary^2 + sigma_stream^2),
the data only constrain the total s^2 = sigma_stream^2 + sigma_binary^2, so the
high-likelihood region in the (sigma_stream, sigma_binary) plane is a thin
quarter-circle of radius s. the split along the arc is set by the prior alone.
this computes the grid posterior on synthetic two-epoch data (N = 679, total
4.66 km/s) under several priors and prints the 16/50/84% of each marginal.

result: equal half-normal scales (10/10 or 30/30), unequal scales (30/10) and a
flat prior all give the same uniform-in-angle answer (median s/sqrt2, 16-84%
~[0.25 s, 0.97 s], identical for both sigmas). only a prior whose *shape* matters
on the scale of s (log-uniform, or a 1 km/s half-normal) changes the result.
"""
import numpy as np
from scipy.stats import norm

# synthetic two-epoch data, independent-epoch likelihood (Nibauer+26 form):
# v_ik ~ N(mu, err_ik^2 + sb^2 + ss^2), every epoch independent
rng = np.random.default_rng(0)
N, s_true = 679, 4.66
err = rng.uniform(1, 6, size=(N, 2))
v = rng.normal(0, np.sqrt(s_true**2 + err**2))

g = np.linspace(1e-3, 12, 1200)
SS, SB = np.meshgrid(g, g, indexing='ij')
# likelihood depends only on s^2 = ss^2 + sb^2: evaluate on a 1D grid in s, interpolate
s1 = np.linspace(0, 17, 4000)
lnL1 = np.array([norm.logpdf(v, 0, np.sqrt(s**2 + err**2)).sum() for s in s1])
lnL = np.interp(np.hypot(SS, SB), s1, lnL1)

def summarize(name, lnprior):
    lp = lnL + lnprior
    p = np.exp(lp - lp.max()); p /= p.sum()
    for lab, m in [('sigma_stream', p.sum(1)), ('sigma_binary', p.sum(0))]:
        c = np.cumsum(m)
        q = [g[np.searchsorted(c, x)] for x in (0.16, 0.5, 0.84)]
        print(f'{name:28s} {lab:13s} median {q[1]:.2f}  16-84% [{q[0]:.2f}, {q[2]:.2f}]')

hn = lambda x, tau: -0.5 * (x / tau)**2
summarize('half-normal 10 / 10', hn(SS, 10) + hn(SB, 10))
summarize('half-normal 30 / 30', hn(SS, 30) + hn(SB, 30))
summarize('half-normal ss 30 / sb 10', hn(SS, 30) + hn(SB, 10))
summarize('flat', 0 * SS)
summarize('log-uniform (1/sigma)', -np.log(SS) - np.log(SB))
summarize('half-normal ss 1 / sb 10', hn(SS, 1) + hn(SB, 10))
s_hat = s1[np.argmax(lnL1)]
print(f'\nbest-fit total width {s_hat:.2f}; s/sqrt2 = {s_hat/np.sqrt(2):.2f}; '
      f'uniform-angle 16-84%: [{s_hat*np.cos(0.84*np.pi/2):.2f}, {s_hat*np.cos(0.16*np.pi/2):.2f}]')
