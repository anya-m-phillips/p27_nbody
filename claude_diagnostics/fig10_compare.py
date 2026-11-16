"""
a one-off diagnostic, not part of the pipeline; reads nothing (no data at all).

    python claude_diagnostics/fig10_compare.py      # a few seconds

three_regions.py with the actual Nibauer+26 bin sizes: 35 / 74 / 51 stars in the
phi1 bins [-60,-40], [-40,-20], [-20,0] deg, x 195/160 measurements per star
(160 stars, 195 RVs). one global sigma_binary, one sigma_stream per bin, their
Eq. 1 (independent per-survey measurements), expected (Asimov) likelihood with
equal errors e, i.e. NO binary information. half-normal priors, scale 10 km/s.
the region totals (3.1, 5.5, 4.5) = sqrt(2.6^2 + {1.7, 4.8, 3.7}^2), from their
medians. e is a guess between their MMT (0.8) and SDSS (6.4 km/s) medians, so it
loops over 2, 3, 4 km/s.

prints sigma_binary and each bin's sigma_stream (median +84 -16), plus each bin's
posterior shape on 0-8 km/s in 0.5 km/s bins (normalized to a peak of 1), to compare
by eye with their Fig. 10 (top panel, |dv| < 30 km/s).

result (e = 4): sigma_binary 2.76 +1.15/-1.60; bins 1.82 +1.38/-1.24 (flat from 0 to
~2.5, then falls), 4.75 +0.91/-1.23, 3.49 +1.16/-1.78 (tail down to 0). theirs:
2.6 +1.2/-1.4; 1.7 +-1.2, 4.8 +0.7/-1.2, 3.7 +1.2/-1.8, with the same shapes in Fig. 10.
"""
import numpy as np

# three bins, global sb, per-bin ss; expected likelihood, no binary information
Nmeas = np.array([35, 74, 51]) * 195/160
s_tot = (3.1, 5.5, 4.5)
for e in (2.0, 3.0, 4.0):
    sb = np.linspace(0, 10, 1001); ss = np.linspace(0, 10, 1001)
    SB, SS = np.meshgrid(sb, ss, indexing='ij')
    lnpr = lambda x: -0.5*(x/10)**2
    lnp_sb = lnpr(sb).copy(); joint = []
    for N, s in zip(Nmeas, s_tot):
        S2 = SS**2 + SB**2
        lp = -0.5*N*(np.log(S2 + e**2) + (s**2 + e**2)/(S2 + e**2)) + lnpr(SS)
        m = lp.max(); lnp_sb += np.log(np.exp(lp - m).sum(1)) + m; joint.append(lp)
    p_sb = np.exp(lnp_sb - lnp_sb.max()); p_sb /= p_sb.sum()
    q = lambda x, p: [x[np.searchsorted(np.cumsum(p), t)] for t in (.16, .5, .84)]
    v = q(sb, p_sb); print(f'e={e}: sigma_binary {v[1]:.2f} +{v[2]-v[1]:.2f} -{v[1]-v[0]:.2f}')
    edges = np.arange(0, 8.5, 0.5)
    for k, lp in enumerate(joint):
        c = np.exp(lp - lp.max(1, keepdims=True)); c /= c.sum(1, keepdims=True)
        p = (p_sb[:, None]*c).sum(0); v = q(ss, p)
        h = np.histogram(ss, edges, weights=p)[0]; h /= h.max()
        print(f'  bin {k}: {v[1]:.2f} +{v[2]-v[1]:.2f} -{v[1]-v[0]:.2f}  shape(0..8 by 0.5):', ' '.join(f'{x:.1f}' for x in h))
