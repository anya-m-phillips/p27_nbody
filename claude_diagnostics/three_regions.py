"""
a one-off diagnostic, not part of the pipeline; reads nothing (no data at all).

    python claude_diagnostics/three_regions.py      # ~10 s

the Nibauer+26 setup: one global sigma_binary + one sigma_stream per region of
GD-1, independent-epoch likelihood. region k's likelihood depends only on its
total s_k^2 = ss_k^2 + sb^2, so one direction stays degenerate (sb up, every ss_k
down in quadrature), capped at sb = min_k s_k. near the cap the coldest region's
arc runs parallel to the sb axis, so the posterior piles up toward it.

uses the 'expected' (Asimov) likelihood for N stars per region with equal errors e,
so there's no draw noise and no binary information:
  lnL_k(s) = -N/2 [ ln(s^2 + e^2) + (s_k^2 + e^2)/(s^2 + e^2) ]   (one epoch per star; extra epochs ~ bigger N)
totals (3.1, 5.5, 4.5) = sqrt(2.6^2 + {1.7, 4.8, 3.7}^2), i.e. what the Nibauer+26
medians imply. N = 50 and e = 2 km/s are guesses.

result: reproduces their sigma_binary = 2.6 +1.2/-1.4 and region dispersions
1.7 +-1.2, 4.8 +0.7/-1.2, 3.7 +1.2/-1.8 km/s (medians by construction, error bars
not), unchanged for tau = 30; going to N = 150 tightens the hot regions but not
sigma_binary or the coldest region (prior-dominated).
"""
import numpy as np

def run(s_tot, N, e=2.0, tau=10.0, label=''):
    sb = np.linspace(0, 10, 2001)
    ss = np.linspace(0, 12, 2401)
    SB, SS = np.meshgrid(sb, ss, indexing='ij')
    lnprior = lambda x: -0.5 * (x / tau)**2          # half-normal, same tau on everything
    lnp_sb = lnprior(sb).copy()
    joint = []
    for s_k in s_tot:
        S2 = SS**2 + SB**2
        lnL = -0.5 * N * (np.log(S2 + e**2) + (s_k**2 + e**2) / (S2 + e**2))
        lp = lnL + lnprior(SS)                       # (sb, ss_k)
        m = lp.max()
        lnp_sb += np.log(np.exp(lp - m).sum(1)) + m  # marginalize ss_k at each sb
        joint.append(lp)
    p_sb = np.exp(lnp_sb - lnp_sb.max()); p_sb /= p_sb.sum()

    def q(x, p):
        c = np.cumsum(p); return [x[np.searchsorted(c, t)] for t in (0.16, 0.5, 0.84)]
    def fmt(v): return f'{v[1]:.2f} +{v[2]-v[1]:.2f} / -{v[1]-v[0]:.2f}'
    print(f'--- {label}  totals {s_tot}, N/region {N}, err {e}, tau {tau}')
    print(f'  sigma_binary        {fmt(q(sb, p_sb))}')
    for k, lp in enumerate(joint):
        # p(ss_k) = sum_sb p(sb) p(ss_k | sb)
        cond = np.exp(lp - lp.max(1, keepdims=True)); cond /= cond.sum(1, keepdims=True)
        p_ss = (p_sb[:, None] * cond).sum(0)
        print(f'  sigma_stream[{k}]     {fmt(q(ss, p_ss))}   (total {s_tot[k]})')

totals = (3.1, 5.5, 4.5)   # sqrt(2.6^2 + {1.7, 4.8, 3.7}^2), i.e. what their medians imply
run(totals, N=50, label='three regions, N=50 each')
run(totals, N=50, tau=30.0, label='same, tau = 30')
run(totals, N=150, label='three regions, N=150 each')
run((3.1,), N=50, label='region 1 alone (one sigma_stream, one sigma_binary)')
