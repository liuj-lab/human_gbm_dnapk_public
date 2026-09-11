# Bootstrap analysis of h matrix to infer confidence in perturbation - program edges #
#
# Part 1  refit_h / bootstrap_h            -- empirical standard error of the response
# Part 2  response_se / compare_conditions -- z statistics, P values, sign-flip calls
#
# Significance model
# ------------------
# D-SPIN reports a *relative* response for each sgRNA within a condition:
#       r_g,c = h_g,c - mean(h_NTC,c)
# so its sampling variance has two independent contributions, the perturbation
# arm and the non-targeting arm of that same condition:
#       Var(r_g,c) = Var(h_g,c) + Var(h_NTC,c)
# The RT-vs-control contrast subtracts two such quantities, estimated from four
# disjoint sets of cells and therefore mutually independent:
#       D_g      = r_g,RT - r_g,control
#       Var(D_g) = Var(h_g,RT) + Var(h_NTC,RT) + Var(h_g,ctrl) + Var(h_NTC,ctrl)
# Each variance term is the bootstrap variance from Part 1. P values use the
# normal approximation z = D / SE(D) rather than the empirical bootstrap
# percentile, because with B resamples a percentile P value cannot fall below
# 1/B (0.005 at B = 200) and would floor the strongest contrasts.

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from dspin.compute import pseudol_gradient


# --------------------------------------------------------------------------
# Part 1 -- bootstrap standard error of the field vector h
# --------------------------------------------------------------------------

def refit_h(J, S_cells, h_init=None, n_iter=300, lr=0.5, tol=1e-7):
    """
    Refit the field vector h for one sample's cells, with J fixed.
    S_cells : (n_cells, num_spin) discretized states.
    Returns h of shape (num_spin,).
    """
    state = np.ascontiguousarray(S_cells.T)          # (num_spin, n_cells)
    num_spin = state.shape[0]
    h = np.zeros((num_spin, 1)) if h_init is None else h_init.reshape(-1, 1).copy()

    m, v = np.zeros_like(h), np.zeros_like(h)
    b1, b2, eps = 0.9, 0.999, 1e-8
    prev = h.copy()
    for t in range(1, n_iter + 1):
        _, hgrad = pseudol_gradient(J, h, state)
        m = b1 * m + (1 - b1) * hgrad
        v = b2 * v + (1 - b2) * hgrad ** 2
        mhat = m / (1 - b1 ** t)
        vhat = v / (1 - b2 ** t)
        h = h - lr * mhat / (np.sqrt(vhat) + eps)
        if t % 25 == 0:
            if np.max(np.abs(h - prev)) < tol:
                break
            prev = h.copy()
    return h.ravel()


def bootstrap_h(J, S_cells, n_boot=200, seed=0, h_init=None, **kw):
    """Resample cells with replacement, refit h each time. Returns (n_boot, num_spin)."""
    rng = np.random.default_rng(seed)
    n = S_cells.shape[0]
    out = np.zeros((n_boot, S_cells.shape[1]))
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        out[b] = refit_h(J, S_cells[idx], h_init=h_init, **kw)
    return out


def bootstrap_se(boot):
    """Standard error of h: SD across bootstrap replicates. boot is (n_boot, num_spin)."""
    return boot.std(axis=0, ddof=1)


# --------------------------------------------------------------------------
# Part 2 -- z statistics, P values, sign-flip calls
# --------------------------------------------------------------------------

def response_se(se_gene, se_ntc):
    """
    SE of a *relative* response r = h_gene - mean(h_NTC) within one condition.
    Both arms are bootstrapped from disjoint cells, so their variances add.
    """
    return np.sqrt(np.asarray(se_gene) ** 2 + np.asarray(se_ntc) ** 2)


def response_z(effect, se_gene, se_ntc):
    """
    Per-condition confidence in a single response (used for the z heat map).
    Returns (z, two-sided P).
    """
    z = np.asarray(effect) / response_se(se_gene, se_ntc)
    return z, 2 * norm.sf(np.abs(z))


def compare_conditions(df, control_gene='non-targeting',
                       cond_a='RT', cond_b='noRT', include_ntc_variance=True):
    """
    RT-vs-control contrast for every sgRNA x gene-program pair, with P values.

    Parameters
    ----------
    df : DataFrame with columns
         condition, gene, program_idx, effect, bootstrap_std
         ('effect' is the D-SPIN relative response; 'bootstrap_std' is the
          bootstrap SE of h from Part 1).
    include_ntc_variance : if True, propagate the non-targeting arm's bootstrap
         variance into the contrast as well (recommended; see header).

    Returns
    -------
    DataFrame with delta, se_delta, z, pval, fdr, sign_flip and the symbol
    class used in the Extended Data figure ('+', '-', 'flip').
    """
    ntc = (df[df.gene == control_gene]
           .set_index(['condition', 'program_idx'])['bootstrap_std'])

    w = df[df.gene != control_gene].pivot_table(
        index=['gene', 'program_idx'], columns='condition',
        values=['effect', 'bootstrap_std']).reset_index()
    w.columns = ['_'.join([c for c in col if c]).strip() for col in w.columns]
    w = w.rename(columns={'gene_': 'gene', 'program_idx_': 'program_idx'})

    ea, eb = w[f'effect_{cond_a}'], w[f'effect_{cond_b}']
    sa, sb = w[f'bootstrap_std_{cond_a}'], w[f'bootstrap_std_{cond_b}']

    if include_ntc_variance:
        na = w.program_idx.map(ntc.loc[cond_a])
        nb = w.program_idx.map(ntc.loc[cond_b])
        var = sa ** 2 + sb ** 2 + na ** 2 + nb ** 2
    else:
        var = sa ** 2 + sb ** 2

    w['delta'] = ea - eb
    w['se_delta'] = np.sqrt(var)
    w['z'] = w['delta'] / w['se_delta']
    w['pval'] = 2 * norm.sf(w['z'].abs())
    w['fdr'] = multipletests(w['pval'], method='fdr_bh')[1]

    w['sign_flip'] = np.sign(ea) != np.sign(eb)
    w['symbol'] = np.where(w['sign_flip'], 'flip', np.where(ea > 0, '+', '-'))
    return w


def summarise(w, alpha=0.05):
    """Counts for the figure legend / Results text."""
    return dict(
        n_pairs=len(w),
        n_sig_p=int((w.pval < alpha).sum()),
        n_sig_fdr=int((w.fdr < alpha).sum()),
        n_flips=int(w.sign_flip.sum()),
        n_flips_sig_p=int((w.sign_flip & (w.pval < alpha)).sum()),
        n_flips_sig_fdr=int((w.sign_flip & (w.fdr < alpha)).sum()),
    )


if __name__ == '__main__':
    df = pd.read_csv('bootstrap_vs_analytic.csv')
    w = compare_conditions(df)
    w.to_csv('signflip_significance.csv', index=False)
    for k, v in summarise(w).items():
        print(f'{k:20s} {v}')
