#!/usr/bin/env python3
"""LARRY clone-level tests of the pre-specified contrasts (see docs/prespecification.md).
Stratified by experiment; permutation tests permute clone labels within experiment."""
import pathlib, sys
import numpy as np, pandas as pd
from scipy import stats
root = pathlib.Path(__file__).resolve().parents[1]
d = pd.read_csv(root / "data/processed/larry/larry_clone_scores.tsv", sep="\t")
gs = pd.read_csv(root / "data/processed/gene_sets.tsv", sep="\t")
tier = gs.drop_duplicates("set").set_index("set").tier
out = root / "results/larry"; out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261116)
NPERM = 20000
PRIMARY = ["ISR_adaptive", "ISR_terminal", "RiBi_nonRP", "RP_cytoplasmic", "mTORC1_hallmark"]

def wz(x):  # within-experiment z
    return (x - x.mean()) / x.std(ddof=1)

def strat_spearman(df, a, b, covars=()):
    """Within-experiment rank transform, optional residualisation on covariates, pooled correlation,
    permutation p (permute b within experiment)."""
    parts = []
    for s, g in df.groupby("sample"):
        if len(g) < 5: continue
        ra, rb = stats.rankdata(g[a]), stats.rankdata(g[b])
        if covars:
            C = np.column_stack([np.ones(len(g))] + [stats.rankdata(g[c]) for c in covars])
            ra = ra - C @ np.linalg.lstsq(C, ra, rcond=None)[0]
            rb = rb - C @ np.linalg.lstsq(C, rb, rcond=None)[0]
        parts.append(((ra - ra.mean()) / ra.std(), (rb - rb.mean()) / rb.std()))
    if not parts: return np.nan, np.nan, 0
    obs = np.mean(np.concatenate([x * y for x, y in parts]))
    null = np.empty(NPERM)
    for k in range(NPERM):
        null[k] = np.mean(np.concatenate([x * rng.permutation(y) for x, y in parts]))
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1), sum(len(x) for x, _ in parts)

def strat_group(df, a, grp):
    """Standardised difference (group 1 minus group 0) of within-experiment z-scores, permutation p."""
    z = df.groupby("sample")[a].transform(wz).to_numpy(); g = df[grp].to_numpy(); s = df["sample"].to_numpy()
    obs = z[g == 1].mean() - z[g == 0].mean()
    null = np.empty(NPERM)
    for k in range(NPERM):
        gp = g.copy()
        for ss in np.unique(s):
            m = s == ss; gp[m] = rng.permutation(g[m])
        null[k] = z[gp == 1].mean() - z[gp == 0].mean()
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)

rows, tert = [], []
for minn in (5, 3, 10):
    df = d[(d.n_hsc_scored >= minn) & d.relout.notna()].copy()
    df["log_relout"] = np.log2(df.relout + 0.05)
    lab = df[(df.inactive == 1) | (df.active == 1)].copy()
    lab["low_output"] = (lab.inactive == 1).astype(int)
    for s in sorted(tier.index):
        if s not in df or df[s].isna().all(): continue
        r, p, n = strat_spearman(df, s, "relout")
        r2, p2, _ = strat_spearman(df, s, "relout", covars=("Stem_MolO", "Stem_dHSC")) if s not in ("Stem_MolO", "Stem_dHSC") else (np.nan, np.nan, 0)
        r3, p3, _ = strat_spearman(df, s, "freqHSC")
        r4, p4, _ = strat_spearman(df, s, "Mkbias")
        dd, pd_ = strat_group(lab, s, "low_output")
        rows.append(dict(min_hsc_cells=minn, set=s, tier=tier[s], n_clones=n, rho_relout=r, p_relout=p,
                         rho_relout_adj_stemness=r2, p_relout_adj_stemness=p2, rho_freqHSC=r3, p_freqHSC=p3,
                         rho_Mkbias=r4, p_Mkbias=p4, n_low=int(lab.low_output.sum()), n_high=int((1 - lab.low_output).sum()),
                         d_low_minus_high=dd, p_low_vs_high=pd_))
        if minn == 5 and s in PRIMARY:
            t = df.copy()
            t["tertile"] = t.groupby("sample")[s].transform(lambda x: pd.qcut(x.rank(method="first"), 3, labels=[1, 2, 3]) if len(x) >= 6 else np.nan)
            t = t.dropna(subset=["tertile"])
            t["relout_rank"] = t.groupby("sample")["relout"].rank(pct=True)
            m = t.groupby("tertile", observed=True).agg(n=("relout", "size"), median_relout=("relout", "median"),
                                         mean_relout_pctrank=("relout_rank", "mean"), frac_low_output=("inactive", "mean")).reset_index()
            m.insert(0, "set", s)
            kw = stats.kruskal(*[t.loc[t.tertile == k, "relout_rank"] for k in (1, 2, 3)])
            mid = t.loc[t.tertile == 2, "relout_rank"].mean() - t.loc[t.tertile != 2, "relout_rank"].mean()
            m["kruskal_p"] = kw.pvalue; m["middle_minus_extremes_pctrank"] = mid
            tert.append(m)
res = pd.DataFrame(rows)
for c in ("p_relout", "p_low_vs_high"):
    res[c.replace("p_", "bh_primary_")] = np.nan
    for minn in res.min_hsc_cells.unique():
        m = (res.min_hsc_cells == minn) & res.set.isin(PRIMARY)
        res.loc[m, c.replace("p_", "bh_primary_")] = stats.false_discovery_control(res.loc[m, c])
res.to_csv(out / "larry_clone_associations.tsv", sep="\t", index=False, float_format="%.4g")
pd.concat(tert).to_csv(out / "larry_tertiles_primary.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
show = res[res.min_hsc_cells == 5].sort_values(["tier", "set"])
print(show[["set", "tier", "n_clones", "rho_relout", "p_relout", "bh_primary_relout", "rho_relout_adj_stemness", "p_relout_adj_stemness", "rho_freqHSC", "p_freqHSC", "n_low", "n_high", "d_low_minus_high", "p_low_vs_high", "bh_primary_low_vs_high"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print(pd.concat(tert).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print(d.groupby("sample").agg(clones=("clone", "size"), ge5=("n_hsc_scored", lambda x: (x >= 5).sum()), median_cells=("n_hsc_scored", "median")))
