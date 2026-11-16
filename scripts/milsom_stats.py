#!/usr/bin/env python3
"""Milsom clonal systems: system-level summaries of pre-specified signature scores in annotated HSCs.
n = 5 clonal systems + 2 polyclonal controls. See docs/prespecification.md (incl. Amendment 1)."""
import itertools, pathlib
import numpy as np, pandas as pd
from scipy import stats
root = pathlib.Path(__file__).resolve().parents[1]
out = root / "results/milsom"; out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261116)
c = pd.read_csv(root / "data/processed/milsom/hspc_cell_scores.tsv.gz", sep="\t")
gs = pd.read_csv(root / "data/processed/gene_sets.tsv", sep="\t")
tier = gs.drop_duplicates("set").set_index("set").tier
SETS = sorted(tier.index)
PRIMARY = ["ISR_adaptive", "ISR_terminal", "RiBi_nonRP", "RP_cytoplasmic", "mTORC1_hallmark"]
ORDER = ["M96", "M61", "M32", "M73", "M51"]                      # authors' slow -> fast
RANK = {m: i + 1 for i, m in enumerate(ORDER)}
RUNCTRL = {"M61": "M14", "M73": "M14", "M51": "M14", "M96": "M16", "M32": "M16"}
PAIRS = [("M61", "M73"), ("M61", "M51"), ("M73", "M51"), ("M96", "M32")]   # (slower, faster) within run
EXT = pd.DataFrame({"system": ORDER, "kinetic_rank": [1, 2, 3, 4, 5],
                    "tyMax_wk": [19.77, 19.44, 18.34, 16.37, 10.35],
                    "HSC_share_primary_pct": [98.84, 59.86, 3.80, 0.88, 0.0],
                    "HSC_share_secondary_pct": [69.81, 49.09, 0.41, 0.0, 0.0],
                    "PB_AUC": [1.171, 3.795, 7.698, 3.961, 1.985]}).set_index("system")
EXT.to_csv(out / "system_external_variables.tsv", sep="\t")

def perm_spearman(x, y):
    obs = stats.spearmanr(x, y)[0]
    null = [stats.spearmanr(x, p)[0] for p in itertools.permutations(y)]
    return obs, float(np.mean(np.abs(null) >= abs(obs) - 1e-12))

def summarise(h, label):
    rows = []
    sd_within = h[h.Mice.isin(ORDER)].groupby("Mice")[SETS].std(ddof=1)
    nn = h[h.Mice.isin(ORDER)].groupby("Mice").size()
    pooled = np.sqrt(((sd_within ** 2).mul(nn - 1, axis=0)).sum() / (nn.sum() - len(nn)))
    mean = h.groupby("Mice")[SETS].mean()
    for s in SETS:
        m = mean[s]
        v = m[ORDER].to_numpy()
        rho, p = perm_spearman(np.arange(1, 6), v)
        loo = [stats.spearmanr(np.delete(np.arange(1, 6), i), np.delete(v, i))[0] for i in range(5)]
        delta = np.array([m[k] - m[RUNCTRL[k]] for k in ORDER])
        rho_d, p_d = perm_spearman(np.arange(1, 6), delta)
        conc = sum(np.sign(m[a] - m[b]) for a, b in PAIRS)   # +4: slower system higher in all within-run pairs
        row = dict(population=label, set=s, tier=tier[s], **{f"mean_{k}": m[k] for k in ORDER + ["M14", "M16"]},
                   pooled_within_sd=pooled[s], rho_kinetic_rank=rho, p_exact=p, loo_min=min(loo), loo_max=max(loo),
                   rho_ctrl_referenced=rho_d, p_exact_ctrl_referenced=p_d,
                   within_run_pairs_slower_higher=int((conc + 4) / 2),
                   d_slow2_vs_rest=(v[:2].mean() - v[2:].mean()) / pooled[s],
                   d_slow2_vs_rest_ctrl_ref=(delta[:2].mean() - delta[2:].mean()) / pooled[s],
                   d_M96_vs_M16ctrl=(m["M96"] - m["M16"]) / pooled[s], d_M61_vs_M14ctrl=(m["M61"] - m["M14"]) / pooled[s],
                   d_ctrl_M16_vs_M14=(m["M16"] - m["M14"]) / pooled[s],
                   rho_tyMax=stats.spearmanr(EXT.loc[ORDER, "tyMax_wk"], v)[0],
                   rho_HSC_share=stats.spearmanr(EXT.loc[ORDER, "HSC_share_primary_pct"], v)[0],
                   rho_PB_AUC=stats.spearmanr(EXT.loc[ORDER, "PB_AUC"], v)[0])
        rows.append(row)
    return pd.DataFrame(rows)

hsc = c[c.Celltypes_Louvain == "HSCs"].copy()
res = [summarise(hsc, "HSCs"), summarise(c[c.Celltype == "HSCs end"], "HSCs_end_subset"),
       summarise(c[c.Celltypes_Louvain == "MPPs"], "MPPs")]
# stemness strata defined on pooled annotated HSCs
hsc["molo_stratum"] = pd.qcut(hsc.Stem_MolO, 3, labels=["low", "mid", "high"])
res.append(summarise(hsc[hsc.molo_stratum == "high"], "HSCs_MolO_top_tertile"))
# G1-only HSCs
res.append(summarise(hsc[hsc.Phase == "G1"], "HSCs_G1_only"))
res = pd.concat(res)
m = (res.population == "HSCs") & res.set.isin(PRIMARY)
res.loc[m, "bh_primary"] = stats.false_discovery_control(res.loc[m, "p_exact"])
res.to_csv(out / "milsom_system_level_results.tsv", sep="\t", index=False, float_format="%.4g")

# bootstrap CIs of system means (cells resampled within system; describes sampling of cells only)
boot = []
for k, g in hsc.groupby("Mice"):
    X = g[SETS].to_numpy(); idx = rng.integers(0, len(X), (2000, len(X)))
    bm = X[idx].mean(axis=1)
    for j, s in enumerate(SETS):
        boot.append(dict(system=k, set=s, n_cells=len(X), mean=X[:, j].mean(), lo=np.percentile(bm[:, j], 2.5), hi=np.percentile(bm[:, j], 97.5)))
pd.DataFrame(boot).to_csv(out / "milsom_hsc_system_means_bootstrap.tsv", sep="\t", index=False, float_format="%.5g")

# within-system cell-level correlations (not between-clone inference)
cov = ["Stem_MolO", "Stem_dHSC", "G2M.Score", "S.Score", "E2F_targets", "IEG_dissociation", "nCount_RNA", "percent.mt"]
w = []
for k, g in hsc.groupby("Mice"):
    for s in PRIMARY:
        for v in cov:
            w.append(dict(system=k, n=len(g), set=s, covariate=v, rho=stats.spearmanr(g[s], g[v])[0]))
w = pd.DataFrame(w); w.to_csv(out / "milsom_within_system_cell_correlations.tsv", sep="\t", index=False, float_format="%.3f")

# composition and counts
comp = c.groupby(["Mice", "Celltypes_Louvain"]).size().unstack(fill_value=0)
comp.to_csv(out / "milsom_cells_per_system_celltype.tsv", sep="\t")
qc = hsc.groupby("Mice").agg(n=("cell", "size"), nCount=("nCount_RNA", "median"), nFeature=("nFeature_RNA", "median"),
                             pct_mt=("percent.mt", "median"), pct_ribo=("percent.ribo", "median"),
                             frac_G1=("Phase", lambda x: (x == "G1").mean()), MolO=("Stem_MolO", "mean"))
qc.to_csv(out / "milsom_hsc_qc_by_system.tsv", sep="\t", float_format="%.4g")

pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
f = lambda v: f"{v:.3f}"
cols = ["set", "tier", "mean_M96", "mean_M61", "mean_M32", "mean_M73", "mean_M51", "mean_M14", "mean_M16", "pooled_within_sd",
        "rho_kinetic_rank", "p_exact", "loo_min", "loo_max", "rho_ctrl_referenced", "within_run_pairs_slower_higher",
        "d_slow2_vs_rest", "d_slow2_vs_rest_ctrl_ref", "d_M96_vs_M16ctrl", "d_M61_vs_M14ctrl", "d_ctrl_M16_vs_M14", "rho_PB_AUC"]
for pop in res.population.unique():
    print("\n#####", pop); print(res[res.population == pop].sort_values(["tier", "set"])[cols].to_string(index=False, float_format=f))
print("\n## HSC QC by system\n", qc.loc[ORDER + ["M14", "M16"]].to_string(float_format=f))
print("\n## within-system cell-level Spearman (annotated HSCs), primary sets\n",
      w.pivot_table(index=["set", "covariate"], columns="system", values="rho")[ORDER + ["M14", "M16"]].to_string(float_format=lambda v: f"{v:.2f}"))
