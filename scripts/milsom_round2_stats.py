#!/usr/bin/env python3
"""Round 2 (exploratory, Amendment 2): HSCs compared between Milsom clones.
Cell-level standardised differences for within-run pairs, and genome-wide pseudobulk contrasts."""
import pathlib, numpy as np, pandas as pd
from scipy import stats
root = pathlib.Path(__file__).resolve().parents[1]
out = root / "results/round2"; out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261116); NBOOT = 2000
c = pd.read_csv(root / "data/processed/milsom/hspc_cell_scores.tsv.gz", sep="\t")
c = c[c.Celltypes_Louvain == "HSCs"].merge(pd.read_csv(root / "data/processed/milsom/round2_hsc_cell_scores.tsv.gz", sep="\t"), on="cell")
c["log10_UMI"] = np.log10(c.nCount_RNA)
SETS = ["ISR_stress_specific", "ISR_stress_specific_noIEG", "ATF4_biosynth", "Growth_noATF4", "RiBi_nonRP", "RP_cytoplasmic",
        "ISR_adaptive", "ISR_terminal", "mTORC1_hallmark", "Stem_MolO", "Stem_dHSC", "IEG_dissociation", "log10_UMI"]
PAIRS = [("M96", "M32", "clone_pair_run2"), ("M61", "M73", "clone_pair_run1"), ("M61", "M51", "clone_pair_run1"),
         ("M73", "M51", "clone_pair_run1"), ("M96", "M16", "vs_control_run2"), ("M32", "M16", "vs_control_run2"),
         ("M61", "M14", "vs_control_run1"), ("M73", "M14", "vs_control_run1"), ("M51", "M14", "vs_control_run1")]
def d(a, b):
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp
rows = []
for s1, s2, kind in PAIRS:
    for s in SETS:
        a, b = c.loc[c.Mice == s1, s].to_numpy(), c.loc[c.Mice == s2, s].to_numpy()
        bs = [d(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(NBOOT)]
        rows.append(dict(first=s1, second=s2, kind=kind, n_first=len(a), n_second=len(b), set=s, d_first_minus_second=d(a, b),
                         ci_lo=np.percentile(bs, 2.5), ci_hi=np.percentile(bs, 97.5)))
res = pd.DataFrame(rows); res.to_csv(out / "milsom_hsc_clone_pairs.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
res["txt"] = res.apply(lambda r: f"{r.d_first_minus_second:+.2f} [{r.ci_lo:+.2f},{r.ci_hi:+.2f}]", axis=1)
res["pair"] = res["first"] + "-" + res["second"]
print(res.pivot(index="set", columns="pair", values="txt").loc[SETS, [f"{a}-{b}" for a, b, _ in PAIRS]].to_string())

# genome-wide pseudobulk contrasts
pb = pd.read_csv(root / "data/processed/milsom/pseudobulk_counts_system_by_celltype.tsv.gz", sep="\t", index_col=0)
h = pd.DataFrame({m: pb[f"{m}|HSCs"] for m in ["M96", "M32", "M61", "M73", "M51", "M16", "M14"]})
h["M73_M51"] = h.M73 + h.M51
lc = np.log2(h / h.sum() * 1e6 + 1)
ok = (h[["M96", "M32", "M61", "M73_M51"]] >= 10).all(axis=1)
g = pd.DataFrame({"A_M96_vs_M32": lc.M96 - lc.M32, "B_M61_vs_M73M51": lc.M61 - lc.M73_M51, "D_M96_vs_M16": lc.M96 - lc.M16,
                  "M61_vs_M14": lc.M61 - lc.M14, "M32_vs_M16": lc.M32 - lc.M16, "pass_filter": ok})
g.to_csv(out / "milsom_hsc_pseudobulk_contrasts.tsv.gz", sep="\t", float_format="%.4g")
gf = g[ok]
print(f"\ngenes passing filter: {ok.sum()}")
for x, y in [("A_M96_vs_M32", "B_M61_vs_M73M51"), ("D_M96_vs_M16", "M61_vs_M14"), ("D_M96_vs_M16", "B_M61_vs_M73M51"), ("A_M96_vs_M32", "M61_vs_M14")]:
    print(f"spearman {x} ~ {y}: {stats.spearmanr(gf[x], gf[y])[0]:+.3f}")
sets = pd.concat([pd.read_csv(root / "data/processed/milsom/gene_sets_resolved.tsv", sep="\t"),
                  pd.read_csv(root / "metadata/gene_sets_round2.tsv", sep="\t")[["set", "gene"]]])
rk = gf[["A_M96_vs_M32", "B_M61_vs_M73M51", "D_M96_vs_M16"]].rank(pct=True)
sr = []
for s, dd in sets.groupby("set"):
    m = gf.index.intersection(dd.gene)
    if len(m) < 5: continue
    sr.append(dict(set=s, n_genes=len(m), **{f"mean_log2fc_{k}": gf.loc[m, k].mean() for k in rk.columns},
                   **{f"mean_pctrank_{k}": rk.loc[m, k].mean() for k in rk.columns}))
sr = pd.DataFrame(sr); sr.to_csv(out / "milsom_hsc_set_foldchanges.tsv", sep="\t", index=False, float_format="%.4g")
print(sr.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
both = gf[(np.sign(gf.A_M96_vs_M32) == np.sign(gf.B_M61_vs_M73M51))].copy()
both["minabs"] = both[["A_M96_vs_M32", "B_M61_vs_M73M51"]].abs().min(axis=1) * np.sign(both.A_M96_vs_M32)
print("\nhigher in slower clone in both contrasts:", ", ".join(both.sort_values("minabs").index[-40:][::-1]))
print("lower in slower clone in both contrasts:", ", ".join(both.sort_values("minabs").index[:40]))
