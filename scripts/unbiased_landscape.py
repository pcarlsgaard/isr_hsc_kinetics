#!/usr/bin/env python3
"""Amendment 5 (exploratory): gene-level statistics for every tested gene and a score for every MSigDB set,
for HSC vs MPP (Milsom, 7 mice, paired) and low- vs high-output clones (LARRY)."""
import pathlib, numpy as np, pandas as pd
from scipy import stats
root = pathlib.Path(__file__).resolve().parents[1]
out = root / "results/round5"; out.mkdir(parents=True, exist_ok=True)
# gene sets
sets = {}
for f, pre in [("mh.all", "H"), ("m2.cp.reactome", "R"), ("m5.go.bp", "GO")]:
    for l in open(root / f"data/raw/msigdb/{f}.v2025.1.Mm.symbols.gmt"):
        p = l.rstrip("\n").split("\t"); sets[p[0]] = set(p[2:])
g3 = pd.read_csv(root / "metadata/gene_sets_round3.tsv", sep="\t")
T21 = ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4"]
for k in T21: sets[k] = set(g3[g3.set == k].gene)
def landscape(t, name):
    z = (t - t.mean()) / t.std(); rows = []
    for k, g in sets.items():
        m = z.index.intersection(list(g))
        if 15 <= len(m) <= 300: rows.append(dict(set=k, n=len(m), score=z[m].mean()))
    r = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    r["rank_from_top"] = np.arange(1, len(r) + 1); r["frac_sets_more_extreme"] = [(r.score.abs() > abs(s)).mean() for s in r.score]
    r.to_csv(out / f"{name}_set_landscape.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"\n{name}: {len(r)} sets scored, {len(z)} genes")
    print(r[r.set.isin(T21)].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("top:", "; ".join(f"{s} ({v:+.2f})" for s, v in zip(r.set[:8], r.score[:8])))
    print("bottom:", "; ".join(f"{s} ({v:+.2f})" for s, v in zip(r.set[-8:], r.score[-8:])))
    h = r[r.set.str.startswith("HALLMARK") | r.set.isin(T21)].copy(); h["rank_in_hallmark_panel"] = np.arange(1, len(h) + 1)
    print(h[h.set.isin(T21)][["set", "score", "rank_in_hallmark_panel"]].to_string(index=False), f"(of {len(h)})")
# (a) HSC vs MPP
pb = pd.read_csv(root / "data/processed/milsom/pseudobulk_counts_system_by_celltype.tsv.gz", sep="\t", index_col=0)
M = ["M96", "M61", "M32", "M73", "M51", "M16", "M14"]; lc = np.log2(pb / pb.sum() * 1e6 + 1)
fc = pd.DataFrame({m: lc[f"{m}|HSCs"] - lc[f"{m}|MPPs"] for m in M})
ok = pd.concat([(pb[f"{m}|HSCs"] + pb[f"{m}|MPPs"]) >= 10 for m in M], axis=1).all(axis=1); fc = fc[ok]
tt = stats.ttest_1samp(fc[M], 0, axis=1)
a = pd.DataFrame({"log2fc": fc[M].mean(axis=1), "t": tt.statistic, "p": tt.pvalue}, index=fc.index).dropna()
a.to_csv(out / "hsc_vs_mpp_genes.tsv.gz", sep="\t", float_format="%.4g"); landscape(a.t, "hsc_vs_mpp")
# (b) LARRY
gm = pd.read_csv(root / "data/processed/larry/round2_clone_gene_means.tsv.gz", sep="\t")
cs = pd.read_csv(root / "data/processed/larry/larry_clone_scores.tsv", sep="\t")
cs = cs[(cs.n_hsc_scored >= 5) & cs.relout.notna() & ((cs.inactive == 1) | (cs.active == 1))]; cs["low"] = (cs.inactive == 1).astype(int)
gm = gm.merge(cs[["sample", "clone", "low"]], on=["sample", "clone"])
gm["z"] = gm.groupby(["sample", "gene"]).mean_log1p.transform(lambda x: (x - x.mean()) / x.std(ddof=1) if x.std(ddof=1) > 0 else np.nan)
gm["key"] = gm["sample"].astype(str) + "|" + gm.clone
Z = gm.pivot(index="gene", columns="key", values="z"); mean_expr = gm.groupby("gene").mean_log1p.mean()
nexp = gm.dropna(subset=["z"]).groupby("gene")["sample"].nunique()
keep = (nexp.reindex(Z.index) == 3) & (mean_expr.reindex(Z.index) >= mean_expr.quantile(0.5)) & Z.notna().all(axis=1)
Z = Z[keep]; lowk = set(cs.loc[cs.low == 1, "sample"].astype(str) + "|" + cs.loc[cs.low == 1, "clone"])
lo = Z[[c for c in Z.columns if c in lowk]]; hi = Z[[c for c in Z.columns if c not in lowk]]
tt = stats.ttest_ind(lo, hi, axis=1, equal_var=False)
b = pd.DataFrame({"diff": lo.mean(axis=1) - hi.mean(axis=1), "t": tt.statistic, "p": tt.pvalue}, index=Z.index)
b.to_csv(out / "larry_low_vs_high_genes.tsv.gz", sep="\t", float_format="%.4g"); print(f"\nLARRY clones low {lo.shape[1]}, high {hi.shape[1]}")
landscape(b.t, "larry_low_vs_high")
