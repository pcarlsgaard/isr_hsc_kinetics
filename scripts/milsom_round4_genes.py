#!/usr/bin/env python3
"""Amendment 4 sensitivity (exploratory): which genes carry the HSC-versus-MPP difference in the ATF4 route sets?
Pseudobulk raw counts per mouse, log2 CPM+1 HSC minus MPP, genes with >=10 counts in HSC+MPP in every mouse."""
import pathlib, numpy as np, pandas as pd
root = pathlib.Path(__file__).resolve().parents[1]
pb = pd.read_csv(root / "data/processed/milsom/pseudobulk_counts_system_by_celltype.tsv.gz", sep="\t", index_col=0)
M = ["M96", "M61", "M32", "M73", "M51", "M16", "M14"]
lc = np.log2(pb / pb.sum() * 1e6 + 1)
fc = pd.DataFrame({m: lc[f"{m}|HSCs"] - lc[f"{m}|MPPs"] for m in M})
ok = pd.concat([(pb[f"{m}|HSCs"] + pb[f"{m}|MPPs"]) >= 10 for m in M], axis=1).all(axis=1)
fc = fc[ok]; fc["median"] = fc[M].median(axis=1); fc["n_up"] = (fc[M] > 0).sum(axis=1); fc["pct"] = fc["median"].rank(pct=True)
g3 = pd.read_csv(root / "metadata/gene_sets_round3.tsv", sep="\t"); g1 = pd.read_csv(root / "data/processed/gene_sets.tsv", sep="\t")
S = lambda d, k: set(d[d.set == k].gene)
excl = S(g1, "IEG_dissociation") | S(g1, "Stem_MolO") | S(g1, "Stem_dHSC") | S(g1, "Stem_LowOutput")
rows = []
for k in ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4", "R_ATF4_ER_stress", "R_PERK_gene_expression"]:
    for lab, gs in [("all", S(g3, k)), ("minus IEG and stemness-set genes", S(g3, k) - excl)]:
        x = fc.loc[fc.index.intersection(list(gs))]
        rows.append(dict(set=k, version=lab, n_expressed=len(x), median_log2fc=x["median"].median(), frac_up_in_HSC=(x["median"] > 0).mean(),
                         frac_up_in_all7=(x.n_up == 7).mean(), mean_pctrank=x.pct.mean()))
r = pd.DataFrame(rows); out = root / "results/round4"; r.to_csv(out / "milsom_hsc_vs_mpp_set_genes.tsv", sep="\t", index=False, float_format="%.4g")
print(f"genes: {len(fc)}; background frac up in HSC {(fc['median'] > 0).mean():.2f}, up in all 7 {(fc.n_up == 7).mean():.2f}")
print(r.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
x = fc.loc[fc.index.intersection(list(S(g3, "T21_ATF4_stress_only")))].sort_values("median", ascending=False)
x.to_csv(out / "milsom_hsc_vs_mpp_stress_only_genes.tsv", sep="\t", float_format="%.3g")
print("stress-only genes highest in HSC:", ", ".join(f"{g}({v:+.1f}{'*' if g in excl else ''})" for g, v in x["median"].head(30).items()))
print("stress-only genes lowest in HSC:", ", ".join(f"{g}({v:+.1f})" for g, v in x["median"].tail(12).items()))
y = fc.loc[fc.index.intersection(list(S(g3, "T21_ATF4_shared")))].sort_values("median", ascending=False)
print("shared ATF4 genes:", ", ".join(f"{g}({v:+.1f})" for g, v in y["median"].items()))
