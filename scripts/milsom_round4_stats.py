#!/usr/bin/env python3
"""Amendment 4 (exploratory): HSC versus MPP (and versus all other progenitors) within each mouse; mice are the units."""
import pathlib, numpy as np, pandas as pd
from scipy import stats
root = pathlib.Path(__file__).resolve().parents[1]
out = root / "results/round4"; out.mkdir(parents=True, exist_ok=True)
old = pd.read_csv(root / "data/processed/milsom/hspc_cell_scores.tsv.gz", sep="\t")
c = old[["cell", "Mice", "Celltypes_Louvain", "ISR_adaptive", "ISR_terminal", "Stem_MolO", "IEG_dissociation", "RiBi_nonRP", "RP_cytoplasmic"]].merge(
    pd.read_csv(root / "data/processed/milsom/round4_hspc_cell_scores.tsv.gz", sep="\t"), on="cell")
SETS = [x for x in c.columns if x not in ("cell", "Mice", "Celltypes_Louvain")]
ORDER = ["M96", "M61", "M32", "M73", "M51", "M16", "M14"]
def d(a, b):
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp
res = []
for ref, lab in [("MPPs", "HSC_minus_MPP"), (None, "HSC_minus_all_other")]:
    for s in SETS:
        v = {}
        for m in ORDER:
            g = c[c.Mice == m]; a = g.loc[g.Celltypes_Louvain == "HSCs", s].to_numpy()
            b = g.loc[(g.Celltypes_Louvain == ref) if ref else (g.Celltypes_Louvain != "HSCs"), s].to_numpy()
            v[m] = d(a, b)
        x = np.array(list(v.values()))
        res.append(dict(contrast=lab, set=s, **v, median=np.median(x), n_positive=int((x > 0).sum()), wilcoxon_p=stats.wilcoxon(x).pvalue))
r = pd.DataFrame(res); r.to_csv(out / "milsom_hsc_vs_progenitors.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(c.groupby(["Mice", "Celltypes_Louvain"]).size().unstack().loc[ORDER].to_string())
for lab, g in r.groupby("contrast", sort=False):
    print("\n" + lab); print(g.drop(columns="contrast").to_string(index=False, float_format=lambda v: f"{v:+.2f}" if abs(v) >= 0.0005 else f"{v:.3f}"))
