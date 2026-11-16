#!/usr/bin/env python3
"""LARRY (GSE134242 / Tcf15_HSCs repo): per-cell rank-based signature scores in sorted HSCs,
aggregated to clone level and joined to the authors' clone table (relout, inactive/active, Mk bias).

Uses the authors' filtered normalised matrices (raw counts are not cell-aligned to the
annotations), so scoring is rank-based only (UCell statistic, maxRank=1500); no count-based DE.
Clone assignment follows the original notebook: exact match of the cell's barcode string
to the clone list.
"""
import csv, gzip, re, sys, pathlib
import numpy as np, pandas as pd
from scipy.stats import rankdata
root = pathlib.Path(__file__).resolve().parents[1]
repo = root / "data/raw/tcf15_github/repo"
mat = root / "data/raw/tcf15_dropbox/filtered_counts_csv"
out = root / "data/processed/larry"; out.mkdir(parents=True, exist_ok=True)
MAXRANK = 1500
RP = re.compile(r"^(Rp[ls][0-9]+[a-z]?|Rplp[0-2]|Rpsa)$")
gs = pd.read_csv(root / "data/processed/gene_sets.tsv", sep="\t")
S = {
 1: dict(m="ARF2019_Sample1_new_counts_norm.csv.gz", g="ARF2019_Sample1_genes.txt", bc="Sample1_new_lin_bcs.txt",
         sid="Sample1_new_samp_id.txt", hsc="Sample1_isHSC.txt", cl="Sample1_clonelist_DEFINITIVE.txt"),
 2: dict(m="ARF2019_Sample2_counts_norm.csv.gz", g="ARF2019_Sample2_genes.txt", bc="ARF2019_Sample2_barcode_list.txt",
         sid="ARF2019_Sample2_samp_id_flat.txt", hsc="Sample2_isHSC.txt", cl="Sample2_clonelist_DEFINITIVE.txt"),
 3: dict(m="ARF2019_Sample3_counts_norm.csv.gz", g="ARF2019_Sample3_genes.txt", bc="ARF2019_Sample3_barcode_list.txt",
         sid="ARF2019_Sample3_samp_id_flat.txt", hsc="Sample3_isHSC.txt", cl="Sample3_clonelist_DEFINITIVE.txt"),
}
def lines(p): return open(p).read().split("\n")
cells_all, cov_rows = [], []
for s, f in S.items():
    genes = [x for x in lines(repo / f["g"]) if x != ""]
    bcs, sid, hsc = lines(repo / f["bc"]), lines(repo / f["sid"]), lines(repo / f["hsc"])
    clones = pd.read_csv(repo / f["cl"], sep="\t").set_index("clone")
    gi = {g: i for i, g in enumerate(genes)}
    sets = {}
    for name, d in gs.groupby("set"):
        g = [x for x in genes if RP.match(x)] if name == "RP_cytoplasmic" else list(d.gene)
        idx = sorted({gi[x] for x in g if x in gi})
        sets[name] = np.array(idx, dtype=int)
        cov_rows.append(dict(sample=s, set=name, n_defined=len(set(g)), n_present=len(idx)))
    nlab = min(len(hsc), len(bcs))
    sel = [i for i in range(nlab) if hsc[i].strip() == "1" and bcs[i] in clones.index]
    selset = set(sel)
    rec, n = [], 0
    def score(i, x):
        r = rankdata(-x, method="average")
        r[r > MAXRANK] = MAXRANK + 1
        row = dict(sample=s, cell=i, library=sid[i], clone=bcs[i], total=float(x.sum()), ngenes=int((x > 0).sum()))
        for name, idx in sets.items():
            k = len(idx)
            row[name] = np.nan if k == 0 else 1 - (r[idx].sum() - k * (k + 1) / 2) / (k * MAXRANK)
        return row
    with gzip.open(mat / f["m"], "rt") as fh:
        if s == 1:  # rows are cells, no header
            for i, line in enumerate(fh):
                n += 1
                if i not in selset: continue
                x = np.array(line.rstrip("\n").split(","), dtype=np.float32)
                assert x.size == len(genes), (x.size, len(genes))
                rec.append(score(i, x))
        else:       # rows are genes (name in first field), columns are cells
            sub = np.zeros((len(genes), len(sel)), dtype=np.float32)
            cols = np.array(sel) + 1
            for gidx, line in enumerate(fh):
                p = line.rstrip("\n").split(",")
                assert p[0] == genes[gidx], (p[0], genes[gidx])
                n = len(p) - 1
                sub[gidx, :] = np.array(p, dtype=object)[cols].astype(np.float32)
            assert gidx == len(genes) - 1
            for j, i in enumerate(sel):
                rec.append(score(i, sub[:, j]))
    print(f"sample {s}: {n} cells x {len(genes)} genes; barcoded sorted-HSC cells in clone table: {len(rec)}", flush=True)
    df = pd.DataFrame(rec)
    cells_all.append(df.merge(clones.reset_index(), on="clone", how="left", suffixes=("", "_clone")))
cells = pd.concat(cells_all, ignore_index=True)
cells.to_csv(out / "larry_hsc_cell_scores.tsv.gz", sep="\t", index=False)
pd.DataFrame(cov_rows).to_csv(out / "larry_gene_set_coverage.tsv", sep="\t", index=False)
setnames = sorted(gs.set.unique())
agg = cells.groupby(["sample", "clone"]).agg(n_hsc_scored=("cell", "size"), **{k: (k, "mean") for k in setnames},
                                             total=("total", "mean"), ngenes=("ngenes", "mean")).reset_index()
keep = ["all", "HSC", "nonHSC", "freqHSC", "freqNon", "relout", "Mk", "Ery", "My", "Ly", "inactive", "active", "Mkbias"]
meta = cells.drop_duplicates(["sample", "clone"])[["sample", "clone"] + keep]
agg.merge(meta, on=["sample", "clone"]).to_csv(out / "larry_clone_scores.tsv", sep="\t", index=False)
print(agg.groupby("sample").size())
