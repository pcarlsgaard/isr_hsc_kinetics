#!/usr/bin/env python3
"""Round 2 (exploratory, Amendment 2), LARRY sorted HSCs: round-2 set scores per cell (UCell statistic,
maxRank 1500, authors' normalised matrices), per-clone gene means of log1p normalised expression,
unadjusted and stemness-adjusted clone tests, HSC-rich/low-output clones, genome-wide low-vs-high
output differences and their agreement with the Milsom HSC contrasts."""
import gzip, pathlib, re
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import rankdata
root = pathlib.Path(__file__).resolve().parents[1]
repo = root / "data/raw/tcf15_github/repo"; mat = root / "data/raw/tcf15_dropbox/filtered_counts_csv"
out = root / "results/round2"; out.mkdir(parents=True, exist_ok=True)
proc = root / "data/processed/larry"
MAXRANK, NPERM, MINN = 1500, 20000, 5
rng = np.random.default_rng(20261116)
gs2 = pd.read_csv(root / "metadata/gene_sets_round2.tsv", sep="\t")
gs1 = pd.read_csv(root / "data/processed/gene_sets.tsv", sep="\t")
S = {1: dict(m="ARF2019_Sample1_new_counts_norm.csv.gz", g="ARF2019_Sample1_genes.txt"),
     2: dict(m="ARF2019_Sample2_counts_norm.csv.gz", g="ARF2019_Sample2_genes.txt"),
     3: dict(m="ARF2019_Sample3_counts_norm.csv.gz", g="ARF2019_Sample3_genes.txt")}
old = pd.read_csv(proc / "larry_hsc_cell_scores.tsv.gz", sep="\t")
cache = proc / "round2_cell_scores.tsv.gz"; gcache = proc / "round2_clone_gene_means.tsv.gz"
if not (cache.exists() and gcache.exists()):
    recs, gm = [], []
    for s, f in S.items():
        genes = [x for x in open(repo / f["g"]).read().split("\n") if x != ""]
        gi = {g: i for i, g in enumerate(genes)}
        sets = {k: np.array(sorted({gi[x] for x in d.gene if x in gi}), dtype=int) for k, d in gs2.groupby("set")}
        o = old[old["sample"] == s]; sel = o.cell.to_numpy(); order = np.argsort(sel); sel_sorted = sel[order]
        X = np.zeros((len(sel), len(genes)), dtype=np.float32)  # cells x genes
        with gzip.open(mat / f["m"], "rt") as fh:
            if s == 1:
                pos = {c: j for j, c in enumerate(sel)}
                for i, line in enumerate(fh):
                    if i in pos: X[pos[i], :] = np.array(line.rstrip("\n").split(","), dtype=np.float32)
            else:
                cols = sel + 1
                for gidx, line in enumerate(fh):
                    p = line.rstrip("\n").split(","); assert p[0] == genes[gidx]
                    X[:, gidx] = np.array(p, dtype=object)[cols].astype(np.float32)
        assert np.allclose(X.sum(1), o.total.to_numpy(), rtol=1e-3), "cell alignment check against round-1 totals failed"
        for j in range(len(sel)):
            r = rankdata(-X[j], method="average"); r[r > MAXRANK] = MAXRANK + 1
            row = dict(sample=s, cell=int(sel[j]))
            for k, idx in sets.items():
                n = len(idx); row[k] = 1 - (r[idx].sum() - n * (n + 1) / 2) / (n * MAXRANK)
            recs.append(row)
        L = pd.DataFrame(np.log1p(X), columns=genes); L.insert(0, "clone", o.clone.to_numpy())
        m = L.groupby("clone").mean(); m.insert(0, "sample", s)
        gm.append(m.reset_index().melt(id_vars=["sample", "clone"], var_name="gene", value_name="mean_log1p"))
        print(f"sample {s}: {len(sel)} cells, set sizes {({k: len(v) for k, v in sets.items()})}", flush=True)
    pd.DataFrame(recs).to_csv(cache, sep="\t", index=False)
    pd.concat(gm).to_csv(gcache, sep="\t", index=False, float_format="%.5g")
cells = old.merge(pd.read_csv(cache, sep="\t"), on=["sample", "cell"])
NEW = sorted(gs2.set.unique())
SETS = ["ISR_stress_specific", "ISR_stress_specific_noIEG", "ATF4_biosynth", "Growth_noATF4", "RiBi_nonRP", "RP_cytoplasmic",
        "ISR_adaptive", "ISR_terminal", "mTORC1_hallmark", "IEG_dissociation", "Stem_MolO", "Stem_dHSC"]
cl = cells.groupby(["sample", "clone"]).agg(n=("cell", "size"), **{k: (k, "mean") for k in SETS},
        relout=("relout", "first"), freqHSC=("freqHSC", "first"), inactive=("inactive", "first"), active=("active", "first")).reset_index()
cl = cl[(cl.n >= MINN) & cl.relout.notna()].copy()

def perm_group(df, a, g):
    z = df.groupby("sample")[a].transform(lambda x: (x - x.mean()) / x.std(ddof=1)).to_numpy(); g = df[g].to_numpy(); s = df["sample"].to_numpy()
    obs = z[g == 1].mean() - z[g == 0].mean(); null = np.empty(NPERM)
    for k in range(NPERM):
        gp = g.copy()
        for ss in np.unique(s):
            m = s == ss; gp[m] = rng.permutation(g[m])
        null[k] = z[gp == 1].mean() - z[gp == 0].mean()
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)
def strat_rho(df, a, b, covars=()):
    parts = []
    for s, g in df.groupby("sample"):
        if len(g) < 5: continue
        ra, rb = stats.rankdata(g[a]), stats.rankdata(g[b])
        if covars:
            C = np.column_stack([np.ones(len(g))] + [stats.rankdata(g[c]) for c in covars])
            ra = ra - C @ np.linalg.lstsq(C, ra, rcond=None)[0]; rb = rb - C @ np.linalg.lstsq(C, rb, rcond=None)[0]
        parts.append(((ra - ra.mean()) / ra.std(), (rb - rb.mean()) / rb.std()))
    obs = np.mean(np.concatenate([x * y for x, y in parts])); null = np.empty(NPERM)
    for k in range(NPERM): null[k] = np.mean(np.concatenate([x * rng.permutation(y) for x, y in parts]))
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)

lab = cl[(cl.inactive == 1) | (cl.active == 1)].copy(); lab["low"] = (lab.inactive == 1).astype(int)
cl["hsc_rich_low_out"] = ((cl.groupby("sample").freqHSC.transform(lambda x: x > x.median())) &
                          (cl.groupby("sample").relout.transform(lambda x: x < x.median()))).astype(int)
rows = []
for s in SETS:
    d1, p1 = perm_group(lab, s, "low"); d2, p2 = perm_group(cl, s, "hsc_rich_low_out")
    r, pr = strat_rho(cl, s, "relout")
    ra, pa = strat_rho(cl, s, "relout", ("Stem_MolO", "Stem_dHSC")) if not s.startswith("Stem_") else (np.nan, np.nan)
    rows.append(dict(set=s, n_low=int(lab.low.sum()), n_high=int((1 - lab.low).sum()), d_low_minus_high=d1, p_low_vs_high=p1,
                     rho_output_per_HSC=r, p_rho=pr, rho_adj_stemness=ra, p_rho_adj=pa,
                     n_hsc_rich_low_out=int(cl.hsc_rich_low_out.sum()), n_other=int((1 - cl.hsc_rich_low_out).sum()),
                     d_hscrich_lowout_minus_other=d2, p_hscrich=p2))
res = pd.DataFrame(rows); res.to_csv(out / "larry_round2_clone_tests.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\nclone-level Spearman among scores (within-experiment ranks pooled):")
rk = cl.groupby("sample")[SETS].rank(pct=True)
print(rk[["ISR_stress_specific", "ATF4_biosynth", "Growth_noATF4", "RiBi_nonRP", "Stem_MolO", "IEG_dissociation"]].corr(method="spearman").round(2).to_string())
print("\noverlap of stemness sets with ISR/translation sets:")
A = {**{k: set(d.gene) for k, d in gs1.groupby("set")}, **{k: set(d.gene) for k, d in gs2.groupby("set")}}
for st in ("Stem_MolO", "Stem_dHSC", "Stem_LowOutput"):
    for k in ("ISR_stress_specific", "ISR_terminal", "ATF4_biosynth", "Growth_noATF4", "RiBi_nonRP"):
        print(f"  {st} & {k}: {sorted(A[st] & A[k])}")

# genome-wide
gmn = pd.read_csv(gcache, sep="\t").merge(lab[["sample", "clone", "low"]], on=["sample", "clone"])
gmn["z"] = gmn.groupby(["sample", "gene"]).mean_log1p.transform(lambda x: (x - x.mean()) / x.std(ddof=1) if x.std(ddof=1) > 0 else np.nan)
det = gmn.groupby("gene").mean_log1p.mean()
nsamp = gmn.groupby("gene")["sample"].nunique()
t = gmn.dropna(subset=["z"]).groupby(["gene", "low"]).z.mean().unstack()
g = pd.DataFrame({"d_low_minus_high": t[1] - t[0], "mean_log1p": det, "n_experiments": nsamp}).dropna()
g = g[(g.n_experiments == 3) & (g.mean_log1p >= g.mean_log1p.quantile(0.5))]
g["pctrank"] = g.d_low_minus_high.rank(pct=True)
g.sort_values("d_low_minus_high").to_csv(out / "larry_gene_low_vs_high_output.tsv.gz", sep="\t", float_format="%.4g")
print(f"\ngenes ranked: {len(g)}")
print("top higher in low-output clones:", ", ".join(g.sort_values("d_low_minus_high").index[-40:][::-1]))
print("top lower in low-output clones:", ", ".join(g.sort_values("d_low_minus_high").index[:40]))
sr = []
for k, m in A.items():
    m = g.index.intersection(list(m))
    if k == "RP_cytoplasmic": m = g.index[g.index.str.match(r"^(Rp[ls][0-9]+[a-z]?|Rplp[0-2]|Rpsa)$")]
    if len(m) >= 5: sr.append(dict(set=k, n_genes=len(m), mean_d=g.loc[m, "d_low_minus_high"].mean(), mean_pctrank=g.loc[m, "pctrank"].mean()))
sr = pd.DataFrame(sr).sort_values("mean_pctrank"); sr.to_csv(out / "larry_set_gene_ranks.tsv", sep="\t", index=False, float_format="%.4g")
print(sr.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
mc = pd.read_csv(out / "milsom_hsc_pseudobulk_contrasts.tsv.gz", sep="\t", index_col=0)
mc = mc[mc.pass_filter].join(g[["d_low_minus_high"]], how="inner")
print(f"\ncross-dataset, {len(mc)} shared genes; Spearman of LARRY low-minus-high with Milsom contrasts:")
cr = {k: stats.spearmanr(mc[k], mc.d_low_minus_high)[0] for k in ["A_M96_vs_M32", "B_M61_vs_M73M51", "D_M96_vs_M16", "M61_vs_M14", "M32_vs_M16"]}
pd.Series(cr, name="spearman").to_csv(out / "cross_dataset_concordance.tsv", sep="\t", float_format="%.4g")
print(pd.Series(cr).round(3).to_string())
