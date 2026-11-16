#!/usr/bin/env python3
"""Round 3 (exploratory, Amendment 3): published ISR-route and biosynthesis sets in HSCs.
Milsom: cell-level standardised differences for within-run clone pairs and clone vs control.
LARRY: scores per sorted-HSC cell (UCell statistic, maxRank 1500), clone-level tests as in round 2."""
import gzip, pathlib
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import rankdata
root = pathlib.Path(__file__).resolve().parents[1]
out = root / "results/round3"; out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261116); NBOOT, NPERM, MAXRANK, MINN = 2000, 20000, 1500, 5
gs = pd.read_csv(root / "metadata/gene_sets_round3.tsv", sep="\t")
SETS = list(dict.fromkeys(gs.set))
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)

# ---------- Milsom
c = pd.read_csv(root / "data/processed/milsom/hspc_cell_scores.tsv.gz", sep="\t")[["cell", "Mice", "Celltypes_Louvain"]]
c = c[c.Celltypes_Louvain == "HSCs"].merge(pd.read_csv(root / "data/processed/milsom/round3_hsc_cell_scores.tsv.gz", sep="\t"), on="cell")
PAIRS = [("M96", "M32"), ("M61", "M73"), ("M61", "M51"), ("M73", "M51"), ("M96", "M16"), ("M32", "M16"), ("M61", "M14"), ("M73", "M14"), ("M51", "M14")]
def d(a, b):
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp
rows = []
for s1, s2 in PAIRS:
    for s in SETS:
        a, b = c.loc[c.Mice == s1, s].to_numpy(), c.loc[c.Mice == s2, s].to_numpy()
        bs = [d(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(NBOOT)]
        rows.append(dict(pair=f"{s1}-{s2}", n_first=len(a), n_second=len(b), set=s, d=d(a, b), ci_lo=np.percentile(bs, 2.5), ci_hi=np.percentile(bs, 97.5)))
m = pd.DataFrame(rows); m.to_csv(out / "milsom_hsc_clone_pairs.tsv", sep="\t", index=False, float_format="%.4g")
m["txt"] = m.apply(lambda r: f"{r.d:+.2f}{'*' if r.ci_lo * r.ci_hi > 0 else ' '}", axis=1)
print("Milsom HSCs, standardised difference first minus second (* = bootstrap interval over cells excludes 0)")
print(m.pivot(index="set", columns="pair", values="txt").loc[SETS, [f"{a}-{b}" for a, b in PAIRS]].to_string())

# ---------- LARRY
proc = root / "data/processed/larry"; cache = proc / "round3_cell_scores.tsv.gz"
repo = root / "data/raw/tcf15_github/repo"; mat = root / "data/raw/tcf15_dropbox/filtered_counts_csv"
old = pd.read_csv(proc / "larry_hsc_cell_scores.tsv.gz", sep="\t")
if not cache.exists():
    S = {1: ("ARF2019_Sample1_new_counts_norm.csv.gz", "ARF2019_Sample1_genes.txt"), 2: ("ARF2019_Sample2_counts_norm.csv.gz", "ARF2019_Sample2_genes.txt"),
         3: ("ARF2019_Sample3_counts_norm.csv.gz", "ARF2019_Sample3_genes.txt")}
    recs, cov = [], []
    for s, (mf, gf) in S.items():
        genes = [x for x in open(repo / gf).read().split("\n") if x != ""]; gi = {g: i for i, g in enumerate(genes)}
        sets = {k: np.array(sorted({gi[x] for x in dd.gene if x in gi}), dtype=int) for k, dd in gs.groupby("set")}
        cov += [dict(sample=s, set=k, n_defined=gs[gs.set == k].gene.nunique(), n_present=len(v)) for k, v in sets.items()]
        o = old[old["sample"] == s]; sel = o.cell.to_numpy()
        X = np.zeros((len(sel), len(genes)), dtype=np.float32)
        with gzip.open(mat / mf, "rt") as fh:
            if s == 1:
                pos = {cc: j for j, cc in enumerate(sel)}
                for i, line in enumerate(fh):
                    if i in pos: X[pos[i], :] = np.array(line.rstrip("\n").split(","), dtype=np.float32)
            else:
                for gidx, line in enumerate(fh):
                    p = line.rstrip("\n").split(","); assert p[0] == genes[gidx]
                    X[:, gidx] = np.array(p, dtype=object)[sel + 1].astype(np.float32)
        assert np.allclose(X.sum(1), o.total.to_numpy(), rtol=1e-3)
        for j in range(len(sel)):
            r = rankdata(-X[j], method="average"); r[r > MAXRANK] = MAXRANK + 1
            row = dict(sample=s, cell=int(sel[j]))
            for k, idx in sets.items():
                n = len(idx); row[k] = 1 - (r[idx].sum() - n * (n + 1) / 2) / (n * MAXRANK)
            recs.append(row)
    pd.DataFrame(recs).to_csv(cache, sep="\t", index=False)
    pd.DataFrame(cov).to_csv(proc / "round3_gene_set_coverage.tsv", sep="\t", index=False)
cells = old.merge(pd.read_csv(cache, sep="\t"), on=["sample", "cell"])
cl = cells.groupby(["sample", "clone"]).agg(n=("cell", "size"), **{k: (k, "mean") for k in SETS + ["Stem_MolO", "Stem_dHSC"]},
        relout=("relout", "first"), freqHSC=("freqHSC", "first"), inactive=("inactive", "first"), active=("active", "first")).reset_index()
cl = cl[(cl.n >= MINN) & cl.relout.notna()].copy()
def perm_group(df, a, g):
    z = df.groupby("sample")[a].transform(lambda x: (x - x.mean()) / x.std(ddof=1)).to_numpy(); g = df[g].to_numpy(); s = df["sample"].to_numpy()
    obs = z[g == 1].mean() - z[g == 0].mean(); null = np.empty(NPERM)
    for k in range(NPERM):
        gp = g.copy()
        for ss in np.unique(s):
            mm = s == ss; gp[mm] = rng.permutation(g[mm])
        null[k] = z[gp == 1].mean() - z[gp == 0].mean()
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)
def strat_rho(df, a, b, covars=()):
    parts = []
    for s, g in df.groupby("sample"):
        if len(g) < 5: continue
        ra, rb = stats.rankdata(g[a]), stats.rankdata(g[b])
        if covars:
            C = np.column_stack([np.ones(len(g))] + [stats.rankdata(g[cc]) for cc in covars])
            ra = ra - C @ np.linalg.lstsq(C, ra, rcond=None)[0]; rb = rb - C @ np.linalg.lstsq(C, rb, rcond=None)[0]
        parts.append(((ra - ra.mean()) / ra.std(), (rb - rb.mean()) / rb.std()))
    return np.mean(np.concatenate([x * y for x, y in parts]))
lab = cl[(cl.inactive == 1) | (cl.active == 1)].copy(); lab["low"] = (lab.inactive == 1).astype(int)
cl["hsc_rich_low_out"] = ((cl.groupby("sample").freqHSC.transform(lambda x: x > x.median())) & (cl.groupby("sample").relout.transform(lambda x: x < x.median()))).astype(int)
rows = []
for s in SETS:
    d1, p1 = perm_group(lab, s, "low"); d2, p2 = perm_group(cl, s, "hsc_rich_low_out")
    rows.append(dict(set=s, d_low_minus_high=d1, p_low_vs_high=p1, rho_output_per_HSC=strat_rho(cl, s, "relout"),
                     rho_adj_stemness=strat_rho(cl, s, "relout", ("Stem_MolO", "Stem_dHSC")), d_hscrich_lowout=d2, p_hscrich=p2))
res = pd.DataFrame(rows); res["bh_low_vs_high"] = stats.false_discovery_control(res.p_low_vs_high)
res.to_csv(out / "larry_clone_tests.tsv", sep="\t", index=False, float_format="%.4g")
print(f"\nLARRY clones: low {int(lab.low.sum())}, high {int((1 - lab.low).sum())}; HSC-rich low-output {int(cl.hsc_rich_low_out.sum())} of {len(cl)}")
print(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
rk = cl.groupby("sample")[SETS].rank(pct=True)
print("\nclone-level Spearman (within-experiment ranks):")
print(rk[["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4", "R_ATF4_ER_stress", "LSU_biogenesis", "Purine_biosynthesis"]].corr(method="spearman").round(2).to_string())
print(pd.read_csv(proc / "round3_gene_set_coverage.tsv", sep="\t").pivot(index="set", columns="sample", values="n_present").to_string())
