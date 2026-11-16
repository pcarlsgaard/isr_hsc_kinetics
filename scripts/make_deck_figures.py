#!/usr/bin/env python3
"""Figures for the slide deck (exploratory results, rounds 2-4). Large labels, minimal text.
Reads only derived tables; writes PNGs to results/deck/. Palette validated with the dataviz validator."""
import pathlib, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
root = pathlib.Path(__file__).resolve().parents[1]
out = root / "results/deck"; out.mkdir(parents=True, exist_ok=True)
BG, INK, MUTE, GRID = "#FBFAF7", "#1B2430", "#5B6675", "#D9D6CE"
BLUE, VERM, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#A9AFB8"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 19, "axes.labelsize": 21, "xtick.labelsize": 19, "ytick.labelsize": 20,
    "axes.edgecolor": MUTE, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.linewidth": 1.2, "legend.frameon": False,
    "legend.fontsize": 19})
W, H = 16.64, 6.4
def save(fig, name): fig.savefig(out / f"{name}.png", dpi=200); plt.close(fig); print(name)
NAME = {"T21_ATF4_stress_only": "ATF4 targets: stress route only", "T21_ATF4_shared": "ATF4 targets: shared with mTORC1",
        "T21_mTORC1_nonATF4": "mTORC1 targets (ATF4-independent)", "RiBi_nonRP": "Ribosome biogenesis", "Cyto_translation_nonRP": "Translation factors",
        "Purine_biosynthesis": "Purine synthesis", "Pyrimidine_biosynthesis": "Pyrimidine synthesis", "OXPHOS": "Oxidative phosphorylation",
        "One_carbon": "One-carbon metabolism", "Serine_biosynthesis": "Serine synthesis", "Mito_translation": "Mitochondrial translation",
        "Glycolysis": "Glycolysis", "PolI_transcription": "rRNA transcription (Pol I)", "SSU_biogenesis": "Small ribosomal subunit biogenesis",
        "LSU_biogenesis": "Large ribosomal subunit biogenesis", "Cholesterol_biosynthesis": "Cholesterol synthesis",
        "Fatty_acid_biosynthesis": "Fatty-acid synthesis", "IEG_dissociation": "Immediate-early genes"}
ISR = {"T21_ATF4_stress_only", "T21_ATF4_shared"}

# ---- background: Milsom systems
ext = pd.read_csv(root / "results/milsom/system_external_variables.tsv", sep="\t")
nh = {"M96": 825, "M61": 147, "M32": 243, "M73": 26, "M51": 25}
lab = {"M96": "M9.6", "M61": "M6.1", "M32": "M3.2", "M73": "M7.3", "M51": "M5.1", "M16": "control", "M14": "control"}
fig, ax = plt.subplots(figsize=(9.6, 6.4))
y = np.arange(5)[::-1]
ax.barh(y, ext.tyMax_wk, color=BLUE, height=0.62)
for yy, (_, r) in zip(y, ext.iterrows()): ax.text(r.tyMax_wk + 0.4, yy, f"{nh[r.system]} HSCs", va="center", fontsize=20, color=MUTE)
ax.set_yticks(y); ax.set_yticklabels([lab[s] for s in ext.system], fontsize=24); ax.set_xlim(0, 27)
ax.set_xlabel("Time to peak blood chimerism (weeks)"); ax.tick_params(axis="y", length=0)
ax.text(-0.02, 1.02, "slowest", transform=ax.transAxes, ha="right", fontsize=18, color=MUTE); ax.text(-0.02, -0.04, "fastest", transform=ax.transAxes, ha="right", va="top", fontsize=18, color=MUTE)
fig.tight_layout(); save(fig, "bg_milsom")

# ---- LARRY clone table
old = pd.read_csv(root / "data/processed/larry/larry_hsc_cell_scores.tsv.gz", sep="\t")
cells = old.merge(pd.read_csv(root / "data/processed/larry/round3_cell_scores.tsv.gz", sep="\t"), on=["sample", "cell"])
S3 = ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4"]
cl = cells.groupby(["sample", "clone"]).agg(n=("cell", "size"), **{k: (k, "mean") for k in S3}, relout=("relout", "first"), freqHSC=("freqHSC", "first"),
                                             inactive=("inactive", "first"), active=("active", "first")).reset_index()
cl = cl[(cl.n >= 5) & cl.relout.notna()].copy()
cl["grp"] = np.where(cl.inactive == 1, "low", np.where(cl.active == 1, "high", "none"))
COL = {"low": BLUE, "high": VERM, "none": GREY}
fig, ax = plt.subplots(figsize=(9.6, 6.4))
c2 = cl.sort_values("relout").reset_index(drop=True)
ax.scatter(np.arange(len(c2)) + 1, c2.relout, c=c2.grp.map(COL), s=110, edgecolor=BG, linewidth=1.2, zorder=3)
ax.set_yscale("symlog", linthresh=0.05); ax.set_xlabel(f"{len(c2)} barcoded HSC clones, ranked"); ax.set_ylabel("Blood output per HSC")
ax.text(4, 0.02, f"low-output ({(c2.grp == 'low').sum()})", color=BLUE, fontsize=22, fontweight="bold")
ax.text(40, c2.relout.max() * 0.5, f"high-output ({(c2.grp == 'high').sum()})", color=VERM, fontsize=22, fontweight="bold")
ax.set_ylim(-0.005, c2.relout.max() * 2); fig.tight_layout(); save(fig, "bg_larry")

# ---- Milsom clone pairs
p3 = pd.read_csv(root / "results/round3/milsom_hsc_clone_pairs.tsv", sep="\t")
p2 = pd.read_csv(root / "results/round2/milsom_hsc_clone_pairs.tsv", sep="\t")
p2["pair"] = p2["first"] + "-" + p2["second"]; p2 = p2.rename(columns={"d_first_minus_second": "d"})
pp = pd.concat([p3[["pair", "set", "d", "ci_lo", "ci_hi"]], p2.loc[p2.set == "RiBi_nonRP", ["pair", "set", "d", "ci_lo", "ci_hi"]]])
sets = ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4", "RiBi_nonRP", "Cyto_translation_nonRP", "Purine_biosynthesis", "OXPHOS"]
pairs = [("M96-M32", "M9.6 − M3.2  (825 vs 243 HSCs)", BLUE, "o"), ("M61-M73", "M6.1 − M7.3  (147 vs 26)", VERM, "s"), ("M61-M51", "M6.1 − M5.1  (147 vs 25)", GREEN, "D")]
fig, ax = plt.subplots(figsize=(W, H))
for i, s in enumerate(sets):
    for j, (pk, pl, col, mk) in enumerate(pairs):
        r = pp[(pp.pair == pk) & (pp.set == s)].iloc[0]; yy = len(sets) - 1 - i + (1 - j) * 0.24
        ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=col, lw=3, solid_capstyle="round"); ax.scatter(r.d, yy, color=col, s=170, marker=mk, zorder=3, edgecolor=BG, linewidth=1.5, label=pl if i == 0 else None)
ax.axvline(0, color=INK, lw=1.4); ax.set_yticks(range(len(sets))); ax.set_yticklabels([NAME[s] for s in sets][::-1]); ax.tick_params(axis="y", length=0)
for k in range(len(sets) - 1): ax.axhline(k + 0.5, color=GRID, lw=1)
ax.set_xlabel("Difference in HSCs, slower − faster clone (SD)"); ax.set_ylim(-0.5, len(sets) - 0.5)
fig.legend(loc="upper center", ncol=3, columnspacing=1.6, handletextpad=0.3); fig.tight_layout(rect=(0, 0, 1, 0.92)); save(fig, "milsom_pairs")

# ---- genome-wide concordance
mc = pd.read_csv(root / "results/round2/milsom_hsc_pseudobulk_contrasts.tsv.gz", sep="\t", index_col=0); mc = mc[mc.pass_filter]
lg = pd.read_csv(root / "results/round2/larry_gene_low_vs_high_output.tsv.gz", sep="\t", index_col=0)
fig, axs = plt.subplots(1, 2, figsize=(W, H))
def panel(ax, x, yv, xl, yl):
    ax.scatter(x, yv, s=7, color=BLUE, alpha=0.18, linewidth=0, rasterized=True); ax.axhline(0, color=MUTE, lw=1); ax.axvline(0, color=MUTE, lw=1)
    rho = stats.spearmanr(x, yv)[0]; ax.text(0.04, 0.9, f"ρ = {rho:+.2f}".replace("-", "−"), transform=ax.transAxes, fontsize=30, fontweight="bold")
    ax.text(0.04, 0.82, f"{len(x):,} genes", transform=ax.transAxes, fontsize=19, color=MUTE); ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_xlim(np.percentile(x, 0.3), np.percentile(x, 99.7)); ax.set_ylim(np.percentile(yv, 0.3), np.percentile(yv, 99.7))
panel(axs[0], mc.A_M96_vs_M32, mc.B_M61_vs_M73M51, "M9.6 vs M3.2 (log₂ fold change)", "M6.1 vs M7.3 + M5.1 (log₂ FC)")
j = mc.join(lg[["d_low_minus_high"]], how="inner")
panel(axs[1], j.A_M96_vs_M32, j.d_low_minus_high, "M9.6 vs M3.2 (log₂ fold change)", "LARRY low − high output (SD)")
fig.tight_layout(w_pad=4); save(fig, "concordance")

# ---- LARRY low vs high bars
t = pd.read_csv(root / "results/round3/larry_clone_tests.tsv", sep="\t")
keep = ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4", "Pyrimidine_biosynthesis", "Mito_translation", "One_carbon", "OXPHOS", "Glycolysis",
        "Cyto_translation_nonRP", "Purine_biosynthesis", "PolI_transcription", "SSU_biogenesis", "LSU_biogenesis", "Fatty_acid_biosynthesis", "Cholesterol_biosynthesis"]
t2 = t[t.set.isin(keep)].sort_values("d_low_minus_high", ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(W, H + 0.5))
yy = np.arange(len(t2))
ax.barh(yy, t2.d_low_minus_high, color=[VERM if s in ISR else BLUE for s in t2.set], height=0.68)
for k, r in t2.iterrows():
    if r.p_low_vs_high >= 0.05: ax.text(r.d_low_minus_high - 0.03 if r.d_low_minus_high < 0 else r.d_low_minus_high + 0.03, k, "n.s.", va="center", ha="right" if r.d_low_minus_high < 0 else "left", fontsize=17, color=MUTE)
ax.set_yticks(yy); ax.set_yticklabels([NAME[s] for s in t2.set], fontsize=18.5); ax.tick_params(axis="y", length=0); ax.axvline(0, color=INK, lw=1.4)
for lbl, s in zip(ax.get_yticklabels(), t2.set):
    if s in ISR: lbl.set_color(VERM); lbl.set_fontweight("bold")
ax.set_xlabel("HSC score, low-output − high-output clones (SD)"); ax.set_xlim(-1.2, 0.35); ax.set_ylim(-0.6, len(t2) - 0.4)
fig.tight_layout(); save(fig, "larry_low_vs_high")

# ---- LARRY route scatter
rk = cl.groupby("sample")[S3].rank(pct=True) * 100
fig, axs = plt.subplots(1, 2, figsize=(W, H))
for ax, s, yl in [(axs[0], "T21_ATF4_shared", "ATF4 targets: shared with mTORC1"), (axs[1], "T21_ATF4_stress_only", "ATF4 targets: stress route only")]:
    for g in ["none", "high", "low"]:
        m = cl.grp == g; ax.scatter(rk.loc[m, "T21_mTORC1_nonATF4"], rk.loc[m, s], color=COL[g], s=150, edgecolor=BG, linewidth=1.3, zorder=3,
                                    label={"low": "low-output clone", "high": "high-output clone", "none": "unclassified"}[g])
    rho = stats.spearmanr(rk["T21_mTORC1_nonATF4"], rk[s])[0]
    ax.text(0.04, 0.9, f"ρ = {rho:+.2f}".replace("-", "−"), transform=ax.transAxes, fontsize=30, fontweight="bold", bbox=dict(facecolor=BG, edgecolor="none", pad=3))
    ax.set_xlabel("mTORC1 targets (percentile)"); ax.set_ylabel(yl + "\n(percentile)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h[::-1], l[::-1], loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0), handletextpad=0.2)
fig.tight_layout(rect=(0, 0, 1, 0.92), w_pad=4); save(fig, "larry_routes")

# ---- HSC-rich low-output clones
cl["x"] = cl.groupby("sample").freqHSC.rank(pct=True) * 100; cl["y"] = cl.groupby("sample").relout.rank(pct=True) * 100
q = (cl.groupby("sample").freqHSC.transform(lambda v: v > v.median())) & (cl.groupby("sample").relout.transform(lambda v: v < v.median()))
fig, axs = plt.subplots(1, 2, figsize=(W, H), gridspec_kw=dict(width_ratios=[1, 1.25]))
ax = axs[0]; ax.scatter(cl.x[~q], cl.y[~q], color=GREY, s=140, edgecolor=BG, linewidth=1.2); ax.scatter(cl.x[q], cl.y[q], color=VERM, s=170, edgecolor=BG, linewidth=1.2, zorder=3)
ax.set_title(f"HSC-rich, low-output clones ({int(q.sum())})", color=VERM, fontsize=19, fontweight="bold", loc="left", pad=14)
ax.set_xlabel("HSCs in clone (percentile)"); ax.set_ylabel("Blood output per HSC (percentile)")
qs = ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4", "Cyto_translation_nonRP", "Purine_biosynthesis", "One_carbon", "Pyrimidine_biosynthesis"]
tq = t.set_index("set").loc[qs]; ax = axs[1]; yy = np.arange(len(qs))[::-1]
ax.barh(yy, tq.d_hscrich_lowout, color=[VERM if s in ISR else BLUE for s in qs], height=0.62); ax.axvline(0, color=INK, lw=1.4)
for k, (s, r) in zip(yy, tq.iterrows()):
    txt = f"p = {r.p_hscrich:.3f}" if r.p_hscrich < 0.05 else "n.s."
    ax.text(0.03, k, txt, va="center", ha="left", fontsize=18, color=INK if r.p_hscrich < 0.05 else MUTE)
ax.set_yticks(yy); ax.set_yticklabels([NAME[s] for s in qs], fontsize=18.5); ax.tick_params(axis="y", length=0); ax.set_xlim(-0.9, 0.3)
ax.set_xlabel("HSC-rich low-output − others (SD)"); fig.tight_layout(w_pad=3); save(fig, "larry_hsc_rich")

# ---- HSC vs MPP per mouse
r4 = pd.read_csv(root / "results/round4/milsom_hsc_vs_progenitors.tsv", sep="\t"); r4 = r4[r4.contrast == "HSC_minus_MPP"].set_index("set")
sets = ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4", "RiBi_nonRP", "Cyto_translation_nonRP", "Purine_biosynthesis", "OXPHOS", "IEG_dissociation"]
fig, ax = plt.subplots(figsize=(W, H)); rng = np.random.default_rng(1)
for i, s in enumerate(sets):
    yy = len(sets) - 1 - i
    cv = r4.loc[s, ["M96", "M61", "M32", "M73", "M51"]].to_numpy(float); kv = r4.loc[s, ["M16", "M14"]].to_numpy(float)
    ax.scatter(cv, yy + np.linspace(-0.2, 0.2, 5), color=BLUE, s=190, edgecolor=BG, linewidth=1.5, zorder=3, label="clonal mouse (5)" if i == 0 else None)
    ax.scatter(kv, yy + np.array([-0.1, 0.1]), color=VERM, s=190, marker="D", edgecolor=BG, linewidth=1.5, zorder=4, label="polyclonal control (2)" if i == 0 else None)
    ax.plot([r4.loc[s, "median"]] * 2, [yy - 0.36, yy + 0.36], color=INK, lw=3.5, zorder=2)
ax.axvline(0, color=INK, lw=1.4); ax.set_yticks(range(len(sets))); ax.set_yticklabels([NAME[s] for s in sets][::-1]); ax.tick_params(axis="y", length=0)
for k in range(len(sets) - 1): ax.axhline(k + 0.5, color=GRID, lw=1)
ax.set_xlabel("HSC − MPP within each mouse (SD)"); ax.set_ylim(-0.5, len(sets) - 0.5)
fig.legend(loc="upper center", ncol=2, handletextpad=0.2); fig.tight_layout(rect=(0, 0, 1, 0.92)); save(fig, "hsc_vs_mpp")

# ---- shared ATF4 genes, HSC vs MPP
pb = pd.read_csv(root / "data/processed/milsom/pseudobulk_counts_system_by_celltype.tsv.gz", sep="\t", index_col=0)
M = ["M96", "M61", "M32", "M73", "M51", "M16", "M14"]; lc = np.log2(pb / pb.sum() * 1e6 + 1)
fc = pd.DataFrame({m: lc[f"{m}|HSCs"] - lc[f"{m}|MPPs"] for m in M})
ok = pd.concat([(pb[f"{m}|HSCs"] + pb[f"{m}|MPPs"]) >= 10 for m in M], axis=1).all(axis=1)
g3 = pd.read_csv(root / "metadata/gene_sets_round3.tsv", sep="\t")
gsh = fc[ok].loc[fc[ok].index.intersection(g3[g3.set == "T21_ATF4_shared"].gene)].median(axis=1).sort_values(ascending=False)
fig, ax = plt.subplots(figsize=(W, H)); x = np.arange(len(gsh))
ax.bar(x, gsh.values, color=[BLUE if v > 0 else VERM for v in gsh.values], width=0.72); ax.axhline(0, color=INK, lw=1.4)
for k in list(range(11)): ax.text(x[k], gsh.iloc[k] + 0.06, gsh.index[k], rotation=90, ha="center", va="bottom", fontsize=19, style="italic")
for k in list(range(len(gsh) - 12, len(gsh))): ax.text(x[k], gsh.iloc[k] - 0.06, gsh.index[k], rotation=90, ha="center", va="top", fontsize=19, style="italic")
ax.set_xticks([]); ax.spines["bottom"].set_visible(False); ax.set_ylabel("HSC vs MPP (log₂ fold change,\nmedian of 7 mice)"); ax.set_ylim(-3.0, 4.2)
ax.text(len(gsh) * 0.33, 3.3, "higher in HSCs", color=BLUE, fontsize=24, fontweight="bold"); ax.text(len(gsh) * 0.33, -2.6, "lower in HSCs", color=VERM, fontsize=24, fontweight="bold")
ax.set_xlabel(f"{len(gsh)} ATF4 targets shared by the stress and mTORC1 routes"); fig.tight_layout(); save(fig, "shared_genes")

# ---- per-clone view: each clone vs its same-run control, ordered slowest to fastest
order = [("M96-M16", "M9.6", 825), ("M61-M14", "M6.1", 147), ("M32-M16", "M3.2", 243), ("M73-M14", "M7.3", 26), ("M51-M14", "M5.1", 25)]
fig, axs = plt.subplots(1, 3, figsize=(W, H), sharey=True)
for ax, s, col in zip(axs, ["T21_ATF4_stress_only", "T21_ATF4_shared", "T21_mTORC1_nonATF4"], [VERM, VERM, BLUE]):
    for i, (pk, nm, n) in enumerate(order):
        r = p3[(p3.pair == pk) & (p3.set == s)].iloc[0]
        ax.plot([i, i], [r.ci_lo, r.ci_hi], color=col, lw=4, solid_capstyle="round"); ax.scatter(i, r.d, color=col, s=260, zorder=3, edgecolor=BG, linewidth=2)
    ax.axhline(0, color=INK, lw=1.4); ax.set_xticks(range(5)); ax.set_xticklabels([f"{nm}\n{n}" for _, nm, n in order], fontsize=19)
    ax.set_title({"T21_ATF4_stress_only": "ATF4: stress route only", "T21_ATF4_shared": "ATF4: shared with mTORC1", "T21_mTORC1_nonATF4": "mTORC1, ATF4-independent"}[s], fontsize=21, fontweight="bold", color=col, pad=16); ax.set_xlim(-0.5, 4.5); ax.tick_params(axis="x", length=0)
axs[0].set_ylabel("HSCs, clone − same-run control (SD)")
axs[1].set_xlabel("slowest  →  fastest clone   (number of HSCs below name)", labelpad=14)
fig.tight_layout(w_pad=2.5); save(fig, "milsom_by_clone")

# ---- Amendment 5: volcano plots and gene-set landscapes
T21 = [("T21_mTORC1_nonATF4", "mTORC1 targets, ATF4-independent", BLUE), ("T21_ATF4_shared", "ATF4: shared with mTORC1", "#8A6FB0"), ("T21_ATF4_stress_only", "ATF4: stress route only", VERM)]
def volcano(tab, xcol, xl, name):
    d = pd.read_csv(root / f"results/round5/{tab}", sep="\t", index_col=0); d["y"] = -np.log10(d.p)
    fig, ax = plt.subplots(figsize=(W, H))
    ax.scatter(d[xcol], d.y, s=9, color=GREY, alpha=0.45, linewidth=0, rasterized=True, label=f"all {len(d):,} genes")
    for k, lbl, col in [T21[0], T21[2]]:
        m = d.index.intersection(g3[g3.set == k].gene); ax.scatter(d.loc[m, xcol], d.loc[m, "y"], s=75, color=col, edgecolor=BG, linewidth=0.8, zorder=3, label=f"{lbl} ({len(m)})")
    thr = d[xcol].abs().quantile(0.98); c = d[d[xcol].abs() >= thr]
    for side, xa, ha, head in ((1, 0.99, "right", "top up"), (-1, 0.01, "left", "top down")):
        gl = list(c[np.sign(c[xcol]) == side].sort_values("t", key=abs, ascending=False).head(5).index)
        ax.text(xa, 0.98, head + "\n" + "\n".join(gl), transform=ax.transAxes, ha=ha, va="top", fontsize=19, style="italic", color=INK, linespacing=1.25, zorder=6, bbox=dict(facecolor=BG, edgecolor="none", alpha=0.8, pad=3))
    ax.axvline(0, color=MUTE, lw=1); ax.set_xlabel(xl); ax.set_ylabel("−log₁₀ p"); lim = d[xcol].abs().quantile(0.9995) * 1.1; ax.set_xlim(-lim, lim)
    h, l = ax.get_legend_handles_labels(); lg = fig.legend(h, l, loc="upper center", ncol=3, handletextpad=0.2, markerscale=1.6)
    lg.legend_handles[0].set_sizes([90]); lg.legend_handles[0].set_alpha(1)
    fig.tight_layout(rect=(0, 0, 1, 0.92)); save(fig, name)
volcano("hsc_vs_mpp_genes.tsv.gz", "log2fc", "HSC vs MPP (log₂ fold change, mean of 7 mice)", "volcano_hsc_mpp")
volcano("larry_low_vs_high_genes.tsv.gz", "diff", "Low-output − high-output clones (SD)", "volcano_larry")
def pretty(s): return s.replace("HALLMARK_", "").replace("_", " ").capitalize().replace("Myc", "MYC").replace("E2f", "E2F").replace("G2m", "G2M").replace("Tnfa", "TNFα").replace("nfkb", "NF-κB").replace("Il2 stat5", "IL2–STAT5").replace("Il6 jak stat3", "IL6–JAK–STAT3").replace("Kras", "KRAS").replace("Mtorc1", "mTORC1").replace("Dna", "DNA")
def landscape_fig(tab, name):
    r = pd.read_csv(root / f"results/round5/{tab}", sep="\t"); n = len(r)
    fig, ax = plt.subplots(figsize=(W, H)); cnt, edges, _ = ax.hist(r.score, bins=70, color=GREY, edgecolor=BG, linewidth=0.6); top = cnt.max()
    hall = r[r.set.str.startswith("HALLMARK")]
    refs = [(x.set, pretty(x.set), INK, x.score, x.rank_from_top) for x in list(hall.head(2).itertuples()) + list(hall.tail(2).itertuples())]
    mine = [(k, lbl, col, r.loc[r.set == k, "score"].iloc[0], int(r.loc[r.set == k, "rank_from_top"].iloc[0])) for k, lbl, col in T21]
    items = sorted(refs + mine, key=lambda v: v[3]); med = np.median([v[3] for v in items])
    left = [v for v in items if v[3] < med]; right = [v for v in items if v[3] >= med]
    plan = list(zip(left, [1.72, 1.52, 1.32, 1.12], ["left"] * 4)) + list(zip(right, [1.62, 1.42, 1.22, 1.02], ["right"] * 4))
    for (k, lbl, col, sc, rk), lv, ha in plan:
        yl = top * lv; strong = k.startswith("T21")
        ax.plot([sc, sc], [0, yl], color=col, lw=3.2 if strong else 1.6, ls="-" if strong else (0, (4, 3)))
        ax.text(sc + (0.012 if ha == "left" else -0.012), yl, f"{lbl} · rank {rk:,}", color=col, fontsize=19 if strong else 18, fontweight="bold" if strong else "normal", ha=ha, va="top", zorder=6, bbox=dict(facecolor=BG, edgecolor="none", pad=1.5))
    ax.set_yticks([t for t in ax.get_yticks() if t <= top])
    ax.set_ylim(0, top * 1.8); ax.set_xlabel(f"Gene-set score: mean shift of member genes (SD) · {n:,} gene sets ranked"); ax.set_ylabel("Number of gene sets", y=0.3)
    fig.tight_layout(); save(fig, name)
landscape_fig("hsc_vs_mpp_set_landscape.tsv", "landscape_hsc_mpp"); landscape_fig("larry_low_vs_high_set_landscape.tsv", "landscape_larry")
