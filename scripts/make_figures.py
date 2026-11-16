#!/usr/bin/env python3
"""Figures 1-3 for the exploratory report. Reads results/ and data/processed/ tables only."""
import pathlib
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
root = pathlib.Path(__file__).resolve().parents[1]
fig_dir = root / "results/figures"; fig_dir.mkdir(parents=True, exist_ok=True)
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF, "font.size": 9,
                     "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "text.color": INK, "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
ORDER = ["M96", "M61", "M32", "M73", "M51"]; CTRL = ["M14", "M16"]
LAB = {"M96": "M9.6", "M61": "M6.1", "M32": "M3.2", "M73": "M7.3", "M51": "M5.1", "M14": "ctrl\nM14", "M16": "ctrl\nM16"}
RUN = {"M61": 1, "M73": 1, "M51": 1, "M14": 1, "M96": 2, "M32": 2, "M16": 2}
PRIMARY = [("ISR_adaptive", "ISR adaptive arm"), ("ISR_terminal", "ISR terminal arm"), ("RiBi_nonRP", "Ribosome biogenesis"),
           ("RP_cytoplasmic", "Ribosomal protein genes"), ("mTORC1_hallmark", "mTORC1 hallmark")]

# ---- Figure 1: primary signatures across systems (annotated HSCs)
b = pd.read_csv(root / "results/milsom/milsom_hsc_system_means_bootstrap.tsv", sep="\t")
fig, axes = plt.subplots(1, 5, figsize=(13.5, 3.5), sharex=True)
for ax, (s, title) in zip(axes, PRIMARY):
    d = b[b.set == s].set_index("system")
    for i, k in enumerate(ORDER + CTRL):
        col = BLUE if k in ORDER else MUTED
        ax.plot([i, i], [d.lo[k], d.hi[k]], color=col, lw=2, solid_capstyle="round")
        ax.plot(i, d["mean"][k], marker="o" if RUN[k] == 1 else "s", ms=8, color=col, mec=SURF, mew=2)
    ax.axvline(4.5, color=GRID, lw=1)
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.set_xticks(range(7)); ax.set_xticklabels([LAB[k] for k in ORDER + CTRL], fontsize=8)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.tick_params(length=0)
axes[0].set_ylabel("Mean UCell score in annotated HSCs")
n = b.drop_duplicates("system").set_index("system").n_cells
fig.text(0.008, 0.02, "Clonal systems left to right: slowest to fastest reconstitution.  HSCs per system: " +
         ", ".join(f"{LAB[k].replace(chr(10), ' ')} {n[k]}" for k in ORDER + CTRL) + ".  Bars: 95% bootstrap interval over cells (not over clones).",
         fontsize=7.5, color=INK2)
fig.legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, ms=7, label="clonal system, run 1"),
                    Line2D([], [], marker="s", ls="", color=BLUE, ms=7, label="clonal system, run 2"),
                    Line2D([], [], marker="o", ls="", color=MUTED, ms=7, label="polyclonal control, run 1"),
                    Line2D([], [], marker="s", ls="", color=MUTED, ms=7, label="polyclonal control, run 2")],
           loc="upper right", ncol=4, frameon=False, fontsize=8, bbox_to_anchor=(0.995, 1.0))
fig.suptitle("Translation-related signatures in HSCs across five clonal systems", x=0.008, ha="left", fontsize=12, y=0.995)
fig.tight_layout(rect=(0, 0.06, 1, 0.93)); fig.savefig(fig_dir / "fig1_signatures_across_systems.png", dpi=200); plt.close(fig)

# ---- Figure 2: key genes, pseudobulk log2 CPM, HSCs and MPPs, centred on the same-run control
pb = pd.read_csv(root / "data/processed/milsom/pseudobulk_counts_system_by_celltype.tsv.gz", sep="\t", index_col=0)
groups = [("ISR adaptive", ["Atf4", "Asns", "Psat1", "Phgdh", "Shmt2", "Mthfd2", "Slc7a5", "Slc3a2", "Slc7a11", "Sesn2", "Eif4ebp1", "Nars", "Gars", "Mars"]),
          ("ISR terminal", ["Ddit3", "Ppp1r15a", "Trib3", "Bbc3", "Atf5", "Gadd45a", "Nupr1", "Txnip"]),
          ("Ribosome biogenesis", ["Ncl", "Fbl", "Npm1", "Polr1a", "Polr1b", "Ubtf", "Nop56", "Nop58", "Dkc1"]),
          ("Translation control", ["Myc", "Mycn", "Eif4e", "Eif4ebp2", "Eif2s1", "Eif2ak1", "Eif2ak2", "Eif2ak4", "Eef2", "Eef1a1", "Rpl24", "Rps6"]),
          ("HSC state", ["Hlf", "Mecom", "Procr", "Vwf", "Tcf15", "Mki67"])]
genes = [g for _, gl in groups for g in gl if g in pb.index]
cmap = LinearSegmentedColormap.from_list("div", ["#184f95", BLUE, "#f0efec", ORANGE, "#a8431c"])
fig, axes = plt.subplots(1, 2, figsize=(8.6, 10.5), sharey=True)
for ax, ct in zip(axes, ["HSCs", "MPPs"]):
    m = pb[[f"{k}|{ct}" for k in ORDER + CTRL]]; m.columns = ORDER + CTRL
    lc = np.log2(m / m.sum() * 1e6 + 1).loc[genes]
    rel = pd.DataFrame({k: lc[k] - lc["M14" if RUN[k] == 1 else "M16"] for k in ORDER})
    im = ax.imshow(rel.to_numpy(), cmap=cmap, vmin=-1.5, vmax=1.5, aspect="auto")
    ax.set_xticks(range(5)); ax.set_xticklabels([LAB[k] for k in ORDER]); ax.xaxis.tick_top()
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=8, style="italic", color=INK2)
    ax.set_title(f"{ct}", fontsize=10, pad=22, loc="left")
    ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    y = 0
    for name, gl in groups:
        k = sum(g in pb.index for g in gl)
        if y: ax.axhline(y - 0.5, color=SURF, lw=3)
        if ax is axes[1]: ax.text(5.0, y + k / 2 - 0.5, name, rotation=270, va="center", fontsize=8.5, color=INK2)
        y += k
    ax.set_xticks(np.arange(-.5, 5, 1), minor=True); ax.set_yticks(np.arange(-.5, len(genes), 1), minor=True)
    ax.grid(which="minor", color=SURF, lw=1.2); ax.tick_params(which="minor", length=0)
cb = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.025, pad=0.035, aspect=40)
cb.set_label("log2 CPM difference from the same-run polyclonal control (pseudobulk of raw counts)", fontsize=8.5); cb.outline.set_visible(False)
fig.suptitle("Key genes by clonal system, slowest (left) to fastest (right)", x=0.02, ha="left", fontsize=12, y=0.995)
fig.text(0.02, 0.972, "M7.3 and M5.1 HSC columns rest on 26 and 25 cells and are noisy.", fontsize=8, color=INK2)
fig.savefig(fig_dir / "fig2_key_genes_heatmap.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---- Figure 3: LARRY clone-level association with output per HSC
r = pd.read_csv(root / "results/larry/larry_clone_associations.tsv", sep="\t")
r = r[r.min_hsc_cells == 5].set_index("set")
names = {"ISR_adaptive": "ISR adaptive arm", "ISR_terminal": "ISR terminal arm", "RiBi_nonRP": "Ribosome biogenesis",
         "RP_cytoplasmic": "Ribosomal protein genes", "mTORC1_hallmark": "mTORC1 hallmark", "MYC_targets_V1": "MYC targets V1",
         "Translation_initiation": "Translation initiation", "TOP_nonRP": "Non-RP 5'TOP", "UPR_hallmark_nonISR": "UPR (ISR genes removed)",
         "E2F_targets": "E2F targets (proliferation)", "Stem_MolO": "MolO stemness", "Stem_dHSC": "Dormant-HSC signature", "IEG_dissociation": "Immediate-early genes"}
keys = list(names)
fig, ax = plt.subplots(figsize=(7.8, 5.2))
yy = np.arange(len(keys))[::-1]
for y, k in zip(yy, keys):
    a, bb = r.rho_relout[k], r.rho_relout_adj_stemness[k]
    if not np.isnan(bb):
        ax.plot([a, bb], [y, y], color=GRID, lw=2, zorder=1)
        ax.plot(bb, y, "o", ms=8, color=ORANGE, mec=SURF, mew=2, zorder=3)
    ax.plot(a, y, "o", ms=8, color=BLUE, mec=SURF, mew=2, zorder=3)
ax.axvline(0, color=AXIS, lw=1); ax.axhline(yy[4] - 0.5, color=GRID, lw=1); ax.axhline(yy[8] - 0.5, color=GRID, lw=1)
ax.set_yticks(yy); ax.set_yticklabels([names[k] for k in keys], color=INK2); ax.tick_params(length=0)
ax.set_xlim(-0.7, 0.7); ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines["left"].set_visible(False)
ax.set_xlabel("Spearman correlation of clone-mean HSC score with output per HSC\n(negative: signature higher in low-output clones)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, ms=7, label="unadjusted"),
                   Line2D([], [], marker="o", ls="", color=ORANGE, ms=7, label="adjusted for MolO and dormant-HSC scores")],
          loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=8.5)
fig.suptitle("LARRY: signatures versus clonal output per HSC", x=0.02, ha="left", fontsize=12, y=0.995)
fig.text(0.02, 0.01, f"{int(r.n_clones.iloc[0])} barcoded clones with at least 5 scored sorted HSCs, three experiments, ranks taken within experiment.", fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 0.95)); fig.savefig(fig_dir / "fig3_larry_output.png", dpi=200); plt.close(fig)
print("figures written to", fig_dir)
