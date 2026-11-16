#!/usr/bin/env python3
"""Round-3 (exploratory) gene sets from published sources, see docs/prespecification.md Amendment 3."""
import pathlib, re, sys, collections
import pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from xlsx_min import read_xlsx
root = pathlib.Path(__file__).resolve().parents[1]
raw = root / "data/raw"
RP = re.compile(r"^(Rp[ls][0-9]+[a-z]?[0-9]?|Rplp[0-2]|Rpsa)$")
rows = []
def add(name, genes, src):
    for g in dict.fromkeys(genes): rows.append(dict(set=name, tier="round3_exploratory", gene=g, source=src))
# Torrence 2021
d = read_xlsx(raw / "torrence2021_elife/elife-63326-fig1-data1-v2.xlsx")
sym = lambda k: [r for r in d[k][1:] if len(r) > 2 and r[2]]
tuni = {r[2]: float(r[4]) for r in sym("Tunicamycin Sensitive 2144 gene")}
adep = sym("ATF4 Dep Tuni Sens 774 genes")
shared = [r[2] for r in adep if len(r) > 6 and r[6].strip().upper() == "X"]
ins = {r[2] for r in sym("Ins Sensitive 1246 genes")}
insdep = {r[2] for r in sym("ATF4 Dep Ins Sens 112 genes")}
rap = [r[2] for r in sym("Rap Sensitive 253 genes")]
assert len(shared) == 61
def alias(gs):
    out = []
    for g in gs:
        out.append(g)
        if re.fullmatch(r"[A-Z][a-z]?ars", g) or g == "Eprs": out.append(g + "1")
        if g == "Ero1l": out.append("Ero1a")
    return out
src = "Torrence 2021 eLife 10:e63326 Fig1-source data 1"
add("T21_ATF4_shared", alias(shared), src)
add("T21_ATF4_stress_only", alias([r[2] for r in adep if r[2] not in set(shared) and r[2] not in ins and tuni.get(r[2], 0) >= 1]), src)
add("T21_mTORC1_nonATF4", alias([g for g in rap if g not in insdep and g not in set(shared)]), src)
# human Reactome -> mouse via MGI one-to-one homology classes
h = pd.read_csv(raw / "mgi/HOM_MouseHumanSequence.rpt", sep="\t", usecols=[0, 2, 3], names=["key", "tax", "symbol"], header=0)
cnt = h.groupby(["key", "tax"]).size().unstack(fill_value=0)
ok = cnt.index[(cnt.get(10090, 0) == 1) & (cnt.get(9606, 0) == 1)]
h = h[h.key.isin(ok)]
h2m = h[h.tax == 9606].set_index("key").symbol.to_frame("hs").join(h[h.tax == 10090].set_index("key").symbol.rename("mm")).set_index("hs").mm.to_dict()
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(raw / "msigdb_hs/c2.cp.reactome.v2025.1.Hs.symbols.gmt")}
for name, key in [("R_ATF4_ER_stress", "REACTOME_ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"),
                  ("R_PERK_gene_expression", "REACTOME_PERK_REGULATES_GENE_EXPRESSION"),
                  ("R_GCN2_aa_deficiency", "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY"),
                  ("R_HRI_heme_deficiency", "REACTOME_RESPONSE_OF_EIF2AK1_HRI_TO_HEME_DEFICIENCY")]:
    mm = [h2m[g] for g in gmt[key] if g in h2m]
    print(f"{name}: {len(gmt[key])} human, {len(mm)} mapped, {sum(bool(RP.match(g)) for g in mm)} RP removed")
    add(name, [g for g in mm if not RP.match(g)], f"MSigDB 2025.1.Hs {key}; MGI homology")
# mouse MSigDB biosynthesis sets
mg = {}
for f in ["m5.go.bp", "m2.cp.reactome", "mh.all"]:
    for l in open(raw / f"msigdb/{f}.v2025.1.Mm.symbols.gmt"):
        p = l.rstrip("\n").split("\t"); mg[p[0]] = p[2:]
for name, key, norp in [("AA_biosynthesis", "GOBP_AMINO_ACID_BIOSYNTHETIC_PROCESS", 0), ("Serine_biosynthesis", "GOBP_SERINE_FAMILY_AMINO_ACID_BIOSYNTHETIC_PROCESS", 0),
        ("One_carbon", "GOBP_ONE_CARBON_METABOLIC_PROCESS", 0), ("Purine_biosynthesis", "GOBP_PURINE_NUCLEOTIDE_BIOSYNTHETIC_PROCESS", 0),
        ("Pyrimidine_biosynthesis", "GOBP_PYRIMIDINE_NUCLEOTIDE_BIOSYNTHETIC_PROCESS", 0), ("PolI_transcription", "REACTOME_RNA_POLYMERASE_I_TRANSCRIPTION", 0),
        ("LSU_biogenesis", "GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS", 1), ("SSU_biogenesis", "GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS", 1),
        ("Cyto_translation_nonRP", "GOBP_CYTOPLASMIC_TRANSLATION", 1), ("Mito_translation", "REACTOME_MITOCHONDRIAL_TRANSLATION", 0),
        ("Cholesterol_biosynthesis", "REACTOME_CHOLESTEROL_BIOSYNTHESIS", 0), ("Fatty_acid_biosynthesis", "GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS", 0),
        ("OXPHOS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 0), ("Glycolysis", "HALLMARK_GLYCOLYSIS", 0)]:
    add(name, [g for g in mg[key] if not (norp and RP.match(g))], f"MSigDB 2025.1.Mm {key}" + (" minus RP genes" if norp else ""))
out = pd.DataFrame(rows); out.to_csv(root / "metadata/gene_sets_round3.tsv", sep="\t", index=False)
print(out.groupby("set").size().to_string())
