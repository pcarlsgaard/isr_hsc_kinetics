#!/usr/bin/env python3
"""Assemble the pre-specified gene sets into data/processed/gene_sets.tsv.

Inputs: metadata/gene_sets_curated.tsv, metadata/gene_sets_msigdb.tsv, MSigDB mouse GMTs
(2025.1.Mm) and the stemness sets shipped in the Zenodo manuscript bundle.
RP_cytoplasmic is defined by regex and resolved against each dataset's features at scoring time.
"""
import csv, re, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
RP = re.compile(r"^(Rp[ls][0-9]+[a-z]?|Rplp[0-2]|Rpsa)$")
gmt = {}
for f in (root / "data/raw/msigdb").glob("*.gmt"):
    coll = f.name.split(".v2025")[0]
    for line in open(f):
        p = line.rstrip("\n").split("\t")
        gmt[(coll, p[0])] = p[2:]
cur = {}
tier = {}
for r in csv.DictReader(open(root / "metadata/gene_sets_curated.tsv"), delimiter="\t"):
    if r["gene"] == "__REGEX__":
        tier[r["set"]] = r["tier"]; cur[r["set"]] = ["__REGEX__"]; continue
    cur.setdefault(r["set"], []).append(r["gene"]); tier[r["set"]] = r["tier"]
isr = set(cur["ISR_adaptive"]) | set(cur["ISR_terminal"])
for r in csv.DictReader(open(root / "metadata/gene_sets_msigdb.tsv"), delimiter="\t"):
    g = gmt[(r["collection"], r["msigdb_name"])]
    if "minus RP" in r["modification"]:
        g = [x for x in g if not RP.match(x)]
    if "minus ISR" in r["modification"]:
        g = [x for x in g if x not in isr]
    cur[r["set"]] = g; tier[r["set"]] = r["tier"]
# author-shipped stemness / quiescence sets (controls)
hs = root / "data/raw/zenodo_19486868/extracted/manuscript/suppFig7/data/Genesets_HSCs.csv"
rows = list(csv.reader(open(hs, encoding="utf-8-sig"), delimiter=";"))
want = {"MolO_Wilson_et_al": "Stem_MolO", "dHSC_vs_aHSC_Cabezas-Wallscheid_et_al": "Stem_dHSC",
        "Low-output_Rodriguez_et_al": "Stem_LowOutput", "High-output_Rodriguez_et_al": "HighOutput"}
for j, name in enumerate(rows[0]):
    if name in want:
        genes = [r[j] for r in rows[1:] if j < len(r) and r[j] and not r[j].startswith("ENSMUSG")]
        if "Rodriguez" in name:  # ranked lists; keep top 100 as in the original notebook
            genes = genes[:100]
        cur[want[name]] = list(dict.fromkeys(genes)); tier[want[name]] = "control"
out = root / "data/processed/gene_sets.tsv"
with open(out, "w") as f:
    f.write("set\ttier\tgene\n")
    for s, g in cur.items():
        for x in dict.fromkeys(g):
            f.write(f"{s}\t{tier[s]}\t{x}\n")
for s, g in cur.items():
    print(f"{s:26s} {tier[s]:9s} n={len(set(g))}")
