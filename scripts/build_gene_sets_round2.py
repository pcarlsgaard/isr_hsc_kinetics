#!/usr/bin/env python3
"""Round-2 (exploratory) gene sets, see docs/prespecification.md Amendment 2."""
import pathlib, pandas as pd
root = pathlib.Path(__file__).resolve().parents[1]
gs = pd.read_csv(root / "data/processed/gene_sets.tsv", sep="\t")
S = {k: list(dict.fromkeys(v.gene)) for k, v in gs.groupby("set")}
stress = ["Ddit3", "Trib3", "Chac1", "Ppp1r15a", "Atf5", "Sesn2", "Ddit4", "Nupr1", "Atf3", "Gadd45a", "Stc2"]
drop = {"Sesn2", "Eif4ebp1", "Cebpg", "Stc2", "Vegfa", "Atf4"}
isr_all = set(S["ISR_adaptive"]) | set(S["ISR_terminal"]) | set(stress)
new = {
    "ISR_stress_specific": stress,
    "ISR_stress_specific_noIEG": [g for g in stress if g not in set(S["IEG_dissociation"])],
    "ATF4_biosynth": [g for g in S["ISR_adaptive"] if g not in drop],
    "Growth_noATF4": [g for g in S["mTORC1_hallmark"] if g not in isr_all],
}
rows = [dict(set=k, tier="round2_exploratory", gene=g) for k, v in new.items() for g in v]
pd.DataFrame(rows).to_csv(root / "metadata/gene_sets_round2.tsv", sep="\t", index=False)
for k, v in new.items(): print(k, len(v))
