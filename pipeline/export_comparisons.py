# for every drug and each of its 5 closest look-alikes, works out:
#   - how many fingerprint "building blocks" they share, and how many they have in total
#     (shared / total = the similarity score)
#   - which atoms and bonds make up the biggest piece both molecules share
#     (so the website can highlight the shared part vs. the different part)

import json
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, rdFMCS

drugs = json.load(open("site/data/drugs.json"))   # same order as the website
clean = pd.read_csv("data/processed/drugs_clean.csv", keep_default_na=False)
smiles_by_id = dict(zip(clean["chembl_id"], clean["smiles"]))

smiles = [smiles_by_id[d["id"]] for d in drugs]
mols = [Chem.MolFromSmiles(s) for s in smiles]

# the same fingerprints used everywhere else, as sets of "on" bits
fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
bits = [set(fp_gen.GetFingerprint(m).GetOnBits()) for m in mols]

index_by_id = {d["id"]: i for i, d in enumerate(drugs)}


def shared_piece(i, j):
    # the biggest piece both molecules contain ("maximum common substructure").
    # timeout=1 stops after 1 second on the few very hard pairs and keeps the best piece found
    result = rdFMCS.FindMCS([mols[i], mols[j]], ringMatchesRingOnly=True, completeRingsOnly=True, timeout=1)
    if result.numAtoms == 0:
        return [], [], [], []
    piece = Chem.MolFromSmarts(result.smartsString)

    def atoms_and_bonds(mol):
        match = mol.GetSubstructMatch(piece)   # which atoms of this molecule are the shared piece
        bonds = []
        for bond in piece.GetBonds():
            b = mol.GetBondBetweenAtoms(match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()])
            if b is not None:
                bonds.append(b.GetIdx())
        return list(match), bonds

    atoms_i, bonds_i = atoms_and_bonds(mols[i])
    atoms_j, bonds_j = atoms_and_bonds(mols[j])
    return atoms_i, bonds_i, atoms_j, bonds_j


cache = {}   # drug A vs drug B is the same work as B vs A, so only do it once
pairs = []
for i, d in enumerate(drugs):
    rows = []
    for la in d["lookalikes"]:
        j = index_by_id[la["id"]]
        if (j, i) in cache:
            aj, bj, ai, bi = cache[(j, i)]      # stored the other way round
        else:
            ai, bi, aj, bj = shared_piece(i, j)
            cache[(i, j)] = (ai, bi, aj, bj)
        rows.append({
            "id": la["id"],
            "shared": len(bits[i] & bits[j]),   # building blocks both have
            "total": len(bits[i] | bits[j]),    # building blocks either has
            "a": ai, "ab": bi,                  # shared atoms and bonds in this drug
            "b": aj, "bb": bj,                  # shared atoms and bonds in the look-alike
        })
    pairs.append(rows)

    if (i + 1) % 200 == 0:
        print(i + 1, "of", len(drugs))

with open("site/data/comparisons.json", "w") as f:
    json.dump({"smiles": smiles, "pairs": pairs}, f, separators=(",", ":"))

print("done:", round(os.path.getsize("site/data/comparisons.json") / 1e6, 2), "MB")