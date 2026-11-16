# packs the 3D shape of every drug into one file for the website -> site/data/molecules.json
# reads the .sdf files from fetch_3d.py (site/molecules/) and keeps only what the galaxy draws:
# each heavy atom's element and position, and which atoms are bonded
# needs: pip install rdkit
# run from the main project folder: python pipeline/export_molecules.py

import json
import os
from rdkit import Chem

drugs = json.load(open("site/data/drugs.json"))   # same order as the galaxy uses

molecules = []
missing = 0
for d in drugs:
    path = f"site/molecules/{d['id']}.sdf"
    mol = None
    if os.path.exists(path):
        # sanitize=False: just read the atoms and bonds, don't let RDKit reject odd files
        mol = next(iter(Chem.SDMolSupplier(path, removeHs=False, sanitize=False)), None)

    if mol is None or mol.GetNumConformers() == 0:
        molecules.append({"e": [], "p": [], "b": []})
        missing += 1
        continue

    pos = mol.GetConformer().GetPositions()   # x, y, z of every atom, in angstroms

    # keep heavy atoms only (skip hydrogens), and remember each one's new number
    keep = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    new_index = {old: new for new, old in enumerate(keep)}

    # center the molecule on 0,0,0 so the galaxy can place it anywhere
    center = pos[keep].mean(axis=0)
    elements = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in keep]
    coords = [round(float(v), 1) for i in keep for v in (pos[i] - center)]   # x,y,z,x,y,z,... rounded to 0.1 A

    bonds = []
    for b in mol.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a in new_index and c in new_index:
            bonds += [new_index[a], new_index[c]]

    molecules.append({"e": elements, "p": coords, "b": bonds})

with open("site/data/molecules.json", "w") as f:
    json.dump(molecules, f, separators=(",", ":"))

atoms = sum(len(m["e"]) for m in molecules)
bond_count = sum(len(m["b"]) // 2 for m in molecules)
print("molecules:", len(molecules), "| missing:", missing)
print("atoms:", atoms, "| bonds:", bond_count)
print("size:", round(os.path.getsize("site/data/molecules.json") / 1e6, 2), "MB")