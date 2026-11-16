# downloads a 3D shape (SDF file) for every drug -> site/molecules/
# tries PubChem first, and if PubChem doesn't have one, builds it with RDKit
# needs: pip install requests rdkit pandas
# run from the main project folder: python pipeline/fetch_3d.py
# takes ~15 min. safe to stop and rerun, it skips drugs it already has

import os
import time
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

OUT_DIR = "site/molecules"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/SDF?record_type=3d"

drugs = pd.read_csv("data/processed/drugs_clean.csv", keep_default_na=False)
drugs = drugs[drugs["heavy_atoms"] >= 6].reset_index(drop=True)  # same filter as the notebooks
os.makedirs(OUT_DIR, exist_ok=True)


def from_pubchem(inchi_key):
    # returns the SDF text, or None if pubchem has no 3D shape for it
    try:
        r = requests.get(PUBCHEM.format(inchi_key), timeout=30)
    except requests.RequestException:
        return None
    if r.status_code != 200:
        return None
    # sometimes several records come back, keep just the first one
    return r.text.split("$$$$")[0] + "$$$$\n"


def from_rdkit(smiles):
    # build a 3D shape ourselves: add hydrogens -> embed in 3D -> relax with MMFF
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:  # 0 means it worked
        # big molecules often fail the normal way, so try again from random starting positions
        if AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True) != 0:
            return None
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol) + "$$$$\n"


sources = []
for i, drug in drugs.iterrows():
    path = os.path.join(OUT_DIR, drug["chembl_id"] + ".sdf")

    if os.path.exists(path):  # already done on an earlier run
        sources.append("already had")
        continue

    sdf = from_pubchem(drug["inchi_key"])
    source = "pubchem"
    time.sleep(0.25)  # pubchem asks for at most 5 requests a second

    if sdf is None:
        sdf = from_rdkit(drug["smiles"])
        source = "rdkit"

    if sdf is None:
        sources.append("failed")
        continue

    with open(path, "w") as f:
        f.write(sdf)
    sources.append(source)

    if (i + 1) % 100 == 0:
        print(i + 1, "of", len(drugs))

drugs["shape_source"] = sources
drugs[["chembl_id", "name", "shape_source"]].to_csv("data/processed/molecule_sources.csv", index=False)

print("done")
print(drugs["shape_source"].value_counts())