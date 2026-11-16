# bundles everything the website needs -> site/data/
#   drugs.json  : one entry per drug (name, groups, 3D positions, cluster, closest look-alikes)
#   scores.json : bar chart data (times more than chance, per group) + axis labels
# needs: data/processed/coords.csv (layout.ipynb) and group_scores.csv (analysis.ipynb)
# run from the main project folder: python pipeline/export.py

import os
import json
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.decomposition import PCA

N_LOOKALIKES = 5
os.makedirs("site/data", exist_ok=True)

drugs = pd.read_csv("data/processed/drugs_clean.csv", keep_default_na=False)
drugs = drugs[drugs["heavy_atoms"] >= 6].reset_index(drop=True)
coords = pd.read_csv("data/processed/coords.csv", keep_default_na=False)
scores = pd.read_csv("data/processed/group_scores.csv")

# coords.csv was made from the same filtered list, so the rows should line up exactly
assert list(coords["chembl_id"]) == list(drugs["chembl_id"]), "coords.csv doesn't match drugs_clean.csv, rerun layout.ipynb"

# fingerprints again, to find each drug's closest look-alikes for the molecule panel
fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
mols = [Chem.MolFromSmiles(s) for s in drugs["smiles"]]
fps = [fp_gen.GetFingerprint(m) for m in mols]

# PCA axis labels (same settings as layout.ipynb, so the numbers match)
X = np.array([fp_gen.GetFingerprintAsNumPy(m) for m in mols])
pca = PCA(n_components=3, random_state=42).fit(X)
pca_labels = [f"PC{i+1} ({share:.1%} of variation)" for i, share in enumerate(pca.explained_variance_ratio_)]

print("finding look-alikes...")
entries = []
for i in range(len(drugs)):
    sims = np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps))
    sims[i] = -1  # skip itself
    top = np.argsort(sims)[::-1][:N_LOOKALIKES]

    c = coords.loc[i]
    groups = drugs.loc[i, "atc_groups"]
    entries.append({
        "id": drugs.loc[i, "chembl_id"],
        "name": drugs.loc[i, "name"].title(),
        "groups": groups.split(";") if groups != "" else [],
        "primary": c["primary_group"],
        "cluster": int(c["cluster"]),
        "pca": [round(float(c["pca_x"]), 3), round(float(c["pca_y"]), 3), round(float(c["pca_z"]), 3)],
        "tsne": [round(float(c["tsne_x"]), 3), round(float(c["tsne_y"]), 3), round(float(c["tsne_z"]), 3)],
        "lookalikes": [{"id": drugs.loc[j, "chembl_id"], "sim": round(float(sims[j]), 2)} for j in top],
        "has3d": os.path.exists(f"site/molecules/{drugs.loc[i, 'chembl_id']}.sdf"),
    })

with open("site/data/drugs.json", "w") as f:
    json.dump(entries, f, separators=(",", ":"))  # no extra spaces = smaller file

scores_out = {
    "pca_axis_labels": pca_labels,
    "groups": [
        {
            "group": row["group"],
            "name": row["name"],
            "n_drugs": int(row["n_drugs"]),
            "real": round(float(row["real"]), 3),
            "chance": round(float(row["chance"]), 3),
            "times_chance": round(float(row["times_chance"]), 2),
        }
        for _, row in scores.iterrows()
    ],
}
with open("site/data/scores.json", "w") as f:
    json.dump(scores_out, f, indent=2)

print("drugs.json:", len(entries), "drugs,", round(os.path.getsize("site/data/drugs.json") / 1e6, 2), "MB")
print("with a 3D shape:", sum(e["has3d"] for e in entries))
print("scores.json:", len(scores_out["groups"]), "groups")
print("PCA axes:", pca_labels)