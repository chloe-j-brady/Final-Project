# cleans the raw chembl file -> data/processed/drugs_clean.csv
# strips salts, drops inorganics, merges duplicates

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

raw = pd.read_csv("data/raw/chembl_approved_drugs.csv")
print("raw drugs:", len(raw))

uncharger = rdMolStandardize.Uncharger()

rows = []
bad_smiles = 0
inorganic = 0

for _, drug in raw.iterrows():
    mol = Chem.MolFromSmiles(drug["smiles"])
    if mol is None:  # rdkit couldn't read it
        bad_smiles += 1
        continue

    is_salt = "." in drug["smiles"]

    # split into separate pieces (drug + counter-ions) and keep the biggest one
    pieces = Chem.GetMolFrags(mol, asMols=True)
    parent = max(pieces, key=lambda p: p.GetNumHeavyAtoms())

    # undo charges left over from the salt, e.g. COO- back to COOH
    parent = uncharger.uncharge(parent)

    # no carbon = inorganic (aluminum chloride etc), drop it
    has_carbon = any(atom.GetSymbol() == "C" for atom in parent.GetAtoms())
    if not has_carbon:
        inorganic += 1
        continue

    rows.append({
        "chembl_id": drug["chembl_id"],
        "name": drug["name"],
        "smiles": Chem.MolToSmiles(parent),
        "inchi_key": Chem.MolToInchiKey(parent),  # new key for the stripped molecule
        "atc_codes": drug["atc_codes"],
        "is_salt": is_salt,
        "heavy_atoms": parent.GetNumHeavyAtoms(),
    })

df = pd.DataFrame(rows)
print("couldn't parse:", bad_smiles)
print("inorganic dropped:", inorganic)

# merge duplicates: after stripping, a salt and its parent drug have the same inchi_key
# put non-salt rows first so their name/id is the one we keep
df = df.sort_values("is_salt")


def combine_codes(codes):
    # all ATC codes from every version of the drug, no repeats
    all_codes = set()
    for c in codes.dropna():
        all_codes.update(c.split(";"))
    return ";".join(sorted(all_codes))


merged = df.groupby("inchi_key").agg(
    chembl_id=("chembl_id", "first"),
    name=("name", "first"),
    smiles=("smiles", "first"),
    heavy_atoms=("heavy_atoms", "first"),
    atc_codes=("atc_codes", combine_codes),
    n_versions=("chembl_id", "count"),
).reset_index()

print("duplicates merged:", len(df) - len(merged))

# top-level groups = first letter of each ATC code, e.g. "N02BA01" -> "N"
def get_groups(codes):
    if codes == "":
        return ""
    return ";".join(sorted(set(c[0] for c in codes.split(";"))))


merged["atc_groups"] = merged["atc_codes"].apply(get_groups)
merged["n_groups"] = merged["atc_groups"].apply(lambda g: 0 if g == "" else len(g.split(";")))

os.makedirs("data/processed", exist_ok=True)
merged.to_csv("data/processed/drugs_clean.csv", index=False)

print("final drugs:", len(merged))
print("  with a class:", (merged["n_groups"] > 0).sum())
print("  unclassified:", (merged["n_groups"] == 0).sum())
print("  in more than one group:", (merged["n_groups"] > 1).sum())