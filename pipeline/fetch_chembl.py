# pulls all approved small-molecule drugs from ChEMBL -> data/raw/
# needs: pip install chembl_webresource_client pandas

import os
import pandas as pd
from chembl_webresource_client.new_client import new_client

molecule = new_client.molecule

# max_phase 4 = approved. small molecules only (no antibodies etc)
approved = molecule.filter(
    max_phase=4,
    molecule_type="Small molecule"
).only([
    "molecule_chembl_id",
    "pref_name",
    "molecule_structures",
    "atc_classifications",
])

print("downloading from chembl, takes a few min...")

rows = []
for drug in approved:  # actual download happens here, client handles paging
    structures = drug["molecule_structures"]

    # a few entries don't have a structure, skip them
    if structures is None:
        continue

    rows.append({
        "chembl_id": drug["molecule_chembl_id"],
        "name": drug["pref_name"],
        "smiles": structures["canonical_smiles"],
        "inchi_key": structures["standard_inchi_key"],  # for matching to pubchem later
        "atc_codes": ";".join(drug["atc_classifications"]),  # can have more than one
    })

    if len(rows) % 500 == 0:
        print(len(rows), "so far")

df = pd.DataFrame(rows)

os.makedirs("data/raw", exist_ok=True)
df.to_csv("data/raw/chembl_approved_drugs.csv", index=False)

print("saved", len(df), "drugs")
