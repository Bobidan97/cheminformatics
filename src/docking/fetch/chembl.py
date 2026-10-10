from pathlib import Path
import pandas as pd
from chembl_webresource_client.new_client import new_client

from paths import RAW_DIR


COLUMNS = [
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "assay_type",
    "assay_chembl_id",
    "assay_description",
    "bao_label",
    "target_organism",
    "target_chembl_id",
]

def fetch_chembl_bioactivity(
    target_chembl_id: str,
    out_path: Path = RAW_DIR/"chembl_bioactivity.csv",
    standard_type: str = "IC50",
    standard_units: str = "nM",
) -> pd.DataFrame:
    """
    Pull bioactivity records for a target from ChEMBL and save as CSV.

    target_chembl_id: e.g. "CHEMBL1862" (ABL1)
    out_path: CSV file to write the raw records to. Defaults to data/raw/chembl_bioactivity.csv
    standard_type: activity type to filter on, e.g. "IC50"
    standard_units: restrict to this unit (default "nM")

    Returns the pulled records as a DataFrame (also written to out_path).
    """
    activity = new_client.activity
    query = activity.filter(
        target_chembl_id = target_chembl_id,
        standard_type = standard_type,
        standard_units = standard_units,
        standard_value__isnull = False,
    )

    records = list(query.only(COLUMNS))
    df = pd.DataFrame.from_records(records)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    return df