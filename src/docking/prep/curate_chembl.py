from pathlib import Path
import pandas as pd

from paths import RAW_DIR
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def curate_chembl_bioactivity(
    in_path: Path = RAW_DIR / "chembl_bioactivity.csv",
) -> pd.DataFrame:
    """
    Load raw ChEMBL bioactivity records and filter down to clean,
    comparable IC50 measurements from single-protein binding assays.

    in_path : CSV written by fetch_chembl_bioactivity

    Filters applied:
      - assay_type == "B"                     binding assay, not functional/ADME
      - bao_label == "single protein format"  excludes cell-based/complex formats
      - standard_relation == "="              exact value, not censored (>, <, >=)
      - standard_type == "IC50"
      - standard_units == "nM"                consistent scale for standard_value
      - canonical_smiles not null/blank       drops rows with missing structures

    Returns the filtered DataFrame.
    """
    df = pd.read_csv(in_path)

    mask = (
        (df["assay_type"] == "B") &
        (df["bao_label"] == "single protein format") &
        (df["standard_relation"] == "=") &
        (df["standard_type"] == "IC50") &
        (df["standard_units"] == "nM") &
        (df["canonical_smiles"].notna()) &
        (df["canonical_smiles"].str.strip() != "")
    )

    return df[mask].reset_index(drop=True)

uncharger = rdMolStandardize.Uncharger()
tautomer_enumerator = rdMolStandardize.TautomerEnumerator()

def strip_salts(mol):
    """Remove counterions/solvates, keep the largest fragment as parent."""
    return rdMolStandardize.FragmentParent(mol)

def normalize(mol):
    """Fix inconsistent depictions (nitro groups, charge placement, etc.)."""
    return rdMolStandardize.Normalize(mol)

def neutralize(mol):
    """Neutralize ionizable groups (e.g. carboxylate -> carboxylic acid)."""
    return uncharger.uncharge(mol)

def check_stereo_flags(mol):
    """Return (has_unspecified_stereocenter, has_unspecified_double_bond)."""
    Chem.AssignStereochemistry(mol, cleanIt = True, force = True)

    centers = Chem.FindMolChiralCenters(mol, includeUnassigned = True, useLegacyImplementation = False)
    has_unspecified_center = any(c[1] == "?" for c in centers)

    has_unspecified_bond = any(
        bond.GetStereo() == Chem.BondStereo.STEREOANY for bond in mol.GetBonds()
    )
    return has_unspecified_center, has_unspecified_bond

def canonical_tautomer(mol):
    """Canonicalize to one tautomer -- use for the dedup key, not for display."""
    return tautomer_enumerator.Canonicalize(mol)

def get_inchikey(mol):
    return Chem.MolToInchiKey(mol)

def standardize(smiles: str) -> dict:
    """Run one SMILES through the full pipeline. Returns a dict of results."""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"original_smiles": smiles, "parse_failed": True}

    parent = strip_salts(mol)
    parent = normalize(parent)
    parent = neutralize(parent)

    has_unspecified_center, has_unspecified_bond = check_stereo_flags(parent)

    tautomer = canonical_tautomer(parent)
    canonical_smiles = Chem.MolToSmiles(tautomer)
    inchikey = get_inchikey(tautomer)

    is_clean = not has_unspecified_center and not has_unspecified_bond

    return {
        "original_smiles": smiles,
        "standardized_smiles": canonical_smiles,
        "inchikey": inchikey,
        "is_clean": is_clean,
        "parse_failed": False,
    }

def finalize_compound_set(bioactivity_df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize every compound, drop anything unparseable or flagged
    unclean, then collapse duplicate measurements of the same
    compound (by InChIKey) into one row per compound.
    """
    results_df = pd.DataFrame(list(bioactivity_df["canonical_smiles"].apply(standardize)))
    combined = pd.concat([bioactivity_df.reset_index(drop=True), results_df], axis=1)
    combined = combined[~combined["parse_failed"] & combined["is_clean"]]

    return (
        combined.groupby("inchikey", as_index=False)
        .agg({
            "molecule_chembl_id": "first",
            "standardized_smiles": "first",
            "standard_value": "median",
            "standard_units": "first",
        })
    )