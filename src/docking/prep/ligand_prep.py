from pathlib import Path

import pandas as pd
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

from paths import PROCESSED_DIR
from paths import RAW_DIR

def protonate_at_ph(smiles: str, ph: float = 7.4) -> str:
    """Return a SMILES string with hydrogens added according to the specified pH."""

    mol = pybel.readstring("smi", smiles)
    mol.OBMol.AddHydrogens(False, True, ph)  # polaronly, correctForPH, pH

    return mol.write("smi").strip()


def prepare_ligand_pdbqt(smiles: str, out_path: Path, ph: float = 7.4) -> Path:
    """Prepare a ligand for docking and write it to a PDBQT file."""

    smiles_at_ph = protonate_at_ph(smiles, ph)
    mol = Chem.MolFromSmiles(smiles_at_ph)
    if mol is None:
        raise ValueError(f"Could not parse protonated SMILES: {smiles_at_ph!r}")

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise ValueError(f"3D embedding failed for {smiles_at_ph!r}")
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(setups[0])

    if not success:
        raise ValueError(f"PDBQT writing failed for {smiles_at_ph!r}: {error_msg}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(pdbqt_string)

    return out_path


def prepare_ligand_set(
        curated_df: pd.DataFrame,
        smiles_col: str,
        out_dir: Path = PROCESSED_DIR / "ligands",
        ph: float = 7.4,
) -> pd.DataFrame:
    """Prepare a set of ligands for docking, writing each to a PDBQT file."""

    out_dir.mkdir(parents=True, exist_ok=True)
    records = []

    for _, row in curated_df.iterrows():
        chembl_id = row["molecule_chembl_id"]
        smiles = row[smiles_col]
        out_path = out_dir / f"{chembl_id}.pdbqt"

        try:
            prepare_ligand_pdbqt(smiles, out_path, ph = ph)
            records.append({
                "molecule_chembl_id": chembl_id,
                "pdbqt_path": out_path,
                "prep_error": None,
            })
        except Exception as e:
            records.append({
                "molecule_chembl_id": chembl_id,
                "pdbqt_path": None,
                "prep_error": str(e),
            })

    return pd.DataFrame(records)

if __name__ == "__main__":
    curated_df = pd.read_csv(RAW_DIR / "chembl_bioactivity_curated.csv")
    results_df = prepare_ligand_set(curated_df, smiles_col = "standardized_smiles")  # match your actual column name
    results_df.to_csv(RAW_DIR / "ligand_prep_results.csv", index = False)

    print(results_df["prep_error"].notna().sum(), "failures out of", len(results_df))