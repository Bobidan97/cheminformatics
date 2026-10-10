from paths import RAW_DIR
from pathlib import Path

from rcsbapi.model import ModelQuery

def fetch_pdb_structures(
    pdb_ids: list[str],
    out_dir: Path = RAW_DIR / "pdb_structures",
    encoding_format: str = "cif",
) -> dict:
    """
    Download full coordinates for several PDB entries in one call,
    via the RCSB Model Server API's batch method.

    pdb_ids: accessions
    out_dir: directory to save into (created if missing). Defaults to data/raw/pdb_structures at the project root
    encoding_format: "cif" (default) or "bcif"

    Returns {pdb_id: local_path}.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_ids = [pid.strip().upper() for pid in pdb_ids]

    query = ModelQuery()
    results = query.get_multiple_structures(
        pdb_ids,
        query_type = "full",
        encoding = encoding_format,
        download = True,
        file_directory = str(out_dir),
    )
    return {pdb_id: Path(path) for pdb_id, path in results.items()}



