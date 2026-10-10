from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent  # src/paths.py -> src -> root
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
STRUCTURES_DIR = DATA_DIR / "structures"
PROCESSED_DIR = DATA_DIR / "processed"