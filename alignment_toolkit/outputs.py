
from .config import PRISM_TABLES


def make_output_folders(output_dir):
    folders = {
        "per_dataset": output_dir / "Per_Dataset_Summaries",
        "combined": output_dir / "Combined",
        "qc": output_dir / "QC",
    }
   
    folders.update({name: output_dir / "Prism" / name for name in PRISM_TABLES})

    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)

    return folders