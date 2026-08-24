
from datetime import datetime
from pathlib import Path

import pandas as pd

from .config import PRISM_TABLES
from .processing import combine_split_files
from .quadrants import assign_quadrants
from .summarise import summarise_frames
from .prism import make_prism_table
from .outputs import make_output_folders


def analyse_folder(input_dir, output_dir=None, theta_units="radians",
                   ci=95.0, field_size=1992.0):
    input_dir = Path(input_dir)
    if not input_dir.exists():
        raise SystemExit(f"Input folder does not exist: {input_dir}")

    csv_files = sorted(input_dir.glob("*.csv"))
    if not csv_files:
        raise SystemExit(f"No CSV files found in: {input_dir}")

    output_dir = (
        Path(output_dir) if output_dir
        else input_dir / f"Endothelial_Alignment_{datetime.now():%d%m%Y}"
    )
    folders = make_output_folders(output_dir)
    field = {"xmid": field_size / 2.0, "ymid": field_size / 2.0}

    all_summaries = []
    for csv_path in csv_files:
        dataset_name = csv_path.stem
        combined = combine_split_files([csv_path], theta_units)
        if combined.empty:
            print(f"  {dataset_name}: no valid data, skipped.")
            continue
        combined, xmid, ymid = assign_quadrants(combined, field)
        parts = dataset_name.split("_", 1)
        oxygen = parts[0]
        condition = parts[1] if len(parts) > 1 else ""
        summary = summarise_frames(combined, dataset_name, oxygen, condition, ci)
        summary.to_csv(
            folders["per_dataset"] / f"{dataset_name}_summary.csv", index=False
        )
        for folder_name, value_col in PRISM_TABLES.items():
            make_prism_table(summary, value_col).to_csv(
                folders[folder_name] / f"{dataset_name}_{folder_name}.csv",
                index=False,
            )
        all_summaries.append(summary)
        print(f"  {dataset_name}: {summary['FRAME'].nunique()} hours summarised.")

    if all_summaries:
        pd.concat(all_summaries, ignore_index=True).to_csv(
            folders["combined"] / "ALL_summaries_long.csv", index=False
        )
    return output_dir