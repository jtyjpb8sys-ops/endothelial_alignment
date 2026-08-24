#!/usr/bin/env python3

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from config import (
    DEFAULT_CI, DEFAULT_FIELD_SIZE, DEFAULT_THETA_UNITS,
    MIN_N_FOR_CI, PRISM_TABLES, REGIONS,
)

from cellpose_extract import folder_to_dataset_csv
from processing import combine_split_files
from quadrants import assign_quadrants
from summarise import summarise_frames
from prism import make_prism_table
from outputs import make_output_folders


def run_analyse(args):
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        raise SystemExit(f"Input folder does not exist: {input_dir}")

    csv_files = sorted(input_dir.glob("*.csv"))
    if not csv_files:
        raise SystemExit(f"No CSV files found in: {input_dir}")

    output_dir = (
        Path(args.output_dir) if args.output_dir
        else input_dir / f"Endothelial_Alignment_{datetime.now():%d%m%Y}"
    )
    folders = make_output_folders(output_dir)

    field = {"xmid": args.field_size / 2.0, "ymid": args.field_size / 2.0}

    all_summaries = []
    for csv_path in csv_files:
        dataset_name = csv_path.stem

        combined = combine_split_files([csv_path], args.theta_units)
        if combined.empty:
            print(f"  {dataset_name}: no valid data, skipped.")
            continue

        combined, xmid, ymid = assign_quadrants(combined, field)
        parts = dataset_name.split("_", 1)
        oxygen = parts[0]
        condition = parts[1] if len(parts) > 1 else ""
        summary = summarise_frames(
            combined, dataset_name, oxygen, condition, args.ci
        )

        summary.to_csv(
            folders["per_dataset"] / f"{dataset_name}_summary.csv", index=False
        )

        for folder_name, value_col in PRISM_TABLES.items():
            prism_table = make_prism_table(summary, value_col)
            prism_table.to_csv(
                folders[folder_name] / f"{dataset_name}_{folder_name}.csv", index=False
            )
        all_summaries.append(summary)
        print(f"  {dataset_name}: {summary['FRAME'].nunique()} hours summarised.")

    if all_summaries:
        pd.concat(all_summaries, ignore_index=True).to_csv(
            folders["combined"] / "ALL_summaries_long.csv", index=False
        )
    print(f"\nDone. Outputs in: {output_dir}")

def run_extract(args):
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        raise SystemExit(f"Input folder does not exist: {input_dir}")

    out_dir = Path(args.output_dir) if args.output_dir else input_dir / "ROI_CSVs"
    reference_angle = np.radians(args.reference_deg)

    csv_path, n_hours, n_rois = folder_to_dataset_csv(
        input_dir, out_dir, name=args.name,
        reference_angle=reference_angle, min_area=args.min_area
    )

    print(f"\nWrote {csv_path}")
    print(f"  {n_hours} hour(s), {n_rois} ROIs total")
    print(f"\nNext:  python main.py analyse --input_dir {out_dir}")



def main():
    print("Hello from endothelial-alignment!")


if __name__ == "__main__":
    main()

