#!/usr/bin/env python3

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from alignment_toolkit.batch import choose_folder_with_dialog
from alignment_toolkit.config import (
    DEFAULT_CI, DEFAULT_FIELD_SIZE, DEFAULT_THETA_UNITS,
    MIN_N_FOR_CI, PRISM_TABLES, REGIONS,
)
from alignment_toolkit.cellpose_extract import folder_to_dataset_csv
from alignment_toolkit.processing import combine_split_files
from alignment_toolkit.quadrants import assign_quadrants
from alignment_toolkit.summarise import summarise_frames
from alignment_toolkit.prism import make_prism_table
from alignment_toolkit.outputs import make_output_folders
from alignment_toolkit.analysis import analyse_folder
from alignment_toolkit.batch import run_batch, choose_folder_with_dialog

def run_analyse(args):
    output_dir = analyse_folder(
        args.input_dir, args.output_dir,
        theta_units=args.theta_units, ci=args.ci, field_size=args.field_size,
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

def run_batch_cmd(args):
    # If no --input_dir given, pop up the folder picker.
    cell_line = Path(args.input_dir) if args.input_dir else choose_folder_with_dialog()
    run_batch(
        cell_line,
        theta_units=args.theta_units, ci=args.ci, field_size=args.field_size,
        reference_angle=np.radians(args.reference_deg), min_area=args.min_area,
    )

def main():
    parser = argparse.ArgumentParser(
        description="Endothelial alignment toolkit."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_ex = sub.add_parser("extract", help="Cellpose masks -> per-ROI CSV.")
    p_ex.add_argument("--input_dir", required=True,
                      help="Folder of *_seg.npy mask files.")
    p_ex.add_argument("--output_dir", default=None,
                      help="Where to write the CSV (default: <input>/ROI_CSVs).")
    p_ex.add_argument("--name", default=None,
                      help="Dataset name, e.g. 21_4dyn_1.")
    p_ex.add_argument("--reference_deg", type=float, default=0.0,
                      help="Flow axis in degrees (default 0 = horizontal).")
    p_ex.add_argument("--min_area", type=int, default=0,
                      help="Drop ROIs smaller than this many pixels.")
    p_ex.set_defaults(func=run_extract)

    p_an = sub.add_parser("analyse", help="Per-ROI CSV -> quadrant summaries.")
    p_an.add_argument("--input_dir", required=True,
                      help="Folder containing the dataset CSV(s).")
    p_an.add_argument("--output_dir", default=None)

    p_ba = sub.add_parser("batch", help="Process a whole Cell line tree.")
    p_ba.add_argument("--input_dir", default=None,
                      help="Cell line folder. Omit to open a folder picker.")
    p_ba.add_argument("--theta_units", choices=["radians", "degrees"],
                      default=DEFAULT_THETA_UNITS)
    p_ba.add_argument("--ci", type=float, default=DEFAULT_CI)
    p_ba.add_argument("--field_size", type=float, default=DEFAULT_FIELD_SIZE)
    p_ba.add_argument("--reference_deg", type=float, default=0.0)
    p_ba.add_argument("--min_area", type=int, default=0)
    p_ba.set_defaults(func=run_batch_cmd)

    p_an.add_argument("--theta_units", choices=["radians", "degrees"],
                      default=DEFAULT_THETA_UNITS)
    p_an.add_argument("--ci", type=float, default=DEFAULT_CI)
    p_an.add_argument("--field_size", type=float, default=DEFAULT_FIELD_SIZE)

    p_an.set_defaults(func=run_analyse)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

