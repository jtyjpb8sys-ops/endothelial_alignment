
from datetime import datetime
from pathlib import Path

import pandas as pd
from .prism import make_prism_table, make_prism_combined_table

from .config import (PRISM_TABLES, MIN_N_FOR_CI, DEFAULT_MATCH_DIST,
                     Z_REFERENCE_CONDITION, Z_REFERENCE_OXYGEN)
from .processing import combine_split_files
from .quadrants import assign_quadrants
from .summarise import summarise_frames
from .outputs import make_output_folders
from .tracking import add_movement
from .zscores import write_zscore_outputs


def analyse_folder(input_dir, output_dir=None, theta_units="radians",
                   ci=95.0, field_size=1992.0, match_dist=DEFAULT_MATCH_DIST,
                   z_condition=Z_REFERENCE_CONDITION, z_oxygen=Z_REFERENCE_OXYGEN):
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
    qc_rows = []                                   
    for csv_path in csv_files:
        dataset_name = csv_path.stem
        try:
            combined = combine_split_files([csv_path], theta_units)
        except ValueError as err:          # not a per-ROI CSV (e.g. a Prism / summary table)
            print(f"  {dataset_name}: not a per-ROI CSV, skipped ({str(err).splitlines()[0]})")
            continue
        if combined.empty:
            print(f"  {dataset_name}: no valid data, skipped.")
            continue
        combined, xmid, ymid = assign_quadrants(combined, field)
        combined, dup_frames = add_movement(combined, max_dist=match_dist)
        parts = dataset_name.split("_", 1)
        oxygen = parts[0]
        condition = parts[1] if len(parts) > 1 else ""
        summary = summarise_frames(combined, dataset_name, oxygen, condition, ci,
                                   duplicate_frames=dup_frames)
        summary.to_csv(
            folders["per_dataset"] / f"{dataset_name}_summary.csv", index=False
        )
        for folder_name, value_col in PRISM_TABLES.items():
            make_prism_table(summary, value_col).to_csv(
                folders[folder_name] / f"{dataset_name}_{folder_name}.csv",
                index=False,
            )

        make_prism_combined_table(
            summary, "AP_median", "AP_q3", "AP_q1"
        ).to_csv(folders["AP_Median"] / f"{dataset_name}_AP_median_IQR.csv",
                 index=False)
        make_prism_combined_table(
            summary, "AP_median", "AP_CI_high", "AP_CI_low"
        ).to_csv(folders["AP_Median"] / f"{dataset_name}_AP_median_CI.csv",
                 index=False)

        quad_counts = summary.loc[summary["REGION"] != "Qtotal", "n_ROIs"]
        total = summary[summary["REGION"] == "Qtotal"]
        frame_counts = total["n_ROIs"]
        qc_rows.append({
            "DATASET": dataset_name,
            "oxygen": oxygen,
            "condition": condition,
            "N_FRAMES": int(summary["FRAME"].nunique()),
            "FRAME_MIN": int(summary["FRAME"].min()),
            "FRAME_MAX": int(summary["FRAME"].max()),
            "MIN_ROIS_PER_QUADRANT_FRAME": int(quad_counts.min()),
            "MEDIAN_ROIS_PER_QUADRANT_FRAME": float(quad_counts.median()),
            "MAX_ROIS_PER_QUADRANT_FRAME": int(quad_counts.max()),
            "N_QUADRANT_FRAMES_BELOW_CI_MIN": int((quad_counts < MIN_N_FOR_CI).sum()),
            "MIN_ROIS_PER_FRAME": int(frame_counts.min()),
            "MEDIAN_ROIS_PER_FRAME": float(frame_counts.median()),
            "MEDIAN_CELLS_MATCHED_PER_FRAME": float(total["n_matched"].median()),
            "DUPLICATED_FRAMES": " ".join(str(f) for f in dup_frames),
        })

        all_summaries.append(summary)
        note = f" (duplicated frames: {dup_frames})" if dup_frames else ""
        print(f"  {dataset_name}: {summary['FRAME'].nunique()} hours summarised.{note}")

    if all_summaries:
        combined_all = pd.concat(all_summaries, ignore_index=True)
        combined_all.to_csv(folders["combined"] / "ALL_summaries_long.csv", index=False)
        # z-scores need several datasets (incl. the reference) in one folder; with a single
        # dataset only its all-frame means are written
        status = write_zscore_outputs(combined_all, output_dir, z_condition, z_oxygen)
        print(f"  {status}")

    if qc_rows:                                    
        pd.DataFrame(qc_rows).to_csv(
            folders["qc"] / "QC_summary.csv", index=False
        )

    return output_dir