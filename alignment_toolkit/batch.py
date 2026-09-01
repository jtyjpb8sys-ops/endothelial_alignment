
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import re

import pandas as pd

from .cellpose_extract import folder_to_dataset_csv
from .analysis import analyse_folder


def choose_folder_with_dialog(title="Select the Cell line folder"):
    root = tk.Tk()
    root.withdraw()                    
    root.attributes("-topmost", True)  
    selected = filedialog.askdirectory(title=title)
    root.destroy()

    if not selected:
        raise SystemExit("No folder selected.")
    return Path(selected)

def find_replicate_folders(cell_line_dir):
    cell_line_dir = Path(cell_line_dir)
    replicates = []
    for oxygen_dir in sorted(cell_line_dir.iterdir()):
        if not oxygen_dir.is_dir():
            continue
        for flow_dir in sorted(oxygen_dir.iterdir()):
            if not flow_dir.is_dir():
                continue
            for rep_dir in sorted(flow_dir.iterdir()):
                if not rep_dir.is_dir():
                    continue
                masks = (list(rep_dir.glob("*_seg.npy"))
                     or list(rep_dir.glob("*_cp_masks.png")))
                if not masks:
                    continue
                replicates.append(
                    (oxygen_dir.name, flow_dir.name, rep_dir.name, rep_dir)
                )
    return replicates

def dataset_name_from_path(oxygen, flow, replicate):
    ox = re.sub(r"(?i)oxygen", "", oxygen)
    ox = re.sub(r"[^0-9.]", "", ox)

    cond = re.sub(r"\s+", "_", flow.strip())

    rep = re.sub(r"[^0-9]", "", replicate)

    return f"{ox}_{cond}_{rep}"

def run_batch(cell_line_dir, theta_units="radians", ci=95.0,
              field_size=1992.0, reference_angle=0.0, min_area=0):
    cell_line_dir = Path(cell_line_dir)
    replicates = find_replicate_folders(cell_line_dir)
    if not replicates:
        raise SystemExit(f"No replicate folders with masks under {cell_line_dir}")
    print(f"Found {len(replicates)} replicate(s) to process.\n")

    qc_all = []
    for oxygen, flow, replicate, rep_dir in replicates:

        name = dataset_name_from_path(oxygen, flow, replicate)
        print(f"=== {name}   ({rep_dir.name}) ===")

        out_dir = rep_dir / "ROI_CSVs"
        csv_path, n_hours, n_rois = folder_to_dataset_csv(
            rep_dir, out_dir, name=name,
            reference_angle=reference_angle, min_area=min_area,
        )
        print(f"  extracted {n_rois} ROIs over {n_hours} hour(s)")

        output_dir = analyse_folder(out_dir, theta_units=theta_units, ci=ci,
                                    field_size=field_size)
        print(f"  analysed -> {out_dir}\n")

        qc_path = output_dir / "QC" / "QC_summary.csv"
        if qc_path.exists():
            qc_all.append(pd.read_csv(qc_path))

    if qc_all:
        combined_qc = pd.concat(qc_all, ignore_index=True)
        combined_qc.to_csv(cell_line_dir / "ALL_QC_summary.csv", index=False)
        print(f"Combined QC for {len(qc_all)} datasets -> "
              f"{cell_line_dir / 'ALL_QC_summary.csv'}")