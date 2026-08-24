import re
from pathlib import Path

import numpy as np
import pandas as pd

from skimage.measure import regionprops_table


def load_masks_from_seg_npy(path):
    """Return the label array from a Cellpose *_seg.npy file."""
    data = np.load(path,allow_pickle=True)
    obj = data.item() if data.dtype == object else data
    if isinstance(obj, dict):
        if "masks" not in obj:
            raise KeyError(
                f"No 'masks' key in {path}. Keys present: {list(obj.keys())}"
            )
        return np.asarray(obj["masks"])
    return np.asarray(obj)

def masks_to_roi_table(masks, frame=1, reference_angle=0.0, min_area=0):
    """Fit an ellipse per label; return a per-ROI table (one row per ROI).

    Columns out: FRAME, ELLIPSE_THETA (radians, x-axis convention),
                 POSITION_X, POSITION_Y, MAJOR, MINOR, AREA, label.
    """
    masks = np.asarray(masks)

    # If the frame has no ROIs (max label 0), return an empty table
    # with the right columns so downstream concat doesn't break.
    if masks.max() == 0:
        return pd.DataFrame(
            columns=["FRAME", "ELLIPSE_THETA", "POSITION_X", "POSITION_Y",
                     "MAJOR", "MINOR", "AREA", "label"]
        )

    # Use skimage.measure.regionprops_table to get the ellipse parameters.
    props = regionprops_table(masks, properties=("label", 
                                                 "centroid", 
                                                 "orientation", 
                                                 "major_axis_length", 
                                                 "minor_axis_length", 
                                                 "area")
                                                 )
    df = pd.DataFrame(props)

    # Size filter
    if min_area > 0:
        df = df[df["area"] >= min_area].copy()

    theta_x = np.pi / 2.0 - df["orientation"].to_numpy()
    theta_out = theta_x - float(reference_angle)
    theta_out = (theta_out + np.pi / 2.0) % np.pi - np.pi / 2.0

    # Build and return the output DataFrame.
    return pd.DataFrame({
        "FRAME": int(frame),
        "ELLIPSE_THETA": theta_out,
        "POSITION_X": df["centroid-1"].to_numpy(),
        "POSITION_Y": df["centroid-0"].to_numpy(),
        "MAJOR": df["major_axis_length"].to_numpy(),
        "MINOR": df["minor_axis_length"].to_numpy(),
        "AREA": df["area"].to_numpy(),
        "label": df["label"].to_numpy(),
    })

def mask_file_to_table(path, frame=None, reference_angle=0.0, min_area=0):
    """Load one Cellpose mask file and return its per-ROI table."""
    masks = load_masks_from_seg_npy(path)
    if frame is None:
        frame = frame_from_name(path)
    return masks_to_roi_table(
        masks, frame=frame, reference_angle=reference_angle, min_area=min_area
    )

def _strip_suffix(path):
    """Remove extension and any cp_masks/seg/masks tag from a filename."""
    stem = Path(path).name
    stem = re.sub(r"\.(png|npy|tif|tiff)$", "", stem, flags=re.I)
    stem = re.sub(r"_(cp_masks|seg|masks)$", "", stem, flags=re.I)
    return stem

def frame_from_name(path):
    """Return the hour from a filename's trailing number ('..._001' -> 1)."""
    stem = _strip_suffix(path)
    match = re.search(r"(\d+)\s*$", stem)
    if match:
        return int(match.group(1))
    return 1

def find_mask_files(folder):
    """Return mask files in `folder`, ordered by hour. Prefer .npy over .png."""
    folder = Path(folder)
    files = sorted(folder.glob("*_seg.npy"))
    if not files:
        files = sorted(folder.glob("*_cp_masks.png"))
    files = sorted(files, key=frame_from_name)
    return files

def folder_to_dataset_csv(folder, out_dir, name=None,
                          reference_angle=0.0, min_area=0):
    """Turn one folder of ordered masks into a single dataset CSV.
    Returns (csv_path, n_hours, n_rois_total).
    """
    folder, out_dir = Path(folder), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = find_mask_files(folder)
    if not files:
        raise SystemExit(f"No Cellpose masks found in: {folder}")

    dataset = name if name else folder.name   

    tables = []
    for f in files:
        hour = frame_from_name(f)
        table = mask_file_to_table(f, frame=hour, reference_angle=reference_angle, min_area=min_area)
        tables.append(table)
        print(f"Processed {f.name} (hour {hour})")

    combined = pd.concat(tables, ignore_index=True)
    combined = combined.sort_values(["FRAME", "label"])
    combined = combined.reset_index(drop=True)

    out_path = out_dir / f"{dataset}.csv"
    combined.to_csv(out_path, index=False)
    return out_path, combined["FRAME"].nunique(), len(combined)