import re
from pathlib import Path

import numpy as np
import pandas as pd

from skimage.measure import regionprops_table


def load_masks_from_seg_npy(path):
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
    masks = np.asarray(masks)

    if masks.max() == 0:
        return pd.DataFrame(
            columns=["FRAME", "ELLIPSE_THETA", "POSITION_X", "POSITION_Y",
                     "MAJOR", "MINOR", "AREA", "label"]
        )

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
    masks = load_masks_from_seg_npy(path)
    if frame is None:
        frame = frame_from_name(path)
    return masks_to_roi_table(
        masks, frame=frame, reference_angle=reference_angle, min_area=min_area
    )

def _strip_suffix(path):
    stem = Path(path).name
    stem = re.sub(r"\.(png|npy|tif|tiff)$", "", stem, flags=re.I)
    stem = re.sub(r"_(cp_masks|seg|masks)$", "", stem, flags=re.I)
    return stem

def frame_from_name(path):
    stem = _strip_suffix(path)
    match = re.search(r"(\d+)\s*$", stem)
    if match:
        return int(match.group(1))
    return 1

def find_mask_files(folder):
    folder = Path(folder)
    files = sorted(folder.glob("*_seg.npy"))
    if not files:
        files = sorted(folder.glob("*_cp_masks.png"))
    files = sorted(files, key=frame_from_name)
    return files

def folder_to_dataset_csv(folder, out_dir, name=None,
                          reference_angle=0.0, min_area=0):
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