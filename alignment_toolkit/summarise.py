
import pandas as pd

from .config import METRICS, REGIONS
from .stats import median_ci, quartiles
from .morphology import shape_metrics
from .tracking import movement_metrics


def region_row(rdf, dataset_name, oxygen, condition, frame, region, ci,
               has_previous=True, duplicate=False):
    row = {
        "DATASET": dataset_name,
        "oxygen": oxygen,
        "condition": condition,
        "FRAME": int(frame),
        "REGION": region,
        "n_ROIs": len(rdf),
        "DUPLICATE_FRAME": bool(duplicate),
    }


    for label, column in METRICS:
        values = rdf[column].to_numpy()
        median, ci_low, ci_high = median_ci(values, ci)
        q1, q3 = quartiles(values)
        row.update({
            f"{label}_median": median,
            f"{label}_CI_low": ci_low,
            f"{label}_CI_high": ci_high,
            f"{label}_q1": q1,
            f"{label}_q3": q3,
            f"{label}_iqr": q3 - q1,
        })
    row.update(shape_metrics(rdf))
    row.update(movement_metrics(rdf, frame_has_previous=has_previous))
    return row

def summarise_frames(combined, dataset_name, oxygen, condition, ci,
                     duplicate_frames=()):
    """Per-frame x region (Q1-Q4 + Qtotal) table: AP median, CI and IQR, plus orientation order,
    shape and movement measures. Movement is empty for frame 0, for frames without a preceding
    frame and for duplicated frames (identical to the frame before)."""
    rows = []
    frames = set(combined["FRAME"].unique())
    dup = set(duplicate_frames)
    for frame, fdf in combined.groupby("FRAME", sort=True):
        has_prev = (frame - 1) in frames and frame not in dup
        for region in REGIONS:
            rdf = fdf if region == "Qtotal" else fdf[fdf["QUADRANT"] == region]
            rows.append(
                region_row(rdf, dataset_name, oxygen, condition, frame, region, ci,
                           has_previous=has_prev, duplicate=frame in dup)
            )
    summary = pd.DataFrame(rows)
    summary["REGION"] = pd.Categorical(
        summary["REGION"], categories=REGIONS, ordered=True
    )
    return summary.sort_values(["FRAME", "REGION"]).reset_index(drop=True)