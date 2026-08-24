
import pandas as pd

from .config import METRICS, REGIONS
from .stats import median_ci, quartiles


def region_row(rdf, dataset_name, oxygen, condition, frame, region, ci):
    row = {
        "DATASET": dataset_name,
        "oxygen": oxygen,
        "condition": condition,
        "FRAME": int(frame),
        "REGION": region,
        "n_ROIs": len(rdf),
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
    return row

def summarise_frames(combined, dataset_name, oxygen, condition, ci):
    """Per-frame x region (Q1-Q4 + Qtotal) median, CI and IQR table."""
    rows = []
    for frame, fdf in combined.groupby("FRAME", sort=True):
        for region in REGIONS:
            rdf = fdf if region == "Qtotal" else fdf[fdf["QUADRANT"] == region]
            rows.append(
                region_row(rdf, dataset_name, oxygen, condition, frame, region, ci)
            )
    summary = pd.DataFrame(rows)
    summary["REGION"] = pd.Categorical(
        summary["REGION"], categories=REGIONS, ordered=True
    )
    return summary.sort_values(["FRAME", "REGION"]).reset_index(drop=True)