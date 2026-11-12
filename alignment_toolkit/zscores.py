"""Per-dataset means and z-scores against a reference condition, written as plain tables.

For every dataset (oxygen_condition_replicate) and region (Q1-Q4, Qtotal), each measure is averaged over
all frames that have a value (movement measures start at frame 1; empty frames are skipped; frames
flagged DUPLICATE_FRAME, i.e. identical copies of the previous frame, are left out for every measure). Then, per
measure and region:

    z = (dataset mean - mean of the reference datasets) / SD of all datasets in the run

The reference datasets are those whose oxygen and condition match Z_REFERENCE_OXYGEN and
Z_REFERENCE_CONDITION (default 21 and Static). z = 0 is the reference average; +1 is one
dataset-to-dataset SD above it. Nothing else is calculated.
"""
import re

import numpy as np
import pandas as pd

from .config import Z_MEASURES

QUADRANTS = ["Q1", "Q2", "Q3", "Q4"]


def _split_name(name):
    """oxygen_condition_replicate -> (oxygen, condition, replicate); replicate None if absent."""
    m = re.match(r"^([^_]+)_(.+?)(?:_(\d+))?$", str(name))
    if not m:
        return str(name), "", None
    return m.group(1), m.group(2), (int(m.group(3)) if m.group(3) else None)


def dataset_means(summary):
    """Long table: one row per dataset x region with the all-frame mean of every Z measure."""
    cols = [c for c in Z_MEASURES if c in summary.columns]
    if "DUPLICATE_FRAME" in summary.columns:
        summary = summary[~summary["DUPLICATE_FRAME"].astype(bool)]
    means = (summary.groupby(["DATASET", "REGION"], observed=True)[cols]
             .mean()                     # pandas mean skips empty frames (NaN)
             .reset_index())
    parts = means["DATASET"].map(_split_name)
    means.insert(1, "oxygen", [p[0] for p in parts])
    means.insert(2, "condition", [p[1] for p in parts])
    means.insert(3, "replicate", [p[2] for p in parts])
    return means


def zscores(means, ref_condition, ref_oxygen):
    """Same layout as dataset_means, values replaced by z-scores. Returns (z, n_reference_datasets)."""
    cols = [c for c in Z_MEASURES if c in means.columns]
    is_ref = ((means["condition"].astype(str).str.lower() == str(ref_condition).lower()) &
              (means["oxygen"].astype(str) == str(ref_oxygen)))
    z = means.copy()
    for region, g in means.groupby("REGION", observed=True):
        ref = g[is_ref.loc[g.index]]
        for c in cols:
            sd = g[c].std(ddof=1)
            z.loc[g.index, c] = (g[c] - ref[c].mean()) / sd if sd and np.isfinite(sd) else np.nan
    n_ref = int(means.loc[is_ref, "DATASET"].nunique())
    return z, n_ref


def _ox_key(o):
    try:
        return -float(o)
    except ValueError:
        return 0.0


def prism_grouped(table, column, regions):
    """Grouped Prism table: rows = condition; columns = oxygen x replicate (x quadrant)."""
    t = table[table["REGION"].isin(regions)]
    conds = sorted(t["condition"].unique())
    oxs = sorted(t["oxygen"].unique(), key=_ox_key)
    reps = sorted(r for r in t["replicate"].unique() if r is not None and not pd.isna(r))
    out = pd.DataFrame(index=conds)
    for o in oxs:
        for r in reps:
            for q in regions:
                sel = t[(t["oxygen"] == o) & (t["replicate"] == r) & (t["REGION"] == q)]
                val = sel.set_index("condition")[column]
                lab = f"{o}%: rep{int(r)}" + (f" {q}" if q != "Qtotal" else "")
                out[lab] = val.reindex(conds).values
    out.index.name = "Condition"
    return out


def write_zscore_outputs(summary, output_dir, ref_condition, ref_oxygen):
    """Write dataset means and z-scores (long CSVs + Prism grouped tables) under output_dir.

    Returns a short status string.
    """
    means = dataset_means(summary)
    z, n_ref = zscores(means, ref_condition, ref_oxygen)
    zdir = output_dir / "Z_Scores"
    zdir.mkdir(parents=True, exist_ok=True)
    means.to_csv(zdir / "ALL_dataset_means_long.csv", index=False)
    if n_ref == 0:
        return (f"no reference datasets ({ref_oxygen}_{ref_condition}) found; "
                f"dataset means written, z-scores skipped")
    z.to_csv(zdir / "ALL_zscores_long.csv", index=False)
    for col, name in Z_MEASURES.items():
        if col not in means.columns:
            continue
        for kind, table in (("Dataset_Means", means), ("Z_Scores", z)):
            d = output_dir / "Prism" / kind / name
            d.mkdir(parents=True, exist_ok=True)
            prism_grouped(table, col, ["Qtotal"]).round(6).to_csv(d / f"{name}_{kind}_whole_frame.csv")
            prism_grouped(table, col, QUADRANTS).round(6).to_csv(d / f"{name}_{kind}_quadrants.csv")
    return f"z-scores against {n_ref} reference dataset(s) ({ref_oxygen}_{ref_condition})"
