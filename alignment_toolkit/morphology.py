"""Orientation order, shared axis, aspect ratio and cell length for one group of ROIs (a quadrant or the whole frame).

All angles are the ROI's long-axis angle THETA_RAD, measured from the reference (flow) axis.
Cell axes have no head or tail (a cell at +80 deg and one at -100 deg lie on the same line), so the
averages use the doubled angle 2*theta.
"""
import numpy as np


def orientation_order(theta_rad):
    """Return (C, order_S, axis_deg).

    C = mean cos(2 theta), Sn = mean sin(2 theta)
    order_S = sqrt(C^2 + Sn^2): 0 = random, 1 = every cell on the same axis (any axis)
    axis_deg = 0.5 * atan2(Sn, C) in degrees: the shared axis, -90..+90 from the flow axis
    """
    t = np.asarray(theta_rad, dtype=float)
    t = t[~np.isnan(t)]
    if t.size == 0:
        return np.nan, np.nan, np.nan
    c = np.cos(2 * t).mean()
    s = np.sin(2 * t).mean()
    return float(c), float(np.hypot(c, s)), float(np.degrees(0.5 * np.arctan2(s, c)))


def shape_metrics(rdf):
    """Measures for one quadrant or frame; rdf needs THETA_RAD and optionally AR (major / minor)."""
    _, order_s, axis = orientation_order(rdf["THETA_RAD"].to_numpy())
    out = {
        "order_S": order_s,
        "axis_vs_flow_deg": abs(axis) if not np.isnan(axis) else np.nan,  # 0 = along flow, 90 = across
    }
    if "AR" in rdf.columns:
        ar = rdf["AR"].to_numpy(dtype=float)
        ar = ar[np.isfinite(ar)]
        out["AR_median"] = float(np.median(ar)) if ar.size else np.nan
        out["AR_mean"] = float(ar.mean()) if ar.size else np.nan
    else:
        out["AR_median"] = out["AR_mean"] = np.nan
    if "LENGTH" in rdf.columns:
        ln = rdf["LENGTH"].to_numpy(dtype=float)
        ln = ln[np.isfinite(ln)]
        out["length_median"] = float(np.median(ln)) if ln.size else np.nan
        out["length_mean"] = float(ln.mean()) if ln.size else np.nan
    else:
        out["length_median"] = out["length_mean"] = np.nan
    return out
