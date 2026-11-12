"""Frame-to-frame cell matching and movement measures, from X/Y centroids only.

For consecutive frames t-1 and t, every cell at t-1 is paired with its nearest cell at t and vice
versa. A pair is kept only if the two cells are each other's nearest neighbour (mutual) and closer
than max_dist pixels. Cells that leave the field, divide, are missed by Cellpose, or are ambiguous
stay unmatched.

The movement is stored on the cell at frame t (the end of the interval), so the value at hour t
describes movement from hour t-1 to hour t, and frame 0 has none. The quadrant a moving cell is
counted in is its position at hour t.

Pairs of identical consecutive frames (the same segmentation saved twice) are detected and left
empty rather than reported as zero movement.
"""
import numpy as np
from scipy.spatial import cKDTree


def duplicate_frames(combined):
    """Frames t whose cell positions are identical to frame t-1 (interval t-1 -> t is not real)."""
    def rows_sorted(xy):
        xy = xy.round(6)
        return xy[np.lexsort((xy[:, 1], xy[:, 0]))]   # sort whole (x, y) rows, not columns separately
    pos = {f: rows_sorted(g[["X", "Y"]].to_numpy()) for f, g in combined.groupby("FRAME")}
    frames = sorted(pos)
    dup = []
    for a, b in zip(frames[:-1], frames[1:]):
        if b == a + 1 and pos[a].shape == pos[b].shape and np.array_equal(pos[a], pos[b]):
            dup.append(b)
    return dup


def add_movement(combined, max_dist=35.0):
    """Add DX, DY, DIST, MATCHED columns to the per-ROI table (NaN / False where not matched).

    Returns (combined, duplicated_frames).
    """
    combined = combined.copy()
    combined["DX"] = np.nan
    combined["DY"] = np.nan
    combined["DIST"] = np.nan
    combined["MATCHED"] = False
    dup = set(duplicate_frames(combined))
    groups = {f: g for f, g in combined.groupby("FRAME")}
    for f, cur in groups.items():
        prev = groups.get(f - 1)
        if prev is None or f in dup or len(prev) == 0 or len(cur) == 0:
            continue
        a = prev[["X", "Y"]].to_numpy()          # hour t-1
        b = cur[["X", "Y"]].to_numpy()           # hour t
        d_ba, i_ba = cKDTree(a).query(b)         # nearest t-1 cell for every t cell
        d_ab, i_ab = cKDTree(b).query(a)         # nearest t cell for every t-1 cell
        mutual = (i_ab[i_ba] == np.arange(len(b))) & (d_ba < max_dist)
        idx = cur.index[mutual]
        delta = b[mutual] - a[i_ba[mutual]]
        combined.loc[idx, "DX"] = delta[:, 0]
        combined.loc[idx, "DY"] = delta[:, 1]
        combined.loc[idx, "DIST"] = d_ba[mutual]
        combined.loc[idx, "MATCHED"] = True
    return combined, sorted(dup)


def movement_metrics(rdf, frame_has_previous=True):
    """Movement measures for one quadrant or frame (rdf = ROIs at hour t with DX/DY/DIST/MATCHED).

    Values are means over matched cells (px per frame interval = px/h for hourly frames):
      speed_mean       mean distance moved
      net_dx_mean      mean x displacement: + = moving left -> right (with flow), - = against flow
      net_dy_mean      mean y displacement across the flow (check; ~0 without drift)
      direction_index  mean of dx / distance: +1 all moving with flow, -1 all moving against it
    """
    out = {"n_matched": np.nan, "speed_mean": np.nan,
           "net_dx_mean": np.nan, "net_dy_mean": np.nan, "direction_index": np.nan}
    if not frame_has_previous or "MATCHED" not in rdf.columns:
        return out
    m = rdf[rdf["MATCHED"]]
    out["n_matched"] = int(len(m))
    if len(m) == 0:
        return out
    dist = m["DIST"].to_numpy()
    moved = dist > 0
    out["speed_mean"] = float(dist.mean())
    out["net_dx_mean"] = float(m["DX"].mean())
    out["net_dy_mean"] = float(m["DY"].mean())
    out["direction_index"] = (float((m["DX"].to_numpy()[moved] / dist[moved]).mean())
                              if moved.any() else np.nan)
    return out
