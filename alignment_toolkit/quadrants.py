
import numpy as np

def assign_quadrants(combined, field):
    x, y = combined["X"], combined["Y"]

    xmid = field["xmid"] if field["xmid"] is not None else (x.min() + x.max()) / 2.0
    ymid = field["ymid"] if field["ymid"] is not None else (y.min() + y.max()) / 2.0

    left = x < xmid
    top = y < ymid

    quad = np.where(top & left, "Q1",
           np.where(top & ~left, "Q2",
           np.where(~top & left, "Q3", "Q4")))

    combined = combined.copy()
    combined["QUADRANT"] = quad
    return combined, xmid, ymid