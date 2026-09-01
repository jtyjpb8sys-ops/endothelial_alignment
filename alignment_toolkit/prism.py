

import pandas as pd

from .config import REGIONS


def make_prism_table(summary, value_col):
    wide = (summary
            .pivot_table(index="FRAME", columns="REGION",
                         values=value_col, aggfunc="first", observed=False)
            .reindex(columns=REGIONS)
            .reset_index()
            .sort_values("FRAME"))
    wide.columns.name = None
    return wide

def make_prism_combined_table(summary, median_col, high_col, low_col, decimals=3):
    s = summary.copy()
    def fmt(row):
        med = row[median_col]
        hi = row[high_col]
        lo = row[low_col]
        if pd.isna(med):
            return ""
        med_str = f"{med:.{decimals}f}"
        hi_str = f"{hi:.{decimals}f}" if not pd.isna(hi) else ""
        lo_str = f"{lo:.{decimals}f}" if not pd.isna(lo) else ""
        return f"{med_str} [{lo_str}, {hi_str}]"

    s["COMBINED"] = s.apply(fmt, axis=1)

    wide = (s
            .pivot_table(index="FRAME", columns="REGION",
                         values="COMBINED", aggfunc="first", observed=False)
            .reindex(columns=REGIONS)
            .reset_index()
            .sort_values("FRAME"))
    wide.columns.name = None
    return wide