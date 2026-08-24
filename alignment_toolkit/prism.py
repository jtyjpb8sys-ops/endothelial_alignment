
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
