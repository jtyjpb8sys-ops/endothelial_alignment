
import numpy as np
import pandas as pd

from columns import read_one_csv
from cellpose_extract import frame_from_name   # for ordering, if needed


def combine_split_files(file_list, theta_units):
    dfs = []
    for file_path in file_list:
        df = read_one_csv(file_path)
        dfs.append(df)

    if not dfs:
        return pd.DataFrame()
    combined = pd.concat(dfs, ignore_index=True)

    if theta_units == "degrees":
        combined["THETA_RAD"] = np.radians(combined["THETA_RAW"])
    else:
        combined["THETA_RAD"] = combined["THETA_RAW"]
    combined["THETA_DEG"] = np.degrees(combined["THETA_RAD"])
    combined["AP"] = np.cos(2 * combined["THETA_RAD"])
    combined["THETA_DEG_FOLDED"] = np.degrees(
        np.arccos(np.sqrt((combined["AP"].clip(-1, 1) + 1) / 2))
    )
    return combined