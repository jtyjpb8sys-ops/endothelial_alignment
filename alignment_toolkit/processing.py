
import numpy as np
import pandas as pd

from .columns import read_one_csv


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
    if "MAJOR" in combined.columns and "MINOR" in combined.columns:
        minor = combined["MINOR"].where(combined["MINOR"] > 0)
        combined["AR"] = combined["MAJOR"] / minor      # aspect ratio, 1 = round
        if "AREA" in combined.columns:
            # cell length from mask area and aspect ratio: 2 * sqrt(A * AR / pi)
            # (the major axis of an ellipse with the measured pixel area and aspect ratio)
            combined["LENGTH"] = 2.0 * np.sqrt(combined["AREA"] * combined["AR"] / np.pi)
    return combined