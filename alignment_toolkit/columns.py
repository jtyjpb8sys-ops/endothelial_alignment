
import pandas as pd

from .config import FRAME_NAMES, THETA_NAMES, X_NAMES, Y_NAMES


def find_column(df, possible_names, required=True):
    normalized = {
        str(col).strip().lower().replace(" ", "_"): col for col in df.columns
    }

    for name in possible_names:
        key = str(name).strip().lower().replace(" ", "_")
        if key in normalized:
            return normalized[key]

    if required:
        raise ValueError(
            f"Could not find any of these columns: {possible_names}\n"
            f"Available columns were: {list(df.columns)}"
        )
    return None


def read_one_csv(file_path):
    df = pd.read_csv(file_path)

    frame_col = find_column(df, FRAME_NAMES)
    theta_col = find_column(df, THETA_NAMES)
    x_col = find_column(df, X_NAMES)
    y_col = find_column(df, Y_NAMES)

    out = df[[frame_col, theta_col, x_col, y_col]].copy()
    out.columns = ["FRAME", "THETA_RAW", "X", "Y"]

    for c in ["FRAME", "THETA_RAW", "X", "Y"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.dropna(subset=["FRAME", "THETA_RAW", "X", "Y"])

    out["FRAME"] = out["FRAME"].astype(int)
    return out