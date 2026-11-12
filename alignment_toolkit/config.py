DEFAULT_CI = 95.0
DEFAULT_THETA_UNITS = "radians"
DEFAULT_FIELD_SIZE = 1992.0 #square image edge length in pixels

MIN_N_FOR_CI = 6

# --- Cell matching between frames (movement measures) --- #
DEFAULT_MATCH_DIST = 35.0   # px; a matched pair must be closer than this (mutual nearest neighbours)

# --- Regions --- #
REGIONS = ["Q1", "Q2", "Q3", "Q4", "Qtotal"]

# --- Metrics --- #
METRICS = [("AP", "AP"), ("theta_deg", "THETA_DEG_FOLDED")]

# --- Column names ---#
FRAME_NAMES = ["FRAME", "frame"]
THETA_NAMES = ["ELLIPSE_THETA", "ellipse_theta"]
X_NAMES = ["POSITION_X", "x"]
Y_NAMES = ["POSITION_Y", "y"]
MAJOR_NAMES = ["MAJOR", "major"]
MINOR_NAMES = ["MINOR", "minor"]

# --- Prism export map ---#
PRISM_TABLES = {
    "AP_Median":       "AP_median",
    "AP_CI_low":       "AP_CI_low",
    "AP_CI_high":      "AP_CI_high",
    "ThetaDeg_Median": "theta_deg_median",
    # orientation order and shape (morphology.py)
    "Order_S":            "order_S",
    "Axis_vs_Flow_Deg":   "axis_vs_flow_deg",
    "AspectRatio_Median": "AR_median",
    "AspectRatio_Mean":   "AR_mean",
    # movement between frames (tracking.py)
    "Cells_Matched":      "n_matched",
    "Speed_Mean":         "speed_mean",
    "NetDX_Mean":         "net_dx_mean",
    "NetDY_Mean":         "net_dy_mean",
    "Direction_Index":    "direction_index",
}

# --- Z-scores (zscores.py) --- #
# Each dataset's mean over all frames, per region, is z-scored against the reference datasets:
#   z = (dataset mean - mean of reference datasets) / SD of all datasets in the run
Z_REFERENCE_CONDITION = "Static"   # condition part of the dataset name (oxygen_condition_replicate)
Z_REFERENCE_OXYGEN = "21"          # oxygen part of the dataset name
Z_MEASURES = {                     # summary column : output name
    "AP_median":        "AP_Median",
    "order_S":          "Order_S",
    "axis_vs_flow_deg": "Axis_vs_Flow_Deg",
    "AR_median":        "AspectRatio_Median",
    "AR_mean":          "AspectRatio_Mean",
    "n_matched":        "Cells_Matched",
    "speed_mean":       "Speed_Mean",
    "net_dx_mean":      "NetDX_Mean",
    "net_dy_mean":      "NetDY_Mean",
    "direction_index":  "Direction_Index",
}