DEFAULT_CI = 95.0
DEFAULT_THETA_UNITS = "radians"
DEFAULT_FIELD_SIZE = 1992.0 #square image edge length in pixels

MIN_N_FOR_CI = 6

# --- Regions --- #
REGIONS = ["Q1", "Q2", "Q3", "Q4", "Qtotal"]

# --- Metrics --- #
METRICS = [("AP", "AP"), ("theta_deg", "THETA_DEG_FOLDED")]

# --- Column names ---#
FRAME_NAMES = ["FRAME", "frame"]
THETA_NAMES = ["ELLIPSE_THETA", "ellipse_theta"]
X_NAMES = ["POSITION_X", "x"]
Y_NAMES = ["POSITION_Y", "y"]

# --- Prism export map ---#
PRISM_TABLES = {
    "AP_Median":       "AP_median",
    "AP_CI_low":       "AP_CI_low",
    "AP_CI_high":      "AP_CI_high",
    "ThetaDeg_Median": "theta_deg_median",
}