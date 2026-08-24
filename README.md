# Endothelial Cell Alignment Workflow

A two-stage pipeline that turns Cellpose segmentation masks into per-quadrant
alignment statistics across a time series (one frame = one hour).

**Stage 1 — `extract`:** read Cellpose masks, fit an ellipse to every ROI, and
write one per-dataset CSV of orientation angles and centroids.

**Stage 2 — `analyse`:** read that CSV, compute the Alignment Parameter (AP),
split each frame into four equal-area quadrants plus the whole frame, and report
per-frame median, confidence interval and IQR.

For every frame you get five sets of outputs:

- `Q1`, `Q2`, `Q3`, `Q4` — the four equal squares of the imaging field
- `Qtotal` — all ROIs in the frame

Each of the five reports the Alignment Parameter (AP) and the folded orientation
angle, with: median, confidence interval, n ROIs, and interquartile range.

## Project layout

```
endothelial_alignment/
├── main.py                     # command-line entry point (extract / analyse)
└── alignment_toolkit/          # the package of modules
    ├── __init__.py
    ├── config.py               # defaults, region/metric definitions, column names
    ├── cellpose_extract.py     # masks -> ellipse angles -> per-ROI CSV
    ├── columns.py              # locate and read the needed CSV columns
    ├── processing.py           # combine CSVs and derive AP + folded angle
    ├── quadrants.py            # assign each ROI to a quadrant
    ├── stats.py                # median CI (order-statistic), quartiles
    ├── summarise.py            # per-frame x region summary table
    ├── prism.py                # reshape into wide GraphPad Prism tables
    └── outputs.py              # create the output folder layout
```

## Requirements

- Python 3.10+
- `numpy`, `pandas`, `scipy`, `scikit-image` (and `Pillow` only if reading
  `*_cp_masks.png` instead of `*_seg.npy`)

Install into your environment:

```
python -m pip install numpy pandas scipy scikit-image
```

## Quadrant definition

Each frame is a fixed square image (default 1992 x 1992 px). Each quadrant is
exactly one quarter of the field (996 x 996 px); `Qtotal` is the whole field.

| Quadrant | Position     | Condition            |
|----------|--------------|----------------------|
| Q1       | top-left     | x < 996,  y < 996    |
| Q2       | top-right    | x >= 996, y < 996    |
| Q3       | bottom-left  | x < 996,  y >= 996   |
| Q4       | bottom-right | x >= 996, y >= 996   |

X/Y coordinates follow image convention (x: left -> right, y: top -> bottom).
Set a different image size with `--field_size`.

## Measurements

The Alignment Parameter (AP) is calculated per ROI:

- `AP = cos(2 * theta) = 2 * cos^2(theta) - 1`
- AP ranges from +1 (parallel to the reference axis) to -1 (perpendicular)
- 0 indicates random / non-specific alignment

The reference axis is the flow direction. By default flow is assumed horizontal
(the image x-axis, `--reference_deg 0`). If flow runs at another angle, pass it
in degrees with `--reference_deg` during extraction and AP is measured relative
to that axis.

A folded orientation angle is also reported:

- `theta_deg` — orientation folded to 0-90 degrees, giving an absolute measure
  of deviation from the reference axis.

## Naming convention

Name each dataset `oxygen_flowcondition_replicate`, e.g. `21_4dyn_1`. During
analysis this is split on the first underscore into `oxygen` (`21`) and
`condition` (`4dyn_1`) for the summary and QC columns.

## Frame ordering (time series)

Each Cellpose mask file is one hour. Ordering comes from the trailing number in
the filename: `..._001` -> hour 1, `..._002` -> hour 2, and so on. Files are read
in ascending order and stacked into a single dataset, so a folder of hourly masks
becomes one CSV with consecutive `FRAME` values.

## Input files

Stage 1 expects a folder of Cellpose outputs:

- `*_seg.npy` (preferred — the exact label array), or
- `*_cp_masks.png` (used only if no `.npy` files are present)

Stage 2 expects the per-dataset CSV produced by stage 1, containing `FRAME`,
`ELLIPSE_THETA`, `POSITION_X`, `POSITION_Y`. Rows with missing values are
dropped automatically.

## Usage

Run from the `endothelial_alignment` folder (the one containing `main.py`).

**Stage 1 — extract angles from a folder of masks:**

```
python main.py extract --input_dir /path/to/masks --name 21_4dyn_1
```

Writes `/path/to/masks/ROI_CSVs/21_4dyn_1.csv`.

Extract options:

- `--input_dir` (required) — folder of `*_seg.npy` mask files
- `--name` — dataset name, e.g. `21_4dyn_1` (defaults to the folder name)
- `--output_dir` — where to write the CSV (default: `<input>/ROI_CSVs`)
- `--reference_deg` — flow axis in degrees (default 0 = horizontal)
- `--min_area` — drop ROIs smaller than this many pixels (default 0)

**Stage 2 — analyse the extracted CSV:**

```
python main.py analyse --input_dir /path/to/masks/ROI_CSVs
```

Analyse options:

- `--input_dir` (required) — folder containing the dataset CSV(s)
- `--output_dir` — where to write results (default: a dated folder beside the input)
- `--theta_units` — `radians` (default) or `degrees`
- `--ci` — confidence level, default 95
- `--field_size` — square image edge length in px (default 1992)

## Outputs

Under the results folder (`Endothelial_Alignment_<date>` by default):

- `Per_Dataset_Summaries/` — the full long-format summary per dataset
- `Prism/AP_Median`, `Prism/AP_CI_low`, `Prism/AP_CI_high`,
  `Prism/ThetaDeg_Median` — wide tables (frames down, quadrants across) for
  GraphPad Prism
- `Combined/ALL_summaries_long.csv` — every dataset stacked into one table
