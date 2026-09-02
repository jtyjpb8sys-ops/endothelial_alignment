# Endothelial Cell Alignment Workflow

A two-stage pipeline that turns Cellpose segmentation masks into per-quadrant
alignment statistics across a time series. Each mask is one hourly timepoint, and
frames are numbered from the baseline: `FRAME` 0 = baseline, `FRAME` 1 = 1 hour,
and so on.

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
├── main.py                     # command-line entry point (extract / analyse / batch)
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
    ├── analysis.py             # shared analyse logic (CSV folder -> summaries)
    ├── outputs.py              # create the output folder layout
    └── batch.py                # folder picker + walk a cell-line tree
```

## Requirements

- Python 3.10+
- `numpy`, `pandas`, `scipy`, `scikit-image` (and `Pillow` only if reading
  `*_cp_masks.png` instead of `*_seg.npy`)
- `tkinter` for the batch folder picker — ships with Python, not installed via
  pip; the picker is only needed when running `batch` without `--input_dir`

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
to that axis. For example, if flow runs vertically in the images, use
`--reference_deg 90`.

The reference angle is applied at extract time (it is baked into the stored
angles), so different orientations are written to separate folders and never
overwrite each other: `--reference_deg 0` writes to `ROI_CSVs/`, while
`--reference_deg 90` writes to `ROI_CSVs_ref90/`. This lets you extract the same
data both ways and compare. Because the choice lives in the folder name, an
`analyse`-only run must be given the same `--reference_deg` so it reads from the
matching folder — at that stage the flag only selects the folder, it does not
change any maths.

A folded orientation angle is also reported:

- `theta_deg` — orientation folded to 0-90 degrees, giving an absolute measure
  of deviation from the reference axis.

## Naming convention

Datasets are named `oxygen_condition_replicate`, e.g. `21_4dyn_laminar_1`. In
batch mode the name is built automatically from the folder path: the oxygen
folder is reduced to its number (`Oxygen 21%` -> `21`), the flow folder name is
kept verbatim with spaces turned into underscores (`4dyn laminar` ->
`4dyn_laminar`, `4dyn 1hz oscilatory` -> `4dyn_1hz_oscilatory`), and the
replicate folder is reduced to its number (`Replicate_1` -> `1`). Keeping the
full flow name means conditions that share a flow rate (`4dyn laminar` vs
`4dyn 1hz oscilatory`) get distinct names rather than colliding.

During analysis the name is split on the first underscore into `oxygen` (`21`)
and `condition` (`4dyn_laminar_1`) for the summary and QC columns.

## Frame ordering (time series)

Each Cellpose mask file is one hourly timepoint. Ordering comes from the trailing
number in the filename, and the numbering is shifted so the first file is the
baseline (hour 0):

- `..._001` -> `FRAME` 0 (baseline / pre-flow)
- `..._002` -> `FRAME` 1 (1 hour)
- `..._003` -> `FRAME` 2 (2 hours), and so on

So `FRAME` is the number of hours since baseline. Files are read in ascending
order and stacked into a single dataset, giving one CSV with consecutive `FRAME`
values starting at 0. This assumes filenames start numbering at `_001`.

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

**Batch — process a whole cell-line tree in one go:**

```
python main.py batch --input_dir /path/to/CellLine --reference_deg 90
```

With no `--input_dir`, a folder picker opens; select the cell-line folder. The
tool walks three levels down (oxygen / flow / replicate), and for every replicate
folder that contains masks it names the dataset from its folder path and writes
outputs inside that replicate folder.

The expected folder structure is:

```
Cell line/                     (e.g. HUtMEC)
    Oxygen 21%/
        4dyn laminar/
            Replicate_1/
                ..._001_seg.npy   (baseline, hour 0)
                ..._002_seg.npy   (1 hour)
                ...
```

Replicate folders with no masks are skipped, so a partially segmented tree is
fine, and the batch is safe to re-run.

The `--stage` flag chooses which of the two stages to run:

- `--stage extract` — the slow pass: read every replicate's masks and write its
  per-ROI CSV. Nothing is analysed. Run this once per reference angle.
- `--stage analyse` — the fast pass: read the CSVs and write summaries, Prism
  tables and QC. Re-runnable (e.g. to change `--ci`). Reads from the folder
  matching `--reference_deg`, so pass the same angle you extracted with.
- `--stage both` (default) — extract then analyse in one go.

Typical two-pass workflow for vertically-running flow:

```
python main.py batch --input_dir /path/to/HUtMEC --stage extract --reference_deg 90
python main.py batch --input_dir /path/to/HUtMEC --stage analyse --reference_deg 90
```

Other batch options mirror the two stages: `--theta_units`, `--ci`,
`--field_size`, `--min_area`.

To analyse a single replicate on its own (for spot-checks or reprocessing one
dataset), point `analyse` straight at that replicate's ROI folder:

```
python main.py analyse --input_dir '/path/to/HUtMEC/Oxygen 21%/4dyn laminar/Replicate_1/ROI_CSVs'
```

## Outputs

Under the results folder (`Endothelial_Alignment_<date>` by default, written
inside the replicate's ROI folder):

- `Per_Dataset_Summaries/` — the full long-format summary per dataset
- `Prism/AP_Median`, `Prism/AP_CI_low`, `Prism/AP_CI_high`,
  `Prism/ThetaDeg_Median` — wide tables (frames down, quadrants across) for
  GraphPad Prism
- `Prism/AP_Median/<name>_AP_median_IQR.csv` and `..._AP_median_CI.csv` —
  combined tables where each cell reads `median [low, high]` (3 decimals), one
  file using the interquartile range as the bounds, the other the confidence
  interval
- `Combined/ALL_summaries_long.csv` — every dataset stacked into one table
- `QC/QC_summary.csv` — per-dataset quality check: frame count, min/median/max
  ROIs per quadrant-frame, and how many quadrant-frames fell below the minimum
  needed for a confidence interval (`N_QUADRANT_FRAMES_BELOW_CI_MIN`)

After a batch run, a combined QC file for every dataset processed is written at
the cell-line root as `ALL_QC_summary.csv` — one row per dataset, so the whole
experiment's data quality can be scanned at a glance (sparse or failed wells show
up as high `N_QUADRANT_FRAMES_BELOW_CI_MIN` values).
