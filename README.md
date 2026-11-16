# Endothelial Cell Alignment Workflow

A two-stage pipeline that turns Cellpose segmentation masks into per-quadrant
alignment statistics across a time series. Each mask is one hourly timepoint, and
frames are numbered from the baseline: `FRAME` 0 = baseline, `FRAME` 1 = 1 hour,
and so on.

**Stage 1 — `extract`:** read Cellpose masks, fit an ellipse to every ROI, and
write one per-dataset CSV of orientation angles and centroids.

**Stage 2 — `analyse`:** read that CSV, compute the Alignment Parameter (AP),
split each frame into four equal-area quadrants plus the whole frame, and report
per-frame median, confidence interval and IQR. Alongside AP it reports
orientation order, the shared axis, aspect ratio and cell movement between
frames (speed and direction), all from the same CSV, and z-scores of each
dataset's mean against a reference condition. It only calculates and tabulates;
no statistics, comparisons or figures are produced.

For every frame you get five sets of outputs:

- `Q1`, `Q2`, `Q3`, `Q4` — the four equal squares of the imaging field
- `Qtotal` — all ROIs in the frame

Each of the five reports the Alignment Parameter (AP) and the folded orientation
angle, with: median, confidence interval, n ROIs, and interquartile range; plus
the orientation-order, shape and movement measures described under
[Additional measures](#additional-measures). Z-scores are described under
[Z-scores](#z-scores).

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
    ├── morphology.py           # orientation order, shared axis, mean AP, aspect ratio
    ├── tracking.py             # frame-to-frame cell matching, speed and direction
    ├── zscores.py              # per-dataset means and z-scores vs a reference condition
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

## Additional measures

All of these are calculated for every frame and every region (Q1–Q4 and
`Qtotal`) from the stage 1 CSV. Angles are each cell's long-axis angle θ from
the reference (flow) axis, after `--reference_deg` has been applied.

### Orientation order and shared axis (`morphology.py`)

A cell's axis has no head or tail (a cell at +80° and one at −100° lie on the
same line), so the averages use the doubled angle 2θ:

- `C = mean cos(2θ)` and `Sn = mean sin(2θ)` over the cells in the region
- `order_S = √(C² + Sn²)` — orientation order: 0 = cells point in random
  directions, 1 = every cell lies on the same axis, whatever that axis is
- `axis_vs_flow_deg = |½·atan2(Sn, C)|` in degrees — the shared axis relative
  to the flow axis: 0° = along the flow, 90° = across it

The shared axis is only meaningful when `order_S` is reasonably high (above
about 0.1); with random orientations it is noise.

### Aspect ratio

- `AR_median`, `AR_mean` — median and mean of `MAJOR / MINOR` (the fitted
  ellipse's long and short axes). 1 = round; higher = more elongated. Needs the
  `MAJOR` and `MINOR` columns written by the extract stage; if they are absent,
  these columns are empty.

### Cell length

- `length_median`, `length_mean` — median and mean cell length per region, where each
  cell's length is calculated from its mask area `A` (`AREA`, pixel count) and
  aspect ratio `AR`:

      length = 2·√(A · AR / π)

  This is the major axis of an ellipse with the cell's measured area and aspect
  ratio. Because A × AR ∝ length², it combines size and elongation in one value:
  a small elongated cell and a large elongated cell with the same aspect ratio get
  different lengths. Needs `AREA`, `MAJOR` and `MINOR`; otherwise empty. Units: px.

### Cell movement between frames (`tracking.py`)

Cellpose labels restart in every frame, so cells are matched between
consecutive frames using their centroids (`POSITION_X`, `POSITION_Y`) only:

1. every cell at hour t−1 is paired with its nearest cell at hour t, and every
   cell at hour t with its nearest cell at hour t−1;
2. a pair is kept only if the two cells are each other's nearest neighbour and
   are less than `--match_dist` pixels apart (default 35 px, about half a cell
   length). Cells that leave the field, divide, are missed by Cellpose or are
   ambiguous stay unmatched;
3. for each kept pair: `dx = x(t) − x(t−1)`, `dy = y(t) − y(t−1)`,
   `d = √(dx² + dy²)`.

Reported per region at hour t (movement from t−1 to t; hour 0 is empty):

| Column | Meaning |
|---|---|
| `n_matched` | cells matched to the previous frame |
| `speed_mean` | mean d over matched cells (px per frame interval = px/h for hourly frames) |
| `net_dx_mean` | mean dx: + = moving left → right (with flow), − = against flow |
| `net_dy_mean` | mean dy, across the flow: a check that should be about 0 without drift |
| `direction_index` | mean of dx / d over matched cells that moved (d > 0): +1 = all moving with flow, −1 = all against it |

Movement is counted in the region where the cell is at hour t. Flow is assumed
to run left → right along the image x axis; for slides where it runs another
way, the sign of `net_dx_mean` and `direction_index` follows the image x axis,
not the flow.

Pairs of identical consecutive frames (the same segmentation saved twice) are
detected: the copy is flagged `DUPLICATE_FRAME = True` in the summary tables,
its movement is left empty rather than reported as zero, and the frames are
listed in the QC file (`DUPLICATED_FRAMES`).

With hourly frames, cells can move a large part of their own length between
frames, so some pairs will be wrong matches; the distance limit caps their
effect. Shorter frame intervals make matching more reliable and allow a
smaller `--match_dist`.

All distances are in pixels; multiply by the pixel size (µm/px) to convert.

## Z-scores

After analysing a folder (or, in `batch`, after every replicate of a cell line
has been analysed), each dataset's measures are averaged over all frames that
have a value, per region (movement measures start at frame 1; frames flagged
`DUPLICATE_FRAME` — identical copies of the previous frame — are left out for
every measure). Each average is then z-scored per measure and
region:

    z = (dataset mean − mean of the reference datasets) / SD of all datasets in the run

The reference datasets are those named `<oxygen>_<condition>_<replicate>` with
oxygen `--z_reference_oxygen` (default `21`) and condition
`--z_reference_condition` (default `Static`). z = 0 is the reference average,
+1 is one dataset-to-dataset SD above it. Measures z-scored: `AP_median`,
`order_S`, `axis_vs_flow_deg`, `AR_median`, `AR_mean`, `length_median`,
`length_mean`, `n_matched`,
`speed_mean`, `net_dx_mean`, `net_dy_mean`, `direction_index`.

If the run has no reference datasets (for example a single replicate folder),
only the dataset means are written.

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
`ELLIPSE_THETA`, `POSITION_X`, `POSITION_Y` (and `MAJOR`, `MINOR` for aspect
ratio; without them the aspect-ratio columns are left empty). Rows with missing
values are dropped automatically. Other CSVs in the folder that lack these
columns (e.g. exported Prism tables) are skipped with a message.

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
- `--match_dist` — maximum distance in px for matching a cell between frames
  (default 35)
- `--z_reference_oxygen`, `--z_reference_condition` — datasets used as the
  z-score reference (default `21` and `Static`)
- `--per_cell` — also write every cell's values (position, quadrant, angle,
  area, axes, aspect ratio, cell length, movement) to `Per_Cell/`

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
`--field_size`, `--min_area`, `--match_dist`, `--z_reference_oxygen`,
`--z_reference_condition`, `--per_cell`.

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
- `Prism/Order_S`, `Prism/Axis_vs_Flow_Deg`, `Prism/AspectRatio_Median`,
  `Prism/AspectRatio_Mean`, `Prism/CellLength_Median`, `Prism/CellLength_Mean`,
  `Prism/Cells_Matched`, `Prism/Speed_Mean`,
  `Prism/NetDX_Mean`, `Prism/NetDY_Mean`, `Prism/Direction_Index` — the same
  layout for each additional measure
- `Prism/Dataset_Means/<measure>/` and `Prism/Z_Scores/<measure>/` — grouped
  tables, one row per condition and one column per oxygen × replicate
  (`..._whole_frame.csv`) or oxygen × replicate × quadrant (`..._quadrants.csv`)
- `Z_Scores/ALL_dataset_means_long.csv`, `Z_Scores/ALL_zscores_long.csv` — the
  same values in long format (one row per dataset × region)
- `Prism/AP_Median/<name>_AP_median_IQR.csv` and `..._AP_median_CI.csv` —
  combined tables where each cell reads `median [low, high]` (3 decimals), one
  file using the interquartile range as the bounds, the other the confidence
  interval
- `Combined/ALL_summaries_long.csv` — every dataset stacked into one table,
  one row per dataset × frame × region with every measure (the table to open
  in Excel)
- `QC/QC_summary.csv` — per-dataset quality check: frame count, min/median/max
  ROIs per quadrant-frame, and how many quadrant-frames fell below the minimum
  needed for a confidence interval (`N_QUADRANT_FRAMES_BELOW_CI_MIN`); min and
  median ROIs per whole frame, the median number of cells matched between
  frames, and any duplicated frames

After a batch run, every dataset of the cell line is also gathered into
`Cell_line_summary_<ROI folder>/` at the cell-line root: `Combined/ALL_summaries_long.csv`
(all datasets), `Z_Scores/` and `Prism/Dataset_Means`, `Prism/Z_Scores` as
above, calculated across the whole cell line. A combined QC file for every
dataset processed is written at the cell-line root as `ALL_QC_summary.csv` — one row per dataset, so the whole
experiment's data quality can be scanned at a glance (sparse or failed wells show
up as high `N_QUADRANT_FRAMES_BELOW_CI_MIN` values).
