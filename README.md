
==================================================
Endothelial Cell Alignment Workflow
==================================================
endothelial cell alignment script
    - per frame median + Confidence Interval (CI) for each of four equal-area quadrants and for the whole frame

For every frame you therefore get five sets of outputs:
    - Q1, Q2, Q3, Q4  (the four equal squares of the imaging field)
    - Qtotal          (all ROIs in the frame)

Each of the five reports the Alignmnet Parameter (AP) and initial orientation angles:
    - median
    - Confidence intervals
    - n ROIs
    - Interquartile Range

==================================================
Quadrants Definition
==================================================
Each frame/field is a fixed square image (default 1992 x 1992 px)
    - Each quadrant is exactly one quarter of the field: 996 x 996 px
    - Q Total = whole 1992 x 1992 field
            Q1 = top-left      (x < 996,  y < 996)
            Q2 = top-right     (x >= 996, y < 996)
            Q3 = bottom-left   (x < 996,  y >= 996)
            Q4 = bottom-right  (x >= 996, y >= 996)

Note: X/Y coordinates are in TrackMate/image convention (left -> right, top -> bottom).       

Note: Set the image size with --field_size (default 1992)

==================================================
Measurements
==================================================
The Alignment Parameter (AP) is calculated per object by:
    - AP = cos(2*theta) = 2*cos^2(theta) -1
    - AP values are between +1 (parallel with reference axis) or -1 (perpendicular with reference axis)
    - A value of 0 is considered random (non-specific alignment)

Degrees are calculated from the radians of each cell against the reference axis:
    - theta_deg = orientation angle folded to 0-90 degrees

NOTE: theta_deg provides an absolute positive value to assess outcomes

==================================================
Usage and File Requirements
==================================================
Expected filenames: oxygen_condition_sequence.csv  (e.g. 21_4dyn_1.csv)
    - If analysing frames individually, sequentially numbers files of the same oxygen + condition are stitched together and offset to avoid overlap
    - e.g. 21_4dyn_1.csv (has 3 frames) + 21_4dyn_2.csv (has 12 frames) ==> final output will have 15 consecutive frames

Required CSV columns: Frame, Elliipse_Theta, and X/Y position
    - Non essential rows are removed

Place script in your preferred directory
    - Script requires ## 'pandas and numpy'
    - Run in an appropriate venv (e.g. conda activate pandas)
 
Options:
    - Default: python alignment.py   # Select folder containing CSV files (Will use all files within the folder)
    - Manual Folder Selection: alignment.py --input_dir ./data  # or point at a folder
  

