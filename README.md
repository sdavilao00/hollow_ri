# Hollow soil transport and recurrence-interval workflow

This repository contains a two-stage workflow for modeling hollow soil-depth evolution and estimating hollow failure recurrence intervals. The workflow was designed for characteristic hollow DEM snippets from the Oregon Coast Range, where smaller DEM subsets (`extX`) were clipped around representative convergent hollows to reduce computational cost while preserving local hillslope and hollow geometry.

## Files

- `config.py`  
  Stores the project folder, extent name, CRS settings, model parameters, cohesion values, saturation values, and output folders. Edit this file first.

- `01_run_soil_transport.py`  
  Runs the soil transport/soil production model, writes time-stepped GeoTIFFs, and reprojects raster and shapefile outputs.

- `02_extract_and_calculate_RI.py`  
  Reads the reprojected soil-depth rasters, extracts soil depth by hollow and candidate buffer size, calculates factor of safety (FS), interpolates the first FS = 1 crossing, and saves optimal-buffer recurrence-interval outputs.

- `environment.yml`  
  Defines the conda environment needed to run the workflow.

## DEM preparation and hollow selection

The workflow is intended for clipped DEM subsets that contain characteristic hollows and adjacent hillslopes. Each DEM subset should be large enough to preserve the local topographic setting around the hollow, but small enough to keep the Landlab soil-transport simulation computationally efficient.

Each model extent requires:

- A clipped DEM, for example `extX.tif`
- - A slope DEM, for example `slope.tif`
- A point shapefile marking hollow center locations, for example `extX.shp`
- A downslope polyline shapefile used to extract representative hollow slope, for example `polylines/extX_lines.shp`

The scripts reproject outputs to a common projected coordinate reference system, typically EPSG:32610.

## Installation

Create the conda environment from the project folder:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate hollow-ri-model
```

## Folder setup

A simple project folder can look like this:

```text
newcodesoil/
├── config.py
├── 01_run_soil_transport.py
├── 02_extract_and_calculate_RI.py
├── environment.yml
├── README_workflow.md
├── ext26.tif
├── ext26.shp
├── ext26.shx
├── ext26.dbf
├── ext26.prj
├── dem_smooth_m_warp.tif
├── slope_smooth_m_warp.tif
├── polylines/
│   ├── ext26_lines.shp
│   ├── ext26_lines.shx
│   ├── ext26_lines.dbf
│   └── ext26_lines.prj
├── simulation_results/
├── reproj_shp/
└── results/
```

## 1. Edit `config.py`

Set the project folder:

```python
PROJECT_DIR = Path(r"C:/Users/sdavilao/Documents/newcodesoil")
```

Set the extent-specific inputs:

```python
basename = "ext26"
input_tiff = "ext26.tif"
points_shp = "ext26.shp"
```

Check the model and stability scenarios:

```python
m_values = (0.85,)
cohesion_values = (760, 1920)
buffer_sizes = (3, 4, 5, 6, 7, 8, 9)
```

## 2. Run soil transport once per extent

From a terminal opened in the project folder, run:

```bash
python 01_run_soil_transport.py
```

In Spyder, open `01_run_soil_transport.py` and press the green Run button or press `F5`.

This is the expensive stage. It generates time-stepped soil-depth rasters and reprojects them for the RI analysis.

## 3. Run extraction and RI analysis

Then run:

```bash
python 02_extract_and_calculate_RI.py
```

In Spyder, open `02_extract_and_calculate_RI.py` and press Run.

The first time, use these settings near the bottom of `02_extract_and_calculate_RI.py`:

```python
RUN_EXTRACTION = True
RUN_FS_ANALYSIS = True
LOAD_EXISTING_EXTRACTION = False
```

This saves an intermediate file like:

```text
results/new/intermediate/ext26_buffer_soil_depths.csv
```

## 4. Fast reruns for new cohesion or saturation values

After the extraction CSV exists, you can change cohesion, saturation, or FS logic without rerunning soil transport or raster extraction.

Use:

```python
RUN_EXTRACTION = False
RUN_FS_ANALYSIS = True
LOAD_EXISTING_EXTRACTION = True
```

Then rerun:

```bash
python 02_extract_and_calculate_RI.py
```

## Outputs

The recurrence-interval outputs are saved in `results_dir`, which is set in `config.py`. Each output CSV contains one row per hollow and scenario, including:

- `Point_ID`
- `Optimal_Buffer_m`
- `Year`, interpreted as modeled recurrence interval
- `FS`
- `Avg_Soil_Depth_m`
- `Avg_Slope_deg`

## Publication logic

The split between scripts makes the workflow clearer and more reproducible. The soil-depth evolution model only needs to be rerun when the DEM, hollow initialization, or transport/production parameters change. The FS and RI analysis can be rerun independently for different cohesion or saturation scenarios.

## Plotting manuscript figures

After running `01_run_soil_transport.py` and `02_extract_and_calculate_RI.py`, the plotting scripts can be run independently. Each script reads the RI output CSVs from the results folder defined in `config.py` and saves figures to a `figures/` folder inside the project directory.

Run them in this order if you want all manuscript figures:

```bash
python 03_plot_ri.py
python 04_plot_soil_depth.py
python 05_plot_volume.py
python 06_plot_normalized_ri.py
```

The plotting scripts are separated so individual figures can be regenerated without rerunning the soil-transport model or the RI calculation. This is useful when changing figure styling, axis limits, legends, or manuscript formatting.

The plotting workflow is:

1. `03_plot_ri.py` plots recurrence interval as a function of hollow slope.
2. `04_plot_soil_depth.py` plots soil depth at failure as a function of hollow slope.
3. `05_plot_volume.py` calculates critical area/failure volume and plots RI with point size scaled by volume.
4. `06_plot_normalized_ri.py` plots recurrence interval normalized by cohesion.

`plot_helpers.py` contains shared helper functions used by the plotting scripts.
