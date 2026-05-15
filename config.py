# -*- coding: utf-8 -*-
"""
Configuration for the hollow soil-transport and recurrence-interval workflow.

Edit this file first when switching extents, folders, or model parameters.
All paths are built from PROJECT_DIR so the workflow can be moved more easily.
"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class WorkflowConfig:
    # -----------------------------
    # Project folder
    # -----------------------------
    # Put the scripts, input DEMs, shapefiles, and output folders inside this folder.
    PROJECT_DIR: Path = Path(r"C:/Users/sdavilao/Documents/newcodesoil")

    # -----------------------------
    # Extent and CRS settings
    # -----------------------------
    basename: str = "ext26"
    target_crs: str = "EPSG:32610"
    source_shp_epsg: int = 6557  # Only assigned when an input shapefile has no CRS

    # -----------------------------
    # Input file names / subfolders
    # -----------------------------
    input_tiff: str = "ext26.tif"
    points_shp: str = "ext26.shp"
    lines_subdir: str = "polylines"
    dem_name: str = "dem_smooth_m_warp.tif"
    slope_name: str = "slope_smooth_m_warp.tif"

    # -----------------------------
    # Output subfolders
    # -----------------------------
    simulation_subdir: str = "simulation_results/new/new"
    reproj_tif_subdir: str = "simulation_results/new/GeoTIFFs/reproj_tif"
    reproj_shp_subdir: str = "reproj_shp"
    results_subdir: str = "results/new/m"
    intermediate_subdir: str = "results/new/intermediate"

    # -----------------------------
    # Soil-transport parameters
    # -----------------------------
    hollow_buffer_distance: float = 16.0
    K: float = 0.0042
    Sc: float = 1.25
    pr: float = 2000.0
    ps: float = 1600.0
    P0: float = 0.0003
    h0: float = 0.5
    dt: int = 50
    target_time: int = 5000
    output_interval: int = 50

    # -----------------------------
    # Candidate landslide-buffer sizes
    # -----------------------------
    buffer_sizes: tuple = (3, 4, 5, 6, 7, 8, 9)

    # -----------------------------
    # Slope-stability parameters
    # -----------------------------
    pw: float = 1000.0
    rho_s: float = 1600.0
    g: float = 9.81
    phi_deg: float = 41.0
    m_values: tuple = (0.85,)
    cohesion_values: tuple = (760, 1920)
    l: float = 10.0
    w: float = 6.7
    j: float = 0.8

    # -----------------------------
    # FS crossing behavior
    # -----------------------------
    fs_threshold: float = 1.0
    skip_failed_at_start: bool = True

    # Internal path aliases created automatically from PROJECT_DIR
    base_dir: Path = field(init=False)
    simulation_out_dir: Path = field(init=False)
    reproj_tif_dir: Path = field(init=False)
    reproj_shp_dir: Path = field(init=False)
    results_dir: Path = field(init=False)
    intermediate_dir: Path = field(init=False)
    downslope_lines: Path = field(init=False)
    dem_path: Path = field(init=False)
    slope_path: Path = field(init=False)

    def __post_init__(self):
        self.base_dir = self.PROJECT_DIR
        self.simulation_out_dir = self.PROJECT_DIR / self.simulation_subdir
        self.reproj_tif_dir = self.PROJECT_DIR / self.reproj_tif_subdir
        self.reproj_shp_dir = self.PROJECT_DIR / self.reproj_shp_subdir
        self.results_dir = self.PROJECT_DIR / self.results_subdir
        self.intermediate_dir = self.PROJECT_DIR / self.intermediate_subdir
        self.downslope_lines = self.PROJECT_DIR / self.lines_subdir / f"{self.basename}_lines.shp"
        self.dem_path = self.PROJECT_DIR / self.dem_name
        self.slope_path = self.PROJECT_DIR / self.slope_name

    @property
    def png_dir(self):
        return self.simulation_out_dir / "PNGs"

    @property
    def tif_dir(self):
        return self.simulation_out_dir / "GeoTIFFs"

    @property
    def asc_dir(self):
        return self.simulation_out_dir / "ASCs"

    @property
    def input_tiff_path(self):
        return self.base_dir / self.input_tiff

    @property
    def points_shp_path(self):
        return self.base_dir / self.points_shp

    @property
    def buffer_shp_path(self):
        return self.base_dir / f"{self.basename}_buff.shp"

    @property
    def soil_depth_pattern(self):
        return str(self.reproj_tif_dir / f"{self.basename}_total_soil_depth_*yrs_32610.tif")

    @property
    def reprojected_points_path(self):
        return self.reproj_shp_dir / f"{self.basename}_32610.shp"

    @property
    def buffer_soil_depth_csv(self):
        return self.intermediate_dir / f"{self.basename}_buffer_soil_depths.csv"

    def make_dirs(self):
        for folder in [
            self.simulation_out_dir,
            self.png_dir,
            self.tif_dir,
            self.asc_dir,
            self.reproj_tif_dir,
            self.reproj_shp_dir,
            self.results_dir,
            self.intermediate_dir,
        ]:
            folder.mkdir(parents=True, exist_ok=True)
