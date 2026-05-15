# -*- coding: utf-8 -*-
"""
02_extract_and_calculate_RI.py

Post-processes modeled soil-depth rasters to estimate recurrence intervals.

Workflow:
1. Extract mean soil depth through time for each hollow and candidate buffer size.
2. Save the extracted buffer-soil-depth table for reuse.
3. Calculate factor of safety for each cohesion/saturation scenario.
4. Interpolate the first FS = 1 crossing.
5. Save optimal-buffer recurrence-interval results.

This script can be rerun many times without rerunning soil transport.
"""

import glob
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.plot import show
from rasterstats import zonal_stats
from scipy.interpolate import interp1d

from config import WorkflowConfig


def ensure_point_id(gdf):
    """Standardize an 'id' field to 'Point_ID' when needed."""
    if "Point_ID" not in gdf.columns and "id" in gdf.columns:
        gdf = gdf.rename(columns={"id": "Point_ID"})
    return gdf


def safe_reproject_gdf(gdf, target_crs, assumed_source_epsg=None):
    """Reproject without overwriting valid CRS metadata."""
    if gdf.crs is None:
        if assumed_source_epsg is None:
            raise ValueError("Input GeoDataFrame has no CRS. Provide assumed_source_epsg.")
        gdf = gdf.set_crs(epsg=assumed_source_epsg)

    if str(gdf.crs).lower() != str(target_crs).lower():
        gdf = gdf.to_crs(target_crs)

    return gdf


def load_hollow_geometries(cfg):
    """Load hollow points and downslope lines in the target CRS."""
    points_path = cfg.reprojected_points_path if cfg.reprojected_points_path.exists() else cfg.points_shp_path
    gdf = gpd.read_file(points_path)
    line_gdf = gpd.read_file(cfg.downslope_lines)

    gdf = ensure_point_id(gdf)
    line_gdf = ensure_point_id(line_gdf)

    gdf = safe_reproject_gdf(gdf, cfg.target_crs, cfg.source_shp_epsg)
    line_gdf = safe_reproject_gdf(line_gdf, cfg.target_crs, cfg.source_shp_epsg)

    return gdf, line_gdf


def preview_buffers(gdf, soil_depth_files, buffer_sizes):
    """Optional visual check of buffer sizes over a soil-depth raster."""
    preview_files = [f for f in soil_depth_files if "100yrs" in f or "200yrs" in f]
    if not preview_files:
        print("No early-year soil-depth raster found for preview.")
        return

    buffer_geoms = []
    for _, row in gdf.iterrows():
        for buffer_distance in buffer_sizes:
            buffer_geoms.append({
                "geometry": row.geometry.buffer(buffer_distance),
                "Point_ID": row["Point_ID"],
                "Buffer_Size": buffer_distance,
            })

    buffer_gdf = gpd.GeoDataFrame(buffer_geoms, crs=gdf.crs)

    with rasterio.open(preview_files[0]) as src:
        fig, ax = plt.subplots(figsize=(10, 10))
        show(src, ax=ax, cmap="terrain", title="Soil depth and candidate buffer zones")
        buffer_gdf.plot(ax=ax, column="Buffer_Size", cmap="viridis", alpha=0.3, legend=True)
        gdf.plot(ax=ax, color="black", markersize=5)
        for x, y, pid in zip(gdf.geometry.x, gdf.geometry.y, gdf["Point_ID"]):
            ax.text(x, y, str(pid), fontsize=7, ha="center", va="center")
        plt.tight_layout()
        plt.show()


def calculate_line_slopes(line_gdf, slope_path):
    """Calculate one average downslope-line slope value per hollow."""
    point_slope = {}
    point_count = {}

    for _, row in line_gdf.iterrows():
        pid = row["Point_ID"]
        try:
            stats = zonal_stats([row.geometry], slope_path, stats=["count", "mean"], nodata=-9999)
            point_slope[pid] = stats[0]["mean"]
            point_count[pid] = stats[0]["count"]
        except Exception as exc:
            point_slope[pid] = np.nan
            point_count[pid] = 0
            print(f"Failed slope extraction for Point {pid}: {exc}")

    return point_slope, point_count


def extract_soil_depth_by_buffer(cfg, make_preview=False):
    """Extract average soil depth for each point, buffer size, and model year."""
    gdf, line_gdf = load_hollow_geometries(cfg)

    with rasterio.open(cfg.dem_path) as dem:
        dem_crs = dem.crs
        if gdf.crs != dem_crs:
            gdf = gdf.to_crs(dem_crs)
        if line_gdf.crs != dem_crs:
            line_gdf = line_gdf.to_crs(dem_crs)

    soil_depth_files = sorted(glob.glob(cfg.soil_depth_pattern))
    if not soil_depth_files:
        raise FileNotFoundError(f"No soil-depth rasters found: {cfg.soil_depth_pattern}")

    if make_preview:
        preview_buffers(gdf, soil_depth_files, cfg.buffer_sizes)

    point_slope, point_count = calculate_line_slopes(line_gdf, cfg.slope_path)
    records = []

    for _, row in gdf.iterrows():
        point_geom = row.geometry
        point_id = row["Point_ID"]
        avg_slope = point_slope.get(point_id, np.nan)
        slope_count = point_count.get(point_id, 0)

        for buffer_distance in cfg.buffer_sizes:
            buffer_geom = point_geom.buffer(buffer_distance)
            buffer_json = [buffer_geom.__geo_interface__]

            for soil_depth_tif in soil_depth_files:
                match = re.search(r"_(\d+)yrs.*\.tif$", soil_depth_tif)
                if not match:
                    continue

                year = int(match.group(1))

                try:
                    with rasterio.open(soil_depth_tif) as soil_raster:
                        soil_image, _ = mask(soil_raster, buffer_json, crop=True)
                        soil_image = soil_image[0]
                        nodata_val = soil_raster.nodata

                        if nodata_val is None:
                            valid = soil_image[np.isfinite(soil_image) & (soil_image > 0)]
                        else:
                            valid = soil_image[
                                np.isfinite(soil_image)
                                & (soil_image != nodata_val)
                                & (soil_image > 0)
                            ]

                        avg_soil_depth = np.nan if valid.size == 0 else float(np.mean(valid))

                except Exception as exc:
                    avg_soil_depth = np.nan
                    print(
                        f"Failed soil extraction: Point {point_id}, "
                        f"Buffer {buffer_distance}, Year {year}: {exc}"
                    )

                records.append({
                    "Point_ID": point_id,
                    "Year": year,
                    "Buffer_Size": buffer_distance,
                    "Avg_Slope": avg_slope,
                    "Slope_Pixel_Count": slope_count,
                    "Avg_Soil_Depth": avg_soil_depth,
                })

    soil_df = pd.DataFrame(records)
    if not soil_df.empty:
        soil_df = soil_df.sort_values(["Point_ID", "Year", "Buffer_Size"]).reset_index(drop=True)

    return soil_df


def calculate_fs_for_row(row, C0, m, cfg):
    """Calculate MD-STAB-style factor of safety for one row."""
    theta = np.radians(row["Avg_Slope"])
    z = row["Avg_Soil_Depth"]

    if np.isnan(z) or np.isnan(theta) or z <= 0:
        return np.nan

    yw = cfg.g * cfg.pw
    ys = cfg.g * cfg.rho_s
    phi = np.deg2rad(cfg.phi_deg)

    Crb = C0 * np.exp(-z * cfg.j)
    Crl = (C0 / (cfg.j * z)) * (1 - np.exp(-z * cfg.j))

    K0 = 1 - np.sin(theta)
    Kp = np.tan(np.deg2rad(45) + phi / 2) ** 2
    Ka = np.tan(np.deg2rad(45) - phi / 2) ** 2

    Frb = (
        Crb
        + (np.cos(theta) ** 2) * z * (ys - yw * m) * np.tan(phi)
    ) * cfg.l * cfg.w

    Frc = (
        Crl
        + K0 * 0.5 * z * (ys - yw * m ** 2) * np.tan(phi)
    ) * (np.cos(theta) * z * cfg.l * 2)

    Frddu = (Kp - Ka) * 0.5 * z ** 2 * (ys - yw * m ** 2) * cfg.w
    Fdc = np.sin(theta) * np.cos(theta) * z * ys * cfg.l * cfg.w

    return (Frb + Frc + Frddu) / Fdc if Fdc != 0 else np.nan


def find_optimal_buffers(soil_df, C0, m, cfg):
    """Find earliest interpolated FS = 1 crossing for each hollow."""
    df = soil_df.copy()
    df["FS"] = df.apply(lambda row: calculate_fs_for_row(row, C0, m, cfg), axis=1)

    optimal_buffers = {}

    for point_id in df["Point_ID"].unique():
        point_data = df[df["Point_ID"] == point_id].copy()
        buffer_crossings = {}

        for buffer_size in sorted(point_data["Buffer_Size"].unique()):
            buffer_data = point_data[point_data["Buffer_Size"] == buffer_size].sort_values("Year")
            years = buffer_data["Year"].values
            fs_values = buffer_data["FS"].values

            if len(years) < 2 or np.all(~np.isfinite(fs_values)):
                continue

            if cfg.skip_failed_at_start and fs_values[0] <= cfg.fs_threshold:
                continue

            diff = fs_values - cfg.fs_threshold
            crossing_idx = np.where(
                np.isfinite(diff[:-1])
                & np.isfinite(diff[1:])
                & (diff[:-1] * diff[1:] <= 0)
            )[0]

            if crossing_idx.size == 0:
                continue

            i = crossing_idx[0]
            fs_i, fs_ip1 = fs_values[i], fs_values[i + 1]
            t_i, t_ip1 = years[i], years[i + 1]

            frac = 0.0 if np.isclose(fs_ip1, fs_i) else (cfg.fs_threshold - fs_i) / (fs_ip1 - fs_i)
            crossing_year = t_i + frac * (t_ip1 - t_i)

            if years.min() < crossing_year < years.max():
                buffer_crossings[buffer_size] = crossing_year

        if not buffer_crossings:
            print(f"No valid FS = 1 crossing for Point_ID {point_id}.")
            continue

        optimal_buffer_size = min(buffer_crossings, key=buffer_crossings.get)
        optimal_year = buffer_crossings[optimal_buffer_size]
        selected = point_data[point_data["Buffer_Size"] == optimal_buffer_size].sort_values("Year")

        soil_depth_interp = interp1d(
            selected["Year"].values,
            selected["Avg_Soil_Depth"].values,
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        slope_interp = interp1d(
            selected["Year"].values,
            selected["Avg_Slope"].values,
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )

        optimal_buffers[point_id] = {
            "Optimal_Buffer_m": optimal_buffer_size,
            "Year": float(optimal_year),
            "FS": cfg.fs_threshold,
            "Avg_Soil_Depth_m": float(soil_depth_interp(optimal_year)),
            "Avg_Slope_deg": float(slope_interp(optimal_year)),
            "Cohesion": int(C0),
            "m": float(m),
        }

    out_df = pd.DataFrame.from_dict(optimal_buffers, orient="index")
    out_df.reset_index(inplace=True)
    out_df.rename(columns={"index": "Point_ID"}, inplace=True)
    return out_df


def run_optimal_area_scenarios(soil_df, cfg):
    """Run all requested cohesion and saturation scenarios."""
    outputs = []

    for m in cfg.m_values:
        for C0 in cfg.cohesion_values:
            out_df = find_optimal_buffers(soil_df, C0=C0, m=m, cfg=cfg)

            m_label = str(m).replace(".", "p")
            out_path = cfg.results_dir / f"optimal_buffer_results_interpolated_{cfg.basename}_{int(C0)}_m{m_label}.csv"
            out_df.to_csv(out_path, index=False)
            outputs.append(out_df)
            print(f"Saved: {out_path}")

    if outputs:
        return pd.concat(outputs, ignore_index=True)
    return pd.DataFrame()


def main(run_extraction=True, run_fs_analysis=True, load_existing_extraction=False, make_preview=False):
    cfg = WorkflowConfig()
    cfg.make_dirs()

    soil_df = None

    if run_extraction:
        soil_df = extract_soil_depth_by_buffer(cfg, make_preview=make_preview)
        soil_df.to_csv(cfg.buffer_soil_depth_csv, index=False)
        print(f"Saved extracted buffer soil depths: {cfg.buffer_soil_depth_csv}")

    if load_existing_extraction:
        if not cfg.buffer_soil_depth_csv.exists():
            raise FileNotFoundError(
                f"No saved extraction CSV found: {cfg.buffer_soil_depth_csv}. "
                "Run with RUN_EXTRACTION = True first."
            )
        soil_df = pd.read_csv(cfg.buffer_soil_depth_csv)
        print(f"Loaded extracted buffer soil depths: {cfg.buffer_soil_depth_csv}")

    if run_fs_analysis:
        if soil_df is None:
            if cfg.buffer_soil_depth_csv.exists():
                soil_df = pd.read_csv(cfg.buffer_soil_depth_csv)
                print(f"Loaded extracted buffer soil depths: {cfg.buffer_soil_depth_csv}")
            else:
                raise FileNotFoundError(
                    f"No extraction dataframe available. Expected: {cfg.buffer_soil_depth_csv}."
                )

        return run_optimal_area_scenarios(soil_df, cfg)

    return soil_df


if __name__ == "__main__":
    # First run after new soil-transport outputs:
    # RUN_EXTRACTION = True
    # RUN_FS_ANALYSIS = True
    # LOAD_EXISTING_EXTRACTION = False

    # Fast rerun after changing only C0, m, or FS logic:
    # RUN_EXTRACTION = False
    # RUN_FS_ANALYSIS = True
    # LOAD_EXISTING_EXTRACTION = True

    RUN_EXTRACTION = True
    RUN_FS_ANALYSIS = True
    LOAD_EXISTING_EXTRACTION = False
    MAKE_PREVIEW = False

    main(
        run_extraction=RUN_EXTRACTION,
        run_fs_analysis=RUN_FS_ANALYSIS,
        load_existing_extraction=LOAD_EXISTING_EXTRACTION,
        make_preview=MAKE_PREVIEW,
    )
