# -*- coding: utf-8 -*-
"""
01_run_soil_transport.py

Runs the soil transport + soil production model and prepares outputs for
recurrence-interval analysis.

Workflow:
1. Initialize hollow soil depth using point buffers.
2. Run Landlab TaylorNonLinearDiffuser with soil production.
3. Save time-stepped soil-depth, elevation-change, and production GeoTIFFs.
4. Reproject output GeoTIFFs and shapefiles to the target CRS.

For publication/reproducibility, this script should be run once per extent.
"""

import glob
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from landlab import imshowhs_grid
from landlab.components import TaylorNonLinearDiffuser
from landlab.io import read_esri_ascii, write_esri_ascii
from osgeo import gdal
from rasterio.features import geometry_mask
from rasterio.warp import Resampling, calculate_default_transform, reproject

from config import WorkflowConfig


FT_TO_M_US = 1200 / 3937
FT_TO_M_INT = 0.3048


def safe_reproject_gdf(gdf, target_crs, assumed_source_epsg=None):
    """Reproject without overwriting valid CRS metadata."""
    if gdf.crs is None:
        if assumed_source_epsg is None:
            raise ValueError("Input GeoDataFrame has no CRS. Provide assumed_source_epsg.")
        gdf = gdf.set_crs(epsg=assumed_source_epsg)

    if str(gdf.crs).lower() != str(target_crs).lower():
        gdf = gdf.to_crs(target_crs)

    return gdf


def create_buffer_from_points(input_points_path, output_buffer_path, buffer_distance, target_crs=None):
    """Create hollow initialization buffers from point shapefile."""
    gdf = gpd.read_file(input_points_path)

    if gdf.crs is None:
        raise ValueError(f"{input_points_path} has no CRS. Define it before creating buffers.")

    if target_crs is not None and gdf.crs != target_crs:
        gdf = gdf.to_crs(target_crs)

    buffered = gdf.copy()
    buffered["geometry"] = buffered.buffer(buffer_distance)
    buffered.to_file(output_buffer_path)

    print(f"Buffer created: {output_buffer_path}")
    return output_buffer_path


def tiff_to_asc(in_path, out_path):
    """Convert GeoTIFF to ESRI ASCII format for Landlab."""
    with rasterio.open(in_path) as src:
        xyz_unit = src.crs.linear_units if src.crs else "meters"
        mean_res = np.mean(src.res)

    gdal.Translate(str(out_path), str(in_path), format="AAIGrid", xRes=mean_res, yRes=mean_res)
    return mean_res, xyz_unit


def asc_to_tiff(asc_path, tiff_path, meta):
    """Convert Landlab ASCII output back to GeoTIFF."""
    data = np.loadtxt(asc_path, skiprows=10)
    out_meta = meta.copy()
    out_meta.update(dtype=rasterio.float32, count=1, compress="deflate")

    with rasterio.open(tiff_path, "w", **out_meta) as dst:
        dst.write(data.astype(rasterio.float32), 1)


def apply_buffer_to_soil_depth(grid, shapefile, buffer_distance, dem_path, show_diagnostics=False):
    """Set initial soil depth to zero inside hollow buffer zones."""
    gdf = gpd.read_file(shapefile)

    with rasterio.open(dem_path) as src:
        dem_crs = src.crs
        transform = src.transform
        out_shape = src.read(1).shape

    if gdf.crs is not None and gdf.crs != dem_crs:
        gdf = gdf.to_crs(dem_crs)

    buffered_geoms = gdf.buffer(buffer_distance)
    unified_geom = buffered_geoms.unary_union
    buffer_mask = geometry_mask([unified_geom], transform=transform, invert=True, out_shape=out_shape)

    soil_depth = grid.at_node["soil__depth"]
    soil_depth[np.flipud(buffer_mask).flatten()] = 0.0
    grid.at_node["soil__depth"] = soil_depth

    print(f"Applied 0 m initial soil depth to {np.sum(soil_depth == 0.0)} cells.")

    if show_diagnostics:
        plt.figure(figsize=(6, 5))
        plt.imshow(buffer_mask, cmap="gray", origin="upper")
        plt.title("Buffer mask")
        plt.colorbar(label="In buffer")
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(6, 5))
        plt.imshow(soil_depth.reshape(grid.shape), cmap="viridis", origin="upper")
        plt.title("Initial soil depth after masking")
        plt.colorbar(label="Soil depth (m)")
        plt.tight_layout()
        plt.show()

    return grid


def init_simulation(asc_file, cfg, xyz_unit, buffer_shapefile):
    """Initialize Landlab grid and nonlinear diffuser."""
    grid, _ = read_esri_ascii(asc_file, name="topographic__elevation")
    grid.set_closed_boundaries_at_grid_edges(False, False, False, False)

    soil_depth = np.full(grid.number_of_nodes, 0.5)
    grid.add_field("soil__depth", soil_depth, at="node", clobber=True)
    grid = apply_buffer_to_soil_depth(
        grid,
        buffer_shapefile,
        cfg.hollow_buffer_distance,
        cfg.input_tiff_path,
    )

    if xyz_unit is None or "meter" in xyz_unit.lower() or "metre" in xyz_unit.lower():
        K_converted = cfg.K
    elif "foot" in xyz_unit.lower():
        K_converted = cfg.K / (FT_TO_M_US ** 2) if "US" in xyz_unit else cfg.K / (FT_TO_M_INT ** 2)
    else:
        raise RuntimeError(f"Unsupported raster unit for K conversion: {xyz_unit}")

    diffuser = TaylorNonLinearDiffuser(
        grid,
        linear_diffusivity=K_converted,
        slope_crit=cfg.Sc,
        dynamic_dt=True,
        nterms=2,
        if_unstable="pass",
    )

    return grid, diffuser


def save_array_as_tiff(data, filename, meta, grid_shape):
    """Save a Landlab node array as a GeoTIFF."""
    if data.ndim == 1:
        data = data.reshape(grid_shape)

    out_meta = meta.copy()
    out_meta.update(dtype=rasterio.float32, count=1, compress="deflate")

    with rasterio.open(filename, "w", **out_meta) as dst:
        dst.write(np.flipud(data.astype(rasterio.float32)), 1)


def save_elevation_outputs(grid, basefilename, time, cfg):
    """Save elevation as PNG and ASCII at a given timestep."""
    plt.figure(figsize=(6, 5.25))
    imshowhs_grid(grid, "topographic__elevation", plot_type="Hillshade")
    plt.title(f"{basefilename}, {time} yr, K = {cfg.K}")
    plt.tight_layout()
    plt.savefig(cfg.png_dir / f"{basefilename}_{time}yrs_K{cfg.K}.png", dpi=150)
    plt.close()

    asc_path = cfg.asc_dir / f"{basefilename}_{time}yrs_K{cfg.K}.asc"
    write_esri_ascii(str(asc_path), grid, names=["topographic__elevation"], clobber=True)
    return asc_path


def run_soil_transport_simulation(cfg):
    """Run soil transport + production and save time-stepped GeoTIFF outputs."""
    cfg.make_dirs()

    with rasterio.open(cfg.input_tiff_path) as src:
        dem_crs = src.crs
        meta = src.meta.copy()

    buffer_shapefile = create_buffer_from_points(
        cfg.points_shp_path,
        cfg.buffer_shp_path,
        cfg.hollow_buffer_distance,
        target_crs=dem_crs,
    )

    basefilename = Path(cfg.input_tiff).stem
    asc_input = cfg.base_dir / f"{basefilename}.asc"
    _, xyz_unit = tiff_to_asc(cfg.input_tiff_path, asc_input)

    grid, diffuser = init_simulation(asc_input, cfg, xyz_unit, buffer_shapefile)

    z_old = grid.at_node["topographic__elevation"].copy()
    total_soil_depth = grid.at_node["soil__depth"].copy()
    initial_soil_depth = total_soil_depth.copy()

    save_elevation_outputs(grid, basefilename, 0, cfg)

    n_steps = int(cfg.target_time / cfg.dt)
    for i in range(n_steps):
        diffuser.run_one_step(cfg.dt)
        time = (i + 1) * cfg.dt

        production_rate = (cfg.pr / cfg.ps) * (cfg.P0 * np.exp(-total_soil_depth / cfg.h0)) * cfg.dt
        z_new = grid.at_node["topographic__elevation"]
        elevation_change = z_new - z_old

        if xyz_unit and "foot" in xyz_unit.lower():
            elevation_change = elevation_change * FT_TO_M_US

        nonzero_soil_mask = initial_soil_depth > 0.0
        erosion_exceeds = (np.abs(elevation_change) > initial_soil_depth) & nonzero_soil_mask

        if np.any(erosion_exceeds):
            z_new[erosion_exceeds] = z_old[erosion_exceeds] - total_soil_depth[erosion_exceeds]
            total_soil_depth[erosion_exceeds] = production_rate[erosion_exceeds]

        change_in_soil_depth = production_rate.copy()
        change_in_soil_depth[~erosion_exceeds] = elevation_change[~erosion_exceeds] + production_rate[~erosion_exceeds]
        total_soil_depth = np.where(erosion_exceeds, production_rate, total_soil_depth + change_in_soil_depth)

        grid.at_node["soil__depth"] = total_soil_depth
        z_old = z_new.copy()

        if time % cfg.output_interval == 0:
            save_array_as_tiff(elevation_change, cfg.tif_dir / f"{basefilename}_change_in_elevation_{time}yrs.tif", meta, grid.shape)
            save_array_as_tiff(change_in_soil_depth, cfg.tif_dir / f"{basefilename}_change_in_soil_depth_{time}yrs.tif", meta, grid.shape)
            save_array_as_tiff(total_soil_depth, cfg.tif_dir / f"{basefilename}_total_soil_depth_{time}yrs.tif", meta, grid.shape)
            save_array_as_tiff(production_rate, cfg.tif_dir / f"{basefilename}_production_rate_{time}yrs.tif", meta, grid.shape)

            asc_path = save_elevation_outputs(grid, basefilename, time, cfg)
            asc_to_tiff(asc_path, cfg.tif_dir / f"{basefilename}_{time}yrs_K{cfg.K}.tif", meta)

    for suffix in [".asc", ".prj"]:
        cleanup = asc_input.with_suffix(suffix)
        if cleanup.exists():
            cleanup.unlink()

    print("Soil transport simulation complete.")


def reproject_all_tifs(input_folder, output_folder, target_crs="EPSG:32610", suffix="_32610"):
    """Reproject all GeoTIFFs in a folder to the target CRS."""
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    tif_files = glob.glob(str(input_folder / "*.tif"))

    for tif_path in tif_files:
        tif_path = Path(tif_path)
        with rasterio.open(tif_path) as src:
            transform, width, height = calculate_default_transform(
                src.crs, target_crs, src.width, src.height, *src.bounds
            )

            kwargs = src.meta.copy()
            kwargs.update({"crs": target_crs, "transform": transform, "width": width, "height": height})

            out_name = tif_path.stem + f"{suffix}.tif"
            out_path = output_folder / out_name

            with rasterio.open(out_path, "w", **kwargs) as dst:
                for band in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, band),
                        destination=rasterio.band(dst, band),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=target_crs,
                        resampling=Resampling.nearest,
                    )

        print(f"Reprojected TIFF: {out_name}")


def reproject_shapefiles_safe(input_folder, output_folder, target_crs="EPSG:32610", assumed_source_epsg=6557):
    """Reproject all shapefiles without overwriting existing CRS metadata."""
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    for shp_path in glob.glob(str(input_folder / "*.shp")):
        shp_path = Path(shp_path)
        try:
            gdf = gpd.read_file(shp_path)
            gdf = safe_reproject_gdf(gdf, target_crs, assumed_source_epsg)
            out_path = output_folder / f"{shp_path.stem}_32610.shp"
            gdf.to_file(out_path)
            print(f"Reprojected shapefile: {shp_path.name} -> {out_path.name}")
        except Exception as exc:
            print(f"Failed to reproject {shp_path.name}: {exc}")


def main(run_transport=True, run_reprojection=True):
    cfg = WorkflowConfig()
    cfg.make_dirs()

    if run_transport:
        run_soil_transport_simulation(cfg)

    if run_reprojection:
        reproject_all_tifs(cfg.tif_dir, cfg.reproj_tif_dir, target_crs=cfg.target_crs)
        reproject_shapefiles_safe(
            cfg.base_dir,
            cfg.reproj_shp_dir,
            target_crs=cfg.target_crs,
            assumed_source_epsg=cfg.source_shp_epsg,
        )


if __name__ == "__main__":
    RUN_TRANSPORT = True
    RUN_REPROJECTION = True

    main(
        run_transport=RUN_TRANSPORT,
        run_reprojection=RUN_REPROJECTION,
    )
