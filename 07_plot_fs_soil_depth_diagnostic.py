# -*- coding: utf-8 -*-
"""
07_plot_fs_soil_depth_diagnostic.py

Create diagnostic figures for a single hollow showing factor of safety (FS),
soil depth, and candidate failure-buffer behavior through time.

Run this after 02_extract_and_calculate_RI.py has generated the intermediate
buffer soil-depth CSV for the selected extent.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from config import WorkflowConfig
from plot_helpers import get_figure_dir


# =============================================================================
# USER SETTINGS
# =============================================================================
TARGET_POINT_ID = 1
DIAGNOSTIC_COHESION = 1920
DIAGNOSTIC_SATURATION = 1.0

# Optional times for the basal-area diagnostic plot. If None, only the failure
# year is used. Change these to values that are meaningful for the selected hollow.
TIMES_TO_PLOT = None  # Example: [1090, 1100]

SAVE_FIGURES = True


plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "legend.title_fontsize": 11,
})


def calculate_fs(df, cfg, cohesion, saturation):
    """Calculate factor of safety for each row of a buffer/year dataframe."""
    out = df.copy()

    theta = np.deg2rad(out["Avg_Slope"].astype(float))
    z = out["Avg_Soil_Depth"].astype(float)

    yw = cfg.g * cfg.pw
    ys = cfg.g * cfg.rho_s
    phi = np.deg2rad(cfg.phi_deg)

    Kp = np.tan(np.deg2rad(45.0) + phi / 2.0) ** 2
    Ka = np.tan(np.deg2rad(45.0) - phi / 2.0) ** 2

    valid = np.isfinite(theta) & np.isfinite(z) & (z > 0)
    fs = np.full(len(out), np.nan, dtype=float)

    theta_v = theta[valid]
    z_v = z[valid]

    Crb = cohesion * np.exp(-z_v * cfg.j)
    Crl = (cohesion / (cfg.j * z_v)) * (1.0 - np.exp(-z_v * cfg.j))
    K0 = 1.0 - np.sin(theta_v)

    Frb = (
        Crb
        + (np.cos(theta_v) ** 2) * z_v * (ys - yw * saturation) * np.tan(phi)
    ) * cfg.l * cfg.w

    Frc = (
        Crl
        + K0 * 0.5 * z_v * (ys - yw * saturation ** 2) * np.tan(phi)
    ) * (np.cos(theta_v) * z_v * cfg.l * 2.0)

    Frddu = (
        (Kp - Ka)
        * 0.5
        * (z_v ** 2)
        * (ys - yw * saturation ** 2)
        * cfg.w
    )

    Fdc = np.sin(theta_v) * np.cos(theta_v) * z_v * ys * cfg.l * cfg.w
    fs[valid] = np.where(Fdc != 0, (Frb + Frc + Frddu) / Fdc, np.nan)

    out["FS"] = fs
    return out


def find_optimal_buffer(point_data, fs_threshold=1.0, skip_failed_at_start=True):
    """Return the first buffer to cross FS threshold and interpolated values."""
    first_crossings = {}

    for buffer_size in sorted(point_data["Buffer_Size"].unique()):
        sub = point_data[point_data["Buffer_Size"] == buffer_size].sort_values("Year")
        years = sub["Year"].to_numpy(dtype=float)
        fs_values = sub["FS"].to_numpy(dtype=float)

        if len(years) < 2 or np.all(~np.isfinite(fs_values)):
            continue

        if skip_failed_at_start and fs_values[0] <= fs_threshold:
            continue

        diff = fs_values - fs_threshold
        idx = np.where(
            np.isfinite(diff[:-1])
            & np.isfinite(diff[1:])
            & (diff[:-1] * diff[1:] <= 0)
        )[0]

        if idx.size == 0:
            continue

        i = idx[0]
        fs_i, fs_ip1 = fs_values[i], fs_values[i + 1]
        t_i, t_ip1 = years[i], years[i + 1]

        frac = 0.0 if np.isclose(fs_ip1, fs_i) else (fs_threshold - fs_i) / (fs_ip1 - fs_i)
        crossing_year = t_i + frac * (t_ip1 - t_i)

        if years.min() < crossing_year < years.max():
            first_crossings[buffer_size] = crossing_year

    if not first_crossings:
        return None

    optimal_buffer = min(first_crossings, key=first_crossings.get)
    optimal_year = first_crossings[optimal_buffer]

    selected = point_data[point_data["Buffer_Size"] == optimal_buffer].sort_values("Year")
    estimated_soil_depth = float(np.interp(optimal_year, selected["Year"], selected["Avg_Soil_Depth"]))
    estimated_slope = float(np.interp(optimal_year, selected["Year"], selected["Avg_Slope"]))

    return {
        "Optimal_Buffer_m": optimal_buffer,
        "Year": optimal_year,
        "FS": fs_threshold,
        "Avg_Soil_Depth_m": estimated_soil_depth,
        "Avg_Slope_deg": estimated_slope,
    }


def savefig(fig, path, save=True):
    if save:
        fig.savefig(path, dpi=450, bbox_inches="tight")


def main():
    cfg = WorkflowConfig()
    fig_dir = get_figure_dir(cfg)

    if not cfg.buffer_soil_depth_csv.exists():
        raise FileNotFoundError(
            f"Missing intermediate extraction CSV: {cfg.buffer_soil_depth_csv}\n"
            "Run 02_extract_and_calculate_RI.py with RUN_EXTRACTION = True first."
        )

    df = pd.read_csv(cfg.buffer_soil_depth_csv)
    point_data = df[df["Point_ID"] == TARGET_POINT_ID].copy()

    if point_data.empty:
        raise ValueError(f"Point_ID {TARGET_POINT_ID} was not found in {cfg.buffer_soil_depth_csv}")

    point_data = calculate_fs(
        point_data,
        cfg=cfg,
        cohesion=DIAGNOSTIC_COHESION,
        saturation=DIAGNOSTIC_SATURATION,
    )

    optimal = find_optimal_buffer(
        point_data,
        fs_threshold=cfg.fs_threshold,
        skip_failed_at_start=cfg.skip_failed_at_start,
    )

    if optimal is None:
        optimal_buffer = sorted(point_data["Buffer_Size"].unique())[0]
        optimal_year = None
        estimated_soil_depth = None
        print("No valid FS=1 crossing found. Plotting the smallest buffer as fallback.")
    else:
        optimal_buffer = optimal["Optimal_Buffer_m"]
        optimal_year = optimal["Year"]
        estimated_soil_depth = optimal["Avg_Soil_Depth_m"]
        print("Optimal buffer result:")
        print(pd.DataFrame([optimal]))

    # -------------------------------------------------------------------------
    # Figure 1: FS and soil depth through time for the optimal/fallback buffer
    # -------------------------------------------------------------------------
    plot_data = point_data[point_data["Buffer_Size"] == optimal_buffer].sort_values("Year")

    fig, ax1 = plt.subplots(figsize=(7, 5))
    ax1.set_facecolor("#f0f0f0")

    ax1.plot(plot_data["Year"], plot_data["FS"], marker="o", lw=1.5, color="C0")
    ax1.set_xlabel("Time (years)")
    ax1.set_ylabel("Factor of safety (FS)", color="C0")
    ax1.tick_params(axis="y", labelcolor="C0")

    ax2 = ax1.twinx()
    ax2.set_facecolor("#f0f0f0")
    ax2.plot(plot_data["Year"], plot_data["Avg_Soil_Depth"], marker="s", lw=1.5, color="C1")
    ax2.set_ylabel("Average soil depth (m)", color="C1")
    ax2.tick_params(axis="y", labelcolor="C1")
    ax2.set_ylim(0, max(0.75, 1.05 * np.nanmax(plot_data["Avg_Soil_Depth"])))

    if optimal_year is not None:
        ax1.set_xlim(0, min(optimal_year * 1.2, plot_data["Year"].max()))
    else:
        ax1.set_xlim(0, plot_data["Year"].max())

    fig.tight_layout()
    savefig(fig, fig_dir / f"{cfg.basename}_point{TARGET_POINT_ID}_fs_soildepth_time.png", SAVE_FIGURES)
    plt.show()

    # -------------------------------------------------------------------------
    # Figure 2: FS by buffer size at the optimal failure year
    # -------------------------------------------------------------------------
    if optimal_year is not None:
        buffers = sorted(point_data["Buffer_Size"].unique())
        fs_at_failure = []
        sd_at_failure = []

        for buffer_size in buffers:
            sub = point_data[point_data["Buffer_Size"] == buffer_size].sort_values("Year")
            fs_at_failure.append(np.interp(optimal_year, sub["Year"], sub["FS"]))
            sd_at_failure.append(np.interp(optimal_year, sub["Year"], sub["Avg_Soil_Depth"]))

        fs_at_failure = np.array(fs_at_failure)
        sd_at_failure = np.array(sd_at_failure)
        basal_areas = np.pi * np.array(buffers) ** 2

        fig, ax1 = plt.subplots(figsize=(9, 6))
        ax1.set_facecolor("#f0f0f0")
        ax1.grid(False)

        norm = colors.Normalize(vmin=min(buffers), vmax=max(buffers))
        cmap = cm.get_cmap("viridis")

        for buffer_size, area, fs_val in zip(buffers, basal_areas, fs_at_failure):
            ax1.plot(area, fs_val, marker="o", markersize=10, color=cmap(norm(buffer_size)))

        ax1.plot(basal_areas, fs_at_failure, lw=1.8, color="black", alpha=0.35)
        ax1.set_xlabel(r"Basal landslide area (m$^2$)")
        ax1.set_ylabel("Factor of safety (FS)")
        ax1.set_xticks(basal_areas)
        ax1.set_xticklabels([f"{area:.0f}" for area in basal_areas])

        ax2 = ax1.twinx()
        ax2.set_facecolor("#f0f0f0")
        for buffer_size, area, sd_val in zip(buffers, basal_areas, sd_at_failure):
            ax2.plot(area, sd_val, marker="s", markersize=10, color=cmap(norm(buffer_size)))

        ax2.plot(basal_areas, sd_at_failure, lw=1.8, color="black", alpha=0.35, linestyle="--")
        ax2.set_ylabel("Average soil depth (m)")

        if optimal_buffer in buffers:
            i_opt = buffers.index(optimal_buffer)
            ax1.plot(
                basal_areas[i_opt], fs_at_failure[i_opt],
                marker="o", markersize=13, markeredgecolor="black",
                markerfacecolor=cmap(norm(optimal_buffer)), zorder=5,
            )
            ax2.plot(
                basal_areas[i_opt], sd_at_failure[i_opt],
                marker="s", markersize=13, markeredgecolor="black",
                markerfacecolor=cmap(norm(optimal_buffer)), zorder=6,
            )

        fig.tight_layout()
        savefig(fig, fig_dir / f"{cfg.basename}_point{TARGET_POINT_ID}_fs_soildepth_basal_area.png", SAVE_FIGURES)
        plt.show()

        # ---------------------------------------------------------------------
        # Figure 3: FS vs basal area at selected times and failure year
        # ---------------------------------------------------------------------
        selected_times = list(TIMES_TO_PLOT) if TIMES_TO_PLOT is not None else []
        selected_times.append(optimal_year)
        selected_labels = [f"{t:.0f} yr" for t in selected_times[:-1]] + [f"Failure year ≈ {optimal_year:.0f} yr"]
        line_styles = ["--", ":", "-"][:len(selected_times)]

        fig, ax = plt.subplots(figsize=(9, 6))
        ax.set_facecolor("#f0f0f0")
        ax.grid(True, which="both", color="#d0d0d0", alpha=0.7)

        buffer_colors = cm.viridis(np.linspace(0, 1, len(buffers)))

        for time, label, line_style in zip(selected_times, selected_labels, line_styles):
            fs_at_time = []
            for buffer_size in buffers:
                sub = point_data[point_data["Buffer_Size"] == buffer_size].sort_values("Year")
                if time < sub["Year"].min() or time > sub["Year"].max():
                    fs_at_time.append(np.nan)
                else:
                    fs_at_time.append(np.interp(time, sub["Year"], sub["FS"]))

            fs_at_time = np.array(fs_at_time)
            valid = np.isfinite(fs_at_time)
            ax.plot(basal_areas[valid], fs_at_time[valid], linestyle=line_style, linewidth=2, color="black", label=label)

            for i, (area, fs_val) in enumerate(zip(basal_areas, fs_at_time)):
                if np.isfinite(fs_val):
                    ax.plot(area, fs_val, marker="o", markersize=8, color=buffer_colors[i])

        if optimal_buffer in buffers:
            i_opt = buffers.index(optimal_buffer)
            fs_opt = np.interp(
                optimal_year,
                point_data[point_data["Buffer_Size"] == optimal_buffer].sort_values("Year")["Year"],
                point_data[point_data["Buffer_Size"] == optimal_buffer].sort_values("Year")["FS"],
            )
            ax.plot(
                basal_areas[i_opt], fs_opt,
                marker="o", markersize=12, markerfacecolor="none",
                markeredgecolor="red", markeredgewidth=2.5, zorder=5,
                label=f"Optimal area = {basal_areas[i_opt]:.0f} m$^2$",
            )

        ax.set_xlabel(r"Basal landslide area (m$^2$)")
        ax.set_ylabel("Factor of safety (FS)")
        ax.set_xticks(basal_areas)
        ax.set_xticklabels([f"{area:.0f}" for area in basal_areas])
        ax.legend(loc="best", frameon=True)

        fig.tight_layout()
        savefig(fig, fig_dir / f"{cfg.basename}_point{TARGET_POINT_ID}_fs_basal_area_times.png", SAVE_FIGURES)
        plt.show()


if __name__ == "__main__":
    main()
