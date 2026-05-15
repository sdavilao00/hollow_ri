# -*- coding: utf-8 -*-
"""
08_plot_erosion.py

Calculate hollow-scale erosion rates from modeled recurrence interval results,
critical failure volume, and contributing-area zonal statistics. Produces the
slope-erosion plot used for manuscript figures.

Run after 02_extract_and_calculate_RI.py and after zonal statistic CSVs are
available.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

from config import WorkflowConfig
from plot_helpers import (
    add_critical_area_and_volume,
    clean_ri_dataframe,
    get_figure_dir,
    read_ri_results,
)


# =============================================================================
# USER SETTINGS
# =============================================================================
TARGET_COHESION = 1920
SATURATION_FOR_VOLUME = 1.0
MIN_SLOPE_DEG = 25.0
DROP_INDICES = [38, 40, 44]

ZONAL_9_CSV = "zonal_9.csv"
ZONAL_15_CSV = "zonal_15.csv"

SAVE_FIGURE = True
SAVE_EROSION_TABLE = True


plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "legend.title_fontsize": 11,
})


def exp_model(x, a, b):
    """Exponential erosion-slope model."""
    return a * np.exp(b * x)


def r_squared(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1 - ss_res / ss_tot


def load_zonal_table(path, required_column):
    """Read zonal-statistics CSV and normalize extent names."""
    df = pd.read_csv(path)
    if "Extent" not in df.columns:
        raise ValueError(f"{path} is missing required column: Extent")
    if required_column not in df.columns:
        raise ValueError(f"{path} is missing required column: {required_column}")

    df["Extent"] = df["Extent"].astype(str).str.replace(r"_\d+$", "", regex=True)
    return df[["Point_ID", "Extent", required_column]]


def main():
    cfg = WorkflowConfig()
    fig_dir = get_figure_dir(cfg)

    ri_df = read_ri_results(cfg.results_dir)
    ri_df = clean_ri_dataframe(ri_df, min_slope=MIN_SLOPE_DEG, drop_indices=DROP_INDICES)

    ri_df = ri_df[ri_df["Cohesion"] == TARGET_COHESION].copy()
    if ri_df.empty:
        raise ValueError(f"No RI results found for Cohesion = {TARGET_COHESION} Pa")

    ri_df = add_critical_area_and_volume(ri_df, cfg, saturation=SATURATION_FOR_VOLUME)

    zonal_9_path = cfg.base_dir / ZONAL_9_CSV
    zonal_15_path = cfg.base_dir / ZONAL_15_CSV

    zonal_9 = load_zonal_table(zonal_9_path, "max_9")
    zonal_15 = load_zonal_table(zonal_15_path, "max_15")

    merged_ero = (
        ri_df
        .merge(zonal_9, on=["Point_ID", "Extent"], how="inner")
        .merge(zonal_15, on=["Point_ID", "Extent"], how="inner")
    )

    merged_ero["Erosion_9"] = (merged_ero["Volume"] / (merged_ero["Year"] * merged_ero["max_9"])) * 1000.0
    merged_ero["Erosion_15"] = (merged_ero["Volume"] / (merged_ero["Year"] * merged_ero["max_15"])) * 1000.0

    # Remove anomalously low erosion-rate point near 42.7 degrees if present.
    plot_ero = merged_ero[
        ~(
            (merged_ero["Avg_Slope_deg"] > 42.5)
            & (merged_ero["Avg_Slope_deg"] < 43.0)
            & (merged_ero["Erosion_9"] < 0.05)
        )
    ].copy()

    plot_ero = plot_ero.dropna(subset=["Avg_Slope_deg", "Erosion_9"])
    plot_ero = plot_ero[plot_ero["Erosion_9"] > 0]

    if SAVE_EROSION_TABLE:
        out_csv = cfg.results_dir / f"erosion_results_C{TARGET_COHESION}.csv"
        merged_ero.to_csv(out_csv, index=False)
        print(f"Saved erosion table: {out_csv}")

    x = plot_ero["Avg_Slope_deg"].to_numpy()
    y = plot_ero["Erosion_9"].to_numpy()

    params, _ = curve_fit(exp_model, x, y, p0=[0.01, 0.05], maxfev=10000)
    r2 = r_squared(y, exp_model(x, *params))

    slope_range = np.linspace(x.min(), x.max(), 500)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_facecolor("#f0f0f0")
    ax.set_ylim(1e-2, 5e-1)

    ax.scatter(
        x,
        y,
        c="black",
        s=70,
        edgecolors="#111",
        linewidths=0.5,
        alpha=0.9,
        zorder=3,
    )

    ax.plot(
        slope_range,
        exp_model(slope_range, *params),
        c="red",
        linestyle="--",
        linewidth=2.2,
        zorder=2,
    )

    ax.set_xlabel(r"Hollow slope, $\theta_H$ (°)", fontweight="bold")
    ax.set_ylabel(r"Erosion rate, $E$ (m yr$^{-1}$)", fontweight="bold")

    xticks = np.arange(
        np.floor(plot_ero["Avg_Slope_deg"].min() / 3) * 3,
        plot_ero["Avg_Slope_deg"].max() + 3,
        3,
    )
    ax.set_xticks(xticks)
    ax.grid(True, which="both", linestyle="--", alpha=0.40)

    eq_label = fr"$E = {params[0]:.5f}e^{{{params[1]:.4f}\theta_H}}$"
    r2_label = fr"$R^2 = {r2:.3f}$"

    legend_handles = [
        Line2D([], [], color="red", linestyle="--", lw=2.2, label=eq_label),
        Line2D([], [], linestyle="None", label=r2_label),
    ]

    legend = ax.legend(
        handles=legend_handles,
        loc="upper left",
        frameon=True,
        facecolor="white",
        edgecolor="black",
    )
    legend.get_texts()[1].set_color("red")

    fig.tight_layout()

    if SAVE_FIGURE:
        fig.savefig(fig_dir / "erosion_slope_final.png", dpi=450, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()
