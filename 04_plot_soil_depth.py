# -*- coding: utf-8 -*-
"""
04_plot_soil_depth.py

Plot soil depth at failure as a function of hollow slope.
Run this after 02_extract_and_calculate_RI.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from config import WorkflowConfig
from plot_helpers import add_critical_area_and_volume, clean_ri_dataframe, get_figure_dir, read_ri_results


plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 10,
    "legend.title_fontsize": 11,
})


def main():
    cfg = WorkflowConfig()
    fig_dir = get_figure_dir(cfg)

    ri_df = read_ri_results(cfg.results_dir)
    ri_df = clean_ri_dataframe(ri_df, min_slope=25.0, drop_indices=[38, 40, 44])
    ri_df = add_critical_area_and_volume(ri_df, cfg, saturation=1.0)

    okabe_ito = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]
    markers = ["o", "s", "^", "D", "v", "P"]
    grouped = sorted(ri_df.groupby("Cohesion"), key=lambda item: item[0], reverse=True)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_facecolor("#f0f0f0")
    ax.grid(True, alpha=0.4)

    for i, (cohesion, group) in enumerate(grouped):
        ax.scatter(
            group["Avg_Slope_deg"],
            group["Avg_Soil_Depth_m"],
            marker=markers[i % len(markers)],
            s=70,
            alpha=0.85,
            color=okabe_ito[i % len(okabe_ito)],
            edgecolors="grey",
            linewidths=0.4,
            label=f"{int(cohesion)} Pa",
        )

    ax.set_xlabel(r"Hollow slope, $\theta_H$ (°)", fontweight="bold")
    ax.set_ylabel(r"Critical soil depth, $h_{crit}$ (m)", fontweight="bold")

    xticks = np.arange(
        np.floor(ri_df["Avg_Slope_deg"].min() / 3) * 3,
        ri_df["Avg_Slope_deg"].max() + 1,
        3,
    )
    ax.set_xticks(xticks)

    legend_handles = [
        Line2D(
            [0], [0], marker=markers[i % len(markers)], color="w",
            markerfacecolor=okabe_ito[i % len(okabe_ito)],
            markeredgecolor="grey", markersize=8,
            label=f"{int(cohesion)} Pa",
        )
        for i, (cohesion, _) in enumerate(grouped)
    ]

    ax.legend(handles=legend_handles, title="Cohesion", loc="upper right", frameon=True)

    fig.tight_layout()
    out_path = fig_dir / "soildepth_no_color.png"
    fig.savefig(out_path, dpi=450, bbox_inches="tight")
    plt.show()

    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
