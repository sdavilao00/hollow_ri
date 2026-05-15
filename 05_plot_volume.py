# -*- coding: utf-8 -*-
"""
05_plot_volume.py

Plot recurrence interval as a function of hollow slope, with marker size scaled by
modeled failure volume.
Run this after 02_extract_and_calculate_RI.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext

from config import WorkflowConfig
from plot_helpers import add_critical_area_and_volume, clean_ri_dataframe, get_figure_dir, read_ri_results


plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 10,
    "legend.title_fontsize": 12,
})


def main():
    cfg = WorkflowConfig()
    fig_dir = get_figure_dir(cfg)

    ri_df = read_ri_results(cfg.results_dir)
    ri_df = clean_ri_dataframe(ri_df, min_slope=25.0, drop_indices=[38, 40, 44])
    plot_df = add_critical_area_and_volume(ri_df, cfg, saturation=1.0)
    plot_df = plot_df.dropna(subset=["Avg_Slope_deg", "Year", "Volume", "Cohesion"])
    plot_df = plot_df[(plot_df["Year"] > 0) & (plot_df["Volume"] > 0)]

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_facecolor("#f0f0f0")
    ax.set_axisbelow(True)
    ax.minorticks_on()

    groups = sorted(plot_df.groupby("Cohesion"), key=lambda item: item[0])
    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0, 1, len(groups)))
    color_map = {coh: color for color, (coh, _) in zip(colors, groups)}
    size_scale = 10

    for cohesion, group in groups:
        ax.scatter(
            group["Avg_Slope_deg"],
            group["Year"],
            s=group["Volume"] * size_scale,
            color=color_map[cohesion],
            alpha=0.8,
            edgecolors="k",
            label=f"{int(cohesion)} Pa",
        )

    ax.set_xlabel(r"Hollow slope, $\theta_H$ (°)", fontweight="bold")
    ax.set_ylabel("Recurrence interval, RI (years)", fontweight="bold")

    ax.set_xlim(plot_df["Avg_Slope_deg"].min() - 1, plot_df["Avg_Slope_deg"].max() + 1.5)
    xticks = np.arange(
        np.floor(plot_df["Avg_Slope_deg"].min() / 3) * 3,
        np.ceil(plot_df["Avg_Slope_deg"].max() / 3) * 3 + 3,
        3,
    )
    ax.set_xticks(xticks)

    ax.set_yscale("log")
    ax.set_ylim(100, 1e4)
    ax.set_yticks([100, 1000, 10000])
    ax.yaxis.set_major_formatter(LogFormatterMathtext())

    ax.grid(True, which="major", axis="y", color="#d0d0d0", alpha=0.7)
    ax.grid(True, which="minor", axis="y", color="#d0d0d0", alpha=0.3)
    ax.grid(True, which="major", axis="x", color="#bdbdbd", linewidth=0.8)
    ax.grid(False, which="minor", axis="x")

    cohesions = sorted(plot_df["Cohesion"].unique(), reverse=True)
    cohesion_handles = [
        Line2D(
            [], [], marker="o", linestyle="None",
            color=color_map[coh], markersize=6,
            label=f"{int(coh)}",
        )
        for coh in cohesions
    ]

    vol_values = np.sort(plot_df["Volume"].to_numpy())
    volume_examples = [vol_values[0], vol_values[len(vol_values) // 2], vol_values[-1]]
    volume_handles = [
        Line2D(
            [], [], marker="o", linestyle="None",
            color="gray", markeredgecolor="k",
            markersize=np.sqrt(vol * size_scale) * 0.8,
            label=f"{int(round(vol)):,}",
        )
        for vol in volume_examples
    ]

    spacer = Line2D([], [], linestyle="None", label="")
    coh_header = Line2D([], [], linestyle="None", label="Cohesion (Pa)")
    vol_header = Line2D([], [], linestyle="None", label=r"Volume (m$^3$)")

    n_rows = max(len(cohesion_handles), len(volume_handles))

    def pad(handles, n):
        return handles + [Line2D([], [], linestyle="None", label="")] * (n - len(handles))

    handles = [coh_header] + pad(cohesion_handles, n_rows) + [vol_header] + pad(volume_handles, n_rows)

    ax.legend(
        handles=handles,
        ncol=2,
        loc="upper right",
        frameon=True,
        columnspacing=2.0,
        handletextpad=0.8,
        labelspacing=0.5,
        borderpad=1.2,
    )

    fig.tight_layout()
    out_path = fig_dir / "volume_final.png"
    fig.savefig(out_path, dpi=450, bbox_inches="tight")
    plt.show()

    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
