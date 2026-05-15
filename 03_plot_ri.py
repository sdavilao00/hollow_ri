# -*- coding: utf-8 -*-
"""
03_plot_ri.py

Plot modeled recurrence interval (RI) as a function of hollow slope.
Run this after 02_extract_and_calculate_RI.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext
from scipy.optimize import curve_fit

from config import WorkflowConfig
from plot_helpers import clean_ri_dataframe, get_figure_dir, inverse_model_log, read_ri_results


plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "legend.title_fontsize": 16,
})


def main():
    cfg = WorkflowConfig()
    fig_dir = get_figure_dir(cfg)

    ri_df = read_ri_results(cfg.results_dir)
    ri_df = clean_ri_dataframe(ri_df, min_slope=25.0, drop_indices=[38, 40, 44])

    cbf_colors = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2"]
    cohesions = sorted(ri_df["Cohesion"].unique(), reverse=True)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_facecolor("#f0f0f0")
    ax.set_axisbelow(True)
    ax.minorticks_on()

    ax.grid(True, which="major", axis="y", alpha=0.4)
    ax.grid(True, which="minor", axis="y", alpha=0.2)
    ax.grid(True, which="major", axis="x", color="#bdbdbd", linewidth=0.8)

    for i, cohesion in enumerate(cohesions):
        group = ri_df[ri_df["Cohesion"] == cohesion].copy()
        x = group["Avg_Slope_deg"].to_numpy()
        y = group["Year"].to_numpy()
        color = cbf_colors[i % len(cbf_colors)]

        ax.scatter(x, y, color=color, s=60, alpha=0.75)

        if len(x) < 3:
            continue

        popt, _ = curve_fit(
            inverse_model_log,
            x,
            np.log10(y),
            p0=[3.0, -25.0],
            maxfev=10000,
        )

        loga_fit, b_fit = popt
        a_fit = 10 ** loga_fit

        x_fit = np.linspace(x.min() - 0.5, x.max() + 0.5, 300)
        y_fit = a_fit / (x_fit + b_fit)

        ax.plot(x_fit, y_fit, color=color, lw=2.5, linestyle="--")

    ax.set_yscale("log")
    ax.set_ylim(100, 1e4)
    ax.set_yticks([100, 1000, 10000])
    ax.yaxis.set_major_formatter(LogFormatterMathtext())

    ax.set_xlabel(r"Hollow slope, $\theta_H$ (°)", fontweight="bold")
    ax.set_ylabel("Recurrence interval, RI (years)", fontweight="bold")

    xticks = np.arange(
        np.floor(ri_df["Avg_Slope_deg"].min() / 3) * 3,
        ri_df["Avg_Slope_deg"].max() + 1,
        3,
    )
    ax.set_xticks(xticks)

    cohesion_handles = [
        Line2D(
            [], [], marker="o", linestyle="None",
            color=cbf_colors[i % len(cbf_colors)],
            markersize=8,
            label=f"{int(cohesion)} Pa",
        )
        for i, cohesion in enumerate(cohesions)
    ]

    ax.legend(handles=cohesion_handles, title="Cohesion", loc="upper right", frameon=True)

    fig.tight_layout()
    out_path = fig_dir / "RI_3D_final.png"
    fig.savefig(out_path, dpi=450, bbox_inches="tight")
    plt.show()

    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
