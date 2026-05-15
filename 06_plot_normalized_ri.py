# -*- coding: utf-8 -*-
"""
06_plot_normalized_ri.py

Plot recurrence interval normalized by cohesion as a function of hollow slope.
Run this after 02_extract_and_calculate_RI.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from config import WorkflowConfig
from plot_helpers import clean_ri_dataframe, get_figure_dir, inverse_model_log, read_ri_results


plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
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
    ri_df["RI_norm"] = ri_df["Year"] / ri_df["Cohesion"]

    x = ri_df["Avg_Slope_deg"].to_numpy()
    y = ri_df["RI_norm"].to_numpy()

    popt, _ = curve_fit(
        inverse_model_log,
        x,
        np.log10(y),
        p0=[-2, -25],
        maxfev=10000,
    )

    loga_fit, b_fit = popt
    a_fit = 10 ** loga_fit

    y_log = np.log10(y)
    y_log_pred = inverse_model_log(x, loga_fit, b_fit)
    ss_res = np.sum((y_log - y_log_pred) ** 2)
    ss_tot = np.sum((y_log - np.mean(y_log)) ** 2)
    r2_log = 1 - ss_res / ss_tot

    x_fit = np.linspace(x.min() - 0.5, x.max() + 0.5, 300)
    y_fit = a_fit / (x_fit + b_fit)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_facecolor("#f0f0f0")
    ax.minorticks_off()
    ax.set_axisbelow(True)
    ax.grid(True, which="major", axis="both", color="#bdbdbd", linewidth=0.8, alpha=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.scatter(x, y, color="#0072B2", s=60, alpha=0.6)
    ax.plot(
        x_fit,
        y_fit,
        color="#D55E00",
        lw=3,
        linestyle="--",
        label=fr"$RI/C = \frac{{{a_fit:.3f}}}{{\theta_H + {b_fit:.2f}}}$" + "\n" + fr"$R^2_{{\log}} = {r2_log:.2f}$",
    )

    ax.set_xlabel(r"Hollow slope, $\theta_H$ (°)", fontweight="bold")
    ax.set_ylabel("RI / Cohesion", fontweight="bold")

    xticks = np.arange(
        np.floor(ri_df["Avg_Slope_deg"].min() / 3) * 3,
        ri_df["Avg_Slope_deg"].max() + 1,
        3,
    )
    ax.set_xticks(xticks)
    ax.legend(frameon=True)

    fig.tight_layout()
    out_path = fig_dir / "RI_normalized.png"
    fig.savefig(out_path, dpi=450, bbox_inches="tight")
    plt.show()

    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
