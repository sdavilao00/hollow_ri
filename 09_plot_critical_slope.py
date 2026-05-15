# -*- coding: utf-8 -*-
"""
09_plot_critical_slope.py

Calculate and plot the critical hollow slope required for FS = 1 across a range
of saturation ratios and cohesion values.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import WorkflowConfig
from plot_helpers import get_figure_dir


# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_CALCULATION = True
LOAD_EXISTING_TABLE_IF_AVAILABLE = True

ZMIN = 0.01
ZMAX = 8.0
NZ = 2200
M_VALUES = np.round(np.linspace(0.5, 1.0, 51), 2)
COHESIONS = [760, 1920, 3600, 6400]

SAVE_FIGURE = True
SAVE_TABLE = True


plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "legend.title_fontsize": 11,
})


def minimum_fs(theta_deg, saturation, cohesion, cfg, zs):
    """Return minimum FS over the supplied soil-depth range."""
    theta = np.deg2rad(theta_deg)
    sin_t = np.sin(theta)
    cos_t = np.cos(theta)

    yw = cfg.g * cfg.pw
    ys = cfg.g * cfg.rho_s
    phi = np.deg2rad(cfg.phi_deg)
    tan_phi = np.tan(phi)

    Kp = np.tan(np.deg2rad(45.0) + phi / 2.0) ** 2
    Ka = np.tan(np.deg2rad(45.0) - phi / 2.0) ** 2

    exp_term = np.exp(-zs * cfg.j)
    Crb = cohesion * exp_term
    Crl = (cohesion / (cfg.j * zs)) * (1.0 - exp_term)
    K0 = 1.0 - sin_t

    Frb = (Crb + cos_t ** 2 * zs * (ys - yw * saturation) * tan_phi) * cfg.l * cfg.w
    Frc = (
        Crl
        + K0 * 0.5 * zs * (ys - yw * saturation ** 2) * tan_phi
    ) * (cos_t * zs * cfg.l * 2.0)
    Frddu = (Kp - Ka) * 0.5 * zs ** 2 * (ys - yw * saturation ** 2) * cfg.w
    Fdc = sin_t * cos_t * zs * ys * cfg.l * cfg.w

    fs = (Frb + Frc + Frddu) / Fdc
    return float(np.nanmin(fs))


def theta_crit_for_m(saturation, cohesion, cfg, zs, theta_min=10.0, theta_max=80.0, coarse_step=1.0, tol=0.05):
    """Find the first slope angle where minimum FS reaches 1."""
    thetas = np.arange(theta_min, theta_max + coarse_step, coarse_step)
    fvals = np.array([minimum_fs(theta, saturation, cohesion, cfg, zs) - 1.0 for theta in thetas])

    if np.nanmin(fvals) > 0:
        return np.nan

    idx = None
    for i in range(len(thetas) - 1):
        if fvals[i] > 0 and fvals[i + 1] <= 0:
            idx = i
            break

    if idx is None:
        return float(theta_min)

    lo = float(thetas[idx])
    hi = float(thetas[idx + 1])

    while (hi - lo) > tol:
        mid = 0.5 * (lo + hi)
        fmid = minimum_fs(mid, saturation, cohesion, cfg, zs) - 1.0
        if fmid <= 0:
            hi = mid
        else:
            lo = mid

    return 0.5 * (lo + hi)


def calculate_critical_slope_table(cfg):
    zs = np.linspace(ZMIN, ZMAX, NZ)
    records = []

    for cohesion in COHESIONS:
        for saturation in M_VALUES:
            theta_crit = theta_crit_for_m(saturation, cohesion, cfg, zs)
            records.append({
                "Cohesion": int(cohesion),
                "m": float(saturation),
                "theta_crit_deg": theta_crit,
            })

    return pd.DataFrame(records)


def main():
    cfg = WorkflowConfig()
    fig_dir = get_figure_dir(cfg)

    table_path = cfg.results_dir / "critical_slopes.csv"

    if LOAD_EXISTING_TABLE_IF_AVAILABLE and table_path.exists():
        crit_df = pd.read_csv(table_path)
    elif RUN_CALCULATION:
        crit_df = calculate_critical_slope_table(cfg)
        if SAVE_TABLE:
            crit_df.to_csv(table_path, index=False)
            print(f"Saved critical slope table: {table_path}")
    else:
        raise FileNotFoundError(
            f"Critical slope table not found: {table_path}. "
            "Set RUN_CALCULATION = True to create it."
        )

    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    for cohesion in COHESIONS:
        sub = crit_df[crit_df["Cohesion"] == int(cohesion)]
        ax.plot(
            sub["m"],
            sub["theta_crit_deg"],
            linewidth=2,
            label=fr"$C_0$ = {int(cohesion)} Pa",
        )

    ax.axvline(x=0.85, color="gray", linestyle="--", linewidth=1.5, label=r"$m$ = 0.85")
    ax.axvline(x=1.00, color="black", linestyle="--", linewidth=1.5, label=r"$m$ = 1.00")

    ax.set_xlabel("Saturation ratio", fontweight="bold")
    ax.set_ylabel(r"Critical slope, $\theta_{\mathrm{crit}}$ (°)", fontweight="bold")

    xticks = list(np.arange(0.5, 1.01, 0.1))
    xticks.append(0.85)
    xticks = sorted(set(np.round(xticks, 2)))
    xtick_labels = [f"{x:.2f}" if np.isclose(x, 0.85) else f"{x:.1f}" for x in xticks]

    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)

    for tick, value in zip(ax.get_xticklabels(), xticks):
        if np.isclose(value, 0.85) or np.isclose(value, 1.0):
            tick.set_fontweight("bold")

    ax.grid(True, alpha=0.3)

    handles, labels = ax.get_legend_handles_labels()
    cohesion_handles = handles[:len(COHESIONS)][::-1]
    cohesion_labels = labels[:len(COHESIONS)][::-1]
    saturation_handles = handles[len(COHESIONS):]
    saturation_labels = labels[len(COHESIONS):]

    ax.legend(
        cohesion_handles + saturation_handles,
        cohesion_labels + saturation_labels,
        frameon=True,
    )

    fig.tight_layout()

    if SAVE_FIGURE:
        fig.savefig(fig_dir / "critical_slope.png", dpi=450, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()
