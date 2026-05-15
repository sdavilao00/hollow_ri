# -*- coding: utf-8 -*-
"""
Shared plotting utilities for recurrence interval manuscript figures.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd


def get_figure_dir(cfg):
    """Create and return the manuscript figure output directory."""
    fig_dir = cfg.base_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    return fig_dir


def read_ri_results(results_dir, pattern="optimal_buffer_results_interpolated_*.csv"):
    """Read optimal-buffer/RI CSVs and add Extent, Cohesion, and m if missing."""
    results_dir = Path(results_dir)
    csv_files = sorted(results_dir.glob(pattern))

    if not csv_files:
        raise FileNotFoundError(f"No RI result CSVs found in: {results_dir}")

    extent_pattern = re.compile(r"(ext\d+)", re.IGNORECASE)
    cohesion_pattern = re.compile(r"_(\d+)(?:_m|_mp|\.csv)", re.IGNORECASE)
    m_pattern = re.compile(r"_m(\d+p?\d*)", re.IGNORECASE)

    frames = []
    for csv_path in csv_files:
        name = csv_path.name
        df = pd.read_csv(csv_path)

        extent_match = extent_pattern.search(name)
        cohesion_match = cohesion_pattern.search(name)
        m_match = m_pattern.search(name)

        if "Extent" not in df.columns and extent_match:
            df["Extent"] = extent_match.group(1)

        if "Cohesion" not in df.columns and cohesion_match:
            df["Cohesion"] = int(cohesion_match.group(1))

        if "m" not in df.columns and m_match:
            df["m"] = float(m_match.group(1).replace("p", "."))

        frames.append(df)

    return pd.concat(frames, ignore_index=True)


def clean_ri_dataframe(df, min_slope=25.0, drop_indices=None):
    """Apply common filtering used for manuscript RI plots."""
    out = df.copy()

    needed = ["Avg_Slope_deg", "Year", "Avg_Soil_Depth_m", "Cohesion"]
    out = out.dropna(subset=[col for col in needed if col in out.columns])

    out = out[out["Avg_Slope_deg"] > min_slope]
    out = out[out["Year"] > 0]

    if drop_indices is not None:
        out = out.drop(drop_indices, errors="ignore")

    return out.reset_index(drop=True)


def inverse_model_log(theta, loga, b):
    """Log-space inverse recurrence interval model."""
    return loga - np.log10(theta + b)


def calculate_log_r2(y_observed, y_predicted):
    """Calculate R² in log10 space."""
    y_log = np.log10(y_observed)
    y_log_pred = np.log10(y_predicted)
    ss_res = np.sum((y_log - y_log_pred) ** 2)
    ss_tot = np.sum((y_log - np.mean(y_log)) ** 2)
    return 1 - ss_res / ss_tot


def add_critical_area_and_volume(df, cfg, saturation=1.0):
    """Calculate MD-STAB critical area and failure volume for each RI result."""
    out = df.copy()

    z = out["Avg_Soil_Depth_m"]
    hollow_rad = np.deg2rad(out["Avg_Slope_deg"])
    C0 = out["Cohesion"]

    yw = cfg.g * cfg.pw
    ys = cfg.g * cfg.rho_s
    phi = np.deg2rad(cfg.phi_deg)

    Crb = C0 * np.exp(-z * cfg.j)
    Crl = (C0 / (cfg.j * z)) * (1 - np.exp(-z * cfg.j))

    K0 = 1 - np.sin(hollow_rad)
    Kp = np.tan(np.deg2rad(45) + phi / 2) ** 2
    Ka = np.tan(np.deg2rad(45) - phi / 2) ** 2

    A = (
        2 * Crl * z
        + K0 * (z ** 2) * (ys - yw * saturation ** 2) * np.tan(phi)
    ) * np.cos(hollow_rad) * (cfg.l / cfg.w) ** 0.5

    B = (
        (Kp - Ka)
        * 0.5
        * (z ** 2)
        * (ys - yw * saturation ** 2)
        * (cfg.l / cfg.w) ** (-0.5)
    )

    C = (
        np.sin(hollow_rad) * np.cos(hollow_rad) * z * ys
        - Crb
        - (np.cos(hollow_rad) ** 2) * z * (ys - yw * saturation) * np.tan(phi)
    )

    out["Ac"] = ((A + B) / C) ** 2
    out["Volume"] = out["Ac"] * out["Avg_Soil_Depth_m"]

    return out
