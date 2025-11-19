#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AIIAC 5–15 Hz band power → Heatmap (R×C)  [Blue colormap, colorbar min=0]
-----------------------------------------------------------------------
• Scans a folder for files named: AIIAC_R###_C##.txt
• Computes 5–15 Hz band power (mV^2) and plots an 11×11 heatmap
• Colorbar uses a blue gradient ("Blues") and starts at 0
• Saves CSV + PNG next to the input folder
"""

from __future__ import annotations
from pathlib import Path
import argparse
import re
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ======================= USER CONFIG (edit once) ============================
ROOT_DIR: str = "results_txt_RB-AII50%_gap200%"                 # "" → use this script's folder; or set absolute path
PATTERN: str = "AIIAC_R*_C*.txt"   # filename glob pattern
RECURSIVE: bool = False            # also search subfolders
# Analysis window & band (ms and Hz)
CROP_MS: Tuple[float, float] = (1000.0, 6000.0)
BAND_HZ: Tuple[float, float] = (5.0, 15.0)
DETREND_LINEAR: bool = False       # optional linear detrend (after mean removal)
USE_HANN: bool = True
# Rod/ Cone grids (fixed to 11×11 as per project)
ROD_COUNTS: List[int] = [1, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400]
CONE_LEVELS: List[int] = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
# Plot controls
INVERT_X: bool = True             # True: X axis 100→0 (to match past figures)
FIGSIZE: Tuple[float, float] = (6.8, 6.2)
DPI: int = 150
# Colorbar controls
CB_CMAP: str = "Blues"            # blue shades
CB_VMIN: float = 0.0              # start colorbar at 0 (band power ≥ 0 by definition)
CB_VMAX: Optional[float] = None   # set to a number to fix the upper bound; None → auto
# Output names
OUT_CSV: str = "AIIAC_bandpower_matrix_5-15Hz.csv"
OUT_PNG: str = "AIIAC_bandpower_heatmap_5-15Hz.png"
# ===========================================================================

FILENAME_RE = re.compile(r"^AIIAC_R(\d{3})_C(\d{2})\.txt$", re.IGNORECASE)


def linear_detrend(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - (A @ coef)


def periodogram_psd_onesided(x: np.ndarray, fs: float, use_hann: bool = True, remove_mean: bool = True):
    N = len(x)
    if remove_mean:
        x = x - np.mean(x)
    w = np.hanning(N) if use_hann else np.ones(N)
    xw = x * w
    U = (w**2).sum() / N
    X = np.fft.rfft(xw)
    Pxx_two = (np.abs(X) ** 2) / (fs * N * U)
    Pxx = Pxx_two.copy()
    if N % 2 == 0:
        Pxx[1:-1] *= 2.0
    else:
        Pxx[1:] *= 2.0
    f = np.fft.rfftfreq(N, d=1 / fs)
    return f, Pxx


def integrate_band(f: np.ndarray, Pxx: np.ndarray, fmin: float, fmax: float) -> float:
    if fmax <= fmin:
        return 0.0
    idx = (f >= fmin) & (f <= fmax)
    if not np.any(idx):
        return 0.0
    return float(np.trapz(Pxx[idx], f[idx]))


def load_aiiac(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim != 2 or data.shape[1] < 2:
        raise RuntimeError(f"File {path} does not have two numeric columns.")
    return data[:, 0], data[:, 1]


def pick_folder(cli_folder: Optional[Path]) -> Path:
    if cli_folder is not None:
        return cli_folder
    if ROOT_DIR.strip():
        return Path(ROOT_DIR)
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path.cwd()


def iter_files(folder: Path, pattern: str, recursive: bool) -> List[Path]:
    it = folder.rglob(pattern) if recursive else folder.glob(pattern)
    return sorted([p for p in it if p.is_file()], key=lambda p: p.name)


def rc_from_name(path: Path) -> Optional[Tuple[int, int]]:
    m = FILENAME_RE.match(path.name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def compute_bandpower_for_file(path: Path, crop_ms: Tuple[float, float], band_hz: Tuple[float, float], detrend_linear: bool, use_hann: bool) -> float:
    t_ms, v_mV = load_aiiac(path)
    tmin, tmax = crop_ms
    mask = (t_ms >= tmin) & (t_ms <= tmax)
    if not np.any(mask):
        raise RuntimeError(f"No samples in crop range {tmin}-{tmax} ms for {path.name}")
    t = t_ms[mask]
    v = v_mV[mask]
    dt_ms = np.median(np.diff(t))
    fs = 1000.0 / dt_ms
    if detrend_linear:
        v = linear_detrend(v, t * 1e-3)
    f, Pxx = periodogram_psd_onesided(v, fs, use_hann=use_hann, remove_mean=True)
    return integrate_band(f, Pxx, band_hz[0], band_hz[1])


def build_matrix(powers: Dict[Tuple[int, int], float], rod_counts: List[int], cone_levels: List[int]) -> np.ndarray:
    Z = np.full((len(rod_counts), len(cone_levels)), np.nan, dtype=float)
    for i, R in enumerate(rod_counts):
        for j, C in enumerate(cone_levels):
            val = powers.get((R, C))
            if val is not None:
                Z[i, j] = val
    return Z


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Compute 5–15 Hz band power and plot a blue heatmap with colorbar min=0.")
    ap.add_argument("folder", nargs="?", type=Path, default=None, help="Folder with AIIAC_R*_C*.txt (default: ROOT_DIR or script folder)")
    ap.add_argument("--pattern", default=PATTERN)
    ap.add_argument("--recursive", action="store_true", default=RECURSIVE)
    ap.add_argument("--crop", nargs=2, type=float, metavar=("MS_MIN", "MS_MAX"), default=None)
    ap.add_argument("--band", nargs=2, type=float, metavar=("HZ_MIN", "HZ_MAX"), default=None)
    ap.add_argument("--detrend-linear", action="store_true", default=DETREND_LINEAR)
    ap.add_argument("--no-hann", action="store_true")
    ap.add_argument("--invert-x", action="store_true", default=INVERT_X)
    ap.add_argument("--out-csv", default=OUT_CSV)
    ap.add_argument("--out-png", default=OUT_PNG)

    args = ap.parse_args(argv)

    folder = pick_folder(args.folder)
    if not folder.exists() or not folder.is_dir():
        print(f"Error: {folder} is not a directory")
        return 2

    pattern = args.pattern or PATTERN
    recursive = bool(args.recursive)
    crop_ms = tuple(args.crop) if args.crop else CROP_MS
    band_hz = tuple(args.band) if args.band else BAND_HZ
    detrend_lin = bool(args.detrend_linear)
    use_hann = not bool(args.no_hann)
    invert_x = bool(args.invert_x)

    files = iter_files(folder, pattern, recursive)
    if not files:
        print(f"No files matched pattern '{pattern}' in {folder}")
        return 0

    # Compute band power per file
    powers: Dict[Tuple[int, int], float] = {}
    for p in files:
        rc = rc_from_name(p)
        if rc is None:
            continue
        try:
            bp = compute_bandpower_for_file(p, crop_ms, band_hz, detrend_lin, use_hann)
        except Exception as e:
            print(f"[SKIP] {p.name}: {e}")
            continue
        powers[rc] = bp

    # Build matrix (rows=rod ascending, cols=cone ascending)
    Z = build_matrix(powers, ROD_COUNTS, CONE_LEVELS)

    # Export CSV (columns as Cone %, rows as Rod %)
    rod_pct = [r / max(ROD_COUNTS) * 100 for r in ROD_COUNTS]
    cone_pct_asc = [c / max(CONE_LEVELS) * 100 if max(CONE_LEVELS) > 0 else 0 for c in CONE_LEVELS]
    df = pd.DataFrame(Z, index=[f"{rp:.2f}" for rp in rod_pct], columns=[f"{cp:.0f}" for cp in cone_pct_asc])
    df.index.name = "Rod_%"
    df.to_csv(folder / args.out_csv)

    # Prepare plotting data
    if invert_x:
        Z_plot = Z[:, ::-1]
        cone_pct_plot = list(reversed(cone_pct_asc))
    else:
        Z_plot = Z
        cone_pct_plot = cone_pct_asc

    Zm = np.ma.masked_invalid(Z_plot)  # ignore NaNs in autoscale

    # Layout with matched-height colorbar
    CB_WIDTH_RATIO = 0.05
    fig = plt.figure(figsize=FIGSIZE, constrained_layout=True, dpi=DPI)
    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=[1, CB_WIDTH_RATIO])
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    im = ax.imshow(
        Zm,
        origin="lower",
        extent=(-5, 105, -5, 105),
        interpolation="nearest",
        cmap=CB_CMAP,
        vmin=CB_VMIN,
        vmax=CB_VMAX,
    )
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Band power 5–15 Hz (mV$^2$)")

    ax.set_xlabel("Cone survival rate (%)")
    ax.set_ylabel("Rod survival rate (%)")
    ax.set_xticks(np.arange(0, 101, 10))
    ax.set_yticks(np.arange(0, 101, 10))
    ax.set_aspect("equal", adjustable="box")

    if invert_x:
        ax.invert_xaxis()

    title = f"AIIAC band power {band_hz[0]:g}–{band_hz[1]:g} Hz (window {crop_ms[0]:g}–{crop_ms[1]:g} ms)"
    ax.set_title(title)

    out_png = folder / args.out_png
    fig.savefig(out_png, bbox_inches="tight")

    print(f"Saved PNG: {out_png}")
    print(f"Saved CSV: {folder / args.out_csv}")

    expected = len(ROD_COUNTS) * len(CONE_LEVELS)
    got = np.isfinite(Z).sum()
    if got < expected:
        print(f"Note: computed {got}/{expected} cells. Missing files will be NaN in the CSV/heatmap.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
