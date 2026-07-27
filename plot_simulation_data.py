#!/usr/bin/env python3
"""
High-quality plotting script for vehicle/trajectory simulation data.

Expects CSV with columns:
t,x,y,z,vx,vy,vz,wx,wy,wz,qw,qx,qy,qz,lat,lon,alt

- Position/velocity in ECEF (meters, m/s)
- wx,wy,wz : body-frame angular velocity (rad/s)
- qw,qx,qy,qz : quaternion rotating body vectors -> ECEF frame (unit quaternion)
- lat,lon (deg), alt (m)

Modes:
  png  – generate high-DPI PNG plots (default)
  live – interactive ECEF attitude replay of the body axes triad

Body-frame convention (initial / local surface frame):
  +x up, +y east, +z north  (right-handed)

Usage:
    python plot_simulation_data.py --mode png --input sim.csv --output-dir plots --dpi 300
    python plot_simulation_data.py --mode live --input sim.csv

Requirements:
    numpy pandas matplotlib scipy
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.spatial.transform import Rotation


REQUIRED_COLS = [
    "t", "x", "y", "z", "vx", "vy", "vz", "wx", "wy", "wz",
    "qw", "qx", "qy", "qz", "lat", "lon", "alt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot simulation CSV data as PNGs, or replay live body attitude "
            "in an ECEF-fixed 3D view."
        )
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["png", "live"],
        default="png",
        help=(
            "png: write high-quality PNG plots to --output-dir (default). "
            "live: interactive real-time attitude replay window."
        ),
    )
    parser.add_argument(
        "--input", "-i", type=str, default="simulation_data.csv",
        help="Path to input CSV file with simulation data",
    )
    parser.add_argument(
        "--output-dir", "-o", type=str, default="plots",
        help="Directory where PNG plots will be saved (created if needed; png mode only)",
    )
    parser.add_argument(
        "--dpi", type=int, default=300,
        help="Resolution (DPI) for saved PNG figures (png mode only)",
    )
    parser.add_argument(
        "--figsize", type=float, nargs=2, default=[12.0, 6.5],
        help="Base figure size in inches (width height). Some plots adjust automatically.",
    )
    parser.add_argument(
        "--target-fps", type=float, default=30.0,
        help=(
            "Target display frame rate for live replay. Frame skip is chosen "
            "automatically so wall-clock time tracks simulation time 1-to-1."
        ),
    )
    return parser.parse_args()


def load_simulation_csv(input_path: Path) -> pd.DataFrame:
    """Load CSV and validate required columns."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print(f"Reading simulation data from: {input_path}")
    df = pd.read_csv(input_path)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    n_points = len(df)
    print(f"Loaded {n_points} data points")

    df = df.copy()
    df["speed"] = np.sqrt(df["vx"] ** 2 + df["vy"] ** 2 + df["vz"] ** 2)
    df["q_norm"] = np.sqrt(df["qw"] ** 2 + df["qx"] ** 2 + df["qy"] ** 2 + df["qz"] ** 2)
    df["omega_norm"] = np.sqrt(df["wx"] ** 2 + df["wy"] ** 2 + df["wz"] ** 2)
    df["radial_dist"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2 + df["z"] ** 2)

    t_min, t_max = df["t"].min(), df["t"].max()
    alt_min, alt_max = df["alt"].min(), df["alt"].max()
    speed_min, speed_max = df["speed"].min(), df["speed"].max()
    qmin, qmax = df["q_norm"].min(), df["q_norm"].max()

    print(f"  Time range:     {t_min:.3f}  ->  {t_max:.3f} s")
    print(f"  Altitude range: {alt_min:.1f}  ->  {alt_max:.1f} m")
    print(f"  Speed range:    {speed_min:.2f}  ->  {speed_max:.2f} m/s")
    print(f"  q_norm range:   {qmin:.6f}  ->  {qmax:.6f}  (should be ~1.0)")
    return df


def compute_frame_skip(t: np.ndarray, target_fps: float = 30.0) -> tuple[int, float, float]:
    """
    Automatically choose sample stride so replay is wall-clock 1-to-1.

    For log step dt and target display rate f, skip ≈ (1/f) / dt samples.
    FuncAnimation interval is set to skip * dt seconds so displayed sim time
    advances at the same rate as real time.

    Returns
    -------
    skip : int
        Number of log rows between consecutive displayed frames (≥ 1).
    interval_ms : float
        Wall-clock delay between animation frames (milliseconds).
    dt : float
        Median simulation time step (seconds).
    """
    if target_fps <= 0:
        raise ValueError("--target-fps must be positive")

    if len(t) < 2:
        return 1, 1000.0 / target_fps, 1.0 / target_fps

    dt = float(np.median(np.diff(t.astype(float))))
    if not np.isfinite(dt) or dt <= 0:
        span = float(t[-1] - t[0])
        dt = span / max(len(t) - 1, 1) if span > 0 else 1.0 / target_fps

    # Desired wall time between UI frames ≈ 1/target_fps; match that in sim time.
    skip = max(1, int(round((1.0 / target_fps) / dt)))
    interval_ms = 1000.0 * skip * dt
    return skip, interval_ms, dt


def quaternions_to_body_axes_ecef(
    qw: np.ndarray, qx: np.ndarray, qy: np.ndarray, qz: np.ndarray
) -> np.ndarray:
    """
    Convert body->ECEF quaternions into body basis vectors expressed in ECEF.

    Parameters
    ----------
    qw, qx, qy, qz : array-like, shape (N,)
        Unit quaternions (scalar-first / Eigen convention in the log).

    Returns
    -------
    axes : ndarray, shape (N, 3, 3)
        axes[i, :, j] is body axis j (0=x, 1=y, 2=z) in ECEF at sample i.
    """
    # SciPy uses [x, y, z, w]
    rot = Rotation.from_quat(np.column_stack([qx, qy, qz, qw]))
    # Columns of the rotation matrix are body axes mapped into ECEF.
    return rot.as_matrix()  # (N, 3, 3)


def live_attitude_replay(df: pd.DataFrame, target_fps: float = 30.0) -> int:
    """
    Interactive ECEF-fixed attitude replay.

    Draws the rocket as a body-fixed triad:
      red   = body +x  (up at launch for surface-aligned vehicles)
      green = body +y  (east at launch)
      blue  = body +z  (north at launch)

    Axes are rotated by the logged body->ECEF quaternion each frame.
    Frame skip is chosen for 1:1 wall-clock timing. Restart restarts the
    sequence; the process blocks until the user closes the window.
    """
    t = df["t"].to_numpy(dtype=float)
    alt = df["alt"].to_numpy(dtype=float)
    axes_ecef = quaternions_to_body_axes_ecef(
        df["qw"].to_numpy(dtype=float),
        df["qx"].to_numpy(dtype=float),
        df["qy"].to_numpy(dtype=float),
        df["qz"].to_numpy(dtype=float),
    )

    skip, interval_ms, dt = compute_frame_skip(t, target_fps=target_fps)
    indices = np.arange(0, len(t), skip, dtype=int)
    if indices[-1] != len(t) - 1:
        indices = np.append(indices, len(t) - 1)

    duration = float(t[-1] - t[0])
    n_frames = len(indices)
    effective_fps = 1000.0 / interval_ms if interval_ms > 0 else target_fps

    print("\nLive attitude replay (ECEF-fixed view)")
    print(f"  Simulation duration: {duration:.3f} s")
    print(f"  Log dt (median):     {dt * 1000:.4f} ms")
    print(f"  Target display FPS:  {target_fps:.1f}")
    print(f"  Auto frame skip:     {skip} sample(s)")
    print(f"  Animation interval:  {interval_ms:.2f} ms  ({effective_fps:.1f} FPS)")
    print(f"  Displayed frames:    {n_frames}")
    print(f"  Expected wall time:  ~{duration:.1f} s (1-to-1)")
    print("  Controls: [Restart] button restarts; close the window to exit.")

    # Axis length for the body triad (unitless ECEF orientation view)
    axis_len = 1.0
    ref_len = 1.15

    plt.rcParams.update({
        "font.size": 10,
        "figure.facecolor": "white",
        "axes.facecolor": "#f7f7f7",
    })

    fig = plt.figure(figsize=(10.5, 8.5))
    # Leave room at the bottom for the restart button
    ax = fig.add_axes([0.06, 0.14, 0.88, 0.80], projection="3d")

    # Fixed ECEF reference triad (faint)
    ecef_colors = ("#999999", "#999999", "#999999")
    ecef_labels = ("ECEF +X", "ECEF +Y", "ECEF +Z")
    ecef_dirs = np.eye(3)
    for i in range(3):
        d = ecef_dirs[:, i] * ref_len
        ax.plot([0, d[0]], [0, d[1]], [0, d[2]],
                color=ecef_colors[i], linewidth=1.2, linestyle="--", alpha=0.55)
        ax.text(d[0] * 1.05, d[1] * 1.05, d[2] * 1.05,
                ecef_labels[i], color="#666666", fontsize=8)

    # Body triad lines (updated each frame)
    body_colors = ("#d62728", "#2ca02c", "#1f77b4")  # x red, y green, z blue
    body_names = ("body +X", "body +Y", "body +Z")
    body_lines = []
    for i, color in enumerate(body_colors):
        (line,) = ax.plot(
            [0, axes_ecef[0, 0, i] * axis_len],
            [0, axes_ecef[0, 1, i] * axis_len],
            [0, axes_ecef[0, 2, i] * axis_len],
            color=color, linewidth=3.0, solid_capstyle="round",
            label=body_names[i],
        )
        body_lines.append(line)

    # Arrowhead markers at tips for readability
    tip_scat = ax.scatter(
        [axes_ecef[0, 0, i] * axis_len for i in range(3)],
        [axes_ecef[0, 1, i] * axis_len for i in range(3)],
        [axes_ecef[0, 2, i] * axis_len for i in range(3)],
        c=list(body_colors), s=40, depthshade=False,
    )

    origin = ax.scatter([0], [0], [0], c="k", s=35, depthshade=False, label="origin")

    lim = 1.35
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    try:
        ax.set_box_aspect([1, 1, 1])
    except AttributeError:
        pass

    ax.set_xlabel("ECEF X")
    ax.set_ylabel("ECEF Y")
    ax.set_zlabel("ECEF Z")
    ax.set_title("Live Attitude Replay  —  body triad in ECEF", pad=12)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.92)
    ax.view_init(elev=22, azim=45)
    ax.grid(True, alpha=0.3)

    status = ax.text2D(
        0.02, 0.02,
        f"t = {t[0]:.3f} s   |   alt = {alt[0]:.2f} m",
        transform=ax.transAxes,
        fontsize=11,
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#cccccc", alpha=0.92),
    )

    def set_body_pose(sample_idx: int, done: bool = False) -> None:
        R_b2e = axes_ecef[sample_idx]
        tips = np.zeros((3, 3))
        for j in range(3):
            tip = R_b2e[:, j] * axis_len
            tips[:, j] = tip
            body_lines[j].set_data_3d([0.0, tip[0]], [0.0, tip[1]], [0.0, tip[2]])
        tip_scat._offsets3d = (tips[0], tips[1], tips[2])
        alt_m = float(alt[sample_idx])
        if done:
            status.set_text(
                f"t = {t[sample_idx]:7.3f} s   |   alt = {alt_m:8.2f} m"
                f"   |   DONE  —  press Restart or close window"
            )
        else:
            status.set_text(
                f"t = {t[sample_idx]:7.3f} s   |   alt = {alt_m:8.2f} m"
                f"   |   sample {sample_idx}/{len(t) - 1}"
                f"   |   skip={skip}"
            )

    def init_anim():
        set_body_pose(int(indices[0]))
        return (*body_lines, tip_scat, status)

    def update_anim(frame_i: int):
        sample_idx = int(indices[frame_i])
        done = frame_i >= n_frames - 1
        set_body_pose(sample_idx, done=done)
        return (*body_lines, tip_scat, status)

    # interval = skip * dt  ->  wall-clock advances in lockstep with sim time (1-to-1)
    anim = FuncAnimation(
        fig,
        update_anim,
        frames=n_frames,
        init_func=init_anim,
        interval=interval_ms,
        blit=False,
        repeat=False,
        cache_frame_data=False,
    )

    # Restart button
    ax_btn = fig.add_axes([0.40, 0.03, 0.20, 0.06])
    btn = Button(ax_btn, "Restart", color="#e8e8e8", hovercolor="#cfe8ff")

    def on_restart(_event) -> None:
        # Rebuild frame sequence and resume the timer
        anim.frame_seq = anim.new_frame_seq()
        set_body_pose(int(indices[0]))
        fig.canvas.draw_idle()
        if anim.event_source is not None:
            anim.event_source.stop()
            anim.event_source.start()

    btn.on_clicked(on_restart)

    # Keep references alive for the duration of the window
    fig._attitude_anim = anim  # type: ignore[attr-defined]
    fig._attitude_btn = btn    # type: ignore[attr-defined]

    print("Opening live replay window (close the window to exit)...")
    plt.show(block=True)
    print("Live replay window closed.")
    return 0


def generate_png_plots(df: pd.DataFrame, args: argparse.Namespace) -> int:
    """Generate the existing high-quality PNG suite."""
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    n_points = len(df)

    # Matplotlib quality settings
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 13,
        "legend.fontsize": 9,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.edgecolor": "none",
        "axes.grid": True,
        "grid.alpha": 0.35,
        "grid.linestyle": "--",
        "lines.linewidth": 1.7,
        "figure.autolayout": False,
    })

    def save_figure(fig: plt.Figure, filename: str):
        """Save figure with high quality and close it."""
        filepath = output_dir / filename
        fig.savefig(
            filepath,
            dpi=args.dpi,
            bbox_inches="tight",
            pad_inches=0.15,
            facecolor="white",
            edgecolor="none",
        )
        plt.close(fig)
        print(f"\tSaved {filename}")

    t = df["t"].values

    # 1. ECEF Position Components vs Time (stacked subplots)
    fig, axs = plt.subplots(
        3, 1, figsize=(args.figsize[0], 9.5), sharex=True, constrained_layout=True
    )
    components = [("x", "X ECEF"), ("y", "Y ECEF"), ("z", "Z ECEF")]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for ax, (col, label), color in zip(axs, components, colors):
        ax.plot(t, df[col], color=color, linewidth=1.8)
        ax.set_ylabel(f"{label} (m)")
        ax.grid(True, alpha=0.35)
        if df[col].abs().max() > 1e5:
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
    axs[-1].set_xlabel("Time (s)")
    fig.suptitle("ECEF Position Components vs Time", fontsize=14, y=0.995)
    save_figure(fig, "01_ecef_position_vs_time.png")

    # 2. ECEF Velocity Components + Speed vs Time
    fig, axs = plt.subplots(
        4, 1, figsize=(args.figsize[0], 9.5), sharex=True, constrained_layout=True
    )
    components = [("vx", "VX ECEF"), ("vy", "VY ECEF"), ("vz", "VZ ECEF"), ("speed", "||V||")]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
    for ax, (col, label), color in zip(axs, components, colors):
        ax.plot(t, df[col], color=color, linewidth=1.8)
        ax.set_ylabel(f"{label} (m/s)")
        ax.grid(True, alpha=0.35)
        if df[col].abs().max() > 1e5:
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
    axs[-1].set_xlabel("Time (s)")
    fig.suptitle("ECEF Velocity Components and Speed vs Time")
    save_figure(fig, "02_ecef_velocity_vs_time.png")

    # 3. Body Angular Velocity (ω) vs Time
    fig, axs = plt.subplots(
        4, 1, figsize=(args.figsize[0], 9.5), sharex=True, constrained_layout=True
    )
    components = [("wx", "ωX Body"), ("wy", "ωY Body"), ("wz", "ωZ Body"), ("omega_norm", "||ω||")]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
    for ax, (col, label), color in zip(axs, components, colors):
        ax.plot(t, df[col], color=color, linewidth=1.8)
        ax.set_ylabel(f"{label} (rad/s)")
        ax.grid(True, alpha=0.35)
        if df[col].abs().max() > 1e5:
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
    axs[-1].set_xlabel("Time (s)")
    fig.suptitle("Body Angular Velocity Components and Speed vs Time")
    save_figure(fig, "03_body_angular_velocity_vs_time.png")

    # 4. Quaternion Components vs Time (with norm diagnostic)
    fig, ax = plt.subplots(figsize=args.figsize, constrained_layout=True)
    ax.plot(t, df["qw"], label="qw", color="#1f77b4", linewidth=1.8)
    ax.plot(t, df["qx"], label="qx", color="#ff7f0e", linewidth=1.8)
    ax.plot(t, df["qy"], label="qy", color="#2ca02c", linewidth=1.8)
    ax.plot(t, df["qz"], label="qz", color="#d62728", linewidth=1.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Quaternion component value")
    ax.set_title("Attitude Quaternion (Rotation: Body -> ECEF) vs Time")
    ax.legend(loc="upper right", ncol=2, framealpha=0.95)
    ax.grid(True, alpha=0.35)
    ax.set_ylim(-1.05, 1.05)

    ax2 = ax.twinx()
    ax2.plot(t, df["q_norm"], color="#8c564b", linewidth=1.6, linestyle="--", label="|q| norm")
    ax2.axhline(y=1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.7)
    ax2.set_ylabel("Quaternion norm |q|", color="#8c564b")
    ax2.tick_params(axis="y", labelcolor="#8c564b")
    ax2.set_ylim(0.98, 1.02)
    ax2.legend(loc="lower right", framealpha=0.95)

    save_figure(fig, "04_quaternion_vs_time.png")

    # 5. Altitude vs Time
    fig, ax = plt.subplots(figsize=args.figsize, constrained_layout=True)
    ax.plot(t, df["alt"], color="#2ca02c", linewidth=2.0)
    ax.fill_between(t, df["alt"], alpha=0.18, color="#2ca02c")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Altitude (m)")
    ax.set_title("Altitude Above Ellipsoid vs Time")
    ax.grid(True, alpha=0.35)
    save_figure(fig, "05_altitude_vs_time.png")

    # 6. Ground Track (Latitude vs Longitude)
    fig, ax = plt.subplots(figsize=(11, 8), constrained_layout=True)
    ax.plot(df["lon"], df["lat"], color="#1f77b4", linewidth=1.8, alpha=0.85, label="Trajectory path")
    ax.scatter(
        df["lon"].iloc[0], df["lat"].iloc[0],
        c="limegreen", s=140, zorder=10, marker="o", edgecolors="black", linewidths=0.8,
        label="Start"
    )
    ax.scatter(
        df["lon"].iloc[-1], df["lat"].iloc[-1],
        c="red", s=140, zorder=10, marker="o", edgecolors="black", linewidths=0.8,
        label="End"
    )
    step = max(1, n_points // 12)
    ax.scatter(
        df["lon"].iloc[::step], df["lat"].iloc[::step],
        c="orange", s=35, zorder=8, alpha=0.75, edgecolors="none", label="Sample points"
    )
    ax.set_xlabel("Longitude (°)", fontsize=12)
    ax.set_ylabel("Latitude (°)", fontsize=12)
    ax.set_title("Ground Track (Latitude vs Longitude)", fontsize=14)
    ax.grid(True, alpha=0.4, linestyle="--")
    ax.legend(loc="best", framealpha=0.95)
    # Geographic domain limits (degrees)
    lon_lo, lon_hi = -180.0, 180.0
    lat_lo, lat_hi = -90.0, 90.0
    # Minimum longitude (width) and latitude (height) spans so a track with
    # little/no motion in either direction is never collapsed to a line/point.
    min_lon_span_deg = 0.05  # minimum plot width
    min_lat_span_deg = 0.05  # minimum plot height
    lon_min, lon_max = float(df["lon"].min()), float(df["lon"].max())
    lat_min, lat_max = float(df["lat"].min()), float(df["lat"].max())
    lon_c = 0.5 * (lon_min + lon_max)
    lat_c = 0.5 * (lat_min + lat_max)
    lon_span = max(lon_max - lon_min, min_lon_span_deg)
    lat_span = max(lat_max - lat_min, min_lat_span_deg)
    # Prefer a square window when it fits inside the geographic domain
    span = max(lon_span, lat_span)
    lon_span = min(span, lon_hi - lon_lo)
    lat_span = min(span, lat_hi - lat_lo)

    def _window_limits(center: float, half: float, lo: float, hi: float) -> tuple[float, float]:
        """Centered window of half-width `half`, shifted/clamped into [lo, hi]."""
        half = min(half, 0.5 * (hi - lo))
        left = center - half
        right = center + half
        if left < lo:
            right += lo - left
            left = lo
        if right > hi:
            left -= right - hi
            right = hi
        return max(left, lo), min(right, hi)

    lon_left, lon_right = _window_limits(lon_c, 0.5 * lon_span, lon_lo, lon_hi)
    lat_bottom, lat_top = _window_limits(lat_c, 0.5 * lat_span, lat_lo, lat_hi)
    ax.set_xlim(lon_left, lon_right)
    ax.set_ylim(lat_bottom, lat_top)
    try:
        ax.set_aspect("equal", adjustable="box")
    except Exception:
        pass
    save_figure(fig, "06_ground_track.png")

    # 7. 3D Trajectory in ECEF (with reference Earth sphere)
    fig = plt.figure(figsize=(10.5, 9.5), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(
        df["x"], df["y"], df["z"],
        color="#1f77b4", linewidth=2.2, label="Trajectory", zorder=6, alpha=0.95
    )
    ax.scatter(
        [df["x"].iloc[0]], [df["y"].iloc[0]], [df["z"].iloc[0]],
        c="limegreen", s=90, label="Start", zorder=10, edgecolors="black", linewidths=0.7
    )
    ax.scatter(
        [df["x"].iloc[-1]], [df["y"].iloc[-1]], [df["z"].iloc[-1]],
        c="red", s=90, label="End", zorder=10, edgecolors="black", linewidths=0.7
    )

    # Reference Earth ellipsoid (WGS84). Low alpha so trajectory remains visible.
    a = 6378137.0                 # semi-major axis (equatorial radius) [m]
    f = 1.0 / 298.257223563       # flattening
    b = a * (1.0 - f)             # semi-minor axis (polar radius) [m] ≈ 6356752.3142 m
    u = np.linspace(0, 2 * np.pi, 55)
    v = np.linspace(0, np.pi, 35)
    x_e = a * np.outer(np.cos(u), np.sin(v))
    y_e = a * np.outer(np.sin(u), np.sin(v))
    z_e = b * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(
        x_e, y_e, z_e,
        color="#707070", alpha=0.12, linewidth=0.65, zorder=1, label="Earth (WGS84 ellipsoid)"
    )

    ax.set_xlabel("X ECEF (m)", fontsize=10, labelpad=8)
    ax.set_ylabel("Y ECEF (m)", fontsize=10, labelpad=8)
    ax.set_zlabel("Z ECEF (m)", fontsize=10, labelpad=8)
    ax.set_title("3D Trajectory in Earth-Centered Earth-Fixed (ECEF) Coordinates", fontsize=13, pad=12)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.92)
    ax.grid(True, alpha=0.25)

    # Equal meter scales on all axes (cubic limits + cubic box) so the WGS84
    # ellipsoid is not stretched into a "football" by the 3D axes frame.
    # Flattening is preserved; only the plot box is prevented from squeezing.
    r_lim = a * 1.05
    ax.set_xlim(-r_lim, r_lim)
    ax.set_ylim(-r_lim, r_lim)
    ax.set_zlim(-r_lim, r_lim)
    try:
        ax.set_box_aspect((1.0, 1.0, 1.0))
    except AttributeError:
        pass

    ax.view_init(elev=22, azim=52)

    save_figure(fig, "07_trajectory_3d_ecef.png")

    # 8. Density, pressure, and temperature vs Time
    fig, axs = plt.subplots(
        3, 1, figsize=(args.figsize[0], 9.5), sharex=True, constrained_layout=True
    )
    components = [("rho", "Density", "kg/m^3"), ("p", "Pressure", "Pa"), ("T", "Temperature", "C")]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for ax, (col, label, unit), color in zip(axs, components, colors):
        ax.plot(t, df[col], color=color, linewidth=1.8)
        ax.set_ylabel(f"{label} ({unit})")
        ax.grid(True, alpha=0.35)
        if df[col].abs().max() > 1e5:
            ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
    axs[-1].set_xlabel("Time (s)")
    fig.suptitle("Atmospheric Density, Pressure, and Temperature vs Time")
    save_figure(fig, "08_density_pressure_temperature_vs_time.png")

    # 9. Dynamic pressure vs Time
    fig, ax = plt.subplots(figsize=args.figsize, constrained_layout=True)
    ax.plot(t, df["q_infinity"], color="#2ca02c", linewidth=2.0)
    idx_max_q = np.argmax(df["q_infinity"])
    max_q = int(df["q_infinity"][idx_max_q])
    t_max_q = df["t"][idx_max_q]
    ax.text(np.mean(t), np.max(df["q_infinity"]), f"Max Q: {max_q:,} Pa @ {t_max_q:,.1f} s")
    ax.fill_between(t, df["q_infinity"], alpha=0.18, color="#2ca02c")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Dynamic pressure (Pa)")
    ax.set_title("Dynamic Pressure vs Time")
    ax.grid(True, alpha=0.35)
    save_figure(fig, "09_dynamic_pressure_vs_time.png")

    print(f"\nAll 9 high-quality plots generated successfully!")
    print(f"   Output directory: {output_dir}")
    print(f"   DPI used: {args.dpi}")
    return 0


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()

    try:
        df = load_simulation_csv(input_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 1

    if args.mode == "live":
        return live_attitude_replay(df, target_fps=args.target_fps)

    return generate_png_plots(df, args)


if __name__ == "__main__":
    raise SystemExit(main())
