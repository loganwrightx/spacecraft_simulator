#!/usr/bin/env python3
"""
Local Plotly + Dash flight-data exploration dashboard for spacesim CSV logs.

Performance model (replay-first)
  - On load: precompute kinematics, ENU body axes, overview series, KPI pack
  - Heavy 2D charts rebuild only when the CSV, time window, or signals change
  - Playback / scrubbing only Patch()-updates playheads + the attitude sphere
    (static sphere/ENU basis are not rebuilt every frame)

Expects columns:
  t,x,y,z,vx,vy,vz,wx,wy,wz,qw,qx,qy,qz,lat,lon,alt[,T,p,rho,q_infinity]

Usage
  python flight_dashboard.py --input sim.csv --port 8050
  open http://127.0.0.1:8050
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.transform import Rotation
from dash import (
    Dash,
    Input,
    Output,
    Patch,
    State,
    callback,
    ctx,
    dcc,
    html,
    no_update,
)


# ---------------------------------------------------------------------------
# Constants / theme  (high-contrast, readable on dark UI)
# ---------------------------------------------------------------------------

REQUIRED_COLS = [
    "t", "x", "y", "z", "vx", "vy", "vz", "wx", "wy", "wz",
    "qw", "qx", "qy", "qz", "lat", "lon", "alt",
]

# Body frame (spacesim): +x up, +y east, +z north at surface-aligned launch.
# Attitude sphere display frame: East-North-Up (ENU).
ENU_AXIS_COLORS = {
    "E": "#ff5c5c",
    "N": "#3dde7f",
    "U": "#4db8ff",
}
BODY_AXIS_COLORS = {
    "X": "#ff8f8f",
    "Y": "#7dffb0",
    "Z": "#8ed0ff",
}

BG = "#0b0f14"
PANEL = "#121820"
PANEL2 = "#1a2330"
TEXT = "#f2f6fb"
MUTED = "#b0becf"
GRID = "#3a4a5c"
ACCENT = "#6eb6ef"
ACCENT2 = "#ffc14d"
TRACE_COLORS = [
    "#6eb6ef", "#ffc14d", "#3dde7f", "#ff7b7b", "#d28bff",
    "#2ec4b6", "#ff9f43", "#8295ff", "#ff7ab6", "#f0c929",
]

MAX_PLOT_POINTS = 2000
# Wireframe density: high sample counts along each arc (smooth circles),
# fewer displayed meridians/parallels (keeps the sphere readable).
SPHERE_MERIDIANS = 16          # number of longitude lines drawn
SPHERE_PARALLELS = 10          # number of latitude lines drawn (excl. poles)
SPHERE_ARC_SAMPLES = 96        # points per arc (smooth, not faceted)
PLAY_INTERVAL_MS = 50  # wall clock between play ticks

# Attitude figure trace order (dynamic traces first for cheap Patch updates)
# 0: body +X, 1: body +Y, 2: body +Z, 3: rocket mesh, 4: nose marker
# then static sphere wires + ENU axes + origin
ATT_BODY_X, ATT_BODY_Y, ATT_BODY_Z, ATT_MESH, ATT_NOSE = 0, 1, 2, 3, 4
ATT_DYNAMIC_COUNT = 5

# Body-frame rocket extents (before uniform scale). Keep overall on-screen size.
ROCKET_SCALE = 0.55
ROCKET_X_BASE = -0.55
ROCKET_X_TIP = 0.80
ROCKET_RADIUS = 0.072


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class FlightData:
    """Server-side flight log with precomputed kinematics for fast replay."""

    df: pd.DataFrame
    path: str
    t: np.ndarray
    v_east: np.ndarray
    v_north: np.ndarray
    v_up: np.ndarray
    speed_ecef: np.ndarray
    ground_speed: np.ndarray
    body_axes_enu: np.ndarray  # (N, 3, 3) columns = body x,y,z in ENU
    overview_idx: np.ndarray
    # Lightweight parallel arrays for KPI strip (avoid df.iloc every tick)
    kpi_alt: np.ndarray
    kpi_lat: np.ndarray
    kpi_lon: np.ndarray
    kpi_omega: np.ndarray
    kpi_qinf: np.ndarray | None = None
    kpi_rho: np.ndarray | None = None
    # Cached rocket mesh in body frame
    rocket_v: np.ndarray = field(repr=False, default=None)  # type: ignore[assignment]
    rocket_i: np.ndarray = field(repr=False, default=None)  # type: ignore[assignment]
    rocket_j: np.ndarray = field(repr=False, default=None)  # type: ignore[assignment]
    rocket_k: np.ndarray = field(repr=False, default=None)  # type: ignore[assignment]

    @property
    def n(self) -> int:
        return len(self.t)

    @property
    def t0(self) -> float:
        return float(self.t[0])

    @property
    def tf(self) -> float:
        return float(self.t[-1])

    @property
    def dt_med(self) -> float:
        if self.n < 2:
            return 0.01
        return float(np.median(np.diff(self.t)))


_DATA: FlightData | None = None


def ecef_to_enu_matrices(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    """
    Batch ECEF->ENU rotation matrices from WGS84 *geodetic* lat/lon (degrees).

    Rows of each 3x3 are East, North, Up unit vectors in ECEF:
        v_enu = R @ v_ecef

    Up is the geodetic ellipsoid normal (not geocentric radius).
    """
    lat = np.deg2rad(np.asarray(lat_deg, dtype=float))
    lon = np.deg2rad(np.asarray(lon_deg, dtype=float))
    s_lat, c_lat = np.sin(lat), np.cos(lat)
    s_lon, c_lon = np.sin(lon), np.cos(lon)
    n = lat.shape[0]
    R = np.empty((n, 3, 3), dtype=float)
    # East
    R[:, 0, 0] = -s_lon
    R[:, 0, 1] = c_lon
    R[:, 0, 2] = 0.0
    # North
    R[:, 1, 0] = -s_lat * c_lon
    R[:, 1, 1] = -s_lat * s_lon
    R[:, 1, 2] = c_lat
    # Up (geodetic)
    R[:, 2, 0] = c_lat * c_lon
    R[:, 2, 1] = c_lat * s_lon
    R[:, 2, 2] = s_lat
    return R


def _rocket_mesh_body_frame() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Tall finless rocket in body coordinates (nose along +x).

    Slender cylinder + hemispherical nose. Body-x span matches the previous
    mesh so ROCKET_SCALE keeps the on-screen size unchanged.
    """
    n_theta = 16
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    r = ROCKET_RADIUS
    x_base = ROCKET_X_BASE
    x_tip = ROCKET_X_TIP
    x_cyl_top = x_tip - r  # hemisphere sits on cylinder

    verts: list[list[float]] = []
    aft_idx = 0
    verts.append([x_base, 0.0, 0.0])

    n_cyl = 7
    cyl_xs = np.linspace(x_base, x_cyl_top, n_cyl)
    ring_start: list[int] = []
    for xs in cyl_xs:
        ring_start.append(len(verts))
        for th in theta:
            verts.append([float(xs), r * np.cos(th), r * np.sin(th)])

    n_nose = 6
    nose_ring_start: list[int] = []
    for elev in np.linspace(0.0, 0.5 * np.pi, n_nose + 1)[1:]:
        # elev=0 at equator (cylinder join), elev=pi/2 at tip
        xs = x_cyl_top + r * np.sin(elev)
        rho = r * np.cos(elev)
        nose_ring_start.append(len(verts))
        for th in theta:
            verts.append([float(xs), rho * np.cos(th), rho * np.sin(th)])

    tip_idx = len(verts)
    verts.append([x_tip, 0.0, 0.0])
    verts_arr = np.asarray(verts, dtype=float)

    i_list: list[int] = []
    j_list: list[int] = []
    k_list: list[int] = []

    def tri(a: int, b: int, c: int) -> None:
        i_list.append(a)
        j_list.append(b)
        k_list.append(c)

    base_ring = ring_start[0]
    for k in range(n_theta):
        k2 = (k + 1) % n_theta
        tri(aft_idx, base_ring + k2, base_ring + k)

    for s in range(len(ring_start) - 1):
        a0, b0 = ring_start[s], ring_start[s + 1]
        for k in range(n_theta):
            k2 = (k + 1) % n_theta
            tri(a0 + k, a0 + k2, b0 + k2)
            tri(a0 + k, b0 + k2, b0 + k)

    # Equator (last cyl ring) through hemispherical nose rings
    all_nose = ring_start[-1:] + nose_ring_start
    for s in range(len(all_nose) - 1):
        a0, b0 = all_nose[s], all_nose[s + 1]
        for k in range(n_theta):
            k2 = (k + 1) % n_theta
            tri(a0 + k, a0 + k2, b0 + k2)
            tri(a0 + k, b0 + k2, b0 + k)

    last = all_nose[-1]
    for k in range(n_theta):
        k2 = (k + 1) % n_theta
        tri(last + k, last + k2, tip_idx)

    return verts_arr, np.array(i_list), np.array(j_list), np.array(k_list)


_ROCKET_V, _ROCKET_I, _ROCKET_J, _ROCKET_K = _rocket_mesh_body_frame()
_ROCKET_V_SCALED = _ROCKET_V * ROCKET_SCALE


def rocket_vertices_enu(R_body_to_enu: np.ndarray) -> np.ndarray:
    return (R_body_to_enu @ _ROCKET_V_SCALED.T).T


def rocket_nose_enu(R_body_to_enu: np.ndarray) -> np.ndarray:
    """Nose tip in ENU (body +x)."""
    return R_body_to_enu @ np.array([ROCKET_X_TIP * ROCKET_SCALE, 0.0, 0.0])


def load_flight_csv(path: str | Path) -> FlightData:
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {path}")

    df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing}. Have: {list(df.columns)}")

    t = df["t"].to_numpy(dtype=float)
    lat = df["lat"].to_numpy(dtype=float)
    lon = df["lon"].to_numpy(dtype=float)
    v_ecef = df[["vx", "vy", "vz"]].to_numpy(dtype=float)

    R_e2n = ecef_to_enu_matrices(lat, lon)
    v_enu = np.einsum("nij,nj->ni", R_e2n, v_ecef)
    v_east, v_north, v_up = v_enu[:, 0], v_enu[:, 1], v_enu[:, 2]
    speed_ecef = np.linalg.norm(v_ecef, axis=1)
    ground_speed = np.hypot(v_east, v_north)

    rot_b2e = Rotation.from_quat(
        np.column_stack(
            [
                df["qx"].to_numpy(dtype=float),
                df["qy"].to_numpy(dtype=float),
                df["qz"].to_numpy(dtype=float),
                df["qw"].to_numpy(dtype=float),
            ]
        )
    )
    R_b2e = rot_b2e.as_matrix()
    body_axes_enu = np.einsum("nij,njk->nik", R_e2n, R_b2e)

    out = df.copy()
    out["speed_ecef"] = speed_ecef
    out["ground_speed"] = ground_speed
    out["v_east"] = v_east
    out["v_north"] = v_north
    out["v_up"] = v_up
    out["omega_norm"] = np.linalg.norm(out[["wx", "wy", "wz"]].to_numpy(dtype=float), axis=1)
    out["q_norm"] = np.linalg.norm(out[["qw", "qx", "qy", "qz"]].to_numpy(dtype=float), axis=1)

    n = len(out)
    if n <= MAX_PLOT_POINTS:
        overview_idx = np.arange(n, dtype=int)
    else:
        overview_idx = np.unique(np.linspace(0, n - 1, MAX_PLOT_POINTS, dtype=int))

    return FlightData(
        df=out,
        path=str(path),
        t=t,
        v_east=v_east,
        v_north=v_north,
        v_up=v_up,
        speed_ecef=speed_ecef,
        ground_speed=ground_speed,
        body_axes_enu=body_axes_enu,
        overview_idx=overview_idx,
        kpi_alt=out["alt"].to_numpy(dtype=float),
        kpi_lat=lat,
        kpi_lon=lon,
        kpi_omega=out["omega_norm"].to_numpy(dtype=float),
        kpi_qinf=out["q_infinity"].to_numpy(dtype=float) if "q_infinity" in out.columns else None,
        kpi_rho=out["rho"].to_numpy(dtype=float) if "rho" in out.columns else None,
        rocket_v=_ROCKET_V_SCALED,
        rocket_i=_ROCKET_I,
        rocket_j=_ROCKET_J,
        rocket_k=_ROCKET_K,
    )


def index_at_time(data: FlightData, t_query: float) -> int:
    i = int(np.searchsorted(data.t, t_query, side="left"))
    if i <= 0:
        return 0
    if i >= data.n:
        return data.n - 1
    # nearest
    if abs(data.t[i] - t_query) < abs(data.t[i - 1] - t_query):
        return i
    return i - 1


def overview_slice(
    data: FlightData, t_lo: float | None = None, t_hi: float | None = None
) -> np.ndarray:
    idx = data.overview_idx
    if t_lo is None and t_hi is None:
        return idx
    t = data.t
    t_lo = data.t0 if t_lo is None else t_lo
    t_hi = data.tf if t_hi is None else t_hi
    mask = (t[idx] >= t_lo) & (t[idx] <= t_hi)
    sel = idx[mask]
    if sel.size < 2:
        i0 = index_at_time(data, t_lo)
        i1 = index_at_time(data, t_hi)
        return np.unique(np.array([i0, i1], dtype=int))
    return sel


# ---------------------------------------------------------------------------
# Figure builders (static content)
# ---------------------------------------------------------------------------

def empty_fig(msg: str = "Load a simulation CSV") -> go.Figure:
    fig = go.Figure()
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor=PANEL,
        plot_bgcolor=PANEL,
        font=dict(color=TEXT, size=12),
        margin=dict(l=40, r=20, t=40, b=40),
        annotations=[
            dict(text=msg, x=0.5, y=0.5, showarrow=False, font=dict(color=MUTED, size=14))
        ],
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        uirevision="empty",
    )
    return fig


def base_layout(title_text: str | None = None, *, has_legend: bool = True, **kwargs: Any) -> dict[str, Any]:
    """
    Shared chart chrome.

    Title sits in the top paper margin (above the plot). Legends sit *below*
    the plot so they never stack under the title.
    """
    # Top margin reserved for title only; bottom margin holds x-label + legend
    top_m = 62 if title_text else 36
    bot_m = 96 if has_legend else 56
    layout: dict[str, Any] = dict(
        template="plotly_dark",
        paper_bgcolor=PANEL,
        plot_bgcolor=PANEL2,
        font=dict(color=TEXT, family="IBM Plex Sans, Segoe UI, sans-serif", size=12),
        margin=dict(l=64, r=28, t=top_m, b=bot_m),
        hovermode="x unified",
    )
    if title_text is not None:
        layout["title"] = dict(
            text=title_text,
            font=dict(size=14, color=TEXT),
            x=0.0,
            xanchor="left",
            # Paper coords: top of figure, clearly above the plotting domain
            y=0.995,
            yanchor="top",
            pad=dict(t=4, b=16, l=4),
        )
    if has_legend:
        layout["legend"] = dict(
            orientation="h",
            # Below the plot so titles never sit on top of legend entries
            y=-0.20,
            yanchor="top",
            x=0.0,
            xanchor="left",
            bgcolor="rgba(18,24,32,0.92)",
            bordercolor=GRID,
            borderwidth=1,
            font=dict(color=TEXT, size=11),
            tracegroupgap=10,
            itemsizing="constant",
            itemwidth=30,
        )
    else:
        layout["showlegend"] = False

    # Caller kwargs can override height/shapes/uirevision/margins; keep title/legend
    # unless explicitly replaced.
    saved_title = layout.get("title")
    saved_legend = layout.get("legend")
    layout.update(kwargs)
    if saved_title is not None and "title" not in kwargs:
        layout["title"] = saved_title
    if has_legend and saved_legend is not None and "legend" not in kwargs:
        layout["legend"] = saved_legend
    return layout


def style_xy(fig: go.Figure, xlabel: str = "Time (s)", ylabel: str = "") -> None:
    fig.update_xaxes(
        title=dict(text=xlabel, font=dict(color=MUTED)),
        gridcolor=GRID,
        zerolinecolor=GRID,
        showline=True,
        linecolor=GRID,
        tickfont=dict(color=MUTED),
        color=MUTED,
    )
    fig.update_yaxes(
        title=dict(text=ylabel, font=dict(color=MUTED)),
        gridcolor=GRID,
        zerolinecolor=GRID,
        showline=True,
        linecolor=GRID,
        tickfont=dict(color=MUTED),
        color=MUTED,
    )


def _playhead_shape(t_now: float, y0: float = 0.0, y1: float = 1.0, yref: str = "y") -> dict:
    return dict(
        type="line",
        x0=t_now,
        x1=t_now,
        y0=y0,
        y1=y1,
        xref="x",
        yref=yref,
        yrefpath=None,
        line=dict(color=ACCENT2, width=2, dash="dot"),
        layer="above",
    )


def _playhead_shape_paper(t_now: float) -> dict:
    """Vertical playhead spanning full subplot height via paper coords."""
    return dict(
        type="line",
        x0=t_now,
        x1=t_now,
        y0=0,
        y1=1,
        xref="x",
        yref="paper",
        line=dict(color=ACCENT2, width=2, dash="dot"),
        layer="above",
    )


def _unit_sphere_wireframe() -> list[go.Scatter3d]:
    """
    Smooth unit-sphere wireframe in ENU display space.

    Uses many samples along each arc (rounded meridians/parallels) while
    drawing only a modest number of grid lines at low opacity.
    """
    traces: list[go.Scatter3d] = []
    n_arc = SPHERE_ARC_SAMPLES
    wire_mer = "rgba(160, 185, 210, 0.20)"
    wire_par = "rgba(160, 185, 210, 0.14)"

    # Meridians (constant longitude): dense polar angle samples
    polar = np.linspace(0.0, np.pi, n_arc)
    for i in range(SPHERE_MERIDIANS):
        lon = 2.0 * np.pi * i / SPHERE_MERIDIANS
        x = np.sin(polar) * np.cos(lon)
        y = np.sin(polar) * np.sin(lon)
        z = np.cos(polar)
        traces.append(
            go.Scatter3d(
                x=x, y=y, z=z, mode="lines",
                line=dict(color=wire_mer, width=1.5),
                hoverinfo="skip", showlegend=False,
            )
        )

    # Parallels (constant latitude): dense azimuth samples, closed loop
    azim = np.linspace(0.0, 2.0 * np.pi, n_arc)
    # Evenly spaced between poles (exclude exact poles where radius is 0)
    for j in range(1, SPHERE_PARALLELS + 1):
        lat = np.pi * j / (SPHERE_PARALLELS + 1)  # (0, pi)
        r_xy = np.sin(lat)
        z = np.cos(lat)
        x = r_xy * np.cos(azim)
        y = r_xy * np.sin(azim)
        traces.append(
            go.Scatter3d(
                x=x, y=y, z=np.full_like(azim, z), mode="lines",
                line=dict(color=wire_par, width=1.5),
                hoverinfo="skip", showlegend=False,
            )
        )
    return traces


def _axis_line(
    direction: np.ndarray,
    length: float,
    color: str,
    name: str,
    width: int = 7,
    dash: str | None = None,
    label: bool = True,
) -> go.Scatter3d:
    tip = direction * length
    line: dict[str, Any] = dict(color=color, width=width)
    if dash:
        line["dash"] = dash
    # Place label slightly past the tip so it doesn't sit on the shaft
    label_pt = direction * (length * 1.12)
    if label:
        return go.Scatter3d(
            x=[0.0, float(tip[0]), float(label_pt[0])],
            y=[0.0, float(tip[1]), float(label_pt[1])],
            z=[0.0, float(tip[2]), float(label_pt[2])],
            mode="lines+markers+text",
            line=line,
            marker=dict(
                size=[0, 5, 0],
                color=[color, color, color],
                symbol=["circle", "circle", "circle"],
            ),
            text=["", "", name],
            textposition="middle center",
            textfont=dict(color=color, size=12, family="IBM Plex Sans, sans-serif"),
            name=name,
            hovertemplate=f"{name}<extra></extra>",
            showlegend=False,
        )
    return go.Scatter3d(
        x=[0.0, float(tip[0])],
        y=[0.0, float(tip[1])],
        z=[0.0, float(tip[2])],
        mode="lines",
        line=line,
        name=name,
        hovertemplate=f"{name}<extra></extra>",
        showlegend=False,
    )


def build_attitude_sphere(data: FlightData, idx: int) -> go.Figure:
    """
    ENU attitude sphere. Display frame is local East-North-Up at the vehicle
    geodetic LLA (WGS84). Dynamic traces (0-4) are patched during playback;
    sphere shell + ENU basis are static.
    """
    R = data.body_axes_enu[idx]
    fig = go.Figure()

    # --- Dynamic (indices 0-4) ---
    # Body triad (thinner solid) - body +X is vehicle nose / stack axis
    for ax_i, (name, color) in enumerate(
        (("+X", BODY_AXIS_COLORS["X"]), ("+Y", BODY_AXIS_COLORS["Y"]), ("+Z", BODY_AXIS_COLORS["Z"]))
    ):
        d = R[:, ax_i]
        fig.add_trace(_axis_line(d, 0.92, color, f"b{name}", width=5, label=True))

    verts = rocket_vertices_enu(R)
    fig.add_trace(
        go.Mesh3d(
            x=verts[:, 0],
            y=verts[:, 1],
            z=verts[:, 2],
            i=data.rocket_i,
            j=data.rocket_j,
            k=data.rocket_k,
            color="#e8edf4",
            opacity=1.0,
            flatshading=True,
            name="Vehicle",
            hovertemplate="Vehicle (body +X = nose)<extra></extra>",
            showscale=False,
            lighting=dict(ambient=0.68, diffuse=0.9, specular=0.45, roughness=0.4),
            lightposition=dict(x=1.2, y=2.0, z=1.8),
        )
    )
    nose = rocket_nose_enu(R)
    fig.add_trace(
        go.Scatter3d(
            x=[float(nose[0])],
            y=[float(nose[1])],
            z=[float(nose[2])],
            mode="markers",
            marker=dict(size=5, color=ACCENT2, symbol="circle"),
            name="Nose",
            hovertemplate="Nose (body +X)<extra></extra>",
            showlegend=False,
        )
    )

    # --- Static: sphere + ENU basis (local geodetic frame) ---
    for tr in _unit_sphere_wireframe():
        fig.add_trace(tr)

    # ENU triad - primary reference (slightly longer, dashed, labeled)
    for name, d, length in (
        ("E", np.array([1.0, 0.0, 0.0]), 1.12),
        ("N", np.array([0.0, 1.0, 0.0]), 1.12),
        ("U", np.array([0.0, 0.0, 1.0]), 1.12),
    ):
        fig.add_trace(
            _axis_line(d, length, ENU_AXIS_COLORS[name], name, width=7, dash="dash", label=True)
        )

    # Origin marker
    fig.add_trace(
        go.Scatter3d(
            x=[0.0], y=[0.0], z=[0.0],
            mode="markers",
            marker=dict(size=4, color="#ffffff", symbol="circle"),
            hovertemplate="Origin (vehicle)<extra></extra>",
            showlegend=False,
        )
    )

    fig.update_layout(
        **base_layout(
            title_text=(
                f"Local ENU attitude - t={data.t[idx]:.2f}s - "
                f"alt={data.kpi_alt[idx]:.1f} m - "
                f"|v|_ECEF={data.speed_ecef[idx]:.2f} m/s"
            ),
            has_legend=False,
            margin=dict(l=12, r=12, t=64, b=16),
            scene=dict(
                xaxis=dict(
                    title=dict(text="East", font=dict(color=ENU_AXIS_COLORS["E"], size=12)),
                    range=[-1.3, 1.3],
                    backgroundcolor=PANEL2,
                    gridcolor="rgba(58,74,92,0.55)",
                    showbackground=True,
                    zerolinecolor=GRID,
                    color=MUTED,
                    showspikes=False,
                    tickfont=dict(size=10, color=MUTED),
                ),
                yaxis=dict(
                    title=dict(text="North", font=dict(color=ENU_AXIS_COLORS["N"], size=12)),
                    range=[-1.3, 1.3],
                    backgroundcolor=PANEL2,
                    gridcolor="rgba(58,74,92,0.55)",
                    showbackground=True,
                    zerolinecolor=GRID,
                    color=MUTED,
                    showspikes=False,
                    tickfont=dict(size=10, color=MUTED),
                ),
                zaxis=dict(
                    title=dict(text="Up", font=dict(color=ENU_AXIS_COLORS["U"], size=12)),
                    range=[-1.3, 1.3],
                    backgroundcolor=PANEL2,
                    gridcolor="rgba(58,74,92,0.55)",
                    showbackground=True,
                    zerolinecolor=GRID,
                    color=MUTED,
                    showspikes=False,
                    tickfont=dict(size=10, color=MUTED),
                ),
                aspectmode="cube",
                bgcolor=PANEL,
                # Default view: Up = screen-up, East = left, North = toward camera
                # (slight E/U bias so the 3D structure is readable)
                camera=dict(
                    eye=dict(x=0.22, y=1.72, z=0.28),
                    up=dict(x=0.0, y=0.0, z=1.0),
                    center=dict(x=0.0, y=0.0, z=0.0),
                    projection=dict(type="perspective"),
                ),
                uirevision="attitude-camera",
            ),
            uirevision="attitude-sphere",
        )
    )
    return fig


def patch_attitude(data: FlightData, idx: int) -> Patch:
    """Minimal payload: only body triad + rocket mesh + nose + title."""
    R = data.body_axes_enu[idx]
    verts = rocket_vertices_enu(R)
    nose = rocket_nose_enu(R)
    p = Patch()
    for ax_i, tr_i in enumerate((ATT_BODY_X, ATT_BODY_Y, ATT_BODY_Z)):
        tip = R[:, ax_i] * 0.92
        label_pt = R[:, ax_i] * (0.92 * 1.12)
        p["data"][tr_i]["x"] = [0.0, float(tip[0]), float(label_pt[0])]
        p["data"][tr_i]["y"] = [0.0, float(tip[1]), float(label_pt[1])]
        p["data"][tr_i]["z"] = [0.0, float(tip[2]), float(label_pt[2])]
    p["data"][ATT_MESH]["x"] = verts[:, 0].tolist()
    p["data"][ATT_MESH]["y"] = verts[:, 1].tolist()
    p["data"][ATT_MESH]["z"] = verts[:, 2].tolist()
    p["data"][ATT_NOSE]["x"] = [float(nose[0])]
    p["data"][ATT_NOSE]["y"] = [float(nose[1])]
    p["data"][ATT_NOSE]["z"] = [float(nose[2])]
    p["layout"]["title"]["text"] = (
        f"Local ENU attitude - t={data.t[idx]:.2f}s - "
        f"alt={data.kpi_alt[idx]:.1f} m - "
        f"|v|_ECEF={data.speed_ecef[idx]:.2f} m/s"
    )
    return p


def patch_playhead_shapes(t_now: float, n_shapes: int = 1) -> Patch:
    p = Patch()
    for i in range(n_shapes):
        p["layout"]["shapes"][i]["x0"] = t_now
        p["layout"]["shapes"][i]["x1"] = t_now
    return p


def patch_marker_xy(trace_idx: int, x: float, y: float) -> Patch:
    p = Patch()
    p["data"][trace_idx]["x"] = [x]
    p["data"][trace_idx]["y"] = [y]
    return p


def build_speed_figure(data: FlightData, idx: int, t_lo: float, t_hi: float) -> go.Figure:
    sel = overview_slice(data, t_lo, t_hi)
    t = data.t[sel]
    t_now = float(data.t[idx])
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.20,
        subplot_titles=(
            "Speed relative to Earth surface (ECEF frame)",
            "Local ENU velocity components",
        ),
        row_heights=[0.5, 0.5],
    )
    fig.add_trace(
        go.Scattergl(x=t, y=data.speed_ecef[sel], name="|v|_ECEF",
                     line=dict(color=ACCENT, width=2.2)),
        row=1, col=1,
    )
    fig.add_trace(
        go.Scattergl(x=t, y=data.ground_speed[sel], name="Ground speed (E-N)",
                     line=dict(color=ACCENT2, width=2.2, dash="dash")),
        row=1, col=1,
    )
    fig.add_trace(
        go.Scattergl(x=t, y=data.v_east[sel], name="v_east",
                     line=dict(color=ENU_AXIS_COLORS["E"], width=2)),
        row=2, col=1,
    )
    fig.add_trace(
        go.Scattergl(x=t, y=data.v_north[sel], name="v_north",
                     line=dict(color=ENU_AXIS_COLORS["N"], width=2)),
        row=2, col=1,
    )
    fig.add_trace(
        go.Scattergl(x=t, y=data.v_up[sel], name="v_up (climb)",
                     line=dict(color=ENU_AXIS_COLORS["U"], width=2)),
        row=2, col=1,
    )
    # Playhead shapes (patched during replay) - one per subplot row
    fig.update_layout(
        **base_layout(
            title_text="Surface-relative kinematics",
            has_legend=True,
            height=500,
            margin=dict(l=64, r=28, t=64, b=100),
            shapes=[
                dict(type="line", x0=t_now, x1=t_now, y0=0, y1=1,
                     xref="x", yref="y domain",
                     line=dict(color=ACCENT2, width=2, dash="dot"), layer="above"),
                dict(type="line", x0=t_now, x1=t_now, y0=0, y1=1,
                     xref="x2", yref="y2 domain",
                     line=dict(color=ACCENT2, width=2, dash="dot"), layer="above"),
            ],
            uirevision="speed-static",
        )
    )
    fig.update_xaxes(title_text="Time (s)", row=2, col=1, gridcolor=GRID, color=MUTED, tickfont=dict(color=MUTED))
    fig.update_xaxes(gridcolor=GRID, color=MUTED, tickfont=dict(color=MUTED), row=1, col=1)
    fig.update_yaxes(title_text="m/s", row=1, col=1, gridcolor=GRID, color=MUTED, tickfont=dict(color=MUTED))
    fig.update_yaxes(title_text="m/s", row=2, col=1, gridcolor=GRID, color=MUTED, tickfont=dict(color=MUTED))
    # Subplot titles sit in the plot band; nudge them up with spacing already set
    fig.update_annotations(font=dict(color=MUTED, size=11))
    return fig


def build_altitude_figure(data: FlightData, idx: int, t_lo: float, t_hi: float) -> go.Figure:
    sel = overview_slice(data, t_lo, t_hi)
    t = data.t[sel]
    alt = data.kpi_alt[sel]
    t_now = float(data.t[idx])
    y_now = float(data.kpi_alt[idx])
    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=t, y=alt, name="Altitude",
            fill="tozeroy",
            line=dict(color="#3dde7f", width=2.2),
            fillcolor="rgba(61, 222, 127, 0.22)",
        )
    )
    # marker trace index 1 - patched on scrub/play
    fig.add_trace(
        go.Scattergl(
            x=[t_now], y=[y_now], mode="markers",
            marker=dict(size=11, color=ACCENT2, line=dict(width=1, color="#ffffff")),
            name="Now", showlegend=False,
        )
    )
    fig.update_layout(
        **base_layout(
            title_text="Altitude (ellipsoid)",
            has_legend=False,
            height=300,
            shapes=[
                dict(type="line", x0=t_now, x1=t_now, y0=0, y1=1,
                     xref="x", yref="y domain",
                     line=dict(color=ACCENT2, width=2, dash="dot"), layer="above"),
            ],
            uirevision="alt-static",
        )
    )
    style_xy(fig, ylabel="Altitude (m)")
    return fig


def build_ground_track(data: FlightData, idx: int, t_lo: float, t_hi: float) -> go.Figure:
    sel = overview_slice(data, t_lo, t_hi)
    lon = data.df["lon"].to_numpy()[sel]
    lat = data.df["lat"].to_numpy()[sel]
    t_lon = float(data.kpi_lon[idx])
    t_lat = float(data.kpi_lat[idx])
    fig = go.Figure()
    fig.add_trace(
        go.Scattergl(
            x=lon, y=lat, mode="lines",
            line=dict(color=ACCENT, width=2.5),
            name="Track",
            hovertemplate="lon=%{x:.5f} deg<br>lat=%{y:.5f} deg<extra></extra>",
        )
    )
    # idx 1 = now marker
    fig.add_trace(
        go.Scattergl(
            x=[t_lon], y=[t_lat], mode="markers",
            marker=dict(size=12, color=ACCENT2, line=dict(width=1, color="#ffffff")),
            name="Now",
        )
    )
    fig.add_trace(
        go.Scattergl(
            x=[float(data.kpi_lon[0])], y=[float(data.kpi_lat[0])],
            mode="markers", marker=dict(size=9, color="#3dde7f"), name="Start",
        )
    )
    fig.add_trace(
        go.Scattergl(
            x=[float(data.kpi_lon[-1])], y=[float(data.kpi_lat[-1])],
            mode="markers", marker=dict(size=9, color="#ff7b7b"), name="End",
        )
    )

    min_span = 0.05
    lon_all = data.kpi_lon
    lat_all = data.kpi_lat
    lon_c = 0.5 * (float(np.min(lon_all)) + float(np.max(lon_all)))
    lat_c = 0.5 * (float(np.min(lat_all)) + float(np.max(lat_all)))
    lon_span = max(float(np.max(lon_all) - np.min(lon_all)), min_span)
    lat_span = max(float(np.max(lat_all) - np.min(lat_all)), min_span)
    span = max(lon_span, lat_span)
    lon_span = min(span, 360.0)
    lat_span = min(span, 180.0)

    def _clamp(c: float, half: float, lo: float, hi: float) -> tuple[float, float]:
        half = min(half, 0.5 * (hi - lo))
        a, b = c - half, c + half
        if a < lo:
            b += lo - a
            a = lo
        if b > hi:
            a -= b - hi
            b = hi
        return max(a, lo), min(b, hi)

    x0, x1 = _clamp(lon_c, 0.5 * lon_span, -180.0, 180.0)
    y0, y1 = _clamp(lat_c, 0.5 * lat_span, -90.0, 90.0)

    fig.update_layout(
        **base_layout(
            title_text="Ground track",
            has_legend=True,
            height=360,
            yaxis=dict(scaleanchor="x", scaleratio=1, title="Latitude (deg)"),
            xaxis=dict(title="Longitude (deg)"),
            uirevision="track-static",
        )
    )
    fig.update_xaxes(range=[x0, x1], gridcolor=GRID, color=MUTED, tickfont=dict(color=MUTED))
    fig.update_yaxes(range=[y0, y1], gridcolor=GRID, color=MUTED, tickfont=dict(color=MUTED))
    return fig


def build_signal_figure(
    data: FlightData, signals: list[str], idx: int, t_lo: float, t_hi: float
) -> go.Figure:
    if not signals:
        return empty_fig("Select one or more signals")
    sel = overview_slice(data, t_lo, t_hi)
    t = data.t[sel]
    t_now = float(data.t[idx])
    fig = go.Figure()
    for i, sig in enumerate(signals):
        if sig not in data.df.columns:
            continue
        fig.add_trace(
            go.Scattergl(
                x=t, y=data.df[sig].to_numpy()[sel], name=sig,
                line=dict(color=TRACE_COLORS[i % len(TRACE_COLORS)], width=2),
            )
        )
    fig.update_layout(
        **base_layout(
            title_text="Signal explorer",
            has_legend=True,
            height=380,
            shapes=[
                dict(type="line", x0=t_now, x1=t_now, y0=0, y1=1,
                     xref="x", yref="y domain",
                     line=dict(color=ACCENT2, width=2, dash="dot"), layer="above"),
            ],
            uirevision="signals-static",
        )
    )
    style_xy(fig, ylabel="Value")
    return fig


def build_rates_figure(data: FlightData, idx: int, t_lo: float, t_hi: float) -> go.Figure:
    sel = overview_slice(data, t_lo, t_hi)
    t = data.t[sel]
    t_now = float(data.t[idx])
    fig = go.Figure()
    for col, color, name in (
        ("wx", BODY_AXIS_COLORS["X"], "w_x"),
        ("wy", BODY_AXIS_COLORS["Y"], "w_y"),
        ("wz", BODY_AXIS_COLORS["Z"], "w_z"),
    ):
        fig.add_trace(
            go.Scattergl(
                x=t, y=data.df[col].to_numpy()[sel], name=name,
                line=dict(color=color, width=2),
            )
        )
    fig.update_layout(
        **base_layout(
            title_text="Body angular rates",
            has_legend=True,
            height=320,
            shapes=[
                dict(type="line", x0=t_now, x1=t_now, y0=0, y1=1,
                     xref="x", yref="y domain",
                     line=dict(color=ACCENT2, width=2, dash="dot"), layer="above"),
            ],
            uirevision="rates-static",
        )
    )
    style_xy(fig, ylabel="rad/s")
    return fig


def build_atmosphere_figure(data: FlightData, idx: int, t_lo: float, t_hi: float) -> go.Figure:
    avail = [c for c in ("rho", "p", "T", "q_infinity") if c in data.df.columns]
    if not avail:
        return empty_fig("No atmosphere columns in this log")
    sel = overview_slice(data, t_lo, t_hi)
    t = data.t[sel]
    t_now = float(data.t[idx])
    n = len(avail)
    fig = make_subplots(
        rows=n, cols=1, shared_xaxes=True, vertical_spacing=0.12,
        subplot_titles=tuple(avail),
    )
    shapes = []
    for i, col in enumerate(avail):
        fig.add_trace(
            go.Scattergl(
                x=t, y=data.df[col].to_numpy()[sel], name=col,
                line=dict(color=TRACE_COLORS[i % len(TRACE_COLORS)], width=2),
                showlegend=False,
            ),
            row=i + 1, col=1,
        )
        xref = "x" if i == 0 else f"x{i + 1}"
        yref = "y domain" if i == 0 else f"y{i + 1} domain"
        shapes.append(
            dict(
                type="line", x0=t_now, x1=t_now, y0=0, y1=1,
                xref=xref, yref=yref,
                line=dict(color=ACCENT2, width=2, dash="dot"), layer="above",
            )
        )
        fig.update_yaxes(gridcolor=GRID, color=MUTED, tickfont=dict(color=MUTED), row=i + 1, col=1)
    fig.update_xaxes(
        title_text="Time (s)", gridcolor=GRID, color=MUTED,
        tickfont=dict(color=MUTED), row=n, col=1,
    )
    fig.update_layout(
        **base_layout(
            title_text="Atmosphere & dynamic pressure",
            has_legend=False,
            height=130 * n + 120,
            shapes=shapes,
            uirevision="atm-static",
        )
    )
    fig.update_annotations(font=dict(color=MUTED, size=11))
    return fig


def kpi_children(data: FlightData, idx: int) -> list:
    cards = [
        ("Time", f"{data.t[idx]:.3f} s"),
        ("Altitude", f"{data.kpi_alt[idx]:.2f} m"),
        ("|v| ECEF", f"{data.speed_ecef[idx]:.2f} m/s"),
        ("Ground spd", f"{data.ground_speed[idx]:.2f} m/s"),
        ("Climb rate", f"{data.v_up[idx]:+.2f} m/s"),
        ("Lat", f"{data.kpi_lat[idx]:.5f} deg"),
        ("Lon", f"{data.kpi_lon[idx]:.5f} deg"),
        ("|omega|", f"{data.kpi_omega[idx]:.4f} rad/s"),
    ]
    if data.kpi_qinf is not None:
        cards.append(("q_inf", f"{data.kpi_qinf[idx]:.1f} Pa"))
    if data.kpi_rho is not None:
        cards.append(("rho", f"{data.kpi_rho[idx]:.4f} kg/m^3"))

    return [
        html.Div(
            [html.Div(label, className="kpi-label"), html.Div(value, className="kpi-value")],
            className="kpi-card",
        )
        for label, value in cards
    ]


def signal_options(data: FlightData) -> list[dict[str, str]]:
    preferred = [
        "alt", "speed_ecef", "ground_speed", "v_east", "v_north", "v_up",
        "vx", "vy", "vz", "wx", "wy", "wz", "omega_norm",
        "qw", "qx", "qy", "qz", "q_norm",
        "lat", "lon", "x", "y", "z",
        "T", "p", "rho", "q_infinity",
    ]
    cols = [c for c in preferred if c in data.df.columns]
    for c in data.df.columns:
        if c not in cols and c != "t":
            cols.append(c)
    return [{"label": c, "value": c} for c in cols]


def _window(data: FlightData, window: list | None) -> tuple[float, float]:
    if window and len(window) == 2:
        t_lo, t_hi = float(window[0]), float(window[1])
    else:
        t_lo, t_hi = data.t0, data.tf
    if t_hi <= t_lo:
        t_hi = t_lo + max(data.dt_med, 1e-3)
    return t_lo, t_hi


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

EXTERNAL_STYLES = [
    {
        "href": (
            "https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600"
            "&family=IBM+Plex+Mono:wght@400;500&display=swap"
        ),
        "rel": "stylesheet",
    }
]

app = Dash(__name__, title="spacesim Flight Explorer", external_stylesheets=EXTERNAL_STYLES)

app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            :root {
                --bg: #0b0f14;
                --panel: #121820;
                --panel2: #1a2330;
                --text: #f2f6fb;
                --muted: #b0becf;
                --accent: #6eb6ef;
                --accent2: #ffc14d;
                --border: #3a4a5c;
            }
            * { box-sizing: border-box; }
            html, body { margin: 0; padding: 0; }
            body {
                background: var(--bg);
                color: var(--text);
                font-family: "IBM Plex Sans", "Segoe UI", sans-serif;
                overflow-x: hidden;
            }
            .app-shell {
                max-width: 1600px;
                margin: 0 auto;
                padding: 16px 18px 48px;
                position: relative;
            }
            .topbar {
                display: flex;
                flex-wrap: wrap;
                align-items: flex-end;
                gap: 14px 20px;
                padding: 14px 16px;
                background: linear-gradient(180deg, #171e28, var(--panel));
                border: 1px solid var(--border);
                border-radius: 12px;
                margin-bottom: 12px;
                position: relative;
                z-index: 5;
            }
            .brand h1 {
                margin: 0;
                font-size: 1.25rem;
                font-weight: 600;
                letter-spacing: 0.02em;
                color: var(--text);
            }
            .brand .sub {
                color: var(--muted);
                font-size: 0.85rem;
                margin-top: 4px;
                word-break: break-all;
            }
            .control-block {
                display: flex;
                flex-direction: column;
                gap: 6px;
                min-width: 140px;
                position: relative;
                z-index: 6;
            }
            .control-block label {
                font-size: 0.72rem;
                text-transform: uppercase;
                letter-spacing: 0.08em;
                color: var(--muted);
            }
            .control-row {
                display: flex;
                flex-wrap: wrap;
                gap: 10px;
                align-items: center;
            }
            button.dash-btn, .dash-btn {
                background: #243244;
                color: var(--text);
                border: 1px solid var(--border);
                border-radius: 8px;
                padding: 8px 14px;
                cursor: pointer;
                font-weight: 500;
                font-family: inherit;
                line-height: 1.2;
            }
            button.dash-btn:hover { border-color: var(--accent); color: #fff; }
            button.dash-btn.primary {
                background: #1f4a6e;
                border-color: var(--accent);
            }
            .kpi-strip {
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(128px, 1fr));
                gap: 8px;
                margin-bottom: 12px;
                position: relative;
                z-index: 1;
            }
            .kpi-card {
                background: var(--panel);
                border: 1px solid var(--border);
                border-radius: 10px;
                padding: 10px 12px;
                min-height: 58px;
            }
            .kpi-label {
                font-size: 0.68rem;
                text-transform: uppercase;
                letter-spacing: 0.07em;
                color: var(--muted);
            }
            .kpi-value {
                font-family: "IBM Plex Mono", monospace;
                font-size: 1.05rem;
                margin-top: 4px;
                color: #ffffff;
                font-weight: 500;
            }
            .grid-main {
                display: grid;
                grid-template-columns: 1.05fr 1fr;
                gap: 12px;
                position: relative;
                z-index: 1;
            }
            @media (max-width: 1100px) {
                .grid-main { grid-template-columns: 1fr; }
            }
            .panel {
                background: var(--panel);
                border: 1px solid var(--border);
                border-radius: 12px;
                padding: 10px;
                overflow: hidden;
                position: relative;
                z-index: 1;
            }
            .panel + .panel-spacer { margin-top: 12px; }
            .panel-title {
                font-size: 0.75rem;
                text-transform: uppercase;
                letter-spacing: 0.08em;
                color: var(--muted);
                padding: 8px 8px 12px;
                margin: 0 0 4px;
                line-height: 1.35;
                border-bottom: 1px solid var(--border);
            }
            .scrubber-panel {
                margin: 0 0 12px;
                padding: 12px 16px 14px;
                background: var(--panel);
                border: 1px solid var(--border);
                border-radius: 12px;
                position: relative;
                z-index: 4;
            }
            .hint {
                color: var(--muted);
                font-size: 0.82rem;
                margin-top: 6px;
                line-height: 1.4;
            }
            /* Dropdown menus above plots */
            .Select-menu-outer, .VirtualizedSelectOption,
            .dash-dropdown .Select-menu-outer {
                z-index: 1000 !important;
            }
            .dash-dropdown { position: relative; z-index: 20; }
            /* Slider contrast */
            .rc-slider-track { background-color: var(--accent) !important; }
            .rc-slider-handle {
                border-color: var(--accent2) !important;
                background-color: #fff !important;
            }
            .rc-slider-rail { background-color: #2a3848 !important; }
            .rc-slider-dot-active { border-color: var(--accent) !important; }
            .rc-slider-mark-text { color: var(--muted) !important; }
            /* Graph containers: no overlap, clear stacking */
            .js-plotly-plot, .plot-container { width: 100% !important; }
            .dash-graph { isolation: isolate; }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""


def build_layout(data: FlightData | None) -> html.Div:
    if data is None:
        t0, tf, marks = 0.0, 1.0, {0: "0", 1: "1"}
        sig_opts: list[dict[str, str]] = []
        path_str = "sim.csv"
        n_pts = 0
        default_sigs: list[str] = []
        step = 0.01
    else:
        t0, tf = data.t0, data.tf
        marks = {
            t0: f"{t0:.0f}s",
            0.5 * (tf - t0) + t0: f"{0.5 * (tf - t0) + t0:.0f}s",
            tf: f"{tf:.0f}s",
        }
        sig_opts = signal_options(data)
        path_str = data.path
        n_pts = data.n
        default_sigs = [
            s for s in ("alt", "speed_ecef", "ground_speed", "q_infinity")
            if s in data.df.columns
        ][:3]
        step = max(data.dt_med, (tf - t0) / max(n_pts, 1))

    return html.Div(
        className="app-shell",
        children=[
            dcc.Store(id="play-state", data={"playing": False}),
            dcc.Interval(
                id="play-clock",
                interval=PLAY_INTERVAL_MS,
                disabled=True,
                n_intervals=0,
            ),
            html.Div(
                className="topbar",
                children=[
                    html.Div(
                        className="brand",
                        children=[
                            html.H1("spacesim Flight Explorer"),
                            html.Div(
                                className="sub",
                                id="dataset-label",
                                children=(
                                    f"{path_str} - {n_pts:,} samples - "
                                    "replay-optimized (patched playhead)"
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        className="control-block",
                        children=[
                            html.Label("CSV path"),
                            dcc.Input(
                                id="csv-path",
                                type="text",
                                value=path_str,
                                style={
                                    "width": "280px",
                                    "background": "#243244",
                                    "color": TEXT,
                                    "border": f"1px solid {GRID}",
                                    "borderRadius": "8px",
                                    "padding": "8px 10px",
                                },
                            ),
                        ],
                    ),
                    html.Button("Reload", id="btn-reload", className="dash-btn", n_clicks=0),
                    html.Div(
                        className="control-block",
                        children=[
                            html.Label("Playback speed"),
                            dcc.Dropdown(
                                id="play-speed",
                                options=[
                                    {"label": "0.25x speed", "value": 0.25},
                                    {"label": "0.5x speed", "value": 0.5},
                                    {"label": "1x realtime", "value": 1.0},
                                    {"label": "2x speed", "value": 2.0},
                                    {"label": "4x speed", "value": 4.0},
                                    {"label": "8x speed", "value": 8.0},
                                ],
                                value=1.0,
                                clearable=False,
                                style={"width": "150px"},
                            ),
                        ],
                    ),
                    html.Div(
                        className="control-row",
                        children=[
                            html.Button("Play", id="btn-play", className="dash-btn primary", n_clicks=0),
                            html.Button("Pause", id="btn-pause", className="dash-btn", n_clicks=0),
                            html.Button("Restart", id="btn-restart", className="dash-btn", n_clicks=0),
                        ],
                    ),
                ],
            ),
            html.Div(
                id="kpi-strip",
                className="kpi-strip",
                children=[] if data is None else kpi_children(data, 0),
            ),
            html.Div(
                className="scrubber-panel",
                children=[
                    html.Div(className="panel-title", children="Mission timeline / playhead"),
                    dcc.Slider(
                        id="time-slider",
                        min=t0,
                        max=tf if tf > t0 else t0 + 1.0,
                        step=step,
                        value=t0,
                        marks=marks,
                        tooltip={"placement": "bottom", "always_visible": False},
                        # drag live for responsive scrub; charts only patch
                        updatemode="drag",
                    ),
                    html.Div(
                        className="hint",
                        children=(
                            "Playback only updates the ENU sphere + playhead markers "
                            "(charts are pre-rendered for the selected window). "
                            "Change the analysis window or signals to rebuild overview plots."
                        ),
                    ),
                    html.Div(
                        className="panel-title",
                        children="Analysis time window (rebuilds overview plots)",
                        style={"marginTop": "12px"},
                    ),
                    dcc.RangeSlider(
                        id="window-slider",
                        min=t0,
                        max=tf if tf > t0 else t0 + 1.0,
                        step=max((tf - t0) / 500.0, 1e-4) if data else 0.01,
                        value=[t0, tf if tf > t0 else t0 + 1.0],
                        marks=marks,
                        tooltip={"placement": "bottom"},
                        allowCross=False,
                        updatemode="mouseup",
                    ),
                ],
            ),
            html.Div(
                className="grid-main",
                children=[
                    html.Div(
                        className="panel",
                        children=[
                            html.Div(
                                className="panel-title",
                                children="ENU attitude sphere (Starship-style)",
                            ),
                            dcc.Graph(
                                id="fig-attitude",
                                figure=empty_fig() if data is None else build_attitude_sphere(data, 0),
                                config={"displayModeBar": True, "scrollZoom": True},
                                style={"height": "480px"},
                            ),
                            html.Div(
                                className="hint",
                                style={"padding": "0 6px 4px"},
                                children=(
                                    "Dashed: fixed East / North / Up. "
                                    "Solid triad: body +X/+Y/+Z. Mesh: vehicle (nose = body +X)."
                                ),
                            ),
                        ],
                    ),
                    html.Div(
                        className="panel",
                        children=[
                            dcc.Graph(
                                id="fig-speed",
                                figure=(
                                    empty_fig()
                                    if data is None
                                    else build_speed_figure(data, 0, t0, tf)
                                ),
                                config={"displaylogo": False},
                            ),
                            dcc.Graph(
                                id="fig-altitude",
                                figure=(
                                    empty_fig()
                                    if data is None
                                    else build_altitude_figure(data, 0, t0, tf)
                                ),
                                config={"displaylogo": False},
                            ),
                        ],
                    ),
                ],
            ),
            html.Div(
                className="grid-main",
                style={"marginTop": "12px"},
                children=[
                    html.Div(
                        className="panel",
                        children=[
                            dcc.Graph(
                                id="fig-track",
                                figure=(
                                    empty_fig()
                                    if data is None
                                    else build_ground_track(data, 0, t0, tf)
                                ),
                                config={"displaylogo": False},
                            ),
                        ],
                    ),
                    html.Div(
                        className="panel",
                        children=[
                            dcc.Graph(
                                id="fig-rates",
                                figure=(
                                    empty_fig()
                                    if data is None
                                    else build_rates_figure(data, 0, t0, tf)
                                ),
                                config={"displaylogo": False},
                            ),
                        ],
                    ),
                ],
            ),
            html.Div(
                className="panel",
                style={"marginTop": "12px"},
                children=[
                    html.Div(
                        style={
                            "display": "flex",
                            "flexWrap": "wrap",
                            "gap": "12px",
                            "alignItems": "center",
                            "padding": "4px 6px 8px",
                            "position": "relative",
                            "zIndex": 30,
                        },
                        children=[
                            html.Div(
                                className="panel-title",
                                children="Signal explorer",
                                style={"padding": 0, "margin": 0},
                            ),
                            dcc.Dropdown(
                                id="signal-select",
                                options=sig_opts,
                                value=default_sigs,
                                multi=True,
                                placeholder="Choose telemetry channels...",
                                style={"minWidth": "420px", "flex": "1"},
                            ),
                        ],
                    ),
                    dcc.Graph(
                        id="fig-signals",
                        figure=(
                            empty_fig()
                            if data is None
                            else build_signal_figure(data, default_sigs, 0, t0, tf)
                        ),
                        config={"displaylogo": False},
                    ),
                ],
            ),
            html.Div(
                className="panel",
                style={"marginTop": "12px"},
                children=[
                    dcc.Graph(
                        id="fig-atmosphere",
                        figure=(
                            empty_fig()
                            if data is None
                            else build_atmosphere_figure(data, 0, t0, tf)
                        ),
                        config={"displaylogo": False},
                    ),
                ],
            ),
            html.Div(
                className="hint",
                style={"marginTop": "14px", "textAlign": "center"},
                children=(
                    "Surface-relative speed uses ECEF velocity (Earth-fixed): "
                    "|v|_ECEF = ||(vx,vy,vz)||. Ground speed / climb rate are that velocity "
                    "in local ENU. Not ECI / inertial."
                ),
            ),
        ],
    )


app.layout = build_layout(None)


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@callback(
    Output("dataset-label", "children"),
    Output("time-slider", "min"),
    Output("time-slider", "max"),
    Output("time-slider", "step"),
    Output("time-slider", "marks"),
    Output("time-slider", "value"),
    Output("window-slider", "min"),
    Output("window-slider", "max"),
    Output("window-slider", "step"),
    Output("window-slider", "marks"),
    Output("window-slider", "value"),
    Output("signal-select", "options"),
    Output("signal-select", "value"),
    Output("play-state", "data", allow_duplicate=True),
    Output("play-clock", "disabled", allow_duplicate=True),
    Input("btn-reload", "n_clicks"),
    State("csv-path", "value"),
    prevent_initial_call=True,
)
def reload_csv(n_clicks: int, path: str):
    global _DATA
    try:
        _DATA = load_flight_csv(path or "sim.csv")
    except Exception as exc:
        return (
            f"ERROR loading {path}: {exc}",
            *(no_update,) * 14,
        )

    data = _DATA
    t0, tf = data.t0, data.tf
    marks = {
        t0: f"{t0:.0f}s",
        0.5 * (tf - t0) + t0: f"{0.5 * (tf - t0) + t0:.0f}s",
        tf: f"{tf:.0f}s",
    }
    step = max(data.dt_med, (tf - t0) / max(data.n, 1))
    default_sigs = [
        s for s in ("alt", "speed_ecef", "ground_speed", "q_infinity")
        if s in data.df.columns
    ][:3]
    label = (
        f"{data.path} - {data.n:,} samples - dt~{data.dt_med * 1000:.2f} ms - "
        "replay-optimized"
    )
    return (
        label,
        t0, tf, step, marks, t0,
        t0, tf, max((tf - t0) / 500.0, 1e-4), marks, [t0, tf],
        signal_options(data),
        default_sigs,
        {"playing": False},
        True,
    )


@callback(
    Output("play-clock", "disabled"),
    Output("play-state", "data"),
    Output("time-slider", "value", allow_duplicate=True),
    Input("btn-play", "n_clicks"),
    Input("btn-pause", "n_clicks"),
    Input("btn-restart", "n_clicks"),
    State("play-state", "data"),
    prevent_initial_call=True,
)
def transport_controls(n_play, n_pause, n_restart, play_state):
    global _DATA
    state = dict(play_state or {"playing": False})
    trigger = ctx.triggered_id
    if trigger == "btn-play":
        state["playing"] = True
        return False, state, no_update
    if trigger == "btn-pause":
        state["playing"] = False
        return True, state, no_update
    # restart
    state["playing"] = True
    t0 = _DATA.t0 if _DATA is not None else 0.0
    return False, state, t0


@callback(
    Output("time-slider", "value", allow_duplicate=True),
    Output("play-state", "data", allow_duplicate=True),
    Output("play-clock", "disabled", allow_duplicate=True),
    Input("play-clock", "n_intervals"),
    State("play-state", "data"),
    State("play-speed", "value"),
    State("time-slider", "value"),
    prevent_initial_call=True,
)
def advance_playhead(n_intervals, play_state, speed, t_slider):
    """Advance mission time only - views patch from the slider callback."""
    global _DATA
    if _DATA is None or not play_state or not play_state.get("playing"):
        return no_update, no_update, no_update

    dt_wall = PLAY_INTERVAL_MS / 1000.0
    speed = float(speed or 1.0)
    t_now = float(t_slider if t_slider is not None else _DATA.t0)
    t_next = t_now + dt_wall * speed

    if t_next >= _DATA.tf:
        state = dict(play_state)
        state["playing"] = False
        return float(_DATA.tf), state, True

    return t_next, no_update, False


@callback(
    Output("fig-speed", "figure"),
    Output("fig-altitude", "figure"),
    Output("fig-track", "figure"),
    Output("fig-rates", "figure"),
    Output("fig-signals", "figure"),
    Output("fig-atmosphere", "figure"),
    Input("window-slider", "value"),
    Input("signal-select", "value"),
    Input("btn-reload", "n_clicks"),
    State("time-slider", "value"),
    prevent_initial_call=False,
)
def rebuild_static_charts(window, signals, _reload_clicks, t_slider):
    """
    Expensive path: decimated overview series.
    Only runs on window / signal / reload changes - not every play tick.
    """
    global _DATA
    if _DATA is None:
        e = empty_fig("Load a CSV (path -> Reload)")
        return e, e, e, e, e, e

    data = _DATA
    t_lo, t_hi = _window(data, window)
    t_now = float(t_slider if t_slider is not None else data.t0)
    idx = index_at_time(data, t_now)
    signals = list(signals or [])

    return (
        build_speed_figure(data, idx, t_lo, t_hi),
        build_altitude_figure(data, idx, t_lo, t_hi),
        build_ground_track(data, idx, t_lo, t_hi),
        build_rates_figure(data, idx, t_lo, t_hi),
        build_signal_figure(data, signals, idx, t_lo, t_hi),
        build_atmosphere_figure(data, idx, t_lo, t_hi),
    )


@callback(
    Output("kpi-strip", "children"),
    Output("fig-attitude", "figure"),
    Output("fig-speed", "figure", allow_duplicate=True),
    Output("fig-altitude", "figure", allow_duplicate=True),
    Output("fig-track", "figure", allow_duplicate=True),
    Output("fig-rates", "figure", allow_duplicate=True),
    Output("fig-signals", "figure", allow_duplicate=True),
    Output("fig-atmosphere", "figure", allow_duplicate=True),
    Input("time-slider", "value"),
    State("fig-attitude", "figure"),
    State("fig-speed", "figure"),
    State("fig-altitude", "figure"),
    State("fig-track", "figure"),
    State("fig-rates", "figure"),
    State("fig-signals", "figure"),
    State("fig-atmosphere", "figure"),
    prevent_initial_call=True,
)
def live_playhead(t_slider, fig_att, fig_spd, fig_alt, fig_trk, fig_rts, fig_sig, fig_atm):
    """
    Cheap path during scrub/play: Patch attitude dynamics + playhead markers.
    Never rebuilds full scatter series.
    """
    global _DATA
    if _DATA is None:
        return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

    data = _DATA
    t_now = float(t_slider if t_slider is not None else data.t0)
    idx = index_at_time(data, t_now)

    # Attitude: full build only if figure is empty / wrong structure; else Patch
    if (
        not fig_att
        or not fig_att.get("data")
        or len(fig_att["data"]) < ATT_DYNAMIC_COUNT
        or fig_att["data"][ATT_MESH].get("type") != "mesh3d"
    ):
        att_out: Any = build_attitude_sphere(data, idx)
    else:
        att_out = patch_attitude(data, idx)

    def _ph_time(fig: dict | None, n_shapes: int) -> Any:
        if not fig or not fig.get("layout") or not fig["layout"].get("shapes"):
            return no_update
        shapes = fig["layout"]["shapes"]
        n = min(n_shapes, len(shapes))
        if n <= 0:
            return no_update
        return patch_playhead_shapes(t_now, n_shapes=n)

    def _ph_alt(fig: dict | None) -> Any:
        if not fig or not fig.get("data") or len(fig["data"]) < 2:
            return no_update
        p = Patch()
        p["layout"]["shapes"][0]["x0"] = t_now
        p["layout"]["shapes"][0]["x1"] = t_now
        p["data"][1]["x"] = [t_now]
        p["data"][1]["y"] = [float(data.kpi_alt[idx])]
        return p

    def _ph_track(fig: dict | None) -> Any:
        if not fig or not fig.get("data") or len(fig["data"]) < 2:
            return no_update
        p = Patch()
        p["data"][1]["x"] = [float(data.kpi_lon[idx])]
        p["data"][1]["y"] = [float(data.kpi_lat[idx])]
        return p

    # Atmosphere shape count = number of subplot rows with shapes
    n_atm = 0
    if fig_atm and fig_atm.get("layout") and fig_atm["layout"].get("shapes"):
        n_atm = len(fig_atm["layout"]["shapes"])

    n_spd = 0
    if fig_spd and fig_spd.get("layout") and fig_spd["layout"].get("shapes"):
        n_spd = len(fig_spd["layout"]["shapes"])

    return (
        kpi_children(data, idx),
        att_out,
        _ph_time(fig_spd, n_spd) if n_spd else no_update,
        _ph_alt(fig_alt),
        _ph_track(fig_trk),
        _ph_time(fig_rts, 1),
        _ph_time(fig_sig, 1),
        _ph_time(fig_atm, n_atm) if n_atm else no_update,
    )


# Full attitude figure once after reload (live_playhead is prevent_initial_call)
@callback(
    Output("fig-attitude", "figure", allow_duplicate=True),
    Output("kpi-strip", "children", allow_duplicate=True),
    Input("btn-reload", "n_clicks"),
    prevent_initial_call=True,
)
def attitude_after_reload(_n):
    global _DATA
    if _DATA is None:
        return empty_fig("Load failed"), []
    return build_attitude_sphere(_DATA, 0), kpi_children(_DATA, 0)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="spacesim Plotly/Dash flight data dashboard")
    p.add_argument("--input", "-i", default="sim.csv", help="Path to simulation CSV")
    p.add_argument("--host", default="127.0.0.1", help="Bind address (default localhost)")
    p.add_argument("--port", type=int, default=8050, help="Port (default 8050)")
    p.add_argument("--debug", action="store_true", help="Enable Dash debug reloader")
    return p.parse_args()


def main() -> int:
    global _DATA
    args = parse_args()
    path = Path(args.input).expanduser()
    try:
        _DATA = load_flight_csv(path)
        print(f"Loaded {_DATA.n} samples from {_DATA.path}")
        print(f"  t in [{_DATA.t0:.3f}, {_DATA.tf:.3f}] s   dt~{_DATA.dt_med * 1000:.3f} ms")
        print(f"  alt in [{_DATA.kpi_alt.min():.2f}, {_DATA.kpi_alt.max():.2f}] m")
        print(
            f"  |v|_ECEF in [{_DATA.speed_ecef.min():.2f}, {_DATA.speed_ecef.max():.2f}] m/s "
            f"(surface-relative / Earth-fixed)"
        )
        print("  Replay mode: overview charts are static; playhead + attitude are patched.")
    except Exception as exc:
        print(f"WARNING: could not preload '{path}': {exc}")
        print("Dashboard will start empty - use path + Reload in the UI.")
        _DATA = None

    app.layout = build_layout(_DATA)
    print(f"\nFlight Explorer -> http://{args.host}:{args.port}\n")
    app.run(host=args.host, port=args.port, debug=args.debug)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
