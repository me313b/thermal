# =====================================================================
#  Immersion Pack Lab
#  Static / stirred immersion-cooled 21700 battery pack with internal
#  water-cooled tube heat exchanger (ICDC architecture).
#
#  Physics anchored to:
#   [1] Wang, Zhao, Wang & Huang (2023) "Heat transfer characteristics and
#       influencing factors of immersion coupled direct cooling for battery
#       thermal management", J. Energy Storage 62, 106821.
#   [2] batterydesign.net, "Mercedes AMG HPB80" (Nov 2023) - production
#       reference for a dielectric-cooled 21700 pack.
#   [3] Coolant property table: coolant_comparison_reviewed.xlsx (mb, 2026).
#  Correlations: Churchill-Chu (vertical plate, horizontal cylinder),
#  Churchill-Bernstein (crossflow), Hausen (laminar entry), Gnielinski
#  (turbulent tube), Schmidt (annular fin efficiency).
#
#  Run:  streamlit run app.py       Smoke test:  SMOKE=1 python app.py
# =====================================================================

import os, math, contextlib, json

APP_VERSION = "v10.23"
from pathlib import Path
_APPDIR = Path(__file__).resolve().parent
import json
import numpy as np
from correlations import water_nu as _water_nu
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.io as pio
import fea_export
import streamlit.components.v1 as components
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from livepack import live_pack_html
from cockpit import cockpit_html
from zonal import (KernelBank, solve_zonal, monte_carlo,
                   derived_layout)
import schematics as _S
import pack3d as _P3

pio.templates["packlab"] = go.layout.Template(layout=dict(
    font=dict(family="Inter, -apple-system, 'Segoe UI', Roboto, sans-serif",
              size=13, color="#334155"),
    title=dict(font=dict(size=15, color="#0F172A"), x=0.01,
               xanchor="left", y=0.985, yanchor="top"),
    legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="right",
                x=1.0, bgcolor="rgba(0,0,0,0)", font=dict(size=11.5)),
    paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
    colorway=["#6366F1", "#06B6D4", "#F59E0B", "#EF4444", "#10B981",
              "#8B5CF6", "#64748B"],
    xaxis=dict(gridcolor="#EEF2F7", zerolinecolor="#E2E8F0", linecolor="#E2E8F0"),
    yaxis=dict(gridcolor="#EEF2F7", zerolinecolor="#E2E8F0", linecolor="#E2E8F0"),
    margin=dict(l=20, r=20, t=64, b=20)))
pio.templates.default = "packlab" 

G = 9.81
T_REF = 25.0          # deg C reference for tabulated properties
KELVIN = 273.15

# ------------------------------------------------------------------ #
#  Fluid properties                                                   #
# ------------------------------------------------------------------ #
FALLBACK_CSV = """name,family,k,rho,cp,nu_cSt,beta,B_visc,bp_C,flash_C,dielectric,bdv_kV,notes,review
Transformer oil,Mineral hydrocarbon,0.13,875,1900,9.8,0.00075,3200,280,150,True,50,Baseline mineral oil,
MIVOLT DF7,Dielectric ester,0.13,900,2000,7.0,0.00075,3200,250,170,True,50,Low-viscosity EV ester,
Novec 7100 (HFE-7100),Hydrofluoroether,0.069,1510,1180,0.38,0.0015,1500,61,True,28,Fluorinated reference,
Deionized water,Water,0.6,997,4180,0.89,0.00026,1900,100,False,Thermal reference only,
"""

def _read_coolants() -> pd.DataFrame:
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "coolants.csv")
    try:
        df = pd.read_csv(path)
    except Exception:
        from io import StringIO
        df = pd.read_csv(StringIO(FALLBACK_CSV))
    df["dielectric"] = df["dielectric"].astype(bool)
    for c in ["k", "rho", "cp", "nu_cSt", "beta", "B_visc", "bp_C", "flash_C"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def fluid_dict(row) -> dict:
    return dict(name=row["name"], family=row["family"], k=float(row["k"]),
                rho=float(row["rho"]), cp=float(row["cp"]),
                nu25=float(row["nu_cSt"]) * 1e-6, beta=float(row["beta"]),
                B=float(row["B_visc"]),
                bp=row.get("bp_C", np.nan), flash=row.get("flash_C", np.nan),
                dielectric=bool(row["dielectric"]))

def nu_of_T(fl: dict, T_C: float) -> float:
    """Andrade-type viscosity-temperature law anchored at 25 °C.
    nu(T) = nu25 * exp(B*(1/T - 1/298.15)), T in K. B set per fluid family
    (mineral/ester ~3200 K fitted to transformer-oil data in [1] Table 3)."""
    T = max(T_C, -30.0) + KELVIN
    return fl["nu25"] * math.exp(fl["B"] * (1.0 / T - 1.0 / (T_REF + KELVIN)))

def film_props(fl: dict, T_film_C: float) -> dict:
    """Constant k, rho, cp from the table; nu evaluated at film temperature."""
    nu = nu_of_T(fl, T_film_C)
    alpha = fl["k"] / (fl["rho"] * fl["cp"])
    return dict(k=fl["k"], rho=fl["rho"], cp=fl["cp"], nu=nu,
                alpha=alpha, Pr=nu / alpha, beta=fl["beta"])

# Water-loop fluids (inside the tubes), properties near 20-25 °C
WATER_LOOP = {
    "Water": dict(rho=998.0, cp=4182.0, k=0.60, mu=1.0e-3),
    "Water-glycol 50/50": dict(rho=1070.0, cp=3300.0, k=0.37, mu=3.8e-3),
}

# ------------------------------------------------------------------ #
#  Convection correlations                                            #
# ------------------------------------------------------------------ #
def rayleigh(p: dict, dT: float, L: float) -> float:
    dT = max(abs(dT), 0.05)
    return G * p["beta"] * dT * L ** 3 / (p["nu"] * p["alpha"])

def nu_vertical_cc(Ra: float, Pr: float) -> float:
    """Churchill-Chu, vertical plate, all Ra (used for the cell wall,
    slightly conservative for a slender cylinder in oil)."""
    f = (1.0 + (0.492 / Pr) ** (9.0 / 16.0)) ** (8.0 / 27.0)
    return (0.825 + 0.387 * Ra ** (1.0 / 6.0) / f) ** 2

def nu_horiz_cyl_cc(Ra: float, Pr: float) -> float:
    """Churchill-Chu, horizontal cylinder (the HX tubes)."""
    f = (1.0 + (0.559 / Pr) ** (9.0 / 16.0)) ** (8.0 / 27.0)
    return (0.60 + 0.387 * Ra ** (1.0 / 6.0) / f) ** 2

def nu_crossflow_cb(Re: float, Pr: float) -> float:
    """Churchill-Bernstein, forced crossflow over a cylinder."""
    if Re < 1e-6:
        return 0.0
    a = 0.62 * Re ** 0.5 * Pr ** (1.0 / 3.0)
    b = (1.0 + (0.4 / Pr) ** (2.0 / 3.0)) ** 0.25
    c = (1.0 + (Re / 282000.0) ** (5.0 / 8.0)) ** (4.0 / 5.0)
    return 0.3 + a / b * c

def blend_mixed(h_nat: float, h_for: float) -> float:
    """Mixed convection blend Nu^3 = Nu_n^3 + Nu_f^3 (transverse flow)."""
    return (h_nat ** 3 + h_for ** 3) ** (1.0 / 3.0)

def gap_factor(gap_mm: float, expo: float = 0.6, floor: float = 0.35) -> float:
    """Confinement penalty on cell-side natural convection.
    Calibrated to [1] Fig. 9: full performance for gap >= 6 mm, degrading
    below (gap velocity fell 1.8 -> 0.5 mm/s from 8 -> 2 mm spacing,
    i.e. h roughly halved at 2 mm)."""
    if gap_mm >= 6.0:
        return 1.0
    return max(floor, (max(gap_mm, 0.3) / 6.0) ** expo)

def h_cell_side(fl, T_s, T_bulk, H_cell, D_cell, gap_mm, u_oil,
                axial=False) -> dict:
    """Oil film on the can wall. Forced term: Churchill-Bernstein
    crossflow for horizontal sweeps (stirred / serpentine channels), or
    a laminar flat-plate boundary layer along the can height for AXIAL
    upward flow (bottom propeller) - at equal velocity axial flow gives
    a thinner-growing but longer boundary layer and a LOWER film than
    crossflow. Assisting buoyancy is captured by the cube-root blend
    with the natural-convection term."""
    p = film_props(fl, 0.5 * (T_s + T_bulk))
    Ra = rayleigh(p, T_s - T_bulk, H_cell)
    Nun = nu_vertical_cc(Ra, p["Pr"]) * gap_factor(gap_mm)
    h_n = Nun * p["k"] / H_cell
    h_f = 0.0
    Re = u_oil * D_cell / p["nu"]
    if u_oil > 1e-6:
        if axial:
            Re_L = u_oil * H_cell / p["nu"]
            h_f = (0.664 * math.sqrt(max(Re_L, 1.0))
                   * p["Pr"] ** (1.0 / 3.0)) * p["k"] / H_cell
        else:
            h_f = nu_crossflow_cb(Re, p["Pr"]) * p["k"] / D_cell
    return dict(h=blend_mixed(h_n, h_f), h_nat=h_n, h_for=h_f,
                Ra=Ra, Re=Re, Pr=p["Pr"])

def h_tube_side(fl, T_bulk, T_wall, D_o, u_oil) -> dict:
    p = film_props(fl, 0.5 * (T_bulk + T_wall))
    Ra = rayleigh(p, T_bulk - T_wall, D_o)
    h_n = nu_horiz_cyl_cc(Ra, p["Pr"]) * p["k"] / D_o
    h_f = 0.0
    Re = u_oil * D_o / p["nu"]
    if u_oil > 1e-6:
        h_f = nu_crossflow_cb(Re, p["Pr"]) * p["k"] / D_o
    return dict(h=blend_mixed(h_n, h_f), h_nat=h_n, h_for=h_f,
                Ra=Ra, Re=Re, Pr=p["Pr"])

def h_water_inside(loop: dict, mdot_tube: float, d_i: float, L: float,
                   lam_nu: float = 3.66, P_wet: float = None) -> dict:
    """Hausen (laminar, entry-corrected) / Gnielinski (turbulent) with a
    linear bridge across the 2300-3000 transition. Delegates to
    correlations.water_nu so the zonal solver uses the identical curve.
    d_i is the hydraulic diameter; for non-round ducts pass the wetted
    perimeter P_wet (Re = 4 mdot / (mu P_wet)) and the shape's laminar
    asymptote lam_nu (Shah & London)."""
    mu, k, cp = loop["mu"], loop["k"], loop["cp"]
    Pr = mu * cp / k
    Pw = P_wet if P_wet else math.pi * d_i
    Re = 4.0 * mdot_tube / (mu * Pw) if mdot_tube > 0 else 0.0
    Nu, regime = _water_nu(Re, Pr, d_i, L, lam_nu)
    return dict(h=Nu * k / d_i, Re=Re, Pr=Pr, Nu=Nu, regime=regime)

# ------------------------------------------------------------------ #
#  Annular fins (Schmidt approximation)                               #
# ------------------------------------------------------------------ #
def fin_pack(d_o, H_f, t_f, p_f, k_fin, h_oil) -> dict:
    """Per metre of finned tube: bare area, fin area, Schmidt efficiency."""
    r1 = d_o / 2.0
    r2 = r1 + H_f
    r2c = r2 + t_f / 2.0
    n_per_m = 1.0 / p_f
    A_bare = math.pi * d_o * max(0.0, 1.0 - t_f / p_f)
    A_fin_each = 2.0 * math.pi * (r2c ** 2 - r1 ** 2)
    A_fin = n_per_m * A_fin_each
    m = math.sqrt(2.0 * max(h_oil, 1.0) / (k_fin * t_f))
    Lc = H_f + t_f / 2.0
    phi = 1.0 + 0.35 * math.log(r2c / r1)
    x = m * Lc * phi
    eta = math.tanh(x) / x if x > 1e-9 else 1.0
    A_eff = A_bare + eta * A_fin
    m_per_m = n_per_m * (math.pi * (r2 ** 2 - r1 ** 2) * t_f)  # fin metal volume/m
    return dict(A_eff_per_m=A_eff, A_bare_per_m=A_bare, A_fin_per_m=A_fin,
                eta=eta, area_gain=A_eff / (math.pi * d_o),
                fin_metal_vol_per_m=m_per_m, fin_gap=p_f - t_f)

# ------------------------------------------------------------------ #
#  Geometry and mass build-up                                         #
# ------------------------------------------------------------------ #
K_TUBE = {"Copper": 385.0, "Aluminium": 205.0, "Stainless steel": 16.0}
RHO_TUBE = {"Copper": 8940.0, "Aluminium": 2700.0, "Stainless steel": 7900.0}


# Shah & London fully-developed laminar constants for rectangular ducts,
# by aspect ratio b/a (short/long side): f*Re and Nu_T (const wall T).
_RECT_ASP = [0.125, 0.25, 1.0 / 3.0, 0.5, 1.0]
_RECT_FRE = [82.34, 72.93, 68.36, 62.19, 56.91]
_RECT_NUT = [5.60, 4.44, 3.96, 3.39, 2.98]


def _interp(x, xs, ys):
    if x <= xs[0]:
        return ys[0]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            w = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + w * (ys[i] - ys[i - 1])
    return ys[-1]


def tube_section(d) -> dict:
    """Water-tube cross-section properties for Round / Square /
    Rectangular tubes. Returns flow area, wetted perimeters, hydraulic
    diameter, metal cross-section, the flat contact width available to
    bond a plate, and the shape's laminar constants (Shah & London:
    f*Re and Nu_T). Round keeps the classical 64 / 3.66."""
    t = d.get("tube_wall", 0.001)
    shape = d.get("tube_shape", "Round")
    if shape == "Square":
        a_o = d.get("tube_w", d["tube_od"])
        a_i = max(a_o - 2 * t, 1e-4)
        A_in, P_in, P_out = a_i ** 2, 4 * a_i, 4 * a_o
        D_h = a_i
        A_metal = a_o ** 2 - a_i ** 2
        A_out_cs = a_o ** 2
        contact_w, asp = a_o, 1.0
        fRe, lam_nu = 56.91, 2.98
    elif shape == "Rectangular":
        w_o = d.get("tube_w", 0.012)
        h_o = d.get("tube_h", 0.008)
        w_i, h_i = max(w_o - 2 * t, 1e-4), max(h_o - 2 * t, 1e-4)
        A_in = w_i * h_i
        P_in = 2 * (w_i + h_i)
        P_out = 2 * (w_o + h_o)
        D_h = 4 * A_in / P_in
        A_metal = w_o * h_o - w_i * h_i
        A_out_cs = w_o * h_o
        contact_w = h_o                    # vertical face bonded to plate
        asp = min(w_i, h_i) / max(w_i, h_i)
        fRe = _interp(asp, _RECT_ASP, _RECT_FRE)
        lam_nu = _interp(asp, _RECT_ASP, _RECT_NUT)
    else:                                   # Round
        od = d["tube_od"]
        a_i = max(od - 2 * t, 1e-4)
        A_in = math.pi * a_i ** 2 / 4
        P_in, P_out = math.pi * a_i, math.pi * od
        D_h = a_i
        A_metal = math.pi / 4 * (od ** 2 - a_i ** 2)
        A_out_cs = math.pi * od ** 2 / 4
        contact_w, asp = 0.0, 1.0           # line contact only
        fRe, lam_nu = 64.0, 3.66
    return dict(shape=shape, A_in=A_in, P_in=P_in, P_out=P_out, D_h=D_h,
                A_metal=A_metal, A_out_cs=A_out_cs, contact_w=contact_w,
                aspect=asp, fRe=fRe, lam_nu=lam_nu, t=t)


def build_geometry(d: dict) -> dict:
    N = d["Ns"] * d["Np"]
    D, H, p = d["d_cell"], d["h_cell"], d["pitch"]
    gap = (p - D) * 1000.0                                   # mm
    n_cols = math.ceil(math.sqrt(N))
    n_rows = math.ceil(N / n_cols)
    row_pitch = p * (math.sqrt(3) / 2 if d["arrangement"] == "Hexagonal" else 1.0)
    edge = d["edge_margin"]
    Lx = n_cols * p + 2 * edge
    Ly = (n_rows - 1) * row_pitch + p + 2 * edge
    Lz = d["bottom_gap"] + H + d["tube_zone"] + d["gas_gap"]
    fill_h = Lz - d["gas_gap"]

    # heat exchanger tubes (horizontal, in the tube zone above the cells)
    ts = tube_section(d)
    d_i = ts["D_h"]                 # hydraulic diameter for all shapes
    L_tube = max(Lx - 2 * d["manifold_margin"], 0.05) * d["passes"]
    A_tube_bare = ts["P_out"] * L_tube * d["n_tubes"]
    A_tube_in = ts["P_in"] * L_tube * d["n_tubes"]

    # areas and volumes
    f_ends = d["end_fraction"]
    A_cells = N * (math.pi * D * H + f_ends * 2 * math.pi * D ** 2 / 4)
    V_box_fill = Lx * Ly * fill_h
    V_cells = N * math.pi * D ** 2 / 4 * H
    V_tubes = d["n_tubes"] * L_tube * ts["A_out_cs"]
    A_box_ext = 2 * (Lx * Ly + Lx * Lz + Ly * Lz)

    free_per_cell = max(p * row_pitch - math.pi * D ** 2 / 4, 1e-6)
    blk = 1.0 - d.get("holder_block", 0.0)           # cell-holder blockage
    A_flow = N * free_per_cell * blk                 # riser flow area, m²
    D_h = 4.0 * free_per_cell / (math.pi * D) * math.sqrt(max(blk, 0.05))

    return dict(N=N, gap_mm=gap, n_cols=n_cols, n_rows=n_rows,
                A_flow=A_flow, D_h=D_h, row_pitch=row_pitch,
                Lx=Lx, Ly=Ly, Lz=Lz, fill_h=fill_h,
                d_i=d_i, L_tube=L_tube, A_tube_bare=A_tube_bare,
                A_tube_in=A_tube_in, A_cells=A_cells,
                V_box_fill=V_box_fill, V_cells=V_cells, V_tubes=V_tubes,
                A_box_ext=A_box_ext, tube_sec=ts)

def enclosure_calc(d, g):
    """Wall thickness of the largest flat panel as a stiffened plate under
    the burst-disc set pressure (dominates over oil static head), then mass
    from total surface area. sigma_allow and the stiffening knock-down are
    sliders; override entirely with struct_mass > 0."""
    p_des = max(d["p_des_bar"] * 1e5, 900.0 * G * g["fill_h"])
    b = min(g["Lx"], g["Ly"])
    t_flat = b * math.sqrt(0.31 * p_des / (d["sigma_MPa"] * 1e6))
    t_eff = max(t_flat * d["stiff"], 0.0015)
    m = g["A_box_ext"] * t_eff * 2700.0 * 1.18      # + fasteners, feedthroughs
    return dict(t_mm=t_eff * 1000, m=m, p_des_bar=p_des / 1e5)

def build_masses(d, g, fl, finres) -> dict:
    V_fins = finres["fin_metal_vol_per_m"] * g["L_tube"] * d["n_tubes"] if d["fins_on"] else 0.0
    V_oil = max(g["V_box_fill"] - g["V_cells"] - g["V_tubes"] - V_fins, 1e-4)
    m_oil = V_oil * fl["rho"]
    m_cells = g["N"] * d["m_cell"]
    rho_t = RHO_TUBE[d["tube_mat"]]
    _ts = g.get("tube_sec") or tube_section(d)
    m_tubes = rho_t * d["n_tubes"] * g["L_tube"] * _ts["A_metal"]
    m_fins = V_fins * (2700.0 if d["fin_mat"] == "Aluminium" else 8940.0)
    enc = enclosure_calc(d, g)
    m_struct = d["struct_mass"] if d["struct_mass"] > 0 else enc["m"]
    m_holders = g["N"] * d["m_holder_g"] / 1000.0
    m_bus = busbar_props(d, g)["m"]
    m_plates = plate_fin_area(d, g, 100.0)[2] if d.get("plate_on") else 0.0
    m_pack = m_cells + m_oil + m_tubes + m_fins + m_struct + m_holders + m_bus \
             + m_plates
    E_kwh = d["Ns"] * d["Np"] * d["v_nom"] * d["cap_Ah"] / 1000.0
    V_outer_L = (g["Lx"] + 2 * enc["t_mm"] / 1000) * (g["Ly"] + 2 * enc["t_mm"] / 1000) \
                * (g["Lz"] + 2 * enc["t_mm"] / 1000) * 1000
    return dict(V_oil_L=V_oil * 1000, m_oil=m_oil, m_cells=m_cells,
                m_tubes=m_tubes, m_fins=m_fins, m_struct=m_struct,
                m_holders=m_holders, m_bus=m_bus, m_plates=m_plates, enc=enc,
                m_pack=m_pack, E_kwh=E_kwh, V_outer_L=V_outer_L,
                whkg_pack=E_kwh * 1000 / m_pack,
                whkg_cells=E_kwh * 1000 / m_cells,
                whl_pack=E_kwh * 1000 / max(V_outer_L, 1.0),
                C_batt=m_cells * d["cp_cell"],
                C_oil=m_oil * fl["cp"])

# ------------------------------------------------------------------ #
#  Steady-state network solver                                        #
# ------------------------------------------------------------------ #
def solve_steady(d, g, fl, Q_total, T_amb, C_rate=None) -> dict:
    """Two-node (battery, oil) network with h(dT) fixed-point iteration.
    Chain: cells -> oil film -> bulk oil -> oil film on tubes (finned) ->
    tube wall -> water film -> water; parallel leak oil -> ambient.
    If C_rate is given, heat generation is recomputed from DCIR(T) each
    iteration; otherwise Q_total is a fixed heater power (benchmark mode).
    Oil-film h values carry the calibration factor d['cal_h']."""
    loop = WATER_LOOP[d["loop_fluid"]]
    mdot_tot = d["flow_lpm"] / 60.0 * loop["rho"] / 1000.0
    mdot_tube = mdot_tot / max(d["n_tubes"], 1)
    ts = g.get("tube_sec") or tube_section(d)
    wat = h_water_inside(loop, mdot_tube, g["d_i"], g["L_tube"],
                         lam_nu=ts["lam_nu"], P_wet=ts["P_in"])
    R_in = 1.0 / max(wat["h"] * g["A_tube_in"], 1e-9)
    if ts["shape"] == "Round":
        R_wall = math.log(d["tube_od"] / g["d_i"]) / (
            2 * math.pi * K_TUBE[d["tube_mat"]] * g["L_tube"]
            * d["n_tubes"])
    else:                                   # flat walls: t / (k A_mean)
        P_m = 0.5 * (ts["P_in"] + ts["P_out"])
        R_wall = ts["t"] / (K_TUBE[d["tube_mat"]] * P_m * g["L_tube"]
                            * d["n_tubes"])
    R_atm = 1.0 / max(d["h_ext"] * g["A_box_ext"], 1e-9)

    # initial guesses
    T_w_in = d["T_water_in"]
    T_il = T_w_in + 8.0
    T_b = T_il + 6.0
    T_wall = T_w_in + 2.0
    cell = tube = fin = None
    ch = d.get("cal_h", 1.0)
    u_ts, dT_loop = 0.0, 0.0
    for _ in range(60):
        if C_rate is not None:
            I = C_rate * d["cap_Ah"]
            Q_total = g["N"] * I * I * r_of_T(d, T_b) * 1e-3
            Q_total += (I * d["Np"]) ** 2 * busbar_props(d, g)["R"]
        u_ts, dT_loop = thermosiphon_u(d, g, fl, max(Q_total, 1.0),
                                       0.5 * (T_b + T_il))
        u_eff = max(d["u_oil"], u_ts)
        cell = h_cell_side(fl, T_b, T_il, d["h_cell"], d["d_cell"], g["gap_mm"], u_eff,
                           axial=d.get("circ", "").startswith("Bottom"))
        tube0 = h_tube_side(fl, T_il, T_wall, d["tube_od"], u_eff)
        for hh in (cell, tube0):
            hh["h"] *= ch; hh["h_nat"] *= ch; hh["h_for"] *= ch
        if d["fins_on"]:
            fin = fin_pack(d["tube_od"], d["fin_h"], d["fin_t"], d["fin_p"],
                           205.0 if d["fin_mat"] == "Aluminium" else 385.0, tube0["h"])
            A_oilside = fin["A_eff_per_m"] * g["L_tube"] * d["n_tubes"]
            A_pl, eta_pl, _mpl = plate_fin_area(d, g, ch * tube0["h"])
            A_oilside += A_pl
        else:
            fin = dict(A_eff_per_m=math.pi * d["tube_od"], eta=1.0, area_gain=1.0,
                       fin_metal_vol_per_m=0.0, fin_gap=1.0, A_fin_per_m=0.0,
                       A_bare_per_m=math.pi * d["tube_od"])
            A_oilside = g["A_tube_bare"]
            A_pl, eta_pl, _mpl = plate_fin_area(d, g, ch * tube0["h"])
            A_oilside += A_pl
        R_b = 1.0 / max(cell["h"] * g["A_cells"], 1e-9)
        R_ot = 1.0 / max(tube0["h"] * A_oilside, 1e-9)
        R_chain = R_ot + R_wall + R_in

        # split heat between water chain and ambient leak
        dT_w_rise = Q_total / max(mdot_tot * loop["cp"], 1e-9)
        T_sink = T_w_in + 0.5 * min(dT_w_rise, 60.0)
        T_il_new = (Q_total + T_sink / R_chain + T_amb / R_atm) / (1.0 / R_chain + 1.0 / R_atm)
        Q_w = (T_il_new - T_sink) / R_chain
        dT_w_rise = max(Q_w, 0.0) / max(mdot_tot * loop["cp"], 1e-9)
        T_wall_new = T_sink + Q_w * (R_in + R_wall)
        T_b_new = T_il_new + Q_total * R_b
        # relax
        T_il += 0.6 * (T_il_new - T_il)
        T_b += 0.6 * (T_b_new - T_b)
        T_wall += 0.6 * (T_wall_new - T_wall)
        tube = tube0
    Q_w = (T_il - (T_w_in + 0.5 * Q_total / max(mdot_tot * loop["cp"], 1e-9))) / R_chain
    Q_atm = (T_il - T_amb) / R_atm
    Q_cell = Q_total / g["N"]
    T_core = T_b + Q_cell * r_core(d)
    dT_water = Q_total / max(mdot_tot * loop["cp"], 1e-9)
    if d.get("tube_plane") == "Interstitial (between rows)" or d.get("plate_on"):
        dT_loop *= 0.35                          # distributed sinks
    dT_pos = 0.5 * dT_water + 0.5 * dT_loop     # worst-position penalty
    return dict(T_b=T_b, T_il=T_il, T_wall=T_wall,
                T_core=T_core, dT_core=Q_cell * r_core(d),
                u_ts=u_ts, dT_loop=dT_loop, Q_eff=Q_total,
                T_worst=T_b + dT_pos, T_best=T_b - dT_pos, spread=2 * dT_pos,
                R_b=R_b, R_ot=R_ot, R_wall=R_wall, R_in=R_in, R_atm=R_atm,
                h_cell=cell["h"], h_cell_nat=cell["h_nat"], h_cell_for=cell["h_for"],
                Ra_cell=cell["Ra"], Re_cell=cell["Re"],
                h_tube=tube["h"], h_tube_nat=tube["h_nat"], h_tube_for=tube["h_for"],
                Ra_tube=tube["Ra"],
                h_water=wat["h"], Re_water=wat["Re"], water_regime=wat["regime"],
                A_oilside=A_oilside, A_plate=A_pl, eta_plate=eta_pl, fin=fin,
                dT_water=Q_total / max(mdot_tot * loop["cp"], 1e-9),
                mdot_tot=mdot_tot, Q_w=Q_w, Q_atm=Q_atm,
                gapf=gap_factor(g["gap_mm"]))

def r_of_T(d, T_cell_C: float) -> float:
    """Cell DCIR in mΩ at temperature T. Exponential fall with T
    (default -1.2 %/K around 25 °C) - the reason AMG run a 45 °C
    set point. Set k_dcir = 0 to disable the coupling."""
    return d["r_dc"] * math.exp(-d.get("k_dcir", 0.0) * (T_cell_C - 25.0))

def r_core(d) -> float:
    """Peak core-to-surface resistance of a cylindrical cell with uniform
    volumetric generation: dT = q''' R^2/(4 k_r) -> R_th = 1/(4 pi k_r H).
    Radial jellyroll k_r ~ 0.8-1.0 W/mK. Mean-to-surface is half this."""
    return 1.0 / (4.0 * math.pi * d.get("k_rad", 0.9) * d["h_cell"])

def thermosiphon_u(d, g, fl, Q: float, T_oil_C: float):
    """Self-circulation velocity in the cell gaps: buoyant head
    rho*beta*g*H_loop*dT_loop against laminar loop friction plus minor
    losses, with dT_loop = Q/(rho u A cp). Solved by bisection.
    H_loop = vertical offset between cell mid-height and the tube plane,
    so tube placement is a live design variable. Order-of-magnitude
    validation: Wang et al. measured 0.5-1.8 mm/s gap velocities."""
    if d.get("tube_plane", "Top of pack") == "Top of pack":
        H_loop = d["h_cell"] / 2 + d["tube_zone"] / 2
    elif d.get("tube_plane") == "Interstitial (between rows)":
        H_loop = 0.008
    elif d.get("tube_plane") == "Mid-height":
        H_loop = 0.01
    else:                                   # below the cells
        return 1e-5, Q / max(fl["rho"] * 1e-5 * g["A_flow"] * fl["cp"], 1e-9)
    p = film_props(fl, T_oil_C)
    A, Dh = g["A_flow"], g["D_h"]
    L_loop = 2.2 * g["fill_h"]
    K = d.get("K_loop", 5.0)
    def resid(u):
        drive = p["beta"] * G * H_loop * Q / (u * A * p["cp"])
        fric = 32.0 * p["rho"] * p["nu"] * L_loop * u / Dh ** 2 \
               + K * 0.5 * p["rho"] * u ** 2
        return drive - fric
    lo, hi = 1e-6, 0.08
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if resid(mid) > 0:
            lo = mid
        else:
            hi = mid
    u = lo
    dT_loop = Q / (p["rho"] * u * A * p["cp"])
    return u, min(dT_loop, 60.0)

def q_gen_per_cell(d, C_rate, T_cell_C: float = 25.0) -> float:
    I = C_rate * d["cap_Ah"]
    return I * I * r_of_T(d, T_cell_C) * 1e-3   # DCIR in mΩ

def max_continuous_C(d, g, fl, T_amb, T_limit) -> float:
    lo, hi = 0.05, 12.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        res = solve_steady(d, g, fl, 1.0, T_amb, C_rate=mid)
        T_check = res["T_core"] if d.get("limit_core", False) else res["T_b"]
        if T_check > T_limit:
            hi = mid
        else:
            lo = mid
    return lo

# ------------------------------------------------------------------ #
#  Transient solver                                                   #
# ------------------------------------------------------------------ #
def duty_profile(kind, dur, C1, t1, C2, t2, csv_tc=None):
    """Return times [s] and C-rate arrays. csv_tc = (t, C) from an upload."""
    t = np.arange(0.0, dur + 1e-9, 2.0)
    if kind == "CSV upload" and csv_tc is not None:
        t = np.arange(0.0, min(dur, csv_tc[0][-1]) + 1e-9, 2.0)
        C = np.interp(t, csv_tc[0], csv_tc[1])
    elif kind == "Constant C":
        C = np.full_like(t, C1)
    elif kind == "Fast charge then rest":
        C = np.where(t < t1, C1, 0.0)
    else:  # pulse train
        period = max(t1 + t2, 1.0)
        C = np.where((t % period) < t1, C1, C2)
    return t, C

def duty_from_csv(file, d) -> tuple:
    """Parse a duty CSV with columns t_s plus either C or P_kW."""
    df = pd.read_csv(file)
    cols = {c.lower().strip(): c for c in df.columns}
    t = df[cols["t_s"]].to_numpy(dtype=float)
    if "c" in cols:
        C = df[cols["c"]].to_numpy(dtype=float)
    elif "p_kw" in cols:
        P = df[cols["p_kw"]].to_numpy(dtype=float) * 1000.0
        C = P / (d["Ns"] * d["Np"] * d["v_nom"] * d["cap_Ah"])
    else:
        raise ValueError("CSV needs columns t_s and C (or P_kW)")
    order = np.argsort(t)
    return np.abs(t[order]), np.abs(C[order])

def solve_transient(d, g, fl, masses, T_amb, t_arr, C_arr) -> dict:
    loop = WATER_LOOP[d["loop_fluid"]]
    mdot = d["flow_lpm"] / 60.0 * loop["rho"] / 1000.0
    C_b, C_il = masses["C_batt"], masses["C_oil"]
    T_b = np.zeros_like(t_arr); T_il = np.zeros_like(t_arr)
    T_b[0] = T_il[0] = d["T_start"]
    # freeze wall/water resistances from a representative solve, update oil films each step
    rep = solve_steady(d, g, fl, 1.0, T_amb, C_rate=max(C_arr.max(), 0.5))
    R_wall, R_in, R_atm = rep["R_wall"], rep["R_in"], rep["R_atm"]
    ch = d.get("cal_h", 1.0)
    T_core = np.zeros_like(t_arr); T_core[0] = d["T_start"]
    Q_tr = np.zeros_like(t_arr)
    for i in range(1, len(t_arr)):
        dt = t_arr[i] - t_arr[i - 1]
        Q = q_gen_per_cell(d, C_arr[i - 1], T_b[i - 1]) * g["N"]
        u_ts, _ = thermosiphon_u(d, g, fl, max(Q, 1.0), 0.5 * (T_b[i-1] + T_il[i-1]))
        u_eff = max(d["u_oil"], u_ts)
        cell = h_cell_side(fl, T_b[i-1], T_il[i-1], d["h_cell"], d["d_cell"], g["gap_mm"], u_eff,
                           axial=d.get("circ", "").startswith("Bottom"))
        T_wall_est = T_il[i-1] - 0.6 * (T_il[i-1] - d["T_water_in"])
        tub = h_tube_side(fl, T_il[i-1], T_wall_est, d["tube_od"], u_eff)
        R_b = 1.0 / max(ch * cell["h"] * g["A_cells"], 1e-9)
        R_ot = 1.0 / max(ch * tub["h"] * rep["A_oilside"], 1e-9)
        Q_bi = (T_b[i-1] - T_il[i-1]) / R_b
        Q_w = (T_il[i-1] - (d["T_water_in"] + 0.5 * max(Q_bi, 0) / max(mdot * loop["cp"], 1e-9))) / (R_ot + R_wall + R_in)
        Q_a = (T_il[i-1] - T_amb) / R_atm
        T_b[i] = T_b[i-1] + dt * (Q - Q_bi) / C_b
        T_il[i] = T_il[i-1] + dt * (Q_bi - Q_w - Q_a) / C_il
        T_core[i] = T_b[i] + (Q / g["N"]) * r_core(d)   # quasi-steady radial
        Q_tr[i] = Q
    return dict(t=t_arr, T_b=T_b, T_il=T_il, T_core=T_core, Q=Q_tr)

# ------------------------------------------------------------------ #
#  Benchmark against Wang et al. (2023)                               #
# ------------------------------------------------------------------ #
def benchmark_wang() -> dict:
    """Rebuild the paper's rig with this app's correlations and compare with
    their measured/derived values (Table 6, Figs 5-8). Prismatic cells, so
    the vertical-plate correlation applies directly."""
    oil = dict(name="Transformer oil (paper Table 3)", family="mineral",
               k=0.13, rho=875.0, cp=1900.0, nu25=17.0e-6, beta=7.5e-4,
               B=3900.0, dielectric=True, bp=280, flash=150)
    A_cells = 6 * (2 * (0.148 * 0.097) + 2 * (0.148 * 0.027) + 2 * (0.097 * 0.027))
    N_t, d_o, d_i, L_t = 4, 0.006, 0.005, 0.222
    A_ot = math.pi * d_o * L_t * N_t          # 0.0167 m2, matches paper
    A_in = math.pi * d_i * L_t * N_t
    loop = dict(rho=1000.0, cp=4200.0, k=0.58, mu=1.3e-3)  # water ~8 °C
    mdot = 17.1e-6 * 1000.0
    wat = h_water_inside(loop, mdot / N_t, d_i, L_t)
    R_in = 1.0 / (wat["h"] * A_in)
    R_wall = math.log(d_o / d_i) / (2 * math.pi * 385.0 * L_t * N_t)
    R_atm_meas = 6.7                           # take the paper's measured value
    C_b, C_il = 7620.0, 9184.0
    Q = 16.4 * 6                               # W, model calibration Table 5 at 2C
    T_amb, T_w = 25.0, 5.0
    t = np.arange(0, 1801.0, 2.0)
    T_b = np.full_like(t, 24.5); T_il = np.full_like(t, 24.5)
    hb_last = ht_last = 0.0
    for i in range(1, len(t)):
        pf = film_props(oil, 0.5 * (T_b[i-1] + T_il[i-1]))
        Ra = rayleigh(pf, T_b[i-1] - T_il[i-1], 0.148)
        hb = nu_vertical_cc(Ra, pf["Pr"]) * pf["k"] / 0.148
        Twall = T_w + 3.0
        pt = film_props(oil, 0.5 * (T_il[i-1] + Twall))
        Rad = rayleigh(pt, T_il[i-1] - Twall, d_o)
        ht = nu_horiz_cyl_cc(Rad, pt["Pr"]) * pt["k"] / d_o
        R_b = 1.0 / (hb * A_cells); R_ot = 1.0 / (ht * A_ot)
        Q_bi = (T_b[i-1] - T_il[i-1]) / R_b
        Q_w = (T_il[i-1] - (T_w + 0.5 * max(Q_bi, 0) / (mdot * loop["cp"]))) / (R_ot + R_wall + R_in)
        Q_a = (T_il[i-1] - T_amb) / R_atm_meas
        T_b[i] = T_b[i-1] + 2.0 * (Q - Q_bi) / C_b
        T_il[i] = T_il[i-1] + 2.0 * (Q_bi - Q_w - Q_a) / C_il
        hb_last, ht_last = hb, ht
    rows = [
        ("h battery-to-oil  [W/m²·K]", f"{hb_last:.0f}", "~80-110 (Fig. 11b, measured)"),
        ("h oil-to-tube  [W/m²·K]", f"{ht_last:.0f}", "~170-270 (Fig. 11b, measured)"),
        ("R battery-oil  [K/W]", f"{1.0/(hb_last*A_cells):.3f}", "0.04 (Table 6)"),
        ("R oil-tube  [K/W]", f"{1.0/(ht_last*A_ot):.2f}", "0.30 (Table 6)"),
        ("R water film  [K/W]", f"{R_in:.3f}", "0.06 (Table 6)"),
        ("R tube wall  [K/W]", f"{R_wall:.5f}", "0.0005 (Table 6)"),
        ("T battery at 1800 s  [°C]", f"{T_b[-1]:.1f}", "~32.3 (Fig. 5/6, measured)"),
        ("T oil at 1800 s  [°C]", f"{T_il[-1]:.1f}", "~28.5 (Fig. 8b)"),
        ("dT battery-oil plateau [K]", f"{T_b[-1]-T_il[-1]:.1f}", "~3.2 (Fig. 8b)"),
    ]
    return dict(rows=rows, t=t, T_b=T_b, T_il=T_il)

# ------------------------------------------------------------------ #
#  Advice engine and sensitivity study                                #
# ------------------------------------------------------------------ #
def diagnose(d, g, fl, masses, res, T_limit, Q) -> list:
    msgs = []
    chain = {"cell-to-oil film": res["R_b"], "oil-to-tube film (finned)": res["R_ot"],
             "tube wall": res["R_wall"], "water film inside tubes": res["R_in"]}
    worst = max(chain, key=chain.get)
    tot = sum(chain.values())
    msgs.append(("info", f"**Bottleneck: {worst}** carries {100*chain[worst]/tot:.0f}% of the "
                 f"cell-to-water resistance ({chain[worst]*1000:.1f} mK/W of {tot*1000:.1f} mK/W)."))
    if worst == "cell-to-oil film":
        msgs.append(("do", "Cell film dominates: add gentle stirring (a few cm/s), widen the "
                     "cell gap towards 6 mm, or switch to a lower-viscosity dielectric. "
                     "Fins on the tubes will NOT help while this film dominates."))
    if worst == "oil-to-tube film (finned)":
        msgs.append(("do", "Tube-side film dominates: add/extend fins, add tubes, or stir. "
                     "This is area-starved, exactly as in Wang et al. (R = 0.3 K/W there)."))
    if worst == "water film inside tubes":
        msgs.append(("do", f"Water film dominates and flow is **{res['water_regime']}** "
                     f"(Re = {res['Re_water']:.0f}). Raise flow past Re 3000, or use more "
                     "smaller tubes in parallel; in laminar flow extra velocity does nothing."))
    msgs.append(("info", f"Predicted self-circulation (thermosiphon): "
                 f"**{res['u_ts']*1000:.1f} mm/s** in the cell gaps (Wang et al. measured "
                 f"0.5-1.8 mm/s), giving a top-to-bottom oil stratification of "
                 f"~{res['dT_loop']:.1f} °C. Tube plane: {d.get('tube_plane','Top of pack')}."))
    if d.get("tube_plane") == "Below the cells":
        msgs.append(("bad", "Tubes below the cells: buoyancy stratifies stably, the "
                     "thermosiphon dies, and hot oil strands at the top. The model's "
                     "well-mixed-oil assumption is optimistic here - expect worse."))
    if res["spread"] > 5.0:
        msgs.append(("warn", f"Estimated best-to-worst cell spread ~{res['spread']:.1f} °C "
                     "exceeds the 5 °C uniformity criterion (water rise + stratification). "
                     "Raise water flow, stir, or split the water loop into counterflowing "
                     "halves."))
    if res["dT_core"] > 3.0:
        msgs.append(("info", f"Core runs ~{res['dT_core']:.1f} °C above the can at this duty "
                     f"(k_r = {d.get('k_rad',0.9):.1f} W/mK). No coolant choice touches this "
                     "term; only lower current or tab/axial extraction do."))
    if res["dT_water"] > 5:
        msgs.append(("warn", f"Water heats by {res['dT_water']:.1f} °C end to end (> 5 K): last "
                     "cells see a warmer sink. Raise flow or split the loop."))
    if g["gap_mm"] < 6:
        msgs.append(("warn", f"Cell gap {g['gap_mm']:.1f} mm < 6 mm: buoyant flow in the gaps "
                     f"is throttled (penalty factor {res['gapf']:.2f}, per Wang et al. Fig. 9)."))
    if res["T_b"] > T_limit:
        msgs.append(("bad", f"Steady cell temperature {res['T_b']:.1f} °C exceeds the "
                     f"{T_limit:.0f} °C limit at this duty. See sensitivity chart for the "
                     "cheapest fix."))
    if not fl["dielectric"]:
        msgs.append(("bad", f"**{fl['name']} is not a dielectric.** Fine as a thermal reference "
                     "in this model, unusable as an immersion fluid in a live pack."))
    if not math.isnan(fl.get("bp", float("nan"))) and res["T_b"] + 10 > fl["bp"]:
        msgs.append(("warn", f"Cell temperature is within 10 K of the fluid boiling point "
                     f"({fl['bp']:.0f} °C): you are entering two-phase territory (pressure "
                     "management needed)."))
    if not math.isnan(fl.get("flash", float("nan"))) and fl["flash"] < 120:
        msgs.append(("warn", f"Flash point {fl['flash']:.0f} °C is low for a lithium pack; "
                     "prefer > 150 °C (ester class)."))
    if d["fins_on"] and res["fin"]["fin_gap"] < 0.004 and d["u_oil"] < 0.005:
        msgs.append(("warn", "Fin gap < 4 mm with still oil: natural-convection boundary "
                     "layers will merge between fins and the Schmidt-efficiency estimate "
                     "becomes optimistic. Open the fin pitch or stir."))
    if masses["m_oil"] / masses["m_pack"] > 0.30:
        msgs.append(("info", f"Oil is {100*masses['m_oil']/masses['m_pack']:.0f}% of pack mass. "
                     "Reduce headspace/tube-zone height, or accept it as buffer thermal mass."))
    exp_L = fl["beta"] * masses["V_oil_L"] * (d.get("T_service_max", 60) - (-10))
    msgs.append(("info", f"Thermal expansion over -10 to {d.get('T_service_max',60):.0f} °C: "
                 f"~{exp_L:.1f} L. Size a bellows/bladder for this; do not leave a free air "
                 "headspace (moisture, tilt)."))
    msgs.append(("info", "Water-in-oil leak is the single-point failure: keep oil static "
                 "pressure above water pressure, or use double-walled tubes with leak "
                 "detection (transformer practice)."))
    return msgs

def sensitivity(d, g, fl, cool_df, T_amb, C_duty, T_limit) -> pd.DataFrame:
    """One-at-a-time perturbations; report change in steady cell temperature."""
    def run(dd, ff):
        gg = build_geometry(dd)
        return solve_steady(dd, gg, ff, 1.0, T_amb, C_rate=C_duty)["T_b"]
    base = run(d, fl)
    cases = []
    dd = dict(d); dd["n_tubes"] = max(1, int(round(d["n_tubes"] * 1.5))); cases.append(("Tubes +50%", run(dd, fl)))
    dd = dict(d)
    if d["fins_on"]:
        dd["fin_p"] = d["fin_p"] * 2 / 3; cases.append(("Fin area +50%", run(dd, fl)))
    else:
        dd["fins_on"] = True; cases.append(("Add fins (8 mm, 4 mm pitch)", run(dd, fl)))
    dd = dict(d); dd["u_oil"] = d["u_oil"] + 0.05; cases.append(("Stir oil +5 cm/s", run(dd, fl)))
    dd = dict(d); dd["flow_lpm"] = d["flow_lpm"] * 2; cases.append(("Water flow x2", run(dd, fl)))
    dd = dict(d); dd["T_water_in"] = d["T_water_in"] - 5; cases.append(("Water inlet -5 K", run(dd, fl)))
    dd = dict(d); dd["pitch"] = d["pitch"] + 0.002; cases.append(("Cell pitch +2 mm", run(dd, fl)))
    diel = cool_df[cool_df["dielectric"]].copy()
    best = diel.loc[diel["nu_cSt"].idxmin()]
    if best["name"] != fl["name"]:
        cases.append((f"Fluid: {best['name']}", run(d, fluid_dict(best))))
    rows = [dict(change=n, T_b=t, dT=t - base) for n, t in cases]
    return pd.DataFrame(rows).sort_values("dT"), base

# ------------------------------------------------------------------ #
#  UI helpers                                                         #
# ------------------------------------------------------------------ #
ACCENT, INK, PAPER = "#F59E0B", "#1F2937", "#FFFFFF"
BRAND_A, BRAND_B = "#7C88F8", "#5BC8E8"
CSS = f"""
<style>
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');
  html, body, .stApp, [class*="css"] {{
      font-family: Inter, -apple-system, 'Segoe UI', Roboto, sans-serif; }}
  .stApp {{ background: {PAPER}; }}
  h1, h2, h3, h4 {{ color: {INK}; letter-spacing: -0.02em; font-weight: 700; }}
  /* hero */
  .hero {{ padding: 12px 16px; border-radius: 16px; margin-bottom: 8px;
      background: linear-gradient(120deg, #EEF2FF 0%, #E6F7FD 100%);
      border: 1px solid #E3E8F4; }}
  .hero-top {{ display:flex; align-items:baseline; gap:12px;
      flex-wrap:wrap; }}
  .hero h1 {{ color:{INK}; margin:0; font-size:1.15rem; display:inline; }}
  .hero .sub {{ color:#64748B; font-size:.82rem; }}
  .hero-stats {{ display:flex; gap:8px; flex-wrap:wrap; margin-top:8px; }}
  .hstat {{ background:rgba(255,255,255,.75); border:1px solid #E3E8F4;
      border-radius:10px; padding:5px 11px; font-size:.72rem;
      color:#64748B; }}
  .hstat b {{ display:block; font-size:.95rem; color:{INK};
      font-weight:700; }}
  .hstat b.ok {{ color:#15803D; }} .hstat b.bad {{ color:#B91C1C; }}
  .chip {{ display:inline-block; padding: 2px 11px; border-radius: 999px;
      font-size: .74rem; font-weight: 600; margin-left:auto;
      background: rgba(124,136,248,.12); color:#4A54D8;
      border:1px solid rgba(124,136,248,.35);}}
  .chip.bad {{ background:#FEE2E2; color:#B91C1C; border-color:#FCA5A5; }}
  .chip.ok  {{ background:#DCFCE7; color:#15803D; border-color:#86EFAC; }}
  /* KPI cards */
  .kpis {{ display:flex; gap: 12px; flex-wrap: wrap; margin: 10px 0 4px 0; }}
  .kpi {{ flex:1 1 140px; background:#FFFFFF; border:1px solid #EEF1F6;
      border-radius: 14px; padding: 12px 14px;
      box-shadow: 0 1px 2px rgba(16,24,40,.04); }}
  .kpi .l {{ font-size:.68rem; font-weight:600; letter-spacing:.06em;
      text-transform: uppercase; color:#64748B; }}
  .kpi .v {{ font-size:1.45rem; font-weight:800; color:{INK};
      margin: 2px 0 0 0; line-height:1.15; }}
  .kpi .v.ok {{ color:#15803D; }} .kpi .v.bad {{ color:#B91C1C; }}
  .kpi .v.brand {{ background: linear-gradient(120deg,{BRAND_A},{BRAND_B});
      -webkit-background-clip:text; background-clip:text; color:transparent; }}
  .kpi .s {{ font-size:.75rem; color:#64748B; margin-top:2px; }}
  /* tabs as pills (selectors covering multiple Streamlit versions) */
  .stTabs [data-baseweb="tab-list"],
  div[data-testid="stTabs"] div[role="tablist"] {{ gap: 6px !important;
      border-bottom: none !important; flex-wrap: wrap; }}
  .stTabs [data-baseweb="tab"], div[data-testid="stTabs"] button[role="tab"],
  .stTabs button[role="tab"] {{ background:#FFFFFF !important;
      border:1px solid #EEF1F6 !important; border-radius: 999px !important;
      padding: 6px 16px !important; color:#475569 !important;
      font-weight:600; font-size:.86rem; }}
  .stTabs [aria-selected="true"],
  div[data-testid="stTabs"] button[role="tab"][aria-selected="true"] {{
      background: linear-gradient(120deg,{BRAND_A},{BRAND_B}) !important;
      color: white !important; border-color: transparent !important; }}
  .stTabs [data-baseweb="tab-highlight"], .stTabs [data-baseweb="tab-border"],
  div[data-testid="stTabs"] div[data-testid="stTabsHighlight"]
      {{ display:none !important; }}
  button[role="tab"] p {{ color: inherit !important; }}
  /* cards (bordered containers) */
  [data-testid="stVerticalBlockBorderWrapper"] {{ background:#FFFFFF;
      border:1px solid #EEF1F6 !important; border-radius:16px !important;
      box-shadow: 0 1px 2px rgba(16,24,40,.04);
      padding: 6px 10px !important; }}
  /* sidebar */
  [data-testid="stSidebar"] {{ background:#FFFFFF;
      border-right:1px solid #E7EAF0; }}
  /* buttons */
  .stButton > button, .stDownloadButton > button {{ border-radius: 10px;
      border:1px solid #EEF1F6; font-weight:600; }}
  .stButton > button[kind="primary"] {{
      background: linear-gradient(120deg,{BRAND_A},{BRAND_B}); border:none; }}
  /* expanders */
  [data-testid="stExpander"] {{ border:1px solid #EEF1F6; border-radius:12px;
      background:#FFFFFF; }}
  /* progress */
  [data-testid="stProgress"] > div > div > div {{
      background: linear-gradient(90deg,{BRAND_A},{BRAND_B}); }}
  .small-note {{ color:#64748B; font-size:0.84rem; }}
  /* silk: transitions + hover lift */
  [data-testid="stVerticalBlockBorderWrapper"], .kpi, .stButton > button,
  .stTabs [data-baseweb="tab"] {{ transition: box-shadow .18s ease,
      transform .18s ease, background .18s ease; }}
  [data-testid="stVerticalBlockBorderWrapper"]:hover, .kpi:hover {{
      box-shadow: 0 6px 18px rgba(16,24,40,.10); transform: translateY(-1px); }}
  .stTabs [data-baseweb="tab"]:hover {{ background:#F1F5F9; }}
  ::-webkit-scrollbar {{ width: 8px; height: 8px; }}
  ::-webkit-scrollbar-thumb {{ background:#CBD5E1; border-radius: 8px; }}
  /* sticky live design panel: several strategies for different builds */
  .st-key-liveview {{ position: sticky; top: 0.8rem; z-index: 3;
      max-height: calc(100vh - 1.6rem); overflow-y: auto; }}
  div[data-testid="stColumn"]:has(.st-key-liveview),
  div[data-testid="column"]:has(.st-key-liveview) {{
      position: sticky; top: 0.8rem; align-self: flex-start;
      height: fit-content; z-index: 3; }}
  div[data-testid="stHorizontalBlock"]:has(.st-key-liveview) {{
      overflow: visible !important; align-items: flex-start; }}
</style>"""

PLOTCFG = dict(displaylogo=False,
               modeBarButtonsToAdd=["drawline", "drawrect", "eraseshape",
                                    "toggleSpikelines"])

def kpi_cards(items):
    """items: list of (label, value, sub, tone in '', 'ok', 'bad', 'brand')."""
    cells = "".join(
        f"<div class='kpi'><div class='l'>{l}</div>"
        f"<div class='v {t}'>{v}</div><div class='s'>{s}</div></div>"
        for l, v, s, t in items)
    st.markdown(f"<div class='kpis'>{cells}</div>", unsafe_allow_html=True)

DEFAULTS = dict(
    Ns=108, Np=10, cap_Ah=5.0, v_nom=3.7, r_dc=25.0, m_cell=0.070, cp_cell=950.0,
    d_cell=0.021, h_cell=0.070, arrangement="Square", pitch=0.027, edge_margin=0.010,
    bottom_gap=0.005, tube_zone=0.035, gas_gap=0.010, end_fraction=0.0,
    coolant="MIVOLT DF7", u_oil=0.0,
    n_tubes=16, tube_od=0.010, tube_wall=0.001, tube_mat="Copper", passes=1,
    tube_shape="Round", tube_w=0.012, tube_h=0.008,
    u_loop=0.05, pipe_id=0.019, pipe_len=2.5,
    manifold_margin=0.020,
    fins_on=True, fin_h=0.008, fin_t=0.0005, fin_p=0.004, fin_mat="Aluminium",
    loop_fluid="Water", flow_lpm=10.0, T_water_in=20.0,
    duty="Constant C", C1=2.0, t1=900.0, C2=0.5, t2=600.0, duration=3600.0,
    T_start=25.0, T_limit=45.0, T_amb=25.0, h_ext=5.0, struct_mass=0.0,
    T_service_max=60.0,
    # v2 additions
    k_dcir=0.012, k_rad=0.9, limit_core=False, tube_plane="Top of pack",
    K_loop=5.0, cal_h=1.0,
    E_tr=55.0, frac_oil=0.6, zone_pitches=1.5, vent_L=5.0,
    # v3 additions
    fmt="21700", v_max=4.2, v_cut=0.05, chg_mult=1.10, entropic=True,
    soc0=0.90, soc_min=0.10, C_chg=1.0, dirn="Discharge", track_soc=False,
    cyc_rest=600.0, n_cyc=3,
    R_bus=0.0, bus_J=5.0, m_holder_g=8.0, holder_block=0.20,
    sigma_MPa=80.0, stiff=0.45, p_des_bar=0.5,
    circ="Thermosiphon only", plate_on=False, plate_t=0.0015,
    plate_mat="Aluminium", plate_contact=0.8, u_guided=0.05,
    veh_m=1900.0, CdA=0.62, Crr=0.009, eta_dt=0.92, eta_rg=0.65,
    P_rg=60.0, P_acc=500.0, cycle="WLTP Class 3b", repeat_cyc=True,
)

def resistance_chart(res):
    items = [("Cell-to-oil film", res["R_b"]), ("Oil-to-tube film (finned area)", res["R_ot"]),
             ("Tube wall", res["R_wall"]), ("Water film in tubes", res["R_in"])]
    tot = sum(v for _, v in items)
    worst = max(items, key=lambda x: x[1])[0]
    fig = go.Figure(go.Bar(
        y=[n for n, _ in items][::-1], x=[v * 1000 for _, v in items][::-1],
        orientation="h",
        marker_color=["#EF4444" if n == worst else "#6366F1" for n, _ in items][::-1],
        text=[f"{v*1000:.2f} mK/W  ({100*v/tot:.0f}%)" for _, v in items][::-1],
        textposition="outside"))
    fig.update_layout(height=260, margin=dict(l=10, r=10, t=30, b=10),
                      title="Where the resistance lives (cell -> water chain)",
                      xaxis_title="Thermal resistance [mK/W]",
                      plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)")
    return fig, worst

def waterfall_chart(res, Q, d):
    steps = [("Water inlet", d["T_water_in"], "absolute"),
             ("Water warm-up (mean)", 0.5 * res["dT_water"], "relative"),
             ("Water film", res["Q_w"] * res["R_in"], "relative"),
             ("Tube wall", res["Q_w"] * res["R_wall"], "relative"),
             ("Oil film on tubes", res["Q_w"] * res["R_ot"], "relative"),
             ("Cell-to-oil film", Q * res["R_b"], "relative"),
             ("Cell surface", None, "total")]
    fig = go.Figure(go.Waterfall(
        x=[s[0] for s in steps], measure=[s[2] for s in steps],
        y=[s[1] if s[1] is not None else 0 for s in steps],
        connector=dict(line=dict(color="#94A3B8")),
        increasing=dict(marker=dict(color=ACCENT)),
        totals=dict(marker=dict(color=INK))))
    fig.update_layout(height=330, margin=dict(l=10, r=10, t=30, b=10),
                      title="Temperature waterfall at this duty [°C]",
                      yaxis_title="°C", plot_bgcolor="rgba(255,255,255,0)",
                      paper_bgcolor="rgba(0,0,0,0)", showlegend=False)
    return fig

def chain_schematic(res, Q):
    boxes = [("CELL", "#D96C4F"), ("oil film", "#F3D9A4"), ("BULK OIL", ACCENT),
             ("oil film", "#F3D9A4"), ("TUBE+FINS", "#D97706"), ("wall", "#8C8C8C"),
             ("water film", "#BFD9EA"), ("WATER", "#0EA5E9")]
    drops = [None, Q * res["R_b"], None, res["Q_w"] * res["R_ot"], None,
             res["Q_w"] * res["R_wall"], res["Q_w"] * res["R_in"], None]
    fig = go.Figure()
    x = 0.0
    for (label, colr), dT in zip(boxes, drops):
        w = 1.4 if label.isupper() else 0.9
        fig.add_shape(type="rect", x0=x, x1=x + w, y0=0, y1=1,
                      fillcolor=colr, line=dict(color=INK, width=1))
        fig.add_annotation(x=x + w / 2, y=0.5, text=label, showarrow=False,
                           font=dict(size=12, color=INK))
        if dT is not None:
            fig.add_annotation(x=x + w / 2, y=1.18, text=f"dT = {dT:.1f} °C",
                               showarrow=False, font=dict(size=11, color="#EF4444"))
        x += w + 0.12
    fig.add_annotation(x=x / 2, y=-0.28, showarrow=False,
                       text="Heat flows left to right. The films are where the kelvins are spent; "
                            "the bulk oil and copper are nearly free.",
                       font=dict(size=11, color="#64748B"))
    fig.update_xaxes(visible=False); fig.update_yaxes(visible=False, range=[-0.5, 1.5])
    fig.update_layout(height=210, margin=dict(l=5, r=5, t=10, b=5),
                      plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)")
    return fig

# ------------------------------------------------------------------ #
#  Main app                                                           #
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
#  Parasitic power                                                    #
# ------------------------------------------------------------------ #
def water_pump_power(d, g, loop) -> dict:
    ts = g.get("tube_sec") or tube_section(d)
    mdot_tot = d["flow_lpm"] / 60.0 * loop["rho"] / 1000.0
    mdot_tube = mdot_tot / max(d["n_tubes"], 1)
    A_i = ts["A_in"]
    v = mdot_tube / (loop["rho"] * A_i)
    Re = loop["rho"] * v * g["d_i"] / loop["mu"]
    f = ts["fRe"] / max(Re, 1.0) if Re < 2300 else 0.316 * Re ** -0.25
    dp = (f * g["L_tube"] / g["d_i"] + 6.0) * 0.5 * loop["rho"] * v ** 2
    return dict(P=dp * (mdot_tot / loop["rho"]) / 0.35, dp=dp, v=v, Re=Re)

def stirrer_power(d, g, fl, u, T_oil=35.0) -> float:
    """Sealed circulator driving the whole free section at u through the
    cell bank: laminar bank friction + minor losses, 30% wire-to-fluid."""
    if u <= 1e-6:
        return 0.0
    p = film_props(fl, T_oil)
    dp = 32.0 * p["rho"] * p["nu"] * (2.2 * g["fill_h"]) * u / g["D_h"] ** 2 \
         + d.get("K_loop", 5.0) * 0.5 * p["rho"] * u ** 2
    return dp * (u * g["A_flow"]) / 0.30

def prop_power(d, g, fl, u, T_oil=35.0):
    """Bottom propeller pushing the oil UP through the row semi-channels
    (aligned with buoyancy). Head = laminar bank friction over the cell
    height on the array's hydraulic diameter, plus turning/grid losses;
    flow = u times the free riser area. Small shrouded axial impellers
    run ~30% wire-to-fluid, like the stirrer."""
    if u <= 1e-6:
        return dict(P=0.0, dp=0.0, Vdot_lpm=0.0)
    p = film_props(fl, T_oil)
    Re = u * g["D_h"] / p["nu"]
    f = 64.0 / max(Re, 1.0) if Re < 2300 else 0.316 * Re ** -0.25
    dp = ((f * d["h_cell"] / g["D_h"] + 2.0) * 0.5 * p["rho"] * u ** 2)
    Vdot = u * g["A_flow"]
    return dict(P=dp * Vdot / 0.30, dp=dp, Vdot_lpm=Vdot * 60000.0)


def plate_fin_area(d, g, h_oil):
    """Effective wetted area of thin plates hung from the tubes between cell
    rows (both sides wetted, conduction only). Fin: eta = tanh(mL)/mL,
    m = sqrt(2h/(k t)); base contact to the tube derated by plate_contact."""
    if not d.get("plate_on"):
        return 0.0, 0.0, 0.0
    k = 205.0 if d.get("plate_mat", "Aluminium") == "Aluminium" else 380.0
    t = d.get("plate_t", 0.0015)
    L = d["h_cell"]
    m = math.sqrt(2.0 * max(h_oil, 5.0) / (k * t))
    eta = math.tanh(m * L) / max(m * L, 1e-9)
    n_pl = max(g["n_rows"] - 1, 1)
    length = max(g["Lx"] - 2 * d["manifold_margin"], 0.1)
    served = min(1.0, d["n_tubes"] / n_pl)   # plates without a tube are
    A_eff = (n_pl * 2.0 * L * length * eta   # passive spreaders, not fins
             * d.get("plate_contact", 0.8) * served)
    m_pl = n_pl * L * length * t * (2700.0 if k < 300 else 8940.0)
    return A_eff, eta, m_pl

def serpentine_pump(d, g, fl, u, T_oil=35.0):
    """Parallel row-channels guided by the plates: laminar slot flow,
    Dp = 12 mu L u / s^2 per channel, all channels manifolded in parallel.
    Returns pump electrical power, oil dT along a channel, and flow."""
    if u <= 1e-6:
        return dict(P=0.0, dT_path=0.0, mdot=0.0, dp=0.0)
    p = film_props(fl, T_oil)
    s = max((d["pitch"] - d["d_cell"] - d.get("plate_t", 0.0015)) / 2, 5e-4)
    Lch = max(g["Lx"] - 2 * d["manifold_margin"], 0.1)
    n_ch = max(g["n_rows"] - 1, 1) * 2            # both sides of each plate
    A_ch = s * d["h_cell"]
    dp = 12.0 * p["rho"] * p["nu"] * Lch * u / s ** 2 + 3.0 * 0.5 * p["rho"] * u ** 2
    mdot = p["rho"] * u * A_ch * n_ch
    P = dp * (mdot / p["rho"]) / 0.35
    return dict(P=P, dT_path=0.0 if mdot < 1e-9 else 0.0, mdot=mdot, dp=dp,
                dT_est=lambda Q: Q / max(mdot * p["cp"], 1e-9))

def ext_loop_pump(d, g, fl, u, T_oil=35.0):
    """Closed dielectric loop driven by an EXTERNAL pump: oil leaves the
    pack, passes the pump, and returns - no external heat exchanger; the
    internal water tubes still remove the heat, so oil and water never
    meet. The pump sets the through-pack velocity u.

    Pressure drop = pack side + external pipework:
      pack, plates ON : laminar slot flow between plate and cells,
                        dp = 12 mu L u / s^2 per channel (parallel)
      pack, plates OFF: flow along the row semi-channels of the bare
                        array, f = 64/Re on the array's D_h over Lx
      pipes           : Darcy over pipe_len at bore pipe_id, plus
                        K_fit = 8 minor losses (bends, entries, volute)
    Electrical power = dp * Vdot / eta (0.35 wire-to-water default).
    """
    if u <= 1e-6:
        return dict(P=0.0, dp=0.0, dp_pack=0.0, dp_pipe=0.0, mdot=0.0,
                    v_pipe=0.0, Re_pipe=0.0, Vdot_lpm=0.0)
    p = film_props(fl, T_oil)
    if d.get("plate_on"):
        # plate-wall channels: defined slots either side of each plate
        s = max((d["pitch"] - d["d_cell"] - d.get("plate_t", 0.0015)) / 2,
                5e-4)
        Lch = max(g["Lx"] - 2 * d["manifold_margin"], 0.1)
        n_ch = max(g["n_rows"] - 1, 1) * 2
        A_pack = s * d["h_cell"] * n_ch
        dp_pack = (12.0 * p["rho"] * p["nu"] * Lch * u / s ** 2
                   + 3.0 * 0.5 * p["rho"] * u ** 2)
    else:
        # row semi-channels of the bare array: one channel per row gap,
        # slot width = clear space between rows, with a x2.5 friction
        # penalty for the flow meandering around the cell cylinders
        s = max(g["row_pitch"] - d["d_cell"], 1e-3)
        Lch = max(g["Lx"] - 2 * d["manifold_margin"], 0.1)
        n_ch = max(g["n_rows"] - 1, 1)
        A_pack = s * d["h_cell"] * n_ch
        dp_pack = (2.5 * 12.0 * p["rho"] * p["nu"] * Lch * u / s ** 2
                   + 3.0 * 0.5 * p["rho"] * u ** 2)
    mdot = p["rho"] * u * A_pack
    Vdot = mdot / p["rho"]
    d_p = max(d.get("pipe_id", 0.019), 3e-3)
    A_p = math.pi * d_p ** 2 / 4
    v_p = Vdot / A_p
    Re_p = v_p * d_p / p["nu"]
    f_p = 64.0 / max(Re_p, 1.0) if Re_p < 2300 else 0.316 * Re_p ** -0.25
    L_p = d.get("pipe_len", 2.5)
    dp_pipe = (f_p * L_p / d_p + 8.0) * 0.5 * p["rho"] * v_p ** 2
    dp = dp_pack + dp_pipe
    P = dp * Vdot / 0.35
    return dict(P=P, dp=dp, dp_pack=dp_pack, dp_pipe=dp_pipe, mdot=mdot,
                v_pipe=v_p, Re_pipe=Re_p, Vdot_lpm=Vdot * 60000.0)



def compare_architectures(d, g, fl, masses, T_amb, C_duty=None) -> pd.DataFrame:
    d = dict(d, C1=(C_duty if C_duty else d["C1"]))
    loop = WATER_LOOP[d["loop_fluid"]]
    Pw = water_pump_power(d, g, loop)["P"]
    imm_mass = masses["m_oil"] + masses["m_tubes"] + masses["m_fins"]
    rows = []

    def imm_case(label, u, extraP, mass, note):
        dd = dict(d); dd["u_oil"] = u
        r = solve_steady(dd, g, fl, 1.0, T_amb, C_rate=d["C1"])
        rows.append(dict(Architecture=label, T_cell=r["T_b"], T_core=r["T_core"],
                         Parasitic_W=Pw + extraP, Thermal_mass_kg=mass, Notes=note))

    imm_case("Static immersion + internal HX", 0.0, 0.0, imm_mass,
             "this design; thermosiphon only, no moving parts in oil")
    imm_case("Stirred immersion + internal HX (5 cm/s)", 0.05,
             stirrer_power(d, g, fl, 0.05), imm_mass + 0.5,
             "sealed magnetically-coupled circulator")

    # bottom cold plate (dry pack): axial cell path + TIM + channel film
    Q_cell = q_gen_per_cell(d, d["C1"], 35.0)
    Q_tot = Q_cell * g["N"]
    R_ax = (d["h_cell"] / 2) / (28.0 * math.pi * d["d_cell"] ** 2 / 4)
    R_cp = R_ax + 0.8 + 1.0 / (3000.0 * d["pitch"] ** 2)
    dTw = Q_tot / max((d["flow_lpm"] / 60 * loop["rho"] / 1000) * loop["cp"], 1e-9)
    T_cp = d["T_water_in"] + 0.5 * dTw + Q_cell * R_cp
    plate_mass = g["Lx"] * g["Ly"] * 0.006 * 2700 + 2.0
    rows.append(dict(Architecture="Bottom cold plate (dry pack)",
                     T_cell=T_cp, T_core=T_cp + Q_cell * r_core(d),
                     Parasitic_W=Pw, Thermal_mass_kg=plate_mass,
                     Notes="axial path ~3.6 K/W + TIM 0.8 + channel film; "
                           "1-2 mm pitch possible (no gap rule)"))

    # pumped dielectric + external HX (AMG HPB80 style)
    p35 = film_props(fl, 35.0)
    u_p = 0.20
    hf = nu_crossflow_cb(u_p * d["d_cell"] / p35["nu"], p35["Pr"]) * p35["k"] \
         / d["d_cell"] * d.get("cal_h", 1.0)
    mdot_oil = Q_tot / (fl["cp"] * 5.0)            # sized for 5 K oil rise
    T_oil_mean = d["T_water_in"] + 3.0 + 2.5       # HX approach + half rise
    T_pp = T_oil_mean + Q_cell / (hf * math.pi * d["d_cell"] * d["h_cell"])
    P_oil = (mdot_oil / fl["rho"]) * 30000.0 / 0.40
    pumped_mass = 0.025 * g["N"] * fl["rho"] / 1000 + 4.5   # ~25 mL/cell (HPB80 ratio) + HX/pump
    rows.append(dict(Architecture="Pumped dielectric + external HX (AMG-style)",
                     T_cell=T_pp, T_core=T_pp + Q_cell * r_core(d),
                     Parasitic_W=Pw + P_oil, Thermal_mass_kg=pumped_mass,
                     Notes=f"oil at {u_p*100:.0f} cm/s past cells, h ~ {hf:.0f}; "
                           "pump, plumbing, filter, de-aeration"))
    return pd.DataFrame(rows)

# ------------------------------------------------------------------ #
#  Thermal-runaway screening (order of magnitude)                     #
# ------------------------------------------------------------------ #
def runaway_screen(d, g, fl, masses) -> dict:
    p_, rp = d["pitch"], g["row_pitch"]
    r_zone = d["zone_pitches"] * p_
    plan = math.pi * r_zone ** 2
    n_in = max(plan / (p_ * rp) - 1.0, 0.0)
    oil_frac = max(1.0 - (math.pi * d["d_cell"] ** 2 / 4) / (p_ * rp), 0.05)
    V_zone = plan * d["h_cell"] * oil_frac + plan * d["tube_zone"]
    C_zone = V_zone * fl["rho"] * fl["cp"] + n_in * d["m_cell"] * d["cp_cell"]
    dT_zone = d["frac_oil"] * d["E_tr"] * 1000.0 / max(C_zone, 1.0)
    dT_bulk = d["E_tr"] * 1000.0 / (masses["C_oil"] + masses["C_batt"])
    V_hs = max(g["Lx"] * g["Ly"] * d["gas_gap"], 1e-5)
    P_final = 101325.0 * (V_hs + d["vent_L"] / 1000.0) / V_hs * (380.0 / 298.0)
    return dict(n_in=n_in, V_zone_L=V_zone * 1000, C_zone=C_zone,
                dT_zone=dT_zone, dT_bulk=dT_bulk,
                P_bar_g=(P_final - 101325.0) / 1e5, V_hs_L=V_hs * 1000)

# ------------------------------------------------------------------ #
#  Plan-view layout figure                                            #
# ------------------------------------------------------------------ #
def layout_figure(d, g, height=560):
    xs, ys, cs = [], [], []
    cnt = 0
    for r in range(g["n_rows"]):
        for c in range(g["n_cols"]):
            if cnt >= g["N"]:
                break
            x = d["edge_margin"] + (c + 0.5) * d["pitch"] \
                + (d["pitch"] / 2 if (d["arrangement"] == "Hexagonal" and r % 2) else 0)
            y = d["edge_margin"] + d["pitch"] / 2 + r * g["row_pitch"]
            xs.append(x); ys.append(y)
            cx, cy = g["Lx"] / 2, g["Ly"] / 2
            cs.append(1.0 - math.hypot(x - cx, y - cy) / math.hypot(cx, cy))
            cnt += 1
    size = max(3.0, d["d_cell"] / max(g["Lx"], 1e-3) * 640)
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, y0=0, x1=g["Lx"], y1=g["Ly"],
                  line=dict(color=INK, width=2), fillcolor="rgba(232,161,58,0.06)")
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", name="cells",
                             marker=dict(size=size, color=cs, colorscale="RdYlBu_r",
                                         showscale=False, line=dict(width=0)),
                             hoverinfo="skip"))
    for j in range(d["n_tubes"]):
        yj = (j + 0.5) * g["Ly"] / d["n_tubes"]
        fig.add_trace(go.Scatter(x=[d["manifold_margin"], g["Lx"] - d["manifold_margin"]],
                                 y=[yj, yj], mode="lines", showlegend=False,
                                 line=dict(color="#0EA5E9", width=3), hoverinfo="skip"))
    fig.update_yaxes(scaleanchor="x", scaleratio=1, visible=False)
    fig.update_xaxes(visible=False)
    fig.update_layout(height=height, margin=dict(l=10, r=10, t=40, b=10),
                      title=f"Plan view: {g['n_cols']} x {g['n_rows']} grid "
                            f"({d['arrangement'].lower()}), {d['n_tubes']} tube runs (blue) "
                            "in the zone above - cell colour hints centre-vs-edge tendency",
                      plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)")
    return fig

# ------------------------------------------------------------------ #
#  Calibration against measured data                                  #
# ------------------------------------------------------------------ #
def fit_calibration(d, g, fl, masses, T_amb, t_arr, C_arr, t_m, T_m):
    """Single multiplier on both oil-film h values, fitted to a measured
    cell-temperature trace - same pattern as the spray app's SS-1.0 factor."""
    best = (1.0, 1e9)
    for c in np.linspace(0.4, 2.2, 19):
        dd = dict(d); dd["cal_h"] = float(c)
        tr = solve_transient(dd, g, fl, masses, T_amb, t_arr, C_arr)
        pred = np.interp(t_m, tr["t"], tr["T_b"])
        rmse = float(np.sqrt(np.mean((pred - np.asarray(T_m)) ** 2)))
        if rmse < best[1]:
            best = (float(c), rmse)
    return best

# ------------------------------------------------------------------ #
#  Keyed widgets (enables save/load of whole designs)                 #
# ------------------------------------------------------------------ #
def _w(fn, label, key, default, **kw):
    k = f"w_{key}"
    if k not in st.session_state:
        st.session_state[k] = default
    return fn(label, key=k, **kw)


def arch_tab(d, g, fl, masses, C_steady):
    if True:
        st.markdown(f"Same pack, same {C_steady:.2f}C-rms duty, same water loop - four ways "
                    "to build the thermal system. This is the 'is it worth doing' slide.")
        adf = compare_architectures(d, g, fl, masses, d["T_amb"], C_steady)
        st.dataframe(adf.round(1), hide_index=True, use_container_width=True)
        cA, cB = st.columns(2)
        with cA:
            figA = go.Figure(go.Bar(x=adf["Architecture"], y=adf["T_cell"],
                                    marker_color=[ACCENT, "#D97706", "#6366F1", "#10B981"],
                                    text=[f"{v:.1f}" for v in adf["T_cell"]],
                                    textposition="outside"))
            figA.add_hline(y=d["T_limit"], line_dash="dash", line_color="#B91C1C")
            figA.update_layout(height=380, title="Steady cell temperature [°C]",
                               plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
                               xaxis_tickangle=-15)
            st.plotly_chart(figA, use_container_width=True)
        with cB:
            figP = go.Figure(go.Bar(x=adf["Architecture"], y=adf["Parasitic_W"],
                                    marker_color="#64748B",
                                    text=[f"{v:.0f} W" for v in adf["Parasitic_W"]],
                                    textposition="outside"))
            figP.update_layout(height=380, title="Parasitic power [W] (chiller excluded)",
                               plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
                               xaxis_tickangle=-15)
            st.plotly_chart(figP, use_container_width=True)
        st.caption("Cold-plate constants: axial cell path ~3.6 K/W, TIM 0.8 K/W, channel film "
                   "h = 3000 W/m²·K - edit in code if you have better numbers. Pumped case: "
                   "20 cm/s past cells, oil sized for 5 K rise, 30 kPa loop at 40% pump "
                   "efficiency, ~25 mL of fluid per cell (the HPB80 ratio). Kelvin per watt: "
                   "work the bottleneck, then buy the cheapest watts.")


def coolant_tab(d, g, cool_df):
    if True:
        st.markdown("Every fluid in the reviewed table, run through the **same pack at the "
                    f"same duty ({d['C1']:.1f}C continuous)**, thermosiphon and DCIR(T) "
                    "included. Non-dielectric fluids are thermal references only.")
        rows = []
        for _, r in cool_df.iterrows():
            f2 = fluid_dict(r)
            try:
                r2 = solve_steady(d, g, f2, 1.0, d["T_amb"], C_rate=d["C1"])
                m2 = build_masses(d, g, f2, r2["fin"])
                flags = []
                if not f2["dielectric"]:
                    flags.append("NOT dielectric")
                if not math.isnan(f2["bp"]) and r2["T_b"] + 10 > f2["bp"]:
                    flags.append("near boiling")
                if not math.isnan(f2["flash"]) and f2["flash"] < 120:
                    flags.append("low flash")
                if any(s in str(r["family"]).lower() for s in ("fluor", "hfo", "hydrofluoro")):
                    flags.append("PFAS")
                rows.append(dict(Fluid=f2["name"], Family=r["family"], nu_cSt=r["nu_cSt"],
                                 k=r["k"], u_ts_mms=r2["u_ts"] * 1000, h_cell=r2["h_cell"],
                                 T_cell=r2["T_b"], oil_kg=m2["m_oil"], Whkg=m2["whkg_pack"],
                                 Flags=", ".join(flags)))
            except Exception:
                pass
        sdf = pd.DataFrame(rows).sort_values("T_cell")
        st.dataframe(sdf.round(1), hide_index=True, use_container_width=True, height=420)
        figS = go.Figure()
        lab_all = thin_labels(sdf["oil_kg"], -sdf["T_cell"], list(sdf["Fluid"]),
                              min_dx=0.06, min_dy=0.10)   # favour the coolest
        lab_map = dict(zip(sdf["Fluid"], lab_all))
        for fam, grp in sdf.groupby("Family"):
            figS.add_trace(go.Scatter(x=grp["oil_kg"], y=grp["T_cell"], mode="markers+text",
                                      text=[lab_map[f] for f in grp["Fluid"]],
                                      hovertext=grp["Fluid"], hoverinfo="text+x+y",
                                      textposition="top center", name=fam,
                                      textfont=dict(size=11),
                                      marker=dict(size=9 + 40 * grp["k"] / sdf["k"].max())))
        figS.add_hline(y=d["T_limit"], line_dash="dash", line_color="#B91C1C")
        figS.update_layout(height=460, title="Cooler is down, lighter is left "
                           "(size ~ conductivity; hover for every fluid)",
                           xaxis_title="Coolant mass on board [kg]",
                           yaxis_title="Steady cell temperature [°C]",
                           plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)")
        st.plotly_chart(figS, use_container_width=True)
        st.caption("The fluorinated fluids win on h despite 5x lower conductivity (Ra ~ 1/nu, "
                   "and note their stronger thermosiphon in the u_ts column) but lose on "
                   "density, boiling point, cost and PFAS status (3M exited PFAS manufacture "
                   "end-2025). The esters are the pragmatic middle. Review notes from the "
                   "source spreadsheet apply.")


def bl_temp_profile_fig(delta_mm, k_oil, dT, T_bulk=35.0):
    """Teaching figure: the temperature across a single stagnant film.
    The whole point is that h is nothing more than k / delta - the fluid's
    conductivity divided by how thick the near-wall stuck layer is."""
    delta = max(delta_mm, 0.05) / 1000.0
    h = k_oil / delta
    q_flux = k_oil * dT / delta                      # W/m2
    T_s = T_bulk + dT
    xmax_mm = delta_mm * 3.2
    fig = go.Figure()
    fig.add_vrect(x0=0, x1=delta_mm, fillcolor="rgba(56,189,248,.18)",
                  line_width=0, layer="below")
    fig.add_trace(go.Scatter(
        x=[0, delta_mm, xmax_mm], y=[T_s, T_bulk, T_bulk], mode="lines",
        line=dict(color="#B91C1C", width=4), hoverinfo="skip"))
    fig.add_vline(x=0, line=dict(color="#334155", width=8))
    fig.add_annotation(x=delta_mm / 2, y=(T_s + T_bulk) / 2 + 0.02 * dT,
                       text=f"stagnant film<br>δ = {delta_mm:.2f} mm<br>"
                            "(conduction only)", showarrow=False,
                       font=dict(size=11, color="#0369A1"))
    fig.add_annotation(x=xmax_mm * 0.72, y=T_bulk + 0.06 * dT,
                       text="moving bulk<br>(well mixed)",
                       showarrow=False,
                       font=dict(size=11, color="#475569"))
    fig.add_annotation(x=0, y=T_s, text=f" wall {T_s:.0f} °C",
                       showarrow=False, xanchor="left",
                       font=dict(size=12, color="#B91C1C"))
    fig.update_layout(
        height=300, showlegend=False,
        title=f"h = k/δ = {k_oil:.2f} / {delta_mm:.2f} mm = "
              f"{h:.0f} W/m²·K    (heat flux {q_flux:.0f} W/m²)",
        xaxis_title="distance from the wall [mm]",
        yaxis_title="temperature [°C]",
        plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=10, r=10, t=48, b=10))
    return fig, h, q_flux


def h_scaling_fig(fl, dT, L_now_mm):
    """h = Nu k / L versus the characteristic length L, to make visible
    that h is a result of the geometry, not a fixed property: taller
    surfaces have a LOWER h (h ~ L^-1/4 in the laminar regime)."""
    pf = film_props(fl, 35.0)
    Ls = np.linspace(0.01, 0.20, 70)
    hs = [nu_vertical_cc(rayleigh(pf, dT, L), pf["Pr"]) * pf["k"] / L
          for L in Ls]
    L0 = max(L_now_mm, 10.0) / 1000.0
    h0 = nu_vertical_cc(rayleigh(pf, dT, L0), pf["Pr"]) * pf["k"] / L0
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Ls * 1000, y=hs, mode="lines",
                             line=dict(color=ACCENT, width=3)))
    fig.add_trace(go.Scatter(x=[L_now_mm], y=[h0], mode="markers",
                             marker=dict(size=13, color="#B91C1C")))
    fig.add_vline(x=L_now_mm, line_dash="dash",
                  annotation_text=f"L = {L_now_mm:.0f} mm")
    fig.update_layout(
        height=300, showlegend=False,
        title="Bigger is worse: h falls as the surface gets taller "
              "(h ~ L⁻¹ᐟ⁴)",
        xaxis_title="characteristic length L [mm]",
        yaxis_title=f"natural-convection h in {fl['name']} [W/m²·K]",
        plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=10, r=10, t=48, b=10))
    return fig, h0


def sandwich_profile_fig(res, Q, h_cell, h_tube, g, d):
    """The temperature climbing from the coolant to the cell core across
    every station, each rise sized by its thermal resistance q*R. The two
    oil films are the tall rungs; the tube wall is a whisker. Recomputes
    the two oil films from the given h so the sliders move the picture."""
    A_cell = g["A_cells"]
    A_tube = g.get("A_tube_bare", g["A_tube_in"])
    R_b = 1.0 / max(h_cell * A_cell, 1e-9)          # cell film
    R_ot = 1.0 / max(h_tube * A_tube, 1e-9)         # tube (oil) film
    R_wall, R_in = res["R_wall"], res["R_in"]
    dT_core, dT_water = res["dT_core"], res["dT_water"]
    Tw = d["T_water_in"] + 0.5 * dT_water           # mean water
    steps = [("water", 0.0, "#0EA5E9"),
             ("water film", Q * R_in, "#38BDF8"),
             ("tube wall", Q * R_wall, "#94A3B8"),
             ("oil film\n(tube side)", Q * R_ot, "#6366F1"),
             ("bulk oil", 0.0, "#38BDF8"),
             ("oil film\n(cell side)", Q * R_b, "#EF4444"),
             ("cell can", 0.0, "#F59E0B"),
             ("jellyroll core", dT_core, "#B91C1C")]
    xs, tops, bots, cols, labs, drops = [], [], [], [], [], []
    T = Tw
    for name, dT, col in steps:
        bots.append(T); T += dT; tops.append(T)
        xs.append(name.replace("\n", "<br>")); cols.append(col)
        labs.append(name); drops.append(dT)
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=xs, y=[t - b for t, b in zip(tops, bots)], base=bots,
        marker_color=cols, width=0.6,
        text=[f"+{dd:.1f}" if dd > 0.05 else "" for dd in drops],
        textposition="outside",
        hovertemplate="%{x}<br>ΔT = %{customdata:.2f} °C<extra></extra>",
        customdata=drops))
    fig.add_trace(go.Scatter(
        x=xs, y=tops, mode="lines+markers",
        line=dict(color="#0F172A", width=2, dash="dot"),
        marker=dict(size=6, color="#0F172A"), hoverinfo="skip"))
    fig.add_hline(y=d["T_limit"], line_dash="dash", line_color="#B91C1C",
                  annotation_text=f"limit {d['T_limit']:.0f} °C")
    fig.update_layout(
        height=380, showlegend=False,
        title=f"Every watt climbs {tops[-1] - bots[0]:.1f} °C from the "
              "water to the core - the two oil films are the tall rungs",
        yaxis_title="temperature [°C]",
        plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=10, r=10, t=48, b=10))
    return fig, dict(R_b=R_b, R_ot=R_ot, T_core=tops[-1], T_can=tops[-2])


def _sci(x, sig=1):
    """Format a number as clean LaTeX scientific notation, e.g.
    7.5e-4 -> '7.5\\times10^{-4}'. Plain formatting near unity."""
    import math as _m
    if x == 0:
        return "0"
    e = int(_m.floor(_m.log10(abs(x))))
    if -3 < e < 4:
        return f"{x:.6g}"
    m_ = x / (10 ** e)
    return f"{m_:.{sig}f}\\times10^{{{e}}}"


def learn_convection_html(k_oil, dT):
    """A self-contained animated boundary-layer explorer: warm oil rises
    past a hot cell wall, the no-slip film sets the temperature drop, and
    a flow slider thins the film and updates h = k/delta live. Renders in
    the browser (no Streamlit reruns), so the motion is smooth."""
    return """
<div id="cvroot" style="font-family:Inter,system-ui,sans-serif;color:#334155">
  <canvas id="cv" width="720" height="360"
          style="width:100%;max-width:720px;border-radius:10px;
                 background:#0b1220"></canvas>
  <div style="display:flex;gap:16px;align-items:center;margin-top:8px;
              flex-wrap:wrap">
    <label style="font-size:13px">Flow vigour
      <input id="sp" type="range" min="0" max="100" value="8"
             style="vertical-align:middle;width:200px">
    </label>
    <span id="ro" style="font-size:13px;font-weight:600;color:#0f172a"></span>
  </div>
  <div style="font-size:12px;color:#64748b;margin-top:4px">
    Drag the slider: more flow thins the stuck film &delta;, and since
    <b>h = k/&delta;</b>, the coefficient climbs &mdash; the oil and the
    temperature difference never change.
  </div>
</div>
<script>
(function(){
  const K = __KOIL__, DT = __DT__;
  const cv = document.getElementById('cv'), ctx = cv.getContext('2d');
  const W = cv.width, H = cv.height, wallX = 90;
  const sp = document.getElementById('sp'), ro = document.getElementById('ro');
  const N = 130, P = [];
  for(let i=0;i<N;i++) P.push({x: wallX + 6 + Math.random()*(W-wallX-16),
                               y: Math.random()*H, r: 1+Math.random()*1.6});
  function lerp(a,b,t){return a+(b-a)*t;}
  function tempColor(f){ // f=0 cold(blue) .. 1 hot(red)
    const r=Math.round(lerp(56,239,f)), g=Math.round(lerp(189,68,f)),
          b=Math.round(lerp(248,68,f)); return 'rgb('+r+','+g+','+b+')';
  }
  function frame(){
    const v = +sp.value/100;
    const delta = (2.0/(1+3.2*v));               // mm, 2.0 -> ~0.2
    const h = K/(delta/1000.0);
    const dpx = Math.max(6, delta*26);           // film width in px
    ctx.clearRect(0,0,W,H);
    // bulk temperature wash (thermal layer near wall)
    for(let xx=wallX; xx<W; xx+=6){
      const f = Math.max(0, 1-(xx-wallX)/dpx);
      ctx.fillStyle = 'rgba('+Math.round(lerp(20,239,f))+','+
        Math.round(lerp(30,68,f))+','+Math.round(lerp(60,68,f))+','+
        (0.10+0.55*f)+')';
      ctx.fillRect(xx,0,6,H);
    }
    // the hot wall (cell)
    const wg = ctx.createLinearGradient(wallX-40,0,wallX,0);
    wg.addColorStop(0,'#7f1d1d'); wg.addColorStop(1,'#ef4444');
    ctx.fillStyle = wg; ctx.fillRect(wallX-40,0,40,H);
    ctx.fillStyle='#fecaca'; ctx.font='12px Inter,sans-serif';
    ctx.save(); ctx.translate(wallX-26,H/2); ctx.rotate(-Math.PI/2);
    ctx.textAlign='center'; ctx.fillText('hot cell wall',0,0); ctx.restore();
    // film boundary line
    ctx.strokeStyle='rgba(56,189,248,.9)'; ctx.setLineDash([5,4]);
    ctx.beginPath(); ctx.moveTo(wallX+dpx,0); ctx.lineTo(wallX+dpx,H);
    ctx.stroke(); ctx.setLineDash([]);
    // rising particles: velocity ~0 at wall (no-slip), grows outward,
    // and the whole field rises faster with flow vigour
    for(const p of P){
      const dx = p.x - wallX;
      const prof = Math.min(1, dx/dpx);          // 0 at wall -> 1 at film edge
      const speed = (0.25 + 3.6*v) * (0.15 + prof);
      p.y -= speed;
      if(p.y < -4){ p.y = H+4; p.x = wallX+6+Math.random()*(W-wallX-16); }
      const f = Math.max(0, 1-dx/dpx);
      ctx.fillStyle = tempColor(f);
      ctx.globalAlpha = 0.5+0.5*prof;
      ctx.beginPath(); ctx.arc(p.x,p.y,p.r,0,6.2832); ctx.fill();
    }
    ctx.globalAlpha=1;
    // labels
    ctx.fillStyle='#93c5fd'; ctx.font='12px Inter,sans-serif';
    ctx.textAlign='left';
    ctx.fillText('stuck film  δ = '+delta.toFixed(2)+' mm',
                 wallX+6, H-14);
    ctx.fillStyle='#cbd5e1';
    ctx.fillText('warm oil rises →', W-150, 22);
    ro.innerHTML = 'δ = '+delta.toFixed(2)+' mm &nbsp;→&nbsp; '+
      'h = k/δ = '+K.toFixed(2)+'/'+delta.toFixed(2)+' mm = <b>'+
      h.toFixed(0)+' W/m²·K</b> &nbsp;|&nbsp; flux '+
      (h*DT).toFixed(0)+' W/m²';
    requestAnimationFrame(frame);
  }
  frame();
})();
</script>
""".replace("__KOIL__", f"{k_oil:.3f}").replace("__DT__", f"{dT:.1f}")



# ------------------------------------------------------------------ #
#  Cases tab - a verification ladder of small, hand-checkable cases  #
#  Every case calls the SAME functions the full model uses.          #
# ------------------------------------------------------------------ #
def case_still_bath(fl, Q, D, H, T_bath, f_ends=0.0, gap_mm=50.0):
    """One heated cylinder in a large still bath: bisect the surface
    temperature until Q = h(T_s) * A * (T_s - T_bath), with h from the
    model's own h_cell_side (Churchill-Chu vertical wall + gap factor,
    zero velocity)."""
    A = math.pi * D * H + f_ends * 2 * math.pi * D ** 2 / 4
    lo, hi = T_bath + 0.01, T_bath + 200.0
    film = None
    for _ in range(80):
        Ts = 0.5 * (lo + hi)
        film = h_cell_side(fl, Ts, T_bath, H, D, gap_mm, 0.0)
        if film["h"] * A * (Ts - T_bath) > Q:
            hi = Ts
        else:
            lo = Ts
    Ts = 0.5 * (lo + hi)
    film = h_cell_side(fl, Ts, T_bath, H, D, gap_mm, 0.0)
    resid = film["h"] * A * (Ts - T_bath) - Q
    return dict(T_s=Ts, h=film["h"], Ra=film["Ra"], Pr=film["Pr"],
                A=A, dT=Ts - T_bath, resid=resid,
                gapf=gap_factor(gap_mm))


def case_mode_temp(fl, Q, D, H, T_bath, u, mode, f_ends=0.0,
                   gap_mm=50.0):
    """Same cylinder with moving oil. mode: 'still' | 'cross' |
    'axial'. Uses h_cell_side exactly as the full solver does."""
    A = math.pi * D * H + f_ends * 2 * math.pi * D ** 2 / 4
    ax = (mode == "axial")
    uu = 0.0 if mode == "still" else u
    lo, hi = T_bath + 0.01, T_bath + 200.0
    for _ in range(80):
        Ts = 0.5 * (lo + hi)
        film = h_cell_side(fl, Ts, T_bath, H, D, gap_mm, uu, axial=ax)
        if film["h"] * A * (Ts - T_bath) > Q:
            hi = Ts
        else:
            lo = Ts
    Ts = 0.5 * (lo + hi)
    film = h_cell_side(fl, Ts, T_bath, H, D, gap_mm, uu, axial=ax)
    return dict(T_s=Ts, dT=Ts - T_bath, **{k: film[k] for k in
                ("h", "h_nat", "h_for", "Ra", "Re")}, A=A)


def case_bath_warmup(fl, N, Q, V_L, U_ext, A_ext, T_amb, t_end_min):
    """Sealed bath with no water sink: m cp dT/dt = N Q - U A (T-Tamb).
    Explicit Euler on the single-node balance the transient solver
    integrates; the steady state and time constant are closed-form."""
    m = V_L / 1000.0 * fl["rho"]
    C = m * fl["cp"]
    UA = U_ext * A_ext
    T_ss = T_amb + N * Q / max(UA, 1e-9)
    tau = C / max(UA, 1e-9)
    n = 400
    dt = t_end_min * 60.0 / n
    T = T_amb
    ts, Ts = [0.0], [T]
    for i in range(n):
        T += dt * (N * Q - UA * (T - T_amb)) / C
        ts.append((i + 1) * dt / 60.0)
        Ts.append(T)
    return dict(t_min=ts, T=Ts, T_ss=T_ss, tau_min=tau / 60.0,
                C_kJK=C / 1000.0, UA=UA, m=m)


def case_chain(fl, loop, Q, D, H, T_w_in, mdot_lpm, td, u=0.0,
               mode="still", fins=False, fin_geo=None, plate=None,
               f_ends=0.0):
    """One cell -> well-mixed oil -> ONE water tube. Solves the oil
    temperature from Q = (T_oil - T_w_mean) / (R_ot + R_wall + R_in)
    using h_tube_side, tube_section and h_water_inside - the exact
    chain inside solve_steady - then the can temperature from the
    cell-side film. Returns every resistance so the ladder can be
    summed by hand."""
    ts_ = tube_section(td)
    L_t = td["L_tube_case"]
    mdot = mdot_lpm / 60.0 * loop["rho"] / 1000.0
    T_w = T_w_in + Q / max(2.0 * mdot * loop["cp"], 1e-9)
    wat = h_water_inside(loop, mdot, ts_["D_h"], L_t,
                         lam_nu=ts_["lam_nu"], P_wet=ts_["P_in"])
    R_in = 1.0 / max(wat["h"] * ts_["P_in"] * L_t, 1e-12)
    if ts_["shape"] == "Round":
        R_wall = math.log(td["tube_od"] / ts_["D_h"]) / (
            2 * math.pi * K_TUBE[td["tube_mat"]] * L_t)
    else:
        P_m = 0.5 * (ts_["P_in"] + ts_["P_out"])
        R_wall = ts_["t"] / (K_TUBE[td["tube_mat"]] * P_m * L_t)
    lo, hi = T_w + 0.01, T_w + 150.0
    tube = None
    A_o = A_fin = A_pl = eta_f = eta_p = 0.0
    for _ in range(80):
        T_oil = 0.5 * (lo + hi)
        T_wall_est = T_w + Q * (R_in + R_wall)
        tube = h_tube_side(fl, T_oil, T_wall_est, td["tube_od"], u)
        A_o = ts_["P_out"] * L_t
        A_fin = eta_f = 0.0
        if fins and fin_geo and ts_["shape"] == "Round":
            fp = fin_pack(td["tube_od"], fin_geo["H"], fin_geo["t"],
                          fin_geo["p"], fin_geo["k"], tube["h"])
            A_fin = fp["A_eff_per_m"] * L_t - A_o
            eta_f = fp["eta"]
        A_pl = eta_p = 0.0
        if plate:
            A_pl, eta_p, _m = plate_fin_area(plate["d"], plate["g"],
                                             tube["h"])
        A_eff = A_o + max(A_fin, 0.0) + A_pl
        R_ot = 1.0 / max(tube["h"] * A_eff, 1e-12)
        if (T_oil - T_w) / (R_ot + R_wall + R_in) > Q:
            hi = T_oil
        else:
            lo = T_oil
    T_oil = 0.5 * (lo + hi)
    R_ot = 1.0 / max(tube["h"] * (A_o + max(A_fin, 0) + A_pl), 1e-12)
    cellc = case_mode_temp(fl, Q, D, H, T_oil, u, mode, f_ends)
    return dict(T_oil=T_oil, T_w_mean=T_w, T_s=cellc["T_s"],
                h_cell=cellc["h"], h_tube=tube["h"], h_water=wat["h"],
                regime=wat["regime"], Re_w=wat["Re"],
                R_ot=R_ot, R_wall=R_wall, R_in=R_in,
                A_o=A_o, A_fin=A_fin, A_pl=A_pl, eta_f=eta_f,
                eta_p=eta_p, ts=ts_, dT_w=Q / max(mdot * loop["cp"],
                                                  1e-9))


def learn_tab(d, g, fl, res, masses, cool_df, loop, Q_duty, chil):
    ACC = "#6366F1"
    st.markdown(
        f"The problem this pack sets is specific and unforgiving. It packs "
        f"21700 cells at high energy density and asks them to reject "
        f"several kilowatts of ohmic heat - **{Q_duty/1000:.2f} kW** right "
        f"now - into a dielectric oil that is a *poor* convector: its "
        f"conductivity is only $k \\approx {fl['k']:.2f}$ W/m·K, some "
        f"{400/fl['k']:.0f}× below copper, and it is viscous enough that, "
        f"left to itself, it barely moves. Every watt a cell makes must "
        f"cross **two nearly-stationary oil films** - one clinging to the "
        f"cell, one to the metal that carries the heat onward to water - "
        f"and those two films, not the choice of fluid and not the tube "
        f"material, set how hard this pack can be driven. Every feature of "
        f"the architecture (the guided serpentine plates, the water tubes, "
        f"any stirring) is really an assault on those two films.")
    st.markdown(
        "This tab derives that claim from the ground up, with each "
        "equation set **inside the sentence that explains it** and every "
        "symbol filled from your live design. The one idea it turns on: "
        "the heat-transfer coefficient $h$ is not a value you look up in a "
        "table - it is a **result** of the geometry, the fluid, and how "
        "hard the fluid is moving. That is why the levers that matter here "
        "are geometric and hydraulic, and the sliders let you feel how "
        "much each one actually buys.")

    # ============================================================= #
    st.divider()
    st.markdown("#### 1 · How heat actually moves, and why a thin film "
                "decides everything")
    st.markdown(
        f"Heat travels in only three ways, and here only two matter.\n\n"
        f"**Conduction** is heat diffusing through a material's own "
        f"structure - hot molecules jostling cooler neighbours. Fourier's "
        f"law fixes the rate: the flux (watts per m²) is "
        f"$q'' = -k\\,\\dfrac{{dT}}{{dx}}$, so across a slab of thickness "
        f"$\\Delta x$ and area $A$ the heat carried is "
        f"$q = kA\\,\\dfrac{{\\Delta T}}{{\\Delta x}}$. The conductivity "
        f"$k$ is what separates a good path from a bad one: copper carries "
        f"heat at $k\\approx 400$ W/m·K, but {fl['name'].split('(')[0].strip()} "
        f"manages only $k \\approx {fl['k']:.2f}$ - roughly "
        f"{400/fl['k']:.0f}× worse. That gap is the whole story of this "
        f"pack.\n\n"
        f"**Convection** is a moving fluid physically carrying heat away. "
        f"It beats conduction handily, *except* for one stubborn fact "
        f"that governs the entire design: right at any solid wall the "
        f"fluid is stuck. The no-slip condition means the oil touching a "
        f"cell cannot move, so inside that thin film heat can only "
        f"*conduct* across, at oil's feeble $k$. **That stagnant film is "
        f"the bottleneck**, and almost everything we do to cool the pack "
        f"is really about making it thinner.\n\n"
        f"**Radiation** ($q = \\varepsilon\\sigma A(T_s^4 - T_\\infty^4)$) "
        f"also moves heat, but at these gentle temperatures it is about a "
        f"watt across the whole pack, so we set it aside.")
    st.markdown(
        "Here is the move that ties convection to a single number. Across "
        "a film of thickness $\\delta$, Fourier's law gives "
        "$q = \\dfrac{k}{\\delta}\\,A\\,\\Delta T$. We bundle the awkward "
        "$k/\\delta$ into one symbol, the **heat-transfer coefficient** "
        "$h \\equiv k/\\delta$, and write **Newton's law of cooling**:")
    st.latex(r"q = h\,A\,\Delta T, \qquad h \equiv \frac{k}{\delta}.")
    st.markdown(
        "So $h$ is nothing more than the fluid's conductivity divided by "
        "how thick the stuck layer is. Stir the fluid and the film thins, "
        "$\\delta$ shrinks, and $h$ climbs - the fluid and the temperature "
        "difference never change. Drag the film thickness below and watch "
        "the temperature line steepen and $h$ rise.")
    st.markdown("**See it move.** Warm oil rises past the hot cell wall; "
                "right at the wall it is stuck (no-slip). More flow thins "
                "the stuck film, and $h=k/\\delta$ climbs:")
    components.html(learn_convection_html(fl["k"], 8.0), height=440)
    st.markdown("And here is the same physics as a temperature graph - "
                "the steeper the drop across the film, the higher $h$:")
    lc1, lc2 = st.columns([1, 2])
    with lc1:
        d_bl = st.slider("Film thickness δ [mm]", 0.2, 3.0, 1.5, 0.1,
                         key="ln_delta")
        dT_bl = st.slider("Wall-to-bulk ΔT [°C]", 2.0, 25.0, 8.0, 1.0,
                          key="ln_dtbl")
    figbl, h_bl, q_bl = bl_temp_profile_fig(d_bl, fl["k"], dT_bl)
    with lc2:
        st.plotly_chart(figbl, width='stretch', key="ln_blfig")
    thin_h = fl["k"] / (0.2e-3)
    st.caption(f"At δ = {d_bl:.1f} mm the film gives h = {h_bl:.0f} "
               f"W/m²·K; thin it to 0.2 mm and the same oil would give "
               f"~{thin_h:.0f} - a {thin_h/h_bl:.0f}× jump with nothing "
               f"but motion. This is why circulation matters more than "
               f"fluid choice.")

    # ============================================================= #
    st.divider()
    st.markdown("#### 2 · Why $h$ is a calculated result, not a number "
                "you can look up")
    st.markdown(
        "We cannot measure $\\delta$ directly and it depends on the flow, "
        "so engineers use a dimensionless ratio instead - the **Nusselt "
        "number** $Nu \\equiv \\dfrac{hL}{k}$, which reads as *how many "
        "times more heat this surface moves than if the same fluid across "
        "the same length $L$ just sat there conducting*. $Nu = 1$ is pure "
        "conduction; $Nu = 10$ means the moving fluid is ten times "
        "better. Rearranged, this **is** the recipe for $h$:")
    st.latex(r"h = \frac{Nu\,k}{L}.")
    st.markdown(
        "Everything now depends on finding $Nu$ - and $Nu$ depends on the "
        "geometry and the flow, which is exactly why $h$ can never be a "
        "fixed value. For a warm cell in still oil, the oil against it "
        "heats, expands, grows lighter and rises: **natural, "
        "buoyancy-driven convection**. How vigorous that rising flow is "
        "captured by the **Rayleigh number**,")
    st.latex(r"Ra_L = \frac{g\,\beta\,\Delta T\,L^{3}}{\nu\,\alpha}.")
    st.markdown(
        "Read it as a tug-of-war: gravity $g$ and thermal expansion "
        "$\\beta$ drive the buoyancy; the temperature difference "
        "$\\Delta T$ and the size $L$ set its scale; viscosity $\\nu$ and "
        "thermal diffusivity $\\alpha$ resist. Note the $L^{3}$ - a taller "
        "surface drives *much* stronger flow. The **Prandtl number** "
        "$Pr = \\nu/\\alpha$ (how fast momentum spreads versus heat) sets "
        "the layer shapes. The Churchill-Chu correlation, fitted to a "
        "century of experiments, turns $Ra$ and $Pr$ into $Nu$:")
    st.latex(r"Nu = \left(0.825 + \frac{0.387\,Ra^{1/6}}"
             r"{\left[1+(0.492/Pr)^{9/16}\right]^{8/27}}\right)^{2}.")
    st.markdown(
        "Put the chain together and you have the answer: **dimensions and "
        "$\\Delta T$ → $Ra$ → $Nu$ → $h$.** Set the three inputs and watch "
        "each number flow through the equations to the $h$ your surface "
        "actually gets.")

    hk1, hk2, hk3 = st.columns(3)
    L_hc = hk1.slider("Length L [mm]", 10.0, 200.0,
                      float(round(d["h_cell"] * 1000)), 5.0, key="ln_L")
    dT_hc = hk2.slider("ΔT [°C]", 2.0, 30.0, 8.0, 1.0, key="ln_dThc")
    fl_names = list(cool_df["name"])
    try:
        fidx = fl_names.index(fl["name"])
    except ValueError:
        fidx = 0
    fname = hk3.selectbox("Fluid", fl_names, index=fidx, key="ln_fluid")
    flc = fluid_dict(cool_df[cool_df["name"] == fname].iloc[0])
    L_m = L_hc / 1000.0
    pf = film_props(flc, 35.0)
    Ra_hc = rayleigh(pf, dT_hc, L_m)
    Nu_hc = nu_vertical_cc(Ra_hc, pf["Pr"])
    h_hc = Nu_hc * pf["k"] / L_m
    A_cell = math.pi * 0.021 * 0.07
    q_cell = h_hc * A_cell * dT_hc
    st.markdown(f"With $L = {L_hc:.0f}$ mm, $\\Delta T = {dT_hc:.0f}$ °C, "
                f"and **{fname.split('(')[0].strip()}** "
                f"($k = {pf['k']:.3f}$, $\\nu = {pf['nu']*1e6:.0f}$ cSt, "
                f"$\\alpha = {pf['alpha']*1e6:.3f}$ mm²/s, "
                f"$Pr = {pf['Pr']:.0f}$):")
    st.latex(r"Ra_L=\frac{g\,\beta\,\Delta T\,L^{3}}{\nu\,\alpha}"
             r"=\frac{9.81\times%s\times%d\times(%.3f)^{3}}"
             r"{%s\times%s}=%s"
             % (_sci(pf["beta"]), int(dT_hc), L_m, _sci(pf["nu"]),
                _sci(pf["alpha"]), _sci(Ra_hc, 2)))
    st.latex(r"Nu=\left(0.825+\frac{0.387\,Ra^{1/6}}{[\,\cdots\,]}"
             r"\right)^{2}=%.1f \qquad\Longrightarrow\qquad "
             r"h=\frac{Nu\,k}{L}=\frac{%.1f\times%.3f}{%.3f}"
             r"=\boxed{%.0f}\ \mathrm{W/m^2K}"
             % (Nu_hc, Nu_hc, pf["k"], L_m, h_hc))
    st.markdown(f"One 21700 cell of that surface would then shed "
                f"$q = hA\\Delta T = {h_hc:.0f}\\times{A_cell*1e3:.2f}"
                f"\\times10^{{-3}}\\times{dT_hc:.0f} = "
                f"\\mathbf{{{q_cell:.1f}}}$ **W**.")
    fighs, _ = h_scaling_fig(flc, dT_hc, L_hc)
    st.plotly_chart(fighs, width='stretch', key="ln_hscale")
    st.info(
        f"**Watch what just moved.** Doubling $L$ drops $h$ by about a "
        f"sixth (bigger is *worse*, because $Nu\\propto Ra^{{1/4}}\\propto "
        f"L^{{3/4}}$, so $h = Nu\\,k/L \\propto L^{{-1/4}}$). Switching to "
        f"a thinner, more conductive fluid raises it. Letting the surface "
        f"run hotter raises $h$ only as $\\Delta T^{{1/4}}$ - you cannot "
        f"cool a pack by letting it get hot. $h$ is a *consequence* of "
        f"these choices, never a dial you set directly.")

    # ============================================================= #
    st.divider()
    st.markdown("#### 3 · Following a single watt: where the temperature "
                "is actually lost")
    st.markdown(
        f"Now assemble the whole path a watt takes to escape, from the "
        f"cell core to the coolant. The pack is making **{Q_duty/1000:.2f} "
        f"kW** right now, and each step of the escape is a **thermal "
        f"resistance**. A convective film is $R = \\dfrac{{1}}{{hA}}$; a "
        f"plane wall is $R = \\dfrac{{\\Delta x}}{{kA}}$; a round tube "
        f"wall is $R = \\dfrac{{\\ln(d_o/d_i)}}{{2\\pi k L}}$. Resistances "
        f"in a chain simply add, and - exactly like a voltage divider - "
        f"the biggest one drops the most temperature, "
        f"$\\Delta T_i = q\\,R_i$.\n\n"
        f"The chain, water to core: the **water film**, the **tube wall** "
        f"(a whisker - copper or aluminium barely matters), the **oil "
        f"film on the tube**, the **bulk oil** (which carries but hardly "
        f"resists), the **oil film on the cell**, then conduction through "
        f"the jellyroll to the **core**. The two oil films are the tall "
        f"rungs. Move the slider - stirring thins both oil films at once - "
        f"and watch the tall rungs shrink and the core drop toward the "
        f"limit.")
    u_ln = st.slider("Oil circulation speed [cm/s] (0 = still, buoyancy "
                     "only)", 0.0, 15.0, 0.0, 0.5, key="ln_stir")
    u_eff = max(u_ln / 100.0, res["u_ts"])
    hc_ln = h_cell_side(fl, res["T_b"], res["T_il"], d["h_cell"],
                        d["d_cell"], g["gap_mm"], u_eff)["h"] * d["cal_h"]
    ht_ln = h_tube_side(fl, res["T_il"], res["T_wall"], d["tube_od"],
                        u_eff)["h"] * d["cal_h"]
    figsw, swinfo = sandwich_profile_fig(res, res["Q_eff"], hc_ln, ht_ln,
                                         g, d)
    st.plotly_chart(figsw, width='stretch', key="ln_sandwich")
    st.caption(
        f"At {u_ln:.1f} cm/s the cell film is h = {hc_ln:.0f} and the tube "
        f"film h = {ht_ln:.0f} W/m²·K, and the core sits at "
        f"{swinfo['T_core']:.1f} °C. So what actually matters, in order: "
        f"(1) the two oil films - move the oil and buy sink area; (2) "
        f"getting the water turbulent (next section); (3) how hard you "
        f"push, since heat grows with the *square* of C-rate; (4) the "
        f"set-point trade. What barely matters: tube material, and fluid "
        f"brand beyond its viscosity class.")

    # ============================================================= #
    st.divider()
    st.markdown("#### 4 · The four levers on $h$ and the temperature")
    st.markdown("##### Lever A - viscosity (weakly, through $Ra$)")
    st.markdown(
        "Since $Ra \\propto 1/\\nu$ and $h \\propto Ra^{1/6\\text{–}1/4}$, "
        "a thinner oil gives a higher $h$ - but only as roughly "
        "$h\\propto\\nu^{-1/4}$. The plot places every candidate fluid: "
        "viscosity (horizontal) matters far more than conductivity "
        "(colour and size), yet the whole $h$ axis spans barely 3×. "
        "**Fluid choice cannot buy fast charge; geometry and stirring "
        "can.**")
    xs, ys, names, ks = [], [], [], []
    for _, r in cool_df.iterrows():
        f3 = fluid_dict(r)
        p3 = film_props(f3, 35.0)
        Ra3 = rayleigh(p3, 8.0, d["h_cell"])
        ys.append(nu_vertical_cc(Ra3, p3["Pr"]) * p3["k"] / d["h_cell"])
        xs.append(r["nu_cSt"]); names.append(r["name"]); ks.append(r["k"])
    figV = go.Figure(go.Scatter(
        x=xs, y=ys, mode="markers+text",
        text=thin_labels(xs, ys, names, logx=True),
        hovertext=names, hoverinfo="text+x+y",
        textposition="top center", textfont=dict(size=11),
        marker=dict(size=8 + 60 * np.array(ks) / max(ks), color=ks,
                    colorscale="YlOrBr", colorbar=dict(title="k [W/mK]"))))
    figV.update_layout(height=380, xaxis_type="log",
                       xaxis_title="kinematic viscosity at 25 °C [cSt] (log)",
                       yaxis_title="natural-convection h on a 21700 [W/m²·K]",
                       title="Viscosity is the strong axis; conductivity "
                             "the weak one",
                       plot_bgcolor="rgba(255,255,255,0)",
                       paper_bgcolor="rgba(0,0,0,0)")
    st.plotly_chart(figV, width='stretch', key="ln_visc")

    st.markdown("##### Lever B - the cell gap (throttling the buoyancy)")
    st.markdown(
        "Squeeze the cells closer and you throttle the buoyant flow in "
        "the gaps. The thin-channel buoyant velocity scales as "
        "$u_\\mathrm{gap} \\sim \\dfrac{g\\beta\\Delta T\\,\\delta^{2}}"
        "{\\nu}$, so $u_\\mathrm{gap}\\propto\\delta^{2}$ - halve the gap "
        "and you quarter the flow. Below about 6 mm the penalty bites "
        "hard (Wang et al. measured gap velocity falling 1.8 → 0.5 mm/s "
        "as spacing shrank 8 → 2 mm). This is the energy-density tax of "
        "*static* immersion; stirring largely removes it.")
    gaps = np.linspace(0.5, 10, 60)
    figG = go.Figure(go.Scatter(x=gaps, y=[gap_factor(x) for x in gaps],
                                line=dict(color=ACC, width=3)))
    figG.add_vline(x=g["gap_mm"], line_dash="dash",
                   annotation_text=f"your gap {g['gap_mm']:.1f} mm")
    figG.update_layout(height=280, xaxis_title="cell-to-cell gap [mm]",
                       yaxis_title="h penalty factor",
                       title="Below ~6 mm the buoyant flow is throttled",
                       plot_bgcolor="rgba(255,255,255,0)",
                       paper_bgcolor="rgba(0,0,0,0)")
    st.plotly_chart(figG, width='stretch', key="ln_gap")

    st.markdown("##### Lever C - stirring (forced convection)")
    st.markdown(
        f"A little motion beats none by a lot. Left alone the pack even "
        f"stirs itself: the warm column is lighter than the cold return, "
        f"and that buoyant head $\\rho\\beta g H\\,\\Delta T$ drives a slow "
        f"thermosiphon loop against friction "
        f"$\\left(K + f\\frac{{L}}{{D_h}}\\right)\\frac{{1}}{{2}}\\rho "
        f"u^{{2}}$ until the two balance - here about "
        f"{res['u_ts']*1000:.1f} mm/s. Force it faster with a small pump "
        f"or guided plates and both films thin together.")
    us = np.linspace(0, 0.2, 50)
    hcs, hts = [], []
    for u in us:
        ue = max(u, res["u_ts"])
        hcs.append(h_cell_side(fl, res["T_b"], res["T_il"], d["h_cell"],
                               d["d_cell"], g["gap_mm"], ue)["h"]
                   * d["cal_h"])
        hts.append(h_tube_side(fl, res["T_il"], res["T_wall"],
                               d["tube_od"], ue)["h"] * d["cal_h"])
    figU = go.Figure()
    figU.add_trace(go.Scatter(x=us * 100, y=hcs, name="cell film",
                              line=dict(color="#EF4444", width=3)))
    figU.add_trace(go.Scatter(x=us * 100, y=hts, name="tube film",
                              line=dict(color="#6366F1", width=3)))
    figU.add_vline(x=res["u_ts"] * 100, line_dash="dot",
                   line_color="#10B981",
                   annotation_text=f"thermosiphon {res['u_ts']*1000:.1f} mm/s")
    figU.update_layout(height=280, xaxis_title="oil velocity [cm/s]",
                       yaxis_title="h [W/m²·K]",
                       title="The pack stirs itself a little; a circulator "
                             "does it properly",
                       plot_bgcolor="rgba(255,255,255,0)",
                       paper_bgcolor="rgba(0,0,0,0)")
    st.plotly_chart(figU, width='stretch', key="ln_stirfig")

    st.markdown("##### Lever D - fins (buying area $A$ when $h$ is stuck)")
    st.markdown(
        "When $h$ is stuck low, buy area $A$ instead. A fin works only if "
        "heat can run out along it before the film pulls it off: the fin "
        "parameter $m = \\sqrt{2h/(k\\,t)}$ sets the reach, and the "
        "efficiency $\\eta = \\tanh(mL)/(mL)$ says how much of the fin is "
        "pulling its weight. Oil's low $h$ is the fin-maker's friend - "
        "even long, thin fins stay ~90% efficient, so the effective area "
        "$A_\\mathrm{eff} = A_\\mathrm{base} + \\eta\\,A_\\mathrm{fin}$ "
        "grows almost for free.")
    hs = np.linspace(0.002, 0.02, 40)
    eff, gain = [], []
    for hh in hs:
        fp = fin_pack(d["tube_od"], hh, d["fin_t"] if d["fins_on"] else 0.0005,
                      d["fin_p"] if d["fins_on"] else 0.004, 205.0,
                      max(res["h_tube"], 30))
        eff.append(fp["eta"]); gain.append(fp["area_gain"])
    figF = go.Figure()
    figF.add_trace(go.Scatter(x=hs * 1000, y=gain, name="area gain ×",
                              line=dict(color=INK, width=3)))
    figF.add_trace(go.Scatter(x=hs * 1000, y=eff, name="fin efficiency",
                              yaxis="y2", line=dict(color=ACC, width=3)))
    figF.update_layout(height=280, xaxis_title="fin height [mm]",
                       yaxis_title="effective area multiplier",
                       yaxis2=dict(title="Schmidt efficiency",
                                   overlaying="y", side="right",
                                   range=[0, 1.05]),
                       title="Long thin fins stay efficient in oil",
                       plot_bgcolor="rgba(255,255,255,0)",
                       paper_bgcolor="rgba(0,0,0,0)")
    st.plotly_chart(figF, width='stretch', key="ln_fins")

    # ============================================================= #
    st.divider()
    st.markdown("#### 5 · The water side: the laminar plateau")
    st.markdown(
        "Inside the tubes it is forced flow, and the regime is "
        "everything. The Reynolds number $Re = \\dfrac{4\\dot m}{\\pi\\mu "
        "d_i}$ decides: below about 2300 the flow is laminar and, "
        "remarkably, $Nu$ is a *constant* (3.66 for fixed wall "
        "temperature, 4.36 for fixed flux) - so **pumping harder does "
        "nothing for $h$**. Cross into turbulence and the Gnielinski "
        "correlation takes over, "
        "$Nu = \\dfrac{(f/8)(Re-1000)Pr}{1+12.7\\sqrt{f/8}\\,"
        "(Pr^{2/3}-1)}$, and $h$ climbs steeply. The lesson: get the "
        "water turbulent, and *then* flow buys cooling; before that, it "
        "does not.")
    fls = np.linspace(0.5, 60, 80)
    hws, res_w = [], []
    for q in fls:
        md = q / 60 * loop["rho"] / 1000 / max(d["n_tubes"], 1)
        w = h_water_inside(loop, md, g["d_i"], g["L_tube"])
        hws.append(w["h"]); res_w.append(w["Re"])
    figW = go.Figure(go.Scatter(x=fls, y=hws, line=dict(color="#0EA5E9",
                                                        width=3)))
    figW.add_vline(x=d["flow_lpm"], line_dash="dash",
                   annotation_text="your flow")
    i2300 = int(np.argmin(np.abs(np.array(res_w) - 2300)))
    figW.add_vline(x=fls[i2300], line_dash="dot", line_color="#B91C1C",
                   annotation_text="Re 2300")
    figW.update_layout(height=280, xaxis_title="total water flow [L/min]",
                       yaxis_title="h inside tube [W/m²·K]",
                       title="Laminar: flat. Turbulent: it climbs.",
                       plot_bgcolor="rgba(255,255,255,0)",
                       paper_bgcolor="rgba(0,0,0,0)")
    st.plotly_chart(figW, width='stretch', key="ln_water")

    # ============================================================= #
    st.divider()
    st.markdown("#### 6 · Heat that fights back, and the core the coolant "
                "cannot reach")
    st.markdown(
        f"The cell's own resistance falls as it warms, "
        f"$R_\\mathrm{{dc}}(T) = R_{{25}}\\,e^{{-k_\\mathrm{{dc}}(T-25)}}$, "
        f"so a hotter cell makes *less* heat at the same current - running "
        f"at {res['T_b']:.0f} °C instead of 25 cuts generation by about "
        f"**{100*(1-r_of_T(d, res['T_b'])/d['r_dc']):.0f}%**. That is why "
        f"a warm set point is a rare double win and why the transient "
        f"self-stabilises near the limit. But one resistance no coolant "
        f"can touch is the spread from the jellyroll core to the can, "
        f"$R_\\mathrm{{core}} = \\dfrac{{1}}{{4\\pi k_r H}} = "
        f"{r_core(d):.2f}$ K/W - pure conduction through the winding, "
        f"worth **{res['dT_core']:.1f} °C** here. That is why the tool "
        f"reports the **core** temperature separately: it is the node "
        f"closest to the limit.")
    Ts = np.linspace(0, 60, 61)
    figD = go.Figure(go.Scatter(x=Ts, y=[r_of_T(d, t) for t in Ts],
                                line=dict(color=INK, width=3)))
    figD.add_vline(x=res["T_b"], line_dash="dash",
                   annotation_text="your cell")
    figD.update_layout(height=280, xaxis_title="cell temperature [°C]",
                       yaxis_title="DCIR [mΩ]",
                       title=f"R(T) = R₂₅·exp(−{d['k_dcir']*100:.1f}%/K × "
                             f"(T−25))",
                       plot_bgcolor="rgba(255,255,255,0)",
                       paper_bgcolor="rgba(0,0,0,0)")
    st.plotly_chart(figD, width='stretch', key="ln_dcir")

    # ============================================================= #
    st.divider()
    st.markdown("#### 7 · Buffering: the oil is a thermal flywheel")
    st.markdown(
        "The oil and cells together are a thermal flywheel of capacity "
        "$C = \\sum_i m_i c_{p,i}$. Faced with a burst of excess heat "
        "$Q_\\mathrm{excess}$ they can absorb, they buy time "
        "$\\tau = \\dfrac{C\\,\\Delta T_\\mathrm{allow}}{Q_\\mathrm{excess}}$ "
        "before the temperature drifts up. Size the steady heat exchanger "
        "for the *continuous* duty and let the flywheel eat the peaks.")
    cQ, cT = st.columns(2)
    Q_ex = cQ.slider("Excess heat beyond removal [kW]", 0.1, 15.0, 3.0,
                     0.1, key="ln_qex")
    dT_h = cT.slider("Allowed temperature drift [°C]", 2.0, 25.0, 10.0,
                     1.0, key="ln_drift")
    Ctot = masses["C_oil"] + masses["C_batt"]
    st.markdown(
        f"Thermal mass = oil {masses['C_oil']/1000:.0f} kJ/K + cells "
        f"{masses['C_batt']/1000:.0f} kJ/K = **{Ctot/1000:.0f} kJ/K**, so "
        f"$\\tau = \\dfrac{{{Ctot/1000:.0f}\\times{dT_h:.0f}}}"
        f"{{{Q_ex:.1f}\\times10^{{3}}}} = "
        f"\\mathbf{{{Ctot*dT_h/(Q_ex*1000)/60:.1f}}}$ **minutes** per "
        f"{dT_h:.0f} °C of drift.")

    # ============================================================= #
    st.divider()
    with st.expander("Reference: moving the oil, safety, the production "
                     "benchmark, and sources"):
        st.markdown(f"""
**How you would actually move the oil**

| Option | How | Power | Reaches tight gaps? | Notes |
|---|---|---|---|---|
| Thermosiphon only | buoyancy loop | 0 W | weakly ({res['u_ts']*1000:.1f} mm/s here) | free; dies if the cold plane is low |
| Magnetic stirrer | sealed impeller in the bulk | ~1-3 W | no - bypasses the gaps | cheapest forced option |
| Pump + jet manifold | nozzles along a wall | ~2-5 W | partially | directional; nozzle fouling |
| **Serpentine plates (guided)** | plates form parallel channels; small pump | **~0.5 W at 5 cm/s** | **yes - every gap** | plates double as fins; needs a manifold |
| **Bottom propeller (axial, up)** | shrouded impeller under the array pushes the oil straight up, with buoyancy | ~0.5-3 W | yes - every riser gap | axial flow gives a weaker film than a crossflow sweep at the same speed; no plate area |
| **External pump loop (closed)** | oil out to a pump and straight back - no external HX; the internal water tubes still reject the heat | ~3-30 W (the pipes dominate) | yes - drives the same channels | pump serviceable without opening the pack; oil and water never meet |
| Full pumped immersion | external HX loop | 20+ W | yes | the AMG HPB80 architecture |

**Safety and practicalities**
- **Water-in-oil leak** is the single-point failure: hold oil pressure above water pressure so leaks go oil-to-water, or use double-walled tubes with leak detection (transformer practice).
- **Expansion**: $\\beta = {fl['beta']:.1e}$ /K on {masses['V_oil_L']:.0f} L means ~{fl['beta']*masses['V_oil_L']*70:.1f} L over a −10 to 60 °C band. Use a bellows or bladder, not free air.
- **Materials**: seal and insulation compatibility, ester moisture uptake, copper-oxidation catalysis (use inhibited fluids or plated tubes).
- **Venting**: a cell venting into a sealed flooded box is a pressure spike (see the Safety tab); the flip side is that oxygen exclusion and the oil's heat absorption suppress propagation.
- **Fluid supply**: 3M exited PFAS manufacture end-2025; anchor the programme on esters.

**Production reference - Mercedes-AMG HPB80** (batterydesign.net): 560 × 21700, 112S5P, 6.1 kWh, 89 kg, 68.5 Wh/kg, 150 kW peak / 70 kW continuous, pumped dielectric (14 L) through an external dielectric-to-water HX, 10 kW cooling, 45 °C set point.

**Sources**: Wang et al. 2023 (J. Energy Storage 62, 106821); Zou et al. 2024 (J. Energy Storage 83, 110634); Roe et al. 2022 (J. Power Sources 525, 231094); batterydesign.net.

*A note on units: temperature **differences** are written in °C here; a difference of 1 K and 1 °C are identical in size, only the zero points of the two scales differ.*""")




# ------------------------------------------------------------------ #
#  FEA export - operating-point builders (shared by tab and smoke)   #
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
#  WP3 basic-module analytics: parse a COMSOL field export, build    #
#  the analytical series solution, and judge the run.                #
# ------------------------------------------------------------------ #
def parse_comsol_field(text):
    """COMSOL spreadsheet Data export -> (x, height, T_2d, info).
    3 numeric columns = 2D (x, y, T). 4 columns = 3D (x, y_depth, z,
    T): the field is sliced at mid-depth for the 2D checker, and
    info["zvar"] reports the worst variation ACROSS the depth at any
    (x, z) - for the z-invariant WP3 case this must be solver noise.
    Regular grids are pivoted exactly; scattered exports fall back to
    bin-averaging."""
    rows = []
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("%"):
            continue
        parts = ln.replace(",", " ").split()
        vals = []
        for p_ in parts:
            try:
                vals.append(float(p_))
            except ValueError:
                break
        if len(vals) >= 3:
            rows.append(vals[:4])
    if len(rows) < 50:
        raise ValueError("not a recognisable COMSOL field export "
                         f"({len(rows)} data rows found)")
    ncol = min(len(r) for r in rows)
    info = dict(kind="2d", zvar=None, depth_n=1)
    if ncol >= 4:
        A = np.asarray([r[:4] for r in rows])
        xa, ya, za, Ta = A[:, 0], A[:, 1], A[:, 2], A[:, 3]
        xu = np.unique(np.round(xa, 9))
        yu = np.unique(np.round(ya, 9))
        zu = np.unique(np.round(za, 9))
        info.update(kind="3d", depth_n=len(yu),
                    depth_span=float(ya.max() - ya.min()))
        if len(xu) * len(yu) * len(zu) == len(Ta):
            ix = np.searchsorted(xu, np.round(xa, 9))
            iy = np.searchsorted(yu, np.round(ya, 9))
            iz = np.searchsorted(zu, np.round(za, 9))
            G = np.full((len(zu), len(yu), len(xu)), np.nan)
            G[iz, iy, ix] = Ta
            info["zvar"] = float(np.nanmax(
                np.nanmax(G, axis=1) - np.nanmin(G, axis=1)))
            mid = len(yu) // 2
            return xu, zu, G[:, mid, :], info
        mid_y = np.median(ya)
        keep = np.abs(ya - mid_y) <= (ya.max() - ya.min()) / (
            2 * max(len(yu), 2))
        x = xa[keep]; y = za[keep]; T = Ta[keep]
    else:
        A = np.asarray([r[:3] for r in rows])
        x, y, T = A[:, 0], A[:, 1], A[:, 2]
    xu = np.unique(np.round(x, 9)); yu = np.unique(np.round(y, 9))
    if len(xu) * len(yu) == len(T):
        ix = np.searchsorted(xu, np.round(x, 9))
        iy = np.searchsorted(yu, np.round(y, 9))
        G = np.full((len(yu), len(xu)), np.nan)
        G[iy, ix] = T
        if not np.isnan(G).any():
            return xu, yu, G, info
    nx, ny = 61, 73
    xe = np.linspace(x.min(), x.max(), nx + 1)
    ye = np.linspace(y.min(), y.max(), ny + 1)
    s = np.zeros((ny, nx)); c = np.zeros((ny, nx))
    ix = np.clip(np.searchsorted(xe, x) - 1, 0, nx - 1)
    iy = np.clip(np.searchsorted(ye, y) - 1, 0, ny - 1)
    np.add.at(s, (iy, ix), T); np.add.at(c, (iy, ix), 1.0)
    c[c == 0] = np.nan
    return (0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]),
            s / c, info)


def wp3_series_field(x, y, W, H, bw, bh, xc, gap, q_v, k, T_top,
                     nmax=90):
    """Analytical (eigenfunction-series) steady field for the basic
    module with UNIFORM conductivity k: Poisson source q_v on the bar
    footprint, adiabatic sides and bottom, T_top on the lid. Modes
    cos(n pi x / W); the y-solve uses the mode Green's function
    G = cosh(l*y_min) sinh(l*(H-y_max)) / (l cosh lH), integrated in
    closed form over the source band. Mode 0 reproduces the exact
    mean profile. Note the caveat: the real aluminium bar is nearly
    isothermal; with uniform k the series overshoots INSIDE the bar,
    so field comparisons mask the bar footprint and the bar estimate
    is the series average over it."""
    x = np.asarray(x); y = np.asarray(y)
    x0, x1 = xc - bw / 2, xc + bw / 2
    y0, y1 = gap, gap + bh
    X, Y = np.meshgrid(x, y)
    c0 = q_v * bw / W
    I = np.where(Y >= y1, bh * (H - Y),
        np.where(Y >= y0,
                 (bh ** 2 - (Y - y0) ** 2) / 2 + bh * (H - y1),
                 bh ** 2 / 2 + bh * (H - y1)))
    T = T_top + (c0 / k) * I
    nmax = min(nmax, max(int(300.0 * W / (math.pi * H)), 8))
    for n in range(1, nmax + 1):
        lam = n * math.pi / W
        cn = 2 * q_v * (math.sin(lam * x1) - math.sin(lam * x0)) / (
            W * lam)
        if abs(cn) < 1e-30:
            continue
        chH = math.cosh(lam * H)
        pref = cn / (k * lam ** 2 * chH)
        fb = np.cosh(lam * Y) * (math.cosh(lam * (H - y0)) -
                                 math.cosh(lam * (H - y1)))
        fa = np.sinh(lam * (H - Y)) * (math.sinh(lam * y1) -
                                       math.sinh(lam * y0))
        fm = (np.sinh(lam * (H - Y)) *
              (np.sinh(lam * Y) - math.sinh(lam * y0)) +
              np.cosh(lam * Y) *
              (np.cosh(lam * (H - Y)) - math.cosh(lam * (H - y1))))
        f = np.where(Y <= y0, fb, np.where(Y >= y1, fa, fm)) * pref
        T = T + f * np.cos(lam * X)
    return T


def wp3_exact_mean(y, W, H, bw, bh, gap, q_v, k, T_top):
    """EXACT horizontal-plane mean temperature (energy integral, no
    approximation): flat below the bar, quadratic across the band,
    linear with slope Q'/(kW) above it, T_top at the lid."""
    y = np.asarray(y)
    c0 = q_v * bw / W
    y0, y1 = gap, gap + bh
    I = np.where(y >= y1, bh * (H - y),
        np.where(y >= y0,
                 (bh ** 2 - (y - y0) ** 2) / 2 + bh * (H - y1),
                 bh ** 2 / 2 + bh * (H - y1)))
    return T_top + (c0 / k) * I


def wp3_check_field(x, y, T, W, H, bw, bh, xc, gap, q_v, k, T_top):
    """Judge an uploaded field: rigorous integral anchors first, then
    agreement with the analytical series outside the bar."""
    order = np.argsort(y)
    y = np.asarray(y)[order]; T = np.asarray(T)[order, :]
    xo = np.argsort(x); x = np.asarray(x)[xo]; T = T[:, xo]
    rm = np.nanmean(T, axis=1)
    y0, y1 = gap, gap + bh
    m = dict()
    m["top_dev"] = float(abs(rm[-1] - T_top))
    below = y < y0 - 0.02 * H
    m["below_flat"] = float(rm[below].max() - rm[below].min()) \
        if below.sum() > 2 else 0.0
    above = y > y1 + 0.05 * H
    slope_exact = -q_v * bw * bh / (k * W)
    if above.sum() > 3:
        A_ = np.vstack([y[above], np.ones(above.sum())]).T
        sl = float(np.linalg.lstsq(A_, rm[above], rcond=None)[0][0])
        m["slope_meas"] = sl
        m["slope_rel"] = float(abs(sl - slope_exact) /
                               abs(slope_exact))
    else:
        m["slope_rel"] = 0.0
    Ts = wp3_series_field(x, y, W, H, bw, bh, xc, gap, q_v, k, T_top)
    X, Y = np.meshgrid(x, y)
    bar = ((X > xc - bw / 2 - 0.02 * W) & (X < xc + bw / 2 + 0.02 * W)
           & (Y > y0 - 0.02 * H) & (Y < y1 + 0.02 * H))
    d = (T - Ts)[~bar & ~np.isnan(T)]
    m["rms_out"] = float(np.sqrt(np.mean(d ** 2)))
    m["max_out"] = float(np.max(np.abs(d)))
    m["T_peak"] = float(np.nanmax(T))
    inbar = ((X >= xc - bw / 2) & (X <= xc + bw / 2) & (Y >= y0)
             & (Y <= y1))
    m["T_bar_series"] = float(np.mean(Ts[inbar])) if inbar.any() \
        else float("nan")
    m["anchors_ok"] = (m["top_dev"] < 0.05 and
                       m["below_flat"] < 0.08 and
                       m["slope_rel"] < 0.05)
    m["Ts"] = Ts; m["mean_num"] = rm
    m["mean_exact"] = wp3_exact_mean(y, W, H, bw, bh, gap, q_v, k, T_top)
    m["x"] = x; m["y"] = y; m["T"] = T
    return m



def fan_mean_exact(y, W, H, bw, bh, gap, q_v, k, rho, cp, u, T_in):
    """EXACT plane-mean profile for the fan rung: upward speed u,
    inlet Dirichlet T_in at the floor, outflow at the lid, sides
    adiabatic. Integrating the 2D equation over x gives the 1D ODE
    k Tbar'' - rho cp u Tbar' = -qbar(y) with qbar = q_v a/W on the
    bar band; solved piecewise in closed form. Exponentials are
    clipped to non-positive arguments so any Peclet number is safe.
    Returns (Tbar(y), info) with the exact energy split: heat
    advected out of the top plus the conductive leak back through
    the inlet equals Q' identically."""
    y = np.asarray(y, float)
    y0, y1 = gap, gap + bh
    P = rho * cp * u / k
    qbar = q_v * bw / W
    m = qbar / (rho * cp * u)
    e0 = math.exp(-P * y0)
    e1 = math.exp(-P * y1)
    A3 = m * bh + (m / P) * (e1 - e0)
    ex = lambda arg: np.exp(np.minimum(arg, 0.0))
    th_lo = (m / P) * ((ex(P * (y - y0)) - ex(P * (y - y1)))
                       - (e0 - e1))
    th_md = (((m / P) * (1 - e0 + e1) - m * y0)
             - (m / P) * ex(P * (y - y1)) + m * y)
    th = np.where(y <= y0, th_lo, np.where(y >= y1, A3, th_md))
    return T_in + th, dict(
        T_out=T_in + A3, A3=A3, Pe=P * H,
        adv_W_per_m=rho * cp * u * W * A3,
        cond_bottom_W_per_m=k * W * m * (e0 - e1))


def fan_check_field(x, y, T, W, H, bw, bh, xc, gap, q_v, k, rho, cp, u,
                    T_in):
    """Judge a fan-rung field: inlet row at T_in, plane-mean profile
    against the exact closed form, outlet mean against the exact
    T_out, and the data-side advected power against the exact split.
    The full-2D wake comparison is a later rung; the mean profile is
    exact and carries the verdict."""
    order = np.argsort(y)
    y = np.asarray(y)[order]; T = np.asarray(T)[order, :]
    xo = np.argsort(x); x = np.asarray(x)[xo]; T = T[:, xo]
    rm = np.nanmean(T, axis=1)
    ex, info = fan_mean_exact(y, W, H, bw, bh, gap, q_v, k, rho, cp, u,
                              T_in)
    m = dict(info)
    m["inlet_dev"] = float(abs(rm[0] - T_in))
    m["outlet_dev"] = float(abs(rm[-1] - info["T_out"]))
    m["rms_mean"] = float(np.sqrt(np.nanmean((rm - ex) ** 2)))
    m["max_mean"] = float(np.nanmax(np.abs(rm - ex)))
    m["adv_data"] = float(rho * cp * u * W * (rm[-1] - T_in))
    m["adv_ratio"] = m["adv_data"] / max(info["adv_W_per_m"], 1e-30)
    tol = max(0.03, 0.15 * abs(info["A3"]))
    m["anchors_ok"] = (m["inlet_dev"] < 0.03 and
                       m["outlet_dev"] < tol and
                       m["rms_mean"] < max(0.05, tol))
    m["T_peak"] = float(np.nanmax(T))
    m["rms_out"] = m["rms_mean"]; m["max_out"] = m["max_mean"]
    m["mean_num"] = rm; m["mean_exact"] = ex
    m["x"] = x; m["y"] = y; m["T"] = T
    return m


def fea_p_basic(fl, W, H, a, x_off, gap, L_z, Q, k_b, rho_b,
                cp_b, T0, t_end, t_step, h_ext, T_amb):
    """Operating point + the exact adiabatic slope for the basic
    module. Oil properties are evaluated at T0 through the app's own
    film_props so the FEA and the app share one property source."""
    p = film_props(fl, T0)
    A_cell = a * a
    A_oil = W * H - A_cell
    C_per_m = rho_b * cp_b * A_cell + p["rho"] * p["cp"] * A_oil
    dTdt = Q / L_z / C_per_m
    return dict(W_tank=W, H_tank=H, a_cell=a, x_off=x_off,
                gap_bot=gap, L_z=L_z, Q_cell=Q, k_bat=k_b,
                rho_bat=rho_b, cp_bat=cp_b, k_oil=p["k"],
                rho_oil=p["rho"], cp_oil=p["cp"], T0=T0,
                t_end=t_end, t_step=t_step, h_ext=h_ext,
                T_amb=T_amb, dTdt_pred=dTdt,
                C_kJK=C_per_m * L_z / 1000.0,
                oil_name=fl["name"], bar_name="heat bar",
                top_fixed=True, T_top=25.0, k_mult=1.0,
                steady=True, cls="ipl_basic_2d")


def cases_tab(d, g, fl, res, cool_df, loop):
    st.markdown("#### The verification ladder")
    st.markdown(
        "Six small cases, each simple enough to check by hand, an "
        "experiment, or a textbook, each adding one physical effect. "
        "Every case calls the **same functions the full model uses** - "
        "`h_cell_side`, `h_tube_side`, `h_water_inside`, `fin_pack`, "
        "`plate_fin_area` - so agreement here IS validation of the "
        "model's building blocks, not of a copy. Adjust everything; "
        "the numbers below update live.")
    cs_fl_name = st.selectbox(
        "Coolant for the cases", list(cool_df["name"]),
        index=int((cool_df["name"] == d["coolant"]).idxmax()),
        key="cs_fluid")
    cfl = fluid_dict(cool_df[cool_df["name"] == cs_fl_name].iloc[0])
    case = st.radio(
        "Case", ["1 · One cell, still bath",
                 "2 · One cell, moving oil",
                 "3 · Sealed bath warm-up",
                 "4 · One cell to one water tube",
                 "5 · Circulation shoot-out",
                 "6 · Reconcile the full model"],
        horizontal=True, key="cs_case")

    # ---------------------------------------------------------- 1
    if case.startswith("1"):
        st.markdown(
            "**A single heated cylinder hangs in a large, still bath.** "
            "This is the base experiment - a cartridge-heated cylinder "
            "and a thermocouple reproduce it on a bench in an "
            "afternoon. The film is natural convection only.")
        c1, c2, c3, c4 = st.columns(4)
        Q = c1.slider("Heat in the cell [W]", 0.5, 40.0, 3.0, 0.5,
                      key="cs1_q")
        Dm = c2.slider("Diameter [mm]", 10.0, 80.0, 21.0, 0.5,
                       key="cs1_d") / 1000
        Hm = c3.slider("Height [mm]", 30.0, 300.0, 70.0, 5.0,
                       key="cs1_h") / 1000
        Tb = c4.slider("Bath temperature [°C]", 15.0, 60.0, 35.0, 1.0,
                       key="cs1_tb")
        r = case_still_bath(cfl, Q, Dm, Hm, Tb)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Surface temperature", f"{r['T_s']:.1f} °C")
        m2.metric("Film ΔT", f"{r['dT']:.1f} °C")
        m3.metric("h (natural)", f"{r['h']:.0f} W/m²·K")
        m4.metric("Rayleigh", f"{r['Ra']:.2e}")
        st.markdown("**The numbers, step by step**")
        st.markdown(
            "Area $A = \\pi D H$ = %.4f m². Churchill-Chu on a "
            "vertical wall gives $Nu(Ra=%.2e, Pr=%.0f)$, so "
            "$h = %.1f$ W/m²·K (gap factor %.2f, i.e. unconfined). "
            "The energy balance $Q = h A (T_s - T_\\infty)$ closes to "
            "%.4f W against the %.1f W input - the bisection residual "
            "is %.2e W. Check any line against Incropera table 9.x by "
            "hand." % (r["A"], r["Ra"], r["Pr"], r["h"], r["gapf"],
                       r["h"] * r["A"] * r["dT"], Q, r["resid"]))
        meas = st.number_input(
            "Measured surface temperature from your rig [°C] "
            "(0 = none)", 0.0, 250.0, 0.0, 0.1, key="cs1_meas")
        if meas > 0:
            err = r["T_s"] - meas
            st.markdown(
                f"Model {r['T_s']:.1f} °C vs measured {meas:.1f} °C: "
                f"**{err:+.1f} °C** "
                f"({100*abs(err)/max(meas-Tb,0.1):.0f}% of the film "
                f"ΔT). Natural-convection correlations carry ±15-20% "
                f"on h, roughly ±%.1f °C here."
                % (0.18 * r["dT"]))

    # ---------------------------------------------------------- 2
    elif case.startswith("2"):
        st.markdown(
            "**The same cylinder, but the oil moves.** Direction "
            "matters: a horizontal sweep is crossflow over the "
            "cylinder (Churchill-Bernstein); an upward push is a "
            "boundary layer ALONG the can (laminar flat plate). Same "
            "velocity, different films - this is the propeller "
            "question in isolation.")
        c1, c2, c3, c4 = st.columns(4)
        Q = c1.slider("Heat [W]", 0.5, 40.0, 3.0, 0.5, key="cs2_q")
        u = c2.slider("Oil velocity [m/s]", 0.0, 0.20, 0.05, 0.005,
                      key="cs2_u")
        Dm = c3.slider("Diameter [mm]", 10.0, 80.0, 21.0, 0.5,
                       key="cs2_d") / 1000
        Tb = c4.slider("Bath temperature [°C]", 15.0, 60.0, 35.0, 1.0,
                       key="cs2_tb")
        Hm = st.slider("Height [mm]", 30.0, 300.0, 70.0, 5.0,
                       key="cs2_h") / 1000
        rc = case_mode_temp(cfl, Q, Dm, Hm, Tb, u, "cross")
        ra = case_mode_temp(cfl, Q, Dm, Hm, Tb, u, "axial")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Crossflow (horizontal sweep)**")
            st.metric("Surface temperature", f"{rc['T_s']:.1f} °C")
            st.markdown(
                "$Re_D=%.0f$, Churchill-Bernstein forced film "
                "$h_f=%.0f$, natural $h_n=%.0f$, blended "
                "$h=(h_n^3+h_f^3)^{1/3}=%.0f$ W/m²·K."
                % (rc["Re"], rc["h_for"], rc["h_nat"], rc["h"]))
        with c2:
            st.markdown("**Axial (pushed upward)**")
            st.metric("Surface temperature", f"{ra['T_s']:.1f} °C")
            st.markdown(
                "$Re_L=%.0f$ on the height, flat-plate "
                "$Nu=0.664\\,Re_L^{1/2}Pr^{1/3}$ gives $h_f=%.0f$, "
                "blended $h=%.0f$ W/m²·K."
                % (u * Hm / film_props(cfl, Tb)["nu"], ra["h_for"],
                   ra["h"]))
        st.markdown(
            "At this speed the crossflow film is **%.2fx** the axial "
            "one, worth **%.1f °C** on the can. Buoyancy alignment "
            "helps the axial case only through the blend with "
            "$h_n=%.0f$." % (rc["h"] / max(ra["h"], 1e-9),
                             ra["T_s"] - rc["T_s"], ra["h_nat"]))
        meas = st.number_input(
            "Measured surface temperature [°C] (0 = none)", 0.0,
            250.0, 0.0, 0.1, key="cs2_meas")
        if meas > 0:
            st.markdown(
                f"Crossflow model {rc['T_s']:.1f}, axial model "
                f"{ra['T_s']:.1f}, measured {meas:.1f} °C.")

    # ---------------------------------------------------------- 3
    elif case.startswith("3"):
        st.markdown(
            "**Seal the bath, remove the sink.** N cells heat V "
            "litres of oil; only the box surface loses to the room. "
            "One line of physics: $m c_p\\,dT/dt = NQ - UA(T-T_{amb})$"
            ". The steady level and the time constant are closed-form "
            "- if the model's flywheel numbers are right, this curve "
            "is right.")
        c1, c2, c3, c4 = st.columns(4)
        N = c1.slider("Cells", 1, 2000, 1080, 1, key="cs3_n")
        Q = c2.slider("Heat per cell [W]", 0.1, 10.0, 2.3, 0.1,
                      key="cs3_q")
        V = c3.slider("Oil volume [L]", 1.0, 150.0, 63.0, 1.0,
                      key="cs3_v")
        Ue = c4.slider("Box U to room [W/m²·K]", 1.0, 20.0, 6.0, 0.5,
                       key="cs3_u")
        c1, c2, c3 = st.columns(3)
        Ae = c1.slider("Box area [m²]", 0.5, 12.0,
                       float(f"{g['A_box_ext']:.2f}"), 0.1,
                       key="cs3_a")
        Ta = c2.slider("Room [°C]", 10.0, 45.0, 25.0, 1.0,
                       key="cs3_ta")
        tend = c3.slider("Simulate [min]", 10.0, 600.0, 120.0, 10.0,
                         key="cs3_t")
        r = case_bath_warmup(cfl, N, Q, V, Ue, Ae, Ta, tend)
        m1, m2, m3 = st.columns(3)
        m1.metric("Steady oil temperature",
                  f"{r['T_ss']:.0f} °C" if r["T_ss"] < 500 else
                  "runaway")
        m2.metric("Time constant τ", f"{r['tau_min']:.0f} min")
        m3.metric("Thermal mass m·cp", f"{r['C_kJK']:.0f} kJ/K")
        fig = go.Figure(go.Scatter(x=r["t_min"], y=r["T"],
                                   line=dict(color="#F59E0B",
                                             width=3)))
        fig.add_hline(y=r["T_ss"], line_dash="dash",
                      line_color="#94A3B8")
        fig.update_layout(height=280, xaxis_title="minutes",
                          yaxis_title="oil °C",
                          margin=dict(l=8, r=8, t=8, b=8))
        st.plotly_chart(fig, width='stretch', key="cs3_fig")
        st.markdown(
            "$T_{ss}=T_{amb}+NQ/UA$ = %.0f + %.0f/%.1f = **%.0f °C**; "
            "$\\tau = mc_p/UA$ = %.0f kJ/K / %.1f W/K = **%.0f min**. "
            "Both are one-line hand checks; the curve is just their "
            "exponential." % (Ta, N * Q, r["UA"], r["T_ss"],
                              r["C_kJK"], r["UA"], r["tau_min"]))

    # ---------------------------------------------------------- 4
    elif case.startswith("4"):
        st.markdown(
            "**Now give the heat somewhere to go: one water tube.** "
            "This is the full model's spine in isolation - oil film "
            "on the tube, wall conduction, water film - with fins and "
            "a plate as optional area multipliers. Every resistance "
            "is printed so the ladder sums by hand.")
        c1, c2, c3, c4 = st.columns(4)
        Q = c1.slider("Heat [W]", 0.5, 60.0, 3.0, 0.5, key="cs4_q")
        Twin = c2.slider("Water in [°C]", 5.0, 40.0, 20.0, 1.0,
                         key="cs4_tw")
        lpm = c3.slider("Water flow [L/min]", 0.05, 6.0, 0.6, 0.05,
                        key="cs4_f")
        u = c4.slider("Oil sweep at the tube [m/s]", 0.0, 0.15, 0.0,
                      0.005, key="cs4_u")
        c1, c2, c3, c4 = st.columns(4)
        shape = c1.selectbox("Tube section",
                             ["Round", "Square", "Rectangular"],
                             key="cs4_shape")
        odm = c2.slider("Tube OD / width [mm]", 4.0, 25.0, 10.0, 0.5,
                        key="cs4_od") / 1000
        Lt = c3.slider("Tube length in oil [m]", 0.1, 2.0, 0.85,
                       0.05, key="cs4_l")
        fins = c4.checkbox("Annular fins", False, key="cs4_fins")
        plate_on = st.checkbox(
            "Bond one plate to the tube (both faces wetted)", False,
            key="cs4_pl")
        td = dict(tube_shape=shape, tube_od=odm, tube_w=odm,
                  tube_h=0.008, tube_wall=0.001, tube_mat="Copper",
                  L_tube_case=Lt)
        fin_geo = dict(H=0.008, t=0.0006, p=0.004, k=205.0)
        plate = None
        if plate_on:
            pd_ = dict(plate_on=True, plate_t=0.0015,
                       plate_mat="Aluminium", plate_contact=0.9,
                       n_tubes=1, manifold_margin=0.0, h_cell=0.07)
            pg_ = dict(n_rows=2, Lx=Lt)
            plate = dict(d=pd_, g=pg_)
        r = case_chain(cfl, loop, Q, 0.021, 0.070, Twin, lpm, td,
                       u=u, mode="cross" if u > 0 else "still",
                       fins=fins, fin_geo=fin_geo, plate=plate)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Oil settles at", f"{r['T_oil']:.1f} °C")
        m2.metric("Can surface", f"{r['T_s']:.1f} °C")
        m3.metric("Water film h", f"{r['h_water']:.0f} W/m²·K "
                                  f"({r['regime']})")
        m4.metric("Water rise", f"{r['dT_w']:.2f} °C")
        rows = [("Oil film on the tube",
                 r["R_ot"], Q * r["R_ot"]),
                ("Tube wall", r["R_wall"], Q * r["R_wall"]),
                ("Water film", r["R_in"], Q * r["R_in"])]
        st.markdown(
            "| Resistance | K/W | ΔT at %.1f W |\n|---|---|---|\n"
            % Q + "\n".join("| %s | %.4f | %.2f °C |" % rr
                            for rr in rows))
        s = sum(x[2] for x in rows)
        st.markdown(
            "Sum of the ladder = **%.2f °C**; the solver's "
            "$T_{oil}-\\bar T_w$ = **%.2f °C** - identical by "
            "construction, which is the point: the full model is this "
            "chain, replicated. Oil-side area: bare %.4f m²%s%s at "
            "$h_{tube}=%.0f$ W/m²·K."
            % (s, r["T_oil"] - r["T_w_mean"], r["A_o"],
               (", fins +%.4f m² (η %.2f)" % (r["A_fin"], r["eta_f"]))
               if r["A_fin"] > 0 else "",
               (", plate +%.4f m² (η %.2f)" % (r["A_pl"], r["eta_p"]))
               if r["A_pl"] > 0 else "", r["h_tube"]))
        meas = st.number_input(
            "Measured oil temperature [°C] (0 = none)", 0.0, 200.0,
            0.0, 0.1, key="cs4_meas")
        if meas > 0:
            st.markdown(f"Model {r['T_oil']:.1f} vs measured "
                        f"{meas:.1f} °C: **{r['T_oil']-meas:+.1f} °C**.")

    # ---------------------------------------------------------- 5
    elif case.startswith("5"):
        st.markdown(
            "**Three ways to move the same oil, one honest table.** "
            "Same cell, same heat, same speed - only the flow "
            "direction and its correlation change. This is the "
            "student's propeller argument settled at case level.")
        c1, c2, c3 = st.columns(3)
        Q = c1.slider("Heat [W]", 0.5, 40.0, 3.0, 0.5, key="cs5_q")
        u = c2.slider("Velocity [m/s]", 0.005, 0.15, 0.05, 0.005,
                      key="cs5_u")
        Tb = c3.slider("Oil bulk [°C]", 15.0, 60.0, 35.0, 1.0,
                       key="cs5_tb")
        rows = []
        for mode, label, note in [
                ("still", "Thermosiphon only",
                 "buoyant film, Churchill-Chu"),
                ("axial", "Bottom propeller (axial up)",
                 "flat plate along the can"),
                ("cross", "Stirred / channel sweep (crossflow)",
                 "Churchill-Bernstein")]:
            rr = case_mode_temp(cfl, Q, 0.021, 0.070, Tb, u, mode)
            rows.append((label, rr["h"], rr["T_s"], note))
        st.markdown(
            "| Circulation | h [W/m²·K] | Can surface [°C] | "
            "correlation |\n|---|---|---|---|\n" +
            "\n".join("| %s | %.0f | %.1f | %s |" % rr for rr in rows))
        st.markdown(
            "Crossflow beats axial at equal speed because the "
            "boundary layer restarts around the cylinder instead of "
            "growing along %.0f mm of height; buoyancy alignment "
            "cannot buy that back (it adds millimetres per second). "
            "Serpentine then wins overall in the full model by adding "
            "plate AREA on top of the crossflow film." % 70)

    # ---------------------------------------------------------- 6
    else:
        st.markdown(
            "**Decompose the full model you have on screen.** The "
            "solver's own resistances, multiplied by the heat "
            "actually flowing to water, must sum to the temperature "
            "gap it reports. Any daylight between the two columns "
            "would be a bug; the residual is the numerical solver "
            "tolerance.")
        Qw = res["Q_w"]
        Tw_mean = d["T_water_in"] + res["dT_water"] / 2.0
        rows = [("Cell can -> bulk oil", res["R_b"],
                 Qw * res["R_b"]),
                ("Bulk oil -> tube surface", res["R_ot"],
                 Qw * res["R_ot"]),
                ("Tube wall", res["R_wall"], Qw * res["R_wall"]),
                ("Water film", res["R_in"], Qw * res["R_in"])]
        st.markdown(
            "| Stage | R [K/W] | ΔT = Q_w·R [°C] |\n|---|---|---|\n" +
            "\n".join("| %s | %.5f | %.2f |" % rr for rr in rows))
        s = sum(x[2] for x in rows)
        actual = res["T_b"] - Tw_mean
        st.markdown(
            "Ladder sum **%.2f °C** vs the solver's "
            "$T_b - \\bar T_w$ = **%.2f °C**: residual **%.3f °C** "
            "(%.2f%% of the drop) - the bisection tolerance, nothing "
            "hidden. Heat to water $Q_w$ = %.0f W of %.0f W generated "
            "(the rest leaves through the box wall). Cases 1-5 are "
            "these rows in isolation; the full model adds only "
            "bookkeeping: N cells, n tubes, DCIR(T), busbars, "
            "parasitics." % (s, actual, s - actual,
                             100 * abs(s - actual) / max(actual, 1e-9),
                             Qw, res["Q_eff"]))
        st.markdown(
            "Films in play right now: cell side %.0f (natural %.0f, "
            "forced %.0f), tube side %.0f, water %.0f W/m²·K; "
            "self-driven thermosiphon %.1f mm/s. Flip any Design "
            "input and watch the same rows move."
            % (res["h_cell"], res["h_cell_nat"], res["h_cell_for"],
               res["h_tube"], res["h_water"], res["u_ts"] * 1000))





def battery_tab():
    import battery_hppc as BH
    import battery_card as BC
    st.markdown("#### Battery model from measurements")
    st.markdown(
        "Three ways in, clearly ranked. **(1) Raw test data**: a "
        "zip of BioLogic condition folders (like the JP50 set - "
        ".mpr streams or Split-file Excel exports) or a plain "
        "time/current/voltage CSV; the app extracts OCV(SOC), "
        "R0(SOC) and the RC branch per condition, proves each fit "
        "by re-simulating the whole measured stream, and builds a "
        "**battery model card**. **(2) A card** exported earlier: "
        "one CSV that carries the whole model - upload it and "
        "everything battery-aware uses it, no raw data needed; "
        "the same file loads straight into COMSOL (its '%' header "
        "lines are COMSOL-native comments) and the analytical "
        "layer. **(3) Nothing**: a simplified default is used and "
        "said so - constant 20 mΩ + 10 mΩ, linear 3.0-4.2 V OCV.")
    ccap, cord, cprr = st.columns(3)
    cap = ccap.number_input("Cell capacity [Ah]", 0.5, 200.0, 5.0,
                            0.1, key="bh_cap")
    order = 1 if cord.selectbox("Model order",
                                ["1 RC branch", "2 RC branches"],
                                key="bh_order").startswith("1") \
        else 2
    prr = cprr.number_input("Min rest before a pulse counts [s]",
                            5.0, 3600.0, 20.0, 5.0, key="bh_prr")
    cA, cB = st.columns(2)
    raw_ups = cA.file_uploader(
        "1) Raw data: zip of condition folders, or CSV(s)",
        type=["zip", "csv", "txt", "xlsx"],
        accept_multiple_files=True, key="bh_raw")
    card_up = cB.file_uploader(
        "2) Battery model card (the *_card_upload_to_app.csv "
        "file)", type=["csv"], key="bh_card")

    card = meta = None
    cond_results = []
    if card_up is not None:
        try:
            card, meta = BC.card_read(card_up)
            st.session_state["bat_src"] = (
                f"card: {meta.get('cell', '?')}")
        except Exception as e:
            st.error(f"Card rejected: {e}")
            card = None
    elif raw_ups:
        import tempfile
        for up in raw_ups:
            try:
                if up.name.lower().endswith(".zip"):
                    wd = tempfile.mkdtemp(prefix="bat_")
                    BC.extract_zip(up.read(), wd)
                    for dcond in BC.find_condition_dirs(wd):
                        nm = os.path.basename(dcond)
                        tC, rC, rep = BC.parse_condition(nm)
                        dfc, srcs, tmeas = BC.load_condition_dir(
                            dcond)
                        out = BH.run_pipeline(
                            dfc, cap_Ah=cap, order=order,
                            min_pre_rest=prr)
                        cond_results.append(dict(
                            name=nm,
                            temp_C=(tmeas if tmeas is not None
                                    else (tC if tC is not None
                                          else 25.0)),
                            rate_C=rC if rC is not None else 1.0,
                            rep=rep, out=out, source=srcs))
                else:
                    if up.name.lower().endswith(".xlsx"):
                        dfc = BC.load_split_xlsx([up])
                    else:
                        dfc = pd.read_csv(up)
                    out = BH.run_pipeline(dfc, cap_Ah=cap,
                                          order=order,
                                          min_pre_rest=prr)
                    cond_results.append(dict(
                        name=up.name, temp_C=25.0, rate_C=1.0,
                        rep=1, out=out, source="upload"))
            except Exception as e:
                st.error(f"{up.name}: {e}")
        if cond_results:
            card, meta = BC.campaign_to_card(
                cond_results, cap_Ah=cap,
                cell_name=os.path.splitext(
                    raw_ups[0].name)[0])
            st.session_state["bat_src"] = (
                f"fitted: {meta['cell']}")

    if card is None:
        st.session_state["bat_src"] = "default (simplified)"
        st.session_state["bat_model"] = BC.default_model()
        st.info(
            "No data or card loaded - the **simplified default** "
            "is active: R0 = 20 mΩ, R1 = 10 mΩ, τ = 30 s, OCV "
            "linear 3.0-4.2 V. Upload measurements or a card "
            "above to replace it everywhere.")
        return

    if cond_results:
        st.success("Campaign fitted. Per-condition validation "
                   "(the whole measured stream re-simulated):")
        vt = pd.DataFrame([
            dict(condition=c["name"],
                 T_meas=f"{c['temp_C']:.1f} °C",
                 rate=f"{c['rate_C']:g}C",
                 pulses=len(c["out"]["pulses"]),
                 rmse_mV=round(c["out"]["rmse_mv"], 2),
                 R0_50=f"{1000*float(c['out']['model']['r0'](50)):.2f} mΩ",
                 q_pulsing=f"{c['out']['heat']['q_mean_active']:.2f} W")
            for c in cond_results])
        st.dataframe(vt, width='stretch', height=230)
    else:
        st.success(
            f"Card loaded: {meta.get('cell', '?')}, "
            f"{len(card)} rows, capacity "
            f"{meta.get('capacity_Ah', '?')} Ah. Validation "
            "recorded at fit time: "
            + ("; ".join(meta.get("validation", []))
               if meta.get("validation") else "(none stored)"))

    conds = BC.conditions_in(card)
    labels = [f"{t:g} °C · {r:g}C · run {int(p)}"
              for t, r, p in conds]
    pick = st.selectbox("Condition for the model in use",
                        labels, key="bh_pick")
    tC, rC, rep = conds[labels.index(pick)]
    model = BC.model_from_card(card, temp_C=tC, rate_C=rC,
                               rep=rep)
    st.session_state["bat_model"] = model
    st.session_state["bat_card"] = card

    cq1, cq2, cq3 = st.columns(3)
    rate_q = cq1.number_input("Design C-rate for heat", 0.1, 10.0,
                              1.0, 0.1, key="bh_rq")
    soc_q = cq2.number_input("SOC for heat [%]", 5.0, 95.0, 50.0,
                             5.0, key="bh_sq")
    qcell = BC.q_steady(model, rate_q * cap, soc_q)
    cq3.metric("Heat per cell, sustained",
               f"{qcell*1000:.0f} mW")
    def _use_q():
        st.session_state["fxb_q"] = float(np.round(
            max(qcell, 0.01), 2))
    st.button("Use this heat as the FEA bar heat",
              on_click=_use_q, key="bh_useq")

    ce1, ce2 = st.columns(2)
    ce1.download_button(
        "Battery model card (CSV - app, COMSOL, analytical)",
        data=BC.card_write(card, meta or {}),
        file_name=f"{(meta or {}).get('cell', 'cell')}"
                  "_card_upload_to_app.csv",
        mime="text/csv", type="primary", key="bh_dlcard")
    qty = ce2.selectbox("2-column COMSOL export",
                        ["ocv_V", "r0_ohm", "r1_ohm", "tau1_s",
                         "q_pulse_W"], key="bh_qty")
    ce2.download_button(
        f"Download {qty} vs SOC at {tC:g} °C {rC:g}C",
        data=BC.comsol_two_col(card, qty, tC, rC),
        file_name=f"{qty}_{tC:g}degC_{rC:g}C_for_comsol.csv",
        mime="text/csv",
        key="bh_dl2col")

    fR = go.Figure()
    for (t_, r_, p_), lb in zip(conds, labels):
        sub = card[(card["temp_C"] == t_)
                   & (card["rate_C"] == r_)
                   & (card["rep"] == p_)
                   & (card["dir"] == "dis")].sort_values(
            "soc_pct")
        fR.add_scatter(x=sub["soc_pct"], y=1000 * sub["r0_ohm"],
                       mode="markers+lines", name=lb)
    fR.update_layout(height=300, xaxis_title="SOC [%]",
                     yaxis_title="R0 [mΩ] (discharge)",
                     legend=dict(orientation="h", y=1.12),
                     margin=dict(l=8, r=8, t=30, b=8),
                     title=dict(text="Ohmic resistance across the "
                                "campaign - the temperature and "
                                "rate story in one plot",
                                font=dict(size=12), x=0.02))
    st.plotly_chart(fR, width='stretch', key="bh_fr2")
    cP1, cP2 = st.columns(2)
    sub0 = card[(card["temp_C"] == tC) & (card["rate_C"] == rC)
                ].sort_values("soc_pct")
    fO = go.Figure()
    fO.add_scatter(x=sub0["soc_pct"], y=sub0["ocv_V"],
                   mode="markers", name="card points",
                   marker=dict(size=6, color="#38BDF8"))
    sg = np.linspace(5, 95, 90)
    fO.add_scatter(x=sg, y=model["ocv"](sg), mode="lines",
                   name="model", line=dict(color="#0EA5E9"))
    fO.update_layout(height=280, xaxis_title="SOC [%]",
                     yaxis_title="OCV [V]",
                     margin=dict(l=8, r=8, t=30, b=8),
                     legend=dict(orientation="h", y=1.12),
                     title=dict(text="OCV (selected condition)",
                                font=dict(size=12), x=0.02))
    cP1.plotly_chart(fO, width='stretch', key="bh_fo2")
    fQ = go.Figure()
    rates = np.linspace(0.2, 5.0, 60)
    for s_ in (20.0, 50.0, 80.0):
        fQ.add_scatter(x=rates,
                       y=[1000 * BC.q_steady(model, rr * cap, s_)
                          for rr in rates],
                       mode="lines", name=f"SOC {s_:.0f}%")
    fQ.update_layout(height=280, xaxis_title="C-rate",
                     yaxis_title="heat per cell [mW]",
                     margin=dict(l=8, r=8, t=30, b=8),
                     legend=dict(orientation="h", y=1.12),
                     title=dict(text="Sustained heat map "
                                "q = I²(R0+R1) from the card",
                                font=dict(size=12), x=0.02))
    cP2.plotly_chart(fQ, width='stretch', key="bh_fq2")
    st.caption(
        "The FEA tab shows which battery source is active. "
        "Charge-direction rows are in the card too (dir column); "
        "the plots above use discharge.")


def fea_tab(d, g, fl, cool_df, loop):
    st.markdown("#### FEA - the basic module (COMSOL, 2D + 3D)")
    st.markdown(
        "**COMSOL** conduction models of the benchmark case from your CFD comparison report: a "
        "long **heat bar** (a 5 mm aluminium bar as the battery "
        "surrogate) inside a liquid tank. Sides and bottom "
        "adiabatic; the **top either held at a fixed temperature** "
        "(25 °C sink, so a steady state exists) **or sealed** "
        "(adiabatic, transient only). Every visit generates BOTH "
        "files from the same settings - `ipl2d.java` (plane "
        "section, per metre) and `ipl3d.java` (the section "
        "extruded the full depth) - so 2D, 3D and the analytical "
        "solution can all be compared below. The liquid is a "
        "conducting solid; buoyant circulation is not solved - "
        "the effective-k multiplier stands in for it, and the "
        "value that reproduces the reference CFD (Ansys Fluent, "
        "from your CFD report) IS the circulation's Nusselt number. "
        "Every input is a named parameter in the exported files.")

    st.markdown(
        "**The test order** - do these in sequence, each step "
        "proves the next one's foundation: **1)** Load the "
        "benchmark, download both files, run them, upload the two "
        "`_upload_to_app.txt` files here - anchors and the "
        "three-way comparison must pass. **2)** Press the "
        "calibrate button on the fixed-top upload to set the "
        "k-multiplier from your own run. **3)** Switch the heat "
        "source to the battery model at your design current "
        "(same Q drives the analytical anchors and the COMSOL "
        "files). **4)** Switch Configuration to Fan upflow, "
        "re-run both files, upload - the outlet-mean anchor must "
        "pass; calibrate the fan speed if you ran a different "
        "one. **5)** Tick the water pipes, set count/size/film, "
        "re-run - judge by the BALANCE row in the summary "
        "table. **6)** Switch Bar geometry to Battery row and "
        "repeat 3-5 on the real cell sizes.")

    def _wp3():
        st.session_state.update(fxb_w=25.0, fxb_h=30.0,
                  fxb_gm="Custom bar",
                  fxb_bw=5.0, fxb_bh=5.0, fxb_lz=300.0,
                  fxb_pipes=False,
                  fxb_gap=10.0, fxb_xo=0.0, fxb_q=0.75, fxb_t0=25.0,
                  fxb_mat="Aluminium (benchmark values)",
                  fx_fluid="Deionized water",
                  fxb_cfg="Fixed-top (benchmark)",
                  fxb_ttop=25.0, fxb_km=1.0)
    st.button("Load the benchmark case (one click)",
              on_click=_wp3, type="primary", key="fxb_preset")
    st.caption("Battery source for heat: **"
               + st.session_state.get("bat_src",
                                      "default (simplified)")
               + "** - set in the Battery tab, where one click "
               "sends its sustained heat into the Q field below.")

    c1, c2, c3, c4 = st.columns(4)
    Wt = c1.number_input("Tank width [mm]", 5.0, 400.0, 25.0,
                         0.1, format="%.1f", key="fxb_w") / 1000
    Ht = c2.number_input("Tank height [mm]", 5.0, 400.0, 30.0,
                         0.1, format="%.1f", key="fxb_h") / 1000
    gmode = c3.selectbox("Bar geometry",
                         ["Custom bar",
                          "Battery row (count × cell size)"],
                         key="fxb_gm")
    if gmode.startswith("Battery"):
        bb1, bb2, bb3, bb4 = st.columns(4)
        n_cells = bb1.number_input("Cells in the row", 1, 200, 10,
                                   1, key="fxb_nc")
        cell_d = bb2.number_input("Cell diameter [mm]", 5.0, 80.0,
                                  21.0, 0.1, format="%.1f",
                                  key="fxb_cd") / 1000
        cell_h = bb3.number_input("Cell height [mm]", 10.0, 200.0,
                                  70.0, 0.1, format="%.1f",
                                  key="fxb_ch") / 1000
        b_w, b_h = cell_d, cell_h
        Lz = n_cells * cell_d
        bb4.metric("Row depth (= n × d)", f"{Lz*1000:.0f} mm")
        st.caption("Battery-row mode: the 2D bar is one cell's "
                   "cross-section (d wide, h tall) and the 3D "
                   "block is EXACTLY the row - d × h × (n·d). "
                   "The depth field is derived, not typed.")
    else:
        b_w = c4.number_input("Bar width [mm]", 1.0, 120.0, 5.0,
                              0.1, format="%.1f",
                              key="fxb_bw") / 1000
        cc0, cc1_, _ = st.columns(3)
        b_h = cc0.number_input("Bar height [mm]", 1.0, 200.0, 5.0,
                               0.1, format="%.1f",
                               key="fxb_bh") / 1000
        Lz = cc1_.number_input("Depth into the plane [mm]", 10.0,
                               4000.0, 300.0, 1.0, format="%.0f",
                               key="fxb_lz") / 1000
    c1, c2, c3, c4 = st.columns(4)
    gap = c1.number_input("Bar bottom above the floor [mm]", 0.0,
                          300.0, 10.0, 0.1, format="%.1f",
                          key="fxb_gap") / 1000
    xoff = c2.number_input("Sideways offset [mm]", -150.0, 150.0,
                           0.0, 0.1, format="%.1f",
                           key="fxb_xo") / 1000
    hsrc = c3.selectbox("Heat source",
                        ["Manual Q [W]",
                         "Battery model at a set current"],
                        key="fxb_hs")
    T0 = c4.slider("Initial temperature [°C]", 5.0, 50.0, 25.0, 1.0,
                   key="fxb_t0")
    bat_mode = hsrc.startswith("Battery")
    if bat_mode:
        import battery_card as _BCm
        _bm = st.session_state.get("bat_model") or             _BCm.default_model()
        hb1, hb2, hb3, hb4 = st.columns(4)
        I_cell = hb1.number_input("Constant current [A]", 0.1,
                                  300.0, 5.0, 0.1, format="%.1f",
                                  key="fxb_ic")
        soc_h = hb2.number_input("SOC [%]", 5.0, 95.0, 50.0, 5.0,
                                 key="fxb_socq")
        R0c = float(_bm["r0"](soc_h))
        R1c = float(_bm["r1"](soc_h))
        Q = I_cell * I_cell * (R0c + R1c)
        hb3.metric("R0 + R1 at SOC",
                   f"{1000*(R0c+R1c):.2f} mΩ")
        hb4.metric("Q = I²(R0+R1)", f"{Q:.3f} W")
        st.caption("Battery source in use: **"
                   + st.session_state.get("bat_src",
                                          "default (simplified)")
                   + "** - the exported COMSOL files carry "
                   "I_cell, R0_cell, R1_cell as named parameters "
                   "with Q_cell defined from them, so the current "
                   "is sweepable inside COMSOL itself. The same Q "
                   "drives every analytical anchor on this page - "
                   "one current, three models.")
    else:
        I_cell, R0c, R1c = 0.0, 0.0, 0.0
        Q = st.number_input("Bar heat, total [W]", 0.01, 200.0,
                            0.75, 0.01, key="fxb_q")
    c1, c2, c3, c4 = st.columns(4)
    mat = c1.selectbox("Bar material", ["Aluminium (benchmark values)",
                                        "Battery jelly-roll",
                                        "Custom"], key="fxb_mat")
    fx_fl = c2.selectbox("Liquid", list(cool_df["name"]),
                         index=int((cool_df["name"] ==
                                    "Deionized water").idxmax()),
                         key="fx_fluid")
    cfg = c3.selectbox("Configuration",
                       ["Fixed-top (benchmark)",
                        "Sealed (adiabatic, transient)",
                        "Fan upflow (open channel)"], key="fxb_cfg")
    top_fixed = cfg.startswith("Fixed")
    fan = cfg.startswith("Fan")
    Ttop, u_fan, T_in = 25.0, 0.0005, 25.0
    if top_fixed:
        Ttop = c4.slider("Top temperature [°C]", 5.0, 60.0, 25.0,
                         1.0, key="fxb_ttop")
        study = st.selectbox("Study", ["Stationary (steady state)",
                                       "Transient"], key="fxb_study")
        steady = study.startswith("Stat")
    elif fan:
        u_fan = c4.number_input("Fan upward speed [mm/s]", 0.02,
                                50.0, 0.5, 0.05, format="%.2f",
                                key="fxb_uf") / 1000.0
        T_in = st.slider("Inlet oil temperature [°C]", 5.0, 60.0,
                         25.0, 1.0, key="fxb_tin")
        steady = True
        st.caption("Fan mode: the floor is the inlet at the set "
                   "temperature, the lid is the outflow, and the "
                   "oil advects upward at the fan speed - the "
                   "open-channel abstraction of bottom fans (the "
                   "return path is external to the modelled "
                   "slice). Validation tip: at high speeds the "
                   "outlet rise shrinks below the mesh's mK "
                   "noise; ~0.5 mm/s keeps the anchor "
                   "well-resolved.")
    else:
        steady = False
        st.caption("Sealed tank: no steady state exists, so the "
                   "study is Transient.")
    if mat.startswith("Alum"):
        kb, rb, cb = 202.4, 2719.0, 871.0
        st.caption("Aluminium at the benchmark values (from your CFD report): k 202.4 "
                   "W/m·K, ρ 2719 kg/m³, cp 871 J/kg·K.")
    elif mat.startswith("Batt"):
        kb, rb, cb = 0.9, 2500.0, 900.0
        st.caption("Jelly-roll-like: k 0.9 (transverse), ρ 2500, "
                   "cp 900 - the later rung where the surrogate "
                   "becomes a cell.")
    else:
        cc1, cc2, cc3 = st.columns(3)
        kb = cc1.slider("Bar k [W/m·K]", 0.3, 400.0, 202.4, 0.1,
                        key="fxb_kb")
        rb = cc2.slider("Bar density [kg/m³]", 500.0, 9000.0, 2719.0,
                        10.0, key="fxb_rb")
        cb = cc3.slider("Bar cp [J/kg·K]", 300.0, 1500.0, 871.0,
                        1.0, key="fxb_cb")
    cc1, cc2, cc3 = st.columns([1.2, 1, 1])
    km = cc1.number_input("Liquid effective-k multiplier (1 = "
                          "pure conduction)", 0.5, 30.0, 1.0,
                          0.05, format="%.2f", key="fxb_km")
    if not steady:
        tend = cc2.slider("Transient time [min]", 5.0, 240.0,
                          30.0, 5.0, key="fxb_te") * 60.0
    else:
        tend = 1800.0
        cc2.caption("Transient time: not applicable to a "
                    "stationary study.")
    build_only = cc3.checkbox(
        "Build-only (java -> mph in seconds; press Compute in the "
        "Desktop)", False, key="fxb_bo")
    use_pipes = st.checkbox(
        "Water pipes at the top (circular channels; walls get the "
        "water-side film h_w to T_w)", False, key="fxb_pipes",
        disabled=(not top_fixed and not fan))
    if (not top_fixed and not fan) and st.session_state.get(
            "fxb_pipes"):
        st.caption("Pipes are unavailable in the sealed "
                   "configuration - a sealed tank with a sink is "
                   "no longer sealed.")
    n_pipes, d_pipe, pipe_drop, h_w, T_w = 0, 0.008, 0.007, \
        1500.0, 20.0
    if use_pipes and (top_fixed or fan):
        pp1, pp2, pp3, pp4, pp5 = st.columns(5)
        n_pipes = pp1.number_input("Number of pipes", 1, 12, 3, 1,
                                   key="fxb_np")
        d_pipe = pp2.number_input("Pipe OD [mm]", 2.0, 30.0, 8.0,
                                  0.1, format="%.1f",
                                  key="fxb_pd") / 1000
        pipe_drop = pp3.number_input("Centre below the lid [mm]",
                                     2.0, 100.0, 7.0, 0.1,
                                     format="%.1f",
                                     key="fxb_pdr") / 1000
        h_w = pp4.number_input("Water film h_w [W/m²K]", 100.0,
                               20000.0, 1500.0, 50.0,
                               key="fxb_hw")
        T_w = pp5.number_input("Water T_w [°C]", 0.0, 60.0, 20.0,
                               0.5, key="fxb_tw")
        st.caption("With pipes the closed-form anchors no longer "
                   "apply exactly (heat leaves mid-height); the "
                   "exported models add pipe-heat and BALANCE "
                   "rows to the summary table - judge by "
                   "those.")
    else:
        use_pipes = False

    cfl = fluid_dict(cool_df[cool_df["name"] == fx_fl].iloc[0])
    P = fea_p_basic(cfl, Wt, Ht, b_w, xoff, gap, Lz, Q, kb, rb, cb,
                    T0, tend, max(tend / 30.0, 10.0), 0.0, 25.0)
    # fixed names: every 2D export is ipl2d, every 3D export is
    # ipl3d, whatever the variant - both are generated on every
    # visit from the same settings.
    P.update(bar_name=mat, top_fixed=top_fixed, T_top=Ttop,
             k_mult=km, steady=steady, build_only=build_only,
             flow_mode="fan" if fan else
                       ("top" if top_fixed else "sealed"),
             u_fan=u_fan, T_in=T_in,
             heat_mode="battery" if bat_mode else "manual",
             I_cell=I_cell, R0_cell=R0c, R1_cell=R1c,
             b_w=b_w, b_h=b_h,
             n_pipes=(n_pipes if use_pipes else 0),
             d_pipe=d_pipe, pipe_drop=pipe_drop, h_w=h_w,
             T_w=T_w)
    qv = Q / (b_w * b_h * Lz)

    cx1, cx2 = st.columns([1.4, 1])
    with cx1:
        fig = go.Figure()
        fig.update_layout(height=340, margin=dict(l=8, r=8, t=30,
                                                  b=8),
                          plot_bgcolor="rgba(0,0,0,0)",
                          paper_bgcolor="rgba(0,0,0,0)",
                          showlegend=False,
                          title=dict(text="What the FEA will solve",
                                     x=0.01, font=dict(size=13)))
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False, scaleanchor="x")
        Wm, Hm = Wt * 1000, Ht * 1000
        bwm, bhm = b_w * 1000, b_h * 1000
        gm, xm = gap * 1000, xoff * 1000
        fig.add_shape(type="rect", x0=0, y0=0, x1=Wm, y1=Hm,
                      fillcolor="rgba(129,140,248,.25)",
                      line=dict(color="#6366F1", width=2))
        if top_fixed:
            fig.add_shape(type="line", x0=0, y0=Hm, x1=Wm, y1=Hm,
                          line=dict(color="#0369A1", width=5))
        if fan:
            for fx0 in (0.18, 0.5, 0.82):
                fig.add_annotation(x=Wm * fx0, y=Hm * 0.30,
                                   ax=Wm * fx0, ay=Hm * 0.06,
                                   xref="x", yref="y", axref="x",
                                   ayref="y", showarrow=True,
                                   arrowhead=3, arrowwidth=2,
                                   arrowcolor="#0EA5E9")
        bx0 = Wm / 2 - bwm / 2 + xm
        fig.add_shape(type="rect", x0=bx0, y0=gm, x1=bx0 + bwm,
                      y1=gm + bhm,
                      fillcolor="rgba(250,204,21,.85)",
                      line=dict(color="#A16207", width=2))
        fig.add_annotation(x=Wm / 2, y=Hm + Hm * 0.07,
                           text=(f"top: T = {Ttop:.0f} °C (fixed)"
                                 if top_fixed else
                                 "top: outflow" if fan else
                                 "top: adiabatic (sealed)"),
                           showarrow=False,
                           font=dict(size=10, color="#0369A1"
                                     if top_fixed else "#64748B"))
        fig.add_annotation(x=Wm / 2, y=-Hm * 0.08,
                           text=f"tank {Wm:.0f} × {Hm:.0f} mm · "
                                f"{fx_fl} (k × {km:.1f}) · " +
                                (f"floor inflow {T_in:.0f} °C · "
                                 f"u = {u_fan*1000:.1f} mm/s up"
                                 if fan else
                                 "sides and bottom adiabatic"),
                           showarrow=False,
                           font=dict(size=10, color="#475569"))
        fig.add_annotation(x=bx0 + bwm / 2, y=gm + bhm / 2,
                           text="bar", showarrow=False,
                           font=dict(size=9, color="#713F12"))
        fig.add_annotation(x=bx0 + bwm + Wm * 0.02,
                           y=gm + bhm / 2,
                           text=f"{bwm:.1f} × {bhm:.1f} mm, "
                                f"{gm:.0f} mm off the floor",
                           showarrow=False, xanchor="left",
                           font=dict(size=9, color="#A16207"))
        if use_pipes and n_pipes:
            rp = d_pipe / 2 * 1000
            zc = Hm - pipe_drop * 1000
            for kpi in range(int(n_pipes)):
                xc_p = Wm * (2 * kpi + 1) / (2 * n_pipes)
                fig.add_shape(type="circle",
                              x0=xc_p - rp, x1=xc_p + rp,
                              y0=zc - rp, y1=zc + rp,
                              fillcolor="rgba(56,189,248,.55)",
                              line=dict(color="#0369A1",
                                        width=1.5))
            fig.add_annotation(x=Wm * 0.02, y=zc,
                               text=f"{int(n_pipes)}× water "
                                    f"pipes, T_w = {T_w:.0f} °C",
                               showarrow=False, xanchor="left",
                               font=dict(size=9,
                                         color="#0369A1"))
        st.plotly_chart(fig, width='stretch', key="fxb_fig")
    with cx2:
        st.metric("Volumetric heat q_v", f"{qv:,.0f} W/m³")
        st.caption("The report uses 100,000 W/m³ - the preset lands "
                   "there exactly (0.75 W over 5×5×300 mm).")
        st.caption("3D note: the section is extruded the full "
                   "depth with adiabatic ends, so the exact "
                   "solution is z-invariant - the 3D run must "
                   "match the 2D field to solver noise, and the "
                   "checker below measures that from the 3D "
                   "export directly.")
        if fan:
            _fm, _fi = fan_mean_exact(
                np.array([Ht]), Wt, Ht, b_w, b_h, gap,
                Q / (b_w * b_h * Lz),
                km * P["k_oil"], P["rho_oil"], P["cp_oil"], u_fan,
                T_in)
            st.metric("Exact outlet mean (the anchor)",
                      f"{_fi['T_out']:.4f} °C")
            _leak = 100 * _fi["cond_bottom_W_per_m"] / (
                Q / Lz)
            st.markdown(
                "**Correctness check:** integrating the 2D "
                "equation across the width gives an exact 1D "
                "advection-diffusion balance for the plane-mean "
                "temperature - closed form at any Péclet number "
                f"(here Pe = {_fi['Pe']:.0f}). The exported models "
                "tabulate the FEA outlet mean against this exact "
                "value, and the energy split is exact too: "
                f"advected out of the top plus the "
                f"{_leak:.1f}% conductive leak back through the "
                "inlet equals the bar heat identically. The full "
                "2D wake comparison joins in a later rung; the "
                "mean profile carries the verdict.")
        if top_fixed:
            st.metric("Anchors at the solution",
                      f"top flux = {Q/Lz:.2f} W/m (2D) · "
                      f"{Q:.2f} W total (3D)")
            st.markdown(
                "**Correctness check:** at steady state every watt "
                "must leave through the fixed top; both exported "
                "models integrate that flux and tabulate its "
                "deviation. **Against the reference CFD** (Ansys "
                "Fluent, from your CFD report): with $k_{mult}=1$ "
                "these conduction-only models read hotter than the "
                "CFD's 25.7 °C, because the CFD resolves the "
                "~0.85 mm/s buoyant plume. Raise $k_{mult}$ until "
                "the fields match - that value is the circulation's "
                "Nusselt number, the quantity this app's "
                "correlations predict.")
        else:
            st.metric("Exact heating slope",
                      f"{P['dTdt_pred']*60:.4f} °C/min")
            st.markdown(
                "**Correctness check:** sealed and adiabatic, the "
                "volume-average temperature must climb this exact "
                "line; the exported model tabulates its own "
                "deviation at every step.")
    jav2 = fea_export.comsol_basic_2d(dict(P, cls="ipl2d"))
    jav3 = fea_export.comsol_basic_3d(dict(P, cls="ipl3d"))
    cdl1, cdl2 = st.columns(2)
    cdl1.download_button("COMSOL 2D model (ipl2d.java)", data=jav2,
                         file_name="ipl2d.java", mime="text/plain",
                         use_container_width=True, type="primary",
                         key="fxb_dl2")
    cdl2.download_button("COMSOL 3D model (ipl3d.java)", data=jav3,
                         file_name="ipl3d.java", mime="text/plain",
                         use_container_width=True, type="primary",
                         key="fxb_dl3")
    st.caption(
        "Both files carry the same settings. Run each in two steps "
        "(COMSOL 6.x batch takes .mph or a compiled .class, not "
        "raw .java): `comsol compile ipl2d.java` then `comsol "
        "batch -inputfile ipl2d.class`, and likewise ipl3d. Each "
        "run writes <name>_summary.txt and <name>_upload_to_app.txt "
        "and saves the .mph. The 3D solve takes minutes and its "
        "field file is a few MB. FEMM's steady solver can twin the "
        "fixed-top case as the next cross-check rung if wanted.")

    st.markdown("---")
    st.markdown("##### Check a COMSOL run against the analytics")
    st.markdown(
        "Each solved run writes its full temperature field "
        "automatically: `ipl2d_upload_to_app.txt` and "
        "`ipl3d_upload_to_app.txt`. **Drop BOTH here for the full "
        "three-way comparison** - 2D vs 3D vs analytical - or "
        "either alone (the `_summary.txt` tables are optional and "
        "just displayed). For every field the app applies the "
        "**rigorous anchors** first - top row equal to $T_{top}$, "
        "plane-mean flat below the bar and exactly linear above it "
        "with slope $Q'/(kW)$; pure energy conservation, no "
        "approximation - then compares against the **analytical "
        "series solution** outside the bar footprint. Checked "
        "against the case configured above, so press the "
        "benchmark preset (or match the inputs) before "
        "uploading.")
    ups = st.file_uploader(
        "Upload the exported file(s)", type=["txt", "dat", "csv"],
        accept_multiple_files=True, key="fx_up")
    st.session_state["fx_fields"] = {}
    _fchecks = []
    for up in ups or []:
        try:
            txt = up.read().decode("utf-8", errors="replace")
        except Exception as e:
            st.error(f"{up.name}: could not read ({e})")
            continue
        if "DEVIATION" in txt or "bar peak" in txt.lower():
            st.markdown(f"**{up.name}** (evaluation table):")
            st.code("\n".join(txt.splitlines()[:14]),
                    language=None)
            continue
        try:
            fx_, fy_, fT_, finf = parse_comsol_field(txt)
        except Exception as e:
            st.error(f"{up.name}: {e}")
            continue
        k_eff = km * P["k_oil"]
        xc_ = Wt / 2 + xoff
        if (abs(fx_.max() - fx_.min() - Wt) > 0.01 * Wt or
                abs(fy_.max() - fy_.min() - Ht) > 0.01 * Ht):
            st.warning(
                f"{up.name}: field extents "
                f"{(fx_.max()-fx_.min())*1000:.1f} x "
                f"{(fy_.max()-fy_.min())*1000:.1f} mm do not match "
                f"the tank set above ({Wt*1000:.0f} x "
                f"{Ht*1000:.0f} mm) - fix the sliders or the "
                "upload before trusting the verdict.")
        if finf["kind"] == "3d":
            st.info(
                f"{up.name}: 3D field ({finf['depth_n']} depth "
                f"planes). Worst variation across the depth at any "
                f"point: **{(finf['zvar'] or 0)*1000:.1f} mK** - "
                "for this z-invariant case that is the 2D-vs-3D "
                "agreement, measured; the report's Fluent pair "
                "differed by ~20 mK. Checking the mid-depth slice "
                "below.")
        if use_pipes:
            st.info(
                f"{up.name}: pipes are present, so the "
                "closed-form anchors do not apply - judge this "
                "run by the pipe-heat and BALANCE rows in its "
                "summary table (upload it too). Field shown for "
                "inspection.")
            fgp = go.Figure(go.Heatmap(
                x=fx_ * 1000, y=fy_ * 1000, z=fT_,
                colorscale="Turbo", colorbar=dict(thickness=10)))
            fgp.update_layout(height=320, margin=dict(l=6, r=6,
                              t=28, b=6),
                              yaxis=dict(scaleanchor="x"),
                              title=dict(text="Numerical field "
                                         "(pipes case)",
                                         font=dict(size=12),
                                         x=0.02))
            st.plotly_chart(fgp, width='stretch',
                            key=f"fx_pipe_{up.name}")
            st.session_state.setdefault("fx_fields", {})[
                up.name] = (fx_, fy_, fT_)
            continue
        if fan:
            m = fan_check_field(fx_, fy_, fT_, Wt, Ht, b_w, b_h,
                                xc_, gap,
                                Q / (b_w * b_h * Lz), k_eff,
                                P["rho_oil"], P["cp_oil"], u_fan,
                                T_in)
        else:
            m = wp3_check_field(fx_, fy_, fT_, Wt, Ht, b_w, b_h,
                                xc_, gap,
                                Q / (b_w * b_h * Lz), k_eff, Ttop)
        st.session_state.setdefault("fx_fields", {})[up.name] = (
            fx_, fy_, fT_)
        _fchecks.append((up.name, finf, m))
        if fan:
            if m["anchors_ok"]:
                st.success(
                    f"{up.name}: fan anchors PASS - inlet row "
                    f"within {m['inlet_dev']*1000:.0f} mK of T_in, "
                    f"outlet mean within "
                    f"{m['outlet_dev']*1000:.0f} mK of the exact "
                    f"closed form, mean-profile RMS "
                    f"{m['rms_mean']*1000:.0f} mK, advected/exact "
                    f"= {m['adv_ratio']:.3f}.")
            else:
                st.error(
                    f"{up.name}: fan anchors FAIL (inlet "
                    f"{m['inlet_dev']*1000:.0f} mK, outlet "
                    f"{m['outlet_dev']*1000:.0f} mK, RMS "
                    f"{m['rms_mean']*1000:.0f} mK) - check u_fan "
                    "and T_in match the run, then convergence.")
            g1, g2, g3 = st.columns(3)
            g1.metric("Bar peak (numerical)",
                      f"{m['T_peak']:.3f} °C")
            g2.metric("Exact outlet mean",
                      f"{m['T_out']:.4f} °C")
            g3.metric("Péclet number", f"{m['Pe']:.0f}")
            fgn = go.Figure(go.Heatmap(
                x=m["x"] * 1000, y=m["y"] * 1000, z=m["T"],
                colorscale="Turbo", colorbar=dict(thickness=10)))
            fgn.update_layout(height=320,
                              title=dict(text="Numerical field "
                                         "(wake model joins in a "
                                         "later rung)",
                                         font=dict(size=12),
                                         x=0.02),
                              margin=dict(l=6, r=6, t=28, b=6),
                              yaxis=dict(scaleanchor="x"))
            st.plotly_chart(fgn, width='stretch',
                            key=f"fx_fan_{up.name}")
            fp = go.Figure()
            fp.add_scatter(x=m["mean_num"], y=m["y"] * 1000,
                           mode="markers",
                           name="numerical row mean",
                           marker=dict(size=5, color="#F59E0B"))
            fp.add_scatter(x=m["mean_exact"], y=m["y"] * 1000,
                           mode="lines",
                           name="exact closed form (any Pe)",
                           line=dict(color="#38BDF8", width=2.5))
            fp.update_layout(height=320,
                             xaxis_title="plane-mean T [°C]",
                             yaxis_title="height y [mm]",
                             margin=dict(l=8, r=8, t=26, b=8),
                             legend=dict(orientation="h", y=1.05),
                             title=dict(text="The rigorous anchor: "
                                        "mean profile vs exact "
                                        "advection-diffusion",
                                        font=dict(size=12),
                                        x=0.02))
            st.plotly_chart(fp, width='stretch',
                            key=f"fx_fanprof_{up.name}")
            _dTo = float(m["mean_num"][-1] - T_in)
            if _dTo > 1e-4:
                _ufit = (Q / Lz) / (P["rho_oil"] * P["cp_oil"]
                                    * Wt * _dTo)
                def _cal_u(v=_ufit):
                    st.session_state["fxb_uf"] = float(
                        np.round(v * 1000.0, 2))
                st.button(
                    f"Calibrate: set the fan speed to what this "
                    f"run implies ({_ufit*1000:.2f} mm/s from "
                    f"the outlet rise)",
                    on_click=_cal_u, key=f"fx_calu_{up.name}")
            continue
        if m["anchors_ok"]:
            st.success(
                f"{up.name}: energy anchors PASS - top row within "
                f"{m['top_dev']:.3f} °C of T_top, sub-bar mean flat "
                f"to {m['below_flat']:.3f} °C, above-bar mean slope "
                f"within {100*m['slope_rel']:.1f}% of the exact "
                f"Q'/(kW). The solve is internally consistent.")
        else:
            st.error(
                f"{up.name}: energy anchors FAIL (top "
                f"{m['top_dev']:.2f} °C, flat {m['below_flat']:.2f} "
                f"°C, slope {100*m['slope_rel']:.0f}%) - wrong "
                "parameters, wrong boundary condition, or an "
                "unconverged solve. Do not compare further; fix "
                "this first.")
        g1, g2, g3, g4 = st.columns(4)
        g1.metric("Bar peak (numerical)", f"{m['T_peak']:.2f} °C")
        g2.metric("Series bar estimate",
                  f"{m['T_bar_series']:.2f} °C")
        g3.metric("Field RMS vs series (outside bar)",
                  f"{m['rms_out']:.3f} °C")
        g4.metric("Worst point (outside bar)",
                  f"{m['max_out']:.3f} °C")
        wp3_like = (abs(Wt - 0.025) < 1e-4 and
                    abs(Ht - 0.030) < 1e-4 and
                    abs(a - 0.005) < 1e-4 and
                    abs(Q - 0.75) < 0.02)
        if wp3_like and abs(km - 1.0) < 1e-6:
            nu = (m["T_peak"] - Ttop) / 0.72
            st.markdown(
                f"**Against the reference CFD (Fluent):** conduction-only peak "
                f"rise {m['T_peak']-Ttop:.2f} °C vs Fluent's 0.72 "
                f"°C - the buoyant circulation is worth a factor "
                f"of **{nu:.1f}** here. Raise the k-multiplier to "
                f"about {nu:.1f} and re-export to emulate it.")
        zmin = float(min(np.nanmin(m["T"]), np.nanmin(m["Ts"])))
        zmax = float(max(np.nanmax(m["T"]), np.nanmax(m["Ts"])))
        h1, h2, h3 = st.columns(3)
        for col, Z, ttl in ((h1, m["T"], "Numerical (COMSOL)"),
                            (h2, m["Ts"], "Analytical (series)"),
                            (h3, m["T"] - m["Ts"],
                             "Difference [°C]")):
            fg = go.Figure(go.Heatmap(
                x=m["x"] * 1000, y=m["y"] * 1000, z=Z,
                colorscale="Turbo" if ttl[0] != "D" else "RdBu",
                zmin=None if ttl[0] == "D" else zmin,
                zmax=None if ttl[0] == "D" else zmax,
                colorbar=dict(thickness=10)))
            fg.update_layout(height=300, title=dict(text=ttl,
                             font=dict(size=12), x=0.02),
                             margin=dict(l=6, r=6, t=28, b=6),
                             yaxis=dict(scaleanchor="x"))
            col.plotly_chart(fg, width='stretch',
                             key=f"fx_hm_{up.name}_{ttl[:4]}")
        fp = go.Figure()
        fp.add_scatter(x=m["mean_num"], y=m["y"] * 1000,
                       mode="markers", name="numerical row mean",
                       marker=dict(size=5, color="#F59E0B"))
        fp.add_scatter(x=m["mean_exact"], y=m["y"] * 1000,
                       mode="lines", name="exact energy integral",
                       line=dict(color="#38BDF8", width=2.5))
        fp.update_layout(height=320, xaxis_title="plane-mean T [°C]",
                         yaxis_title="height y [mm]",
                         margin=dict(l=8, r=8, t=26, b=8),
                         legend=dict(orientation="h", y=1.05),
                         title=dict(text="The rigorous anchor: mean "
                                    "profile, numerical vs exact",
                                    font=dict(size=12), x=0.02))
        st.plotly_chart(fp, width='stretch',
                        key=f"fx_prof_{up.name}")
        if m.get("slope_meas") is not None and \
                m["slope_meas"] < 0:
            _keff = -(Q / Lz) / (m["slope_meas"] * Wt)
            _kfit = _keff / P["k_oil"]
            def _cal_k(v=_kfit):
                st.session_state["fxb_km"] = float(np.round(v, 2))
            st.button(
                f"Calibrate: set the k-multiplier to what this "
                f"run implies ({_kfit:.2f} from the exact "
                f"mean-slope relation)",
                on_click=_cal_k, key=f"fx_calk_{up.name}")
        st.caption(
            "Caveat stated once and honestly: the series assumes "
            "uniform liquid conductivity, so INSIDE the aluminium "
            "bar it overshoots (the real bar is nearly isothermal); "
            "that is why the field comparison masks the bar "
            "footprint and the bar estimate above is the series "
            "average over it. The mean-profile anchor carries no "
            "such caveat - it is exact.")

    if len(_fchecks) >= 2:
        st.markdown("##### Three-way comparison: 2D · 3D · "
                    "analytical")
        rows = []
        for nm_, fi_, mm_ in _fchecks:
            rows.append((nm_,
                         "3D" if fi_["kind"] == "3d" else "2D",
                         f"{mm_['T_peak']:.3f}",
                         f"{mm_['rms_out']*1000:.0f}",
                         f"{mm_['max_out']*1000:.0f}",
                         "PASS" if mm_["anchors_ok"] else "FAIL",
                         (f"{fi_['zvar']*1000:.0f}"
                          if fi_.get("zvar") is not None else "-")))
        st.markdown(
            "| file | kind | bar peak [°C] | RMS vs analytical "
            "[mK] | worst vs analytical [mK] | energy anchors | "
            "depth variation [mK] |\n|---|---|---|---|---|---|---|"
            "\n" + "\n".join("| %s | %s | %s | %s | %s | %s | %s |"
                              % r for r in rows))
        _tbs = _fchecks[0][2].get("T_bar_series")
        if _tbs is not None:
            st.markdown(
                "Analytical bar estimate (series average over the "
                f"footprint): **{_tbs:.3f} °C** - remember it "
                "overshoots the near-isothermal aluminium bar by "
                "construction.")
        _flds = list(st.session_state.get("fx_fields", {}).items())
        for i_ in range(len(_flds)):
            for j_ in range(i_ + 1, len(_flds)):
                nA, (xA, yA, TA) = _flds[i_]
                nB, (xB, yB, TB) = _flds[j_]
                if TA.shape != TB.shape:
                    continue
                dd = TA - TB
                st.markdown(
                    f"**{nA} vs {nB}** (FEA vs FEA, same grid): "
                    f"RMS **{float(np.sqrt(np.nanmean(dd**2)))*1000:.0f} mK**, "
                    f"worst **{float(np.nanmax(np.abs(dd)))*1000:.0f} mK**.")
                fgd = go.Figure(go.Heatmap(
                    x=xA * 1000, y=yA * 1000, z=dd,
                    colorscale="RdBu",
                    colorbar=dict(thickness=10)))
                fgd.update_layout(
                    height=300, margin=dict(l=6, r=6, t=28, b=6),
                    title=dict(text=f"{nA} − {nB} [°C]",
                               font=dict(size=12), x=0.02),
                    yaxis=dict(scaleanchor="x"))
                st.plotly_chart(fgd, width='stretch',
                                key=f"fx_x_{i_}_{j_}")


def system_tab(d, g, fl, res, masses, loop, chil, Q_duty, C_steady):
    import pandas as _pd
    st.markdown("#### What the system physically needs: sizing, bill of "
                "materials, two-phase, and how it compares")
    st.caption("Everything here is computed from your live design and a "
               "set of editable assumptions. It turns the thermal answer "
               "into hardware: what to buy, how big, how heavy, how much "
               "it costs, and which architecture is right for the job.")

    Q_w = res.get("Q_w", Q_duty)
    dTw = res.get("dT_water", 2.0)
    Ecap = masses["E_kwh"]

    # ---- editable assumptions ----
    with st.expander("Cost, sizing and duty assumptions (edit these)",
                     expanded=False):
        a1, a2, a3, a4 = st.columns(4)
        oil_L = a1.number_input("Dielectric ester [£/L]", 3.0, 60.0, 9.0,
                                0.5, key="sy_oil")
        cu_kg = a2.number_input("Copper [£/kg]", 4.0, 20.0, 8.5, 0.5,
                                key="sy_cu")
        al_kg = a3.number_input("Aluminium [£/kg]", 1.5, 8.0, 3.2, 0.1,
                                key="sy_al")
        fab = a4.number_input("Fabrication ×", 1.0, 3.0, 1.8, 0.1,
                              key="sy_fab")
        b1, b2, b3, b4 = st.columns(4)
        pump_gbp = b1.number_input("Pump/circulator [£]", 10.0, 400.0,
                                   28.0, 2.0, key="sy_pump")
        chil_kW = b2.number_input("Chiller [£/kW electrical]", 60.0,
                                  400.0, 140.0, 5.0, key="sy_chil")
        hx_U = b3.number_input("HX overall U [W/m²·K]", 500.0, 6000.0,
                               3000.0, 100.0, key="sy_hxU",
                               help="Plate dielectric-to-water HX; "
                                    "3000 is typical for a brazed plate "
                                    "unit.")
        hx_dT = b4.number_input("HX approach ΔT [°C]", 2.0, 15.0, 5.0,
                                0.5, key="sy_hxdt")
        c1_, c2_ = st.columns(2)
        Trange = c1_.slider("Service temperature band [°C] (for "
                            "expansion vessel)", 40.0, 120.0, 70.0, 5.0,
                            key="sy_trange")
        pump_eta = c2_.slider("Pump total efficiency", 0.15, 0.7, 0.35,
                              0.05, key="sy_peta")

    # =============================================================== #
    st.markdown("##### 1 · The duty the hardware must serve")
    dc = st.columns(4)
    dc[0].metric("Continuous heat", f"{Q_duty/1000:.2f} kW",
                 f"at {C_steady:.2f}C RMS")
    dc[1].metric("Into the water", f"{Q_w/1000:.2f} kW",
                 f"ΔT {dTw:.1f} °C rise")
    dc[2].metric("Pack energy", f"{Ecap:.1f} kWh")
    dc[3].metric("Cell limit", f"{d['T_limit']:.0f} °C",
                 "core" if d.get("limit_core") else "can")
    st.caption(f"The loop carries {Q_w/1000:.2f} kW from the dielectric "
               f"to the water and the chiller lifts it to ambient. "
               f"Everything downstream is sized to move that heat with "
               f"margin, hold the cells at or below {d['T_limit']:.0f} °C, "
               f"and keep cell-to-cell spread under ~5 °C.")

    # =============================================================== #
    st.markdown("##### 2 · Component sizing (from your design)")
    # -- water pump --
    wpp = water_pump_power(d, g, loop)
    mdot_w = d["flow_lpm"] / 60.0 * loop["rho"] / 1000.0
    head_m = wpp["dp"] / (loop["rho"] * 9.81)
    P_hyd_w = wpp["dp"] * mdot_w / loop["rho"]
    P_pump_w = P_hyd_w / pump_eta
    # -- oil circulation --
    ext = d["circ"].startswith("External")
    xp = None
    if ext:
        u_oil = d["u_oil"]
        xp = ext_loop_pump(d, g, fl, u_oil)
        P_oil = xp["P"]
        oil_mode = ("external pump loop"
                    + (" + serpentine plates" if d.get("plate_on") else
                       " (row channels)"))
    elif d["circ"].startswith("Bottom"):
        u_oil = d["u_oil"]
        _pp = prop_power(d, g, fl, u_oil)
        P_oil = _pp["P"]
        oil_mode = (f"bottom propeller, axial up "
                    f"({_pp['Vdot_lpm']:.0f} L/min swept)")
    elif d.get("plate_on"):
        u_oil = d.get("u_guided", 0.05)
        P_oil = serpentine_pump(d, g, fl, u_oil)["P"]
        oil_mode = "serpentine plates"
    elif d["u_oil"] > 1e-6:
        u_oil = d["u_oil"]; P_oil = stirrer_power(d, g, fl, u_oil)
        oil_mode = "magnetic stirrer"
    else:
        u_oil = 0.0; P_oil = 0.0; oil_mode = "thermosiphon (passive)"
    V_pipe_L = (math.pi * d.get("pipe_id", 0.019) ** 2 / 4
                * d.get("pipe_len", 2.5) * 1000 + 0.3) if ext else 0.0
    # -- internal dielectric->water HX --
    lm = max(hx_dT, 1.0)
    UA_hx = Q_w / lm
    A_hx = UA_hx / hx_U
    # -- chiller --
    P_chil = chil["P_el"]; COP = chil["COP"]
    Q_reject = Q_w + P_pump_w + P_oil
    # -- coolant --
    V_oil_L = masses["V_oil_L"]; m_oil = masses["m_oil"]
    # -- tubes / plates --
    m_tubes = masses["m_tubes"]; m_plates = masses.get("m_plates", 0.0)
    tube_mat = d.get("tube_mat", "Copper")
    tube_L = d["n_tubes"] * g["L_tube"]
    # -- expansion vessel --
    dV_exp = fl["beta"] * ((V_oil_L + (0.0 if not d["circ"].startswith(
        "External") else math.pi * d.get("pipe_id", 0.019) ** 2 / 4
        * d.get("pipe_len", 2.5) * 1000 + 0.3)) / 1000.0) * Trange
    V_vessel_L = 1.3 * dV_exp * 1000.0

    sc = st.columns(3)
    with sc[0]:
        st.markdown("**Water pump / circulator**")
        st.markdown(
            f"- Flow **{d['flow_lpm']:.0f} L/min**, head "
            f"**{head_m:.2f} m** ({wpp['dp']/1000:.2f} kPa)\n"
            f"- Hydraulic {P_hyd_w:.1f} W → electrical "
            f"**{P_pump_w:.1f} W** at η={pump_eta:.2f}\n"
            f"- Tube velocity {wpp['v']:.2f} m/s, Re "
            f"{wpp['Re']:.0f} ({res.get('water_regime','')})\n"
            f"- Class: {'small brushless circulator' if P_pump_w < 30 else 'automotive coolant pump'}")
        st.markdown("**Oil circulation**")
        if ext and xp:
            st.markdown(
                f"- Mode: **{oil_mode}** - closed loop, pump outside "
                f"the pack, no external HX; oil and water never meet\n"
                f"- Loop flow **{xp['Vdot_lpm']:.1f} L/min** at "
                f"{u_oil*100:.0f} cm/s through the pack\n"
                f"- Δp: pack **{xp['dp_pack']/1000:.2f} kPa** + pipes "
                f"**{xp['dp_pipe']/1000:.2f} kPa** = "
                f"**{xp['dp']/1000:.2f} kPa** "
                f"({xp['dp']/(fl['rho']*9.81):.2f} m head)\n"
                f"- Pipe velocity {xp['v_pipe']:.2f} m/s "
                f"({'OK' if xp['v_pipe']<=1.5 else 'high - use a larger bore'}), "
                f"Re {xp['Re_pipe']:.0f}\n"
                f"- Motor: **{P_oil:.1f} W electrical** at η=0.35 "
                f"wire-to-fluid - a "
                f"{'small sealed BLDC circulator' if P_oil<40 else 'automotive-class oil pump' if P_oil<250 else 'industrial gear/centrifugal pump'}")
        else:
            st.markdown(f"- Mode: **{oil_mode}**\n"
                        f"- {('~%.2f W at %.0f mm/s' % (P_oil, u_oil*1000)) if u_oil>0 else 'no pump - buoyancy only'}")
    with sc[1]:
        st.markdown("**Dielectric→water heat exchanger**")
        st.markdown(
            f"- Duty **{Q_w/1000:.2f} kW**, approach {hx_dT:.0f} °C\n"
            f"- UA = Q/ΔT = **{UA_hx:.0f} W/K**\n"
            f"- Area ≈ UA/U = **{A_hx*1e4:.0f} cm²** "
            f"({A_hx:.3f} m²) at U={hx_U:.0f}\n"
            f"- Type: brazed-plate or a coil in the reservoir")
        st.markdown("**Chiller / dry-cooler**")
        st.markdown(
            f"- Reject **{Q_reject/1000:.2f} kW** to ambient\n"
            f"- Electrical **{P_chil/1000:.2f} kW** at COP {COP:.1f}\n"
            f"- A dry-cooler (no compressor) works if a warm water set "
            f"point is acceptable")
    with sc[2]:
        st.markdown("**Coolant (dielectric)**")
        _vloop = f" (+{V_pipe_L:.1f} L in the external loop)" if ext else ""
        st.markdown(
            f"- Volume **{V_oil_L:.1f} L**{_vloop}, mass "
            f"**{m_oil:.1f} kg**\n"
            f"- β = {fl['beta']:.1e} /K → expands "
            f"{dV_exp*1000:.2f} L over {Trange:.0f} °C\n"
            f"- Expansion vessel **≈ {V_vessel_L:.2f} L** (bladder)")
        st.markdown("**Tubes, plates, fittings**")
        st.markdown(
            f"- {d['n_tubes']} × {d['tube_od']*1000:.0f} mm {tube_mat} "
            f"tube, total **{tube_L:.1f} m** ({m_tubes:.1f} kg)\n"
            f"- Plates **{m_plates:.1f} kg** aluminium\n"
            f"- Manifold, quick-connects, seals, sensors")

    if ext:
        st.plotly_chart(_S.schematic_oil_loop(d, g, xp),
                        width='stretch', key="sy_oilloop")
        st.caption(
            "**Impact of taking the pump outside.** Thermally it is the "
            "same machine - the pump sets the through-pack velocity, and "
            "that velocity thins the two oil films exactly as stirring or "
            "guided flow does, so h and the temperatures match the "
            "equivalent internal option at the same u. What changes is "
            "practical: the pump is serviceable without opening the "
            "flooded box (a failed internal stirrer means draining "
            f"{V_oil_L:.0f} L; a failed external pump is two valves and "
            "four bolts), the loop self-purges air to the reservoir, and "
            "speed control is a simple external drive. The price is the "
            "pipework: it adds Δp (usually more than the pack itself - "
            "size the bore so the pipe velocity stays near 1 m/s), a "
            "little coolant volume, two bulkhead penetrations that must "
            "seal for life, and priming/de-aeration at commissioning.")

    # =============================================================== #
    st.markdown("##### 3 · Bill of materials (this build)")
    def gpair(mat_kg, cu, al):
        return cu if mat_kg else al
    rho_cost = cu_kg if tube_mat == "Copper" else al_kg
    V_cool_tot = V_oil_L + V_pipe_L
    rows = [
        ("Dielectric ester coolant", f"{V_cool_tot:.1f} L"
         + (" (incl. loop)" if ext else ""), 1,
         oil_L * V_cool_tot, m_oil + V_pipe_L * fl["rho"] / 1000,
         "esters (rapeseed/synthetic); post-PFAS default"),
        (f"{tube_mat} tubes", f"{tube_L:.1f} m × "
         f"{d['tube_od']*1000:.0f} mm", d["n_tubes"],
         rho_cost * m_tubes * fab, m_tubes, "drawn tube, brazed"),
        ("Aluminium cooling plates",
         f"{m_plates:.1f} kg" if m_plates else "none (no serpentine)",
         (g["n_rows"] - 1) if m_plates else 0,
         al_kg * m_plates * fab, m_plates,
         "conduction fin + tube carrier"),
        ("Water circulator / pump", f"{P_pump_w:.0f} W, {head_m:.1f} m",
         1, pump_gbp, 0.4, "brushless, sealed"),
        ("Axial propeller + shroud + motor",
         (f"{P_oil:.1f} W el" if d["circ"].startswith("Bottom") else
          "n/a"), 1 if d["circ"].startswith("Bottom") else 0,
         (28.0 + P_oil * 1.2) if d["circ"].startswith("Bottom") else 0.0,
         0.6 if d["circ"].startswith("Bottom") else 0.0,
         "under the array, pushes the oil up with buoyancy"),
        ("External oil pump", (f"{P_oil:.0f} W, "
         f"{xp['dp']/1000:.1f} kPa, {xp['Vdot_lpm']:.0f} L/min"
         if xp else "n/a"), 1 if ext else 0,
         (pump_gbp * (1.5 + P_oil / 60.0)) if ext else 0.0,
         0.9 if ext else 0.0,
         "outside the pack - serviceable without opening it"),
        ("Oil pipework + bulkheads", (f"{d.get('pipe_len',2.5):.1f} m x "
         f"{d.get('pipe_id',0.019)*1000:.0f} mm + 2 fittings"
         if ext else "n/a"), 1 if ext else 0,
         (7.0 * d.get("pipe_len", 2.5) + 36.0) if ext else 0.0,
         (0.5 * d.get("pipe_len", 2.5)) if ext else 0.0,
         "dielectric-rated hose or tube, FKM seals"),
        ("Dielectric→water HX", f"{UA_hx:.0f} W/K, {A_hx:.3f} m²", 1,
         max(60.0, A_hx * 900.0), 0.8 + A_hx * 5, "brazed plate"),
        ("Chiller / dry-cooler", f"{P_chil/1000:.2f} kW el", 1,
         chil_kW * P_chil / 1000.0, 3.0 + Q_reject / 1000.0,
         "compressor or fan-coil"),
        ("Expansion vessel", f"{V_vessel_L:.2f} L bladder", 1,
         25.0 + 6 * V_vessel_L, 0.5, "accommodates thermal swing"),
        ("Reservoir + filter", "~%.0f%% of fill" % 15, 1, 45.0, 0.8,
         "de-aeration, particulate filter"),
        ("Manifold + quick-connects", "parallel inlet/outlet", 1,
         40.0, 0.6, "dielectric-rated seals (FKM)"),
        ("Sensors + controller", "3× RTD, 1× flow, 1× level", 1,
         85.0, 0.3, "BMS thermal interlocks"),
        ("Seals, hose, misc.", "FKM/EPDM as compatible", 1, 35.0, 0.5,
         "fluid-compatibility critical"),
    ]
    rows = [r for r in rows if r[2] > 0]
    bom = _pd.DataFrame(rows, columns=[
        "Item", "Spec / sizing", "Qty", "Cost £", "Mass kg", "Notes"])
    tot_cost = bom["Cost £"].sum(); tot_mass = bom["Mass kg"].sum()
    cell_cost = Ecap * 79.0 / 0.79   # £ at 79 USD/kWh
    st.dataframe(bom.style.format({"Cost £": "£{:.0f}",
                                   "Mass kg": "{:.1f}"}),
                 width='stretch', hide_index=True)
    mtot = st.columns(3)
    mtot[0].metric("Thermal-system cost", f"£{tot_cost:.0f}",
                   f"{100*tot_cost/max(cell_cost,1):.0f}% of cell cost")
    mtot[1].metric("Thermal-system mass", f"{tot_mass:.1f} kg",
                   f"{100*tot_mass/masses['m_pack']:.0f}% of pack")
    mtot[2].metric("For reference, cells", f"£{cell_cost:,.0f}",
                   f"{Ecap:.1f} kWh @ 79 USD/kWh")
    oil_cost = oil_L * V_oil_L
    st.caption(
        f"Costs are order-of-magnitude engineering estimates from the "
        f"editable assumptions, not quotes. The honest headline the "
        f"numbers give: **the coolant volume dominates** - at "
        f"£{oil_L:.0f}/L the dielectric alone is £{oil_cost:.0f} "
        f"(~{100*oil_cost/max(cell_cost,1):.0f}% of cell cost), and the "
        f"full thermal system is ~{100*tot_cost/max(cell_cost,1):.0f}%. "
        f"That flooded-volume cost is single-phase immersion's real "
        f"economic penalty against a cold plate, and it is far worse for "
        f"two-phase (the fluoroketone is 7-13× the ester price). It is "
        f"also why reduced-fill and guided-channel designs matter: less "
        f"trapped fluid, smaller bill.")

    # =============================================================== #
    st.markdown("##### 4 · How the approaches compare")
    comp = _pd.DataFrame([
        ["Forced air", "10-30", "air (free)", "very low", "fan + ducting",
         "low", "highest (weak h, big ΔT)", "mature",
         "cannot hold fast-charge; large gradients"],
        ["Cold plate (indirect)", "80-150", "water-glycol", "low-med",
         "plates, TIM, pump, HX", "medium", "low", "mature (most EVs)",
         "contact resistance; cells cooled on one face"],
        ["Single-phase immersion (this)", "150-400", "dielectric ester",
         "low-med", "tubes/plates, pump, HX, vessel", "medium-high",
         "low-med", "emerging (AMG HPB80)",
         "oil films dominate; needs circulation"],
        ["Two-phase immersion", "1000-3000+", "fluoroketone / low-GWP "
         "refrigerant", "very low (passive)",
         "condenser, vapor space, sealed vessel", "high",
         "med-high (fluid + condenser)", "early / niche",
         "fluid cost & GWP; vapor management"],
    ], columns=["Architecture", "Surface h [W/m²·K]", "Coolant",
                "Pumping power", "Key parts", "Complexity",
                "Relative cost", "Maturity", "Main limitation"])
    st.dataframe(comp, width='stretch', hide_index=True)
    st.caption("The surface h column is the whole story: air struggles to "
               "reach 30, this single-phase immersion pack works in the "
               "low hundreds (and the two oil films still gate it), while "
               "boiling in two-phase reaches thousands because latent heat "
               "carries the load near-isothermally. Capability rises left "
               "to right; so does cost and complexity.")

    # =============================================================== #
    st.markdown("##### 5 · Two-phase immersion: what it takes")
    st.markdown(
        "Two-phase cooling replaces the two stagnant oil films with "
        "**boiling**. Where single-phase convection gives "
        "$q = hA\\Delta T$ with $h$ in the hundreds, nucleate boiling "
        "gives $q = h_b A\\Delta T$ with $h_b$ in the thousands, because "
        "each bubble carries away the latent heat "
        "$q = \\dot m_\\mathrm{vap}\\,h_{fg}$ at a nearly constant "
        "saturation temperature. The cells sit in a bath that boils at "
        "their target temperature, vapour rises to a condenser, and "
        "liquid returns - often with **no pump at all** (a passive "
        "thermosiphon loop). That is the appeal: the highest heat flux "
        "and the lowest parasitic power, at the same time.")
    tp1, tp2 = st.columns(2)
    with tp1:
        st.markdown("**Coolant - the whole design hinges on it**")
        st.markdown(
            "The saturation (boiling) temperature *is* the cell "
            "temperature, so the fluid is chosen to boil near 34-50 °C at "
            "roughly atmospheric pressure, be a strong dielectric, and be "
            "materials-benign:\n"
            "- **Fluoroketone FK-5-1-12** (3M Novec 649/1230): boils "
            "~49 °C, non-flammable, low GWP (~1), the current front-"
            "runner; costly.\n"
            "- **Hydrofluoroethers** (Novec 7000 ~34 °C, 7100 ~61 °C): "
            "good dielectrics but **PFAS** - 3M is exiting production, so "
            "avoid for new programmes.\n"
            "- **Low-GWP refrigerants** R-1233zd(E) (~18 °C) and "
            "R-1336mzz(Z) (~33 °C): excellent boiling, but volatile and "
            "pressure-managed.\n"
            "- Hydrocarbons/esters can two-phase too but are flammable or "
            "boil too high.")
        st.markdown("**Materials compatibility**")
        st.markdown(
            "- Seals: **FKM/FFKM** (Viton) preferred; EPDM and some "
            "elastomers swell.\n"
            "- Avoid polycarbonate and some engineering plastics with "
            "fluoroketones; PTFE/PEEK are safe.\n"
            "- Metals (Al, Cu, steel) are fine; the enemy is water "
            "ingress (hydrolysis) and non-condensable gases.")
    with tp2:
        st.markdown("**Extra items versus single-phase**")
        st.markdown(
            "- A **condenser** (vapour → liquid), water- or air-cooled, "
            "sized on the full duty.\n"
            "- A sealed, **vapour-tight enclosure** with a vapour plenum "
            "above the liquid.\n"
            "- **Pressure management**: relief valve, and a "
            "non-condensable-gas purge; the box runs near the fluid's "
            "saturation pressure.\n"
            "- A **sight glass / level** and fill port; the fluid is "
            "expensive, so leaks matter.\n"
            "- Often you **delete the pump** - the phase change drives "
            "the loop.")
        st.markdown("**When two-phase wins, and when it does not**")
        st.markdown(
            "- **Wins**: fast charge (>4C), aerospace/eVTOL and other "
            "power- and weight-critical duties, and anywhere near-"
            "isothermal cells or minimal parasitic power matter most.\n"
            "- **Loses**: cost-sensitive volume production, GWP/PFAS "
            "regulation exposure, and programmes that need proven fluid "
            "supply chains today - where single-phase ester immersion or "
            "a cold plate is the pragmatic choice.")
    tp_bom = _pd.DataFrame([
        ["Fluoroketone coolant", f"~{V_oil_L:.0f} L", "£40-120/L",
         f"£{V_oil_L*70:,.0f} (≈7-13× the ester)"],
        ["Condenser", f"{Q_reject/1000:.1f} kW", "-",
         "adds mass + cost, offsets pump saving"],
        ["Vapour-tight enclosure", "sealed + plenum", "premium",
         "heavier than a vented box"],
        ["Pump", "often none", "-", "passive thermosiphon: −£/−W"],
    ], columns=["Two-phase item", "Sizing", "Unit", "Cost impact"])
    st.dataframe(tp_bom, width='stretch', hide_index=True)

    # =============================================================== #
    st.markdown("##### 6 · Where each approach is applied")
    st.markdown(
        "- **Mainstream EVs** - cold plates with water-glycol: cheapest "
        "compliant path at pack scale; most production packs today.\n"
        "- **Performance / fast-charge EVs and motorsport** - single-"
        "phase dielectric immersion (Mercedes-AMG HPB80, this class): "
        "uniform temperatures and high sustained C-rate, weight "
        "acceptable.\n"
        "- **eVTOL / aerospace** - two-phase or single-phase immersion: "
        "power density and near-isothermal cells outweigh fluid cost; "
        "weight and safety are paramount.\n"
        "- **Grid / BESS** - forced air or cold plate: cost and "
        "simplicity dominate; energy, not power, is the driver.\n"
        "- **Data-centre servers** (the analogue that de-risked the "
        "supply chain) - two-phase immersion baths for very high heat "
        "flux; the fluids and hardware come straight from that "
        "industry.\n"
        "- **Defence / directed-energy and pulsed loads** - two-phase, "
        "for the transient flux and isothermality.\n\n"
        "For this pack, the model puts you in the single-phase immersion "
        "band: the serpentine plate-channel architecture buys most of the "
        "uniformity of full immersion at a fraction of the pumping power, "
        "and the two-phase step is the reserve you reach for only if the "
        "duty climbs past what circulation and the water film can "
        "carry.")


def improve_core(d, g, fl, masses, res, Q_duty, C_steady, cool_df):
    if True:
        if True:
            st.markdown("Auto-diagnosis of **this** design at **this** duty, then the "
                        "cheapest fixes ranked by what they actually buy.")
            for level, msg in diagnose(d, g, fl, masses, res, d["T_limit"], Q_duty):
                {"info": st.info, "do": st.success, "warn": st.warning,
                 "bad": st.error}[level](msg)
            st.markdown("---")
            st.markdown("**Sensitivity: change one thing, what happens to steady cell "
                        "temperature?**")
            with st.spinner("Re-solving perturbed designs..."):
                sens, base_T = sensitivity(d, g, fl, cool_df, d["T_amb"], C_steady, d["T_limit"])
            figSe = go.Figure(go.Bar(
                y=sens["change"], x=sens["dT"], orientation="h",
                marker_color=np.where(sens["dT"] < 0, "#10B981", "#EF4444"),
                text=[f"{v:+.1f} °C" for v in sens["dT"]], textposition="outside"))
            figSe.add_vline(x=0, line_color=INK)
            figSe.update_layout(height=380, title=f"Change in steady cell T from "
                                f"{base_T:.1f} °C at {C_steady:.2f}C rms (left = cooler)",
                                xaxis_title="ΔT_cell [°C]", plot_bgcolor="rgba(255,255,255,0)",
                                paper_bgcolor="rgba(0,0,0,0)",
                                margin=dict(l=10, r=60, t=40, b=10))
            st.plotly_chart(figSe, use_container_width=True)

            with st.expander("Goal-seek: cheapest single lever to hit "
                             "the limit", expanded=True):
                gl = st.selectbox("Lever", list(GOAL_LEVERS), key="gs_lever")
                if st.button("Solve for the limit", key="gs_go"):
                    key, lo_, hi_, isint = GOAL_LEVERS[gl]
                    inv = key == "T_water_in"
                    with st.spinner("Bisecting..."):
                        xstar, Tst = goal_seek(d, fl, key, lo_, hi_, isint,
                                               d["T_amb"], C_steady,
                                               d["T_limit"], inv)
                    st.session_state["gs_result"] = (gl, xstar, Tst,
                                                     float(C_steady))
                if "gs_result" in st.session_state:
                    gl_, x_, T_, C_ = st.session_state["gs_result"]
                    key_ = GOAL_LEVERS[gl_][0]
                    if x_ is None:
                        st.error(f"{gl_} alone cannot reach "
                                 f"{d['T_limit']:.0f} °C (best "
                                 f"~{T_:.1f} °C at {C_:.2f}C). Combine "
                                 "levers - see the scale-to-C table above.")
                    else:
                        shown = x_ * 1000 if key_ == "pitch" else x_
                        st.success(f"{gl_} = **{shown:.2f}** gives "
                                   f"{T_:.1f} °C at {C_:.2f}C rms. Set it "
                                   "in Design to adopt.")
                else:
                    st.caption("Pick a lever and press solve - the result "
                               "stays here until you solve another.")

            with st.expander("Two-lever sweep (heatmap)", expanded=True):
                AXES = {"Water flow [L/min]": ("flow_lpm", 2.0, 40.0),
                        "Number of tubes": ("n_tubes", 4, 40),
                        "Stirring [m/s]": ("u_oil", 0.0, 0.12),
                        "Water inlet [°C]": ("T_water_in", 5.0, 35.0),
                        "Cell pitch [mm]": ("pitch", 0.023, 0.033)}
                cx, cy = st.columns(2)
                ax_x = cx.selectbox("X axis", list(AXES), 0, key="sw_x")
                ax_y = cy.selectbox("Y axis", list(AXES), 1, key="sw_y")
                if ax_x == ax_y:
                    st.info("Pick two different levers to see the map.")
                else:
                    kx, x0, x1 = AXES[ax_x]; ky, y0, y1 = AXES[ax_y]
                    with st.spinner("Sweeping 7 x 7 designs (cached)..."):
                        xs, ys, Z = _sweep_cached(
                            json.dumps(d, sort_keys=True, default=float),
                            kx, float(x0), float(x1), ky, float(y0),
                            float(y1), float(round(C_steady, 3)))
                    xs_d = xs * 1000 if kx == "pitch" else xs
                    ys_d = ys * 1000 if ky == "pitch" else ys
                    figH = go.Figure(go.Heatmap(
                        x=xs_d, y=ys_d, z=Z, colorscale="RdYlBu_r",
                        colorbar=dict(title="T_cell [°C]"),
                        hovertemplate=(ax_x + " %{x:.1f}<br>" + ax_y
                                       + " %{y:.1f}<br>T = %{z:.1f} °C"
                                       "<extra></extra>")))
                    figH.add_contour(
                        x=xs_d, y=ys_d, z=Z, showscale=False,
                        contours=dict(start=d["T_limit"], end=d["T_limit"],
                                      coloring="lines"),
                        line=dict(color="black", width=3))
                    figH.add_trace(go.Scatter(
                        x=[d[kx] * (1000 if kx == "pitch" else 1)],
                        y=[d[ky] * (1000 if ky == "pitch" else 1)],
                        mode="markers+text", text=["current"],
                        textposition="top center", showlegend=False,
                        marker=dict(symbol="star", size=16, color="white",
                                    line=dict(color="black", width=1.5))))
                    figH.update_layout(
                        height=460, xaxis_title=ax_x, yaxis_title=ax_y,
                        title=f"Steady cell T at {C_steady:.2f}C rms - "
                              f"black contour = {d['T_limit']:.0f} °C "
                              "limit; star = this design")
                    st.plotly_chart(figH, use_container_width=True,
                                    key="sw_heat", config=PLOTCFG)


def runaway_ui(d, g, fl, masses, res):
    if True:
        if True:
            with st.expander("Thermal-runaway screening (order of magnitude)"):
                cr1, cr2 = st.columns(2)
                with cr1:
                    d["E_tr"] = _w(st.slider, "Heat released per cell [kJ]", "etr", 55.0,
                                   min_value=20.0, max_value=120.0, step=5.0,
                                   help="21700 NMC total ~30-80 kJ depending on SoC")
                    d["frac_oil"] = _w(st.slider, "Fraction into local oil zone", "ftr", 0.6,
                                       min_value=0.2, max_value=1.0, step=0.05)
                with cr2:
                    d["zone_pitches"] = _w(st.slider, "Local zone radius [pitches]", "ztr", 1.5,
                                           min_value=1.0, max_value=3.0, step=0.25)
                    d["vent_L"] = _w(st.slider, "Vent gas at STP [L]", "vtr", 5.0,
                                     min_value=1.0, max_value=15.0, step=0.5)
                rw = runaway_screen(d, g, fl, masses)
                margin = 170.0 - (res["T_b"] + rw["dT_zone"])
                verdict = ("looks containable" if margin > 30
                           else "MARGINAL - add spacing, oil, or interstitial barriers")
                st.markdown(f"""
One cell lets go at {res['T_b']:.0f} °C operating temperature:

* Local zone ({d['zone_pitches']:.1f} pitches): **{rw['V_zone_L']:.1f} L of oil +
  {rw['n_in']:.0f} neighbour cells** -> zone rise **{rw['dT_zone']:.0f} °C**, i.e.
  neighbours reach ~**{res['T_b']+rw['dT_zone']:.0f} °C** vs a ~170-200 °C trigger
  (margin {margin:+.0f} °C, before venting jets - {verdict}).
* Spread over the whole pack it is only **{rw['dT_bulk']:.1f} °C** - the flooded pack's
  big argument.
* {d['vent_L']:.0f} L of vent gas into the {rw['V_hs_L']:.0f} L headspace at ~380 K:
  **~{rw['P_bar_g']:.1f} bar gauge** - size the burst disc well below the lid's rating and
  expect oil ejection through it.

Screening numbers only: vent jets, ejecta and local boiling are not modelled.""")
    else:
        st.sidebar.info("Student version: Decide tab hidden.")

def bench_wang_tab():
    if True:
        st.markdown("The same correlations and two-node network, applied to the **exact rig "
                    "of Wang et al. (2023)**: six prismatic dummy cells (148 x 97 x 27 mm) in "
                    "transformer oil, four 6 mm copper tubes at the top, 5 °C water at "
                    "17.1 mL/s, 2C for 1800 s, 25 °C ambient. Fixed heater power, no "
                    "DCIR(T), calibration factor 1.0.")
        bm = benchmark_wang()
        st.dataframe(pd.DataFrame(bm["rows"], columns=["Quantity", "This app", "Paper"]),
                     hide_index=True, use_container_width=True)
        figB = go.Figure()
        figB.add_trace(go.Scatter(x=bm["t"] / 60, y=bm["T_b"], name="Cell (this app)",
                                  line=dict(color="#EF4444", width=3)))
        figB.add_trace(go.Scatter(x=bm["t"] / 60, y=bm["T_il"], name="Oil (this app)",
                                  line=dict(color=ACCENT, width=3)))
        figB.add_trace(go.Scatter(x=[30], y=[32.3], mode="markers",
                                  name="Paper: cell at 1800 s",
                                  marker=dict(color="#B91C1C", size=12, symbol="x")))
        figB.add_trace(go.Scatter(x=[30], y=[28.5], mode="markers",
                                  name="Paper: oil at 1800 s",
                                  marker=dict(color="#D97706", size=12, symbol="x")))
        figB.update_layout(height=380, xaxis_title="Time [min]",
                           yaxis_title="Temperature [°C]",
                           title="Transient rebuild of the paper's 2C experiment",
                           plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",)
        st.plotly_chart(figB, use_container_width=True)
        st.markdown("""
**How to read the agreement.** Resistances land within ~30-50% and end temperatures within
~1-2 K, with zero tuning. Known residuals: the paper's oil conductivity is internally
inconsistent (Table 3 in kelvin gives 0.30 W/mK, matching their Table 7, vs ~0.13 for real
oil); their h values are back-calculated through a lumped model; and the vertical-plate
correlation is slightly conservative. Treat this app's oil-side h as honest to ~+/-30% -
or pin it with the calibration fit in Results once your own rig data exists.""")

def bench_prod_tab(masses, Cmax):
    if True:

        st.markdown("---")
        st.markdown("**Production packs: where this design sits.** Teardown and "
                    "certification data (OEMs do not publish pack masses); peak power is "
                    "the vehicle rating, a proxy for the battery-side limit.")
        prows = []
        for p_ in PRODUCTION_PACKS:
            r_ = dict(p_)
            if r_["Pack"] == "This design":
                r_.update(kWh=round(masses["E_kwh"], 1), kg=round(masses["m_pack"]),
                          Whkg=round(masses["whkg_pack"]), WhL=round(masses["whl_pack"]),
                          kWpk=round(Cmax * masses["E_kwh"], 1))
                r_["Cooling"] += f" (continuous {Cmax:.1f}C; peak not rated)"
            prows.append(r_)
        pdf_ = pd.DataFrame(prows)
        pdf_["kW/kg"] = (pdf_["kWpk"] / pdf_["kg"]).round(2)
        st.dataframe(pdf_[["Pack", "kWh", "kg", "Whkg", "kW/kg", "Cooling"]],
                     hide_index=True, use_container_width=True)
        figPk = go.Figure()
        for _, r_ in pdf_.iterrows():
            figPk.add_trace(go.Scatter(
                x=[r_["Whkg"]], y=[r_["kW/kg"]], mode="markers+text", text=[r_["Pack"]],
                textposition="top center", showlegend=False,
                marker=dict(size=16 if r_["Pack"] == "This design" else 11,
                            color=ACCENT if r_["Pack"] == "This design" else "#64748B")))
        figPk.update_layout(height=420, xaxis_title="Pack energy density [Wh/kg]",
                            yaxis_title="Peak power density [kW/kg] (vehicle rating)",
                            title="Energy vs power density - immersion trades energy density "
                                  "for simplicity and abuse tolerance",
                            plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)")
        st.plotly_chart(figPk, use_container_width=True)
        st.caption("Sources, retrieved July 2026: Model 3 2170 = Rickard teardown "
                   "(Teslarati/EVANNEX, 478 kg, 4416 cells, ~75 kWh usable, glycol ribbon "
                   "side cooling). Model Y 4680 = EU certification 447 kg at ~79 kWh gross "
                   "(Electrek, May 2026); Munro teardown confirms side ribbons retained and "
                   "cells structurally bonded (batterydesign.net). Plaid = Munro 181.5 Wh/kg "
                   "at ~99 kWh (mass derived, ~545 kg), 760 kW vehicle peak (EVKX), "
                   "7920 x 18650 in 110S72P (Ingineerix). Model 3 LFP = batterydesign.net "
                   "(438 kg, 55 kWh, bottom cold plate). AMG HPB80 = batterydesign.net. "
                   "This design's peak column shows thermal continuous rating, which "
                   "understates a 30 s peak.")


# ------------------------------------------------------------------ #
#  Smoke test                                                         #
# ------------------------------------------------------------------ #
def smoke():
    cool_df = _read_coolants()
    d = dict(DEFAULTS)
    fl = fluid_dict(cool_df[cool_df["name"] == d["coolant"]].iloc[0])
    g = build_geometry(d)
    res = solve_steady(d, g, fl, 1.0, d["T_amb"], C_rate=d["C1"])
    Q = res["Q_eff"]
    masses = build_masses(d, g, fl, res["fin"])
    Cmax = max_continuous_C(d, g, fl, d["T_amb"], d["T_limit"])
    t, C = duty_profile(d["duty"], d["duration"], d["C1"], d["t1"], d["C2"], d["t2"])
    tr = solve_transient(d, g, fl, masses, d["T_amb"], t, C)
    loop = WATER_LOOP[d["loop_fluid"]]
    print(f"Q@2C(DCIR-T)={Q:.0f} W  T_b={res['T_b']:.1f}  T_core={res['T_core']:.1f}  "
          f"u_ts={res['u_ts']*1000:.2f} mm/s  dT_loop={res['dT_loop']:.1f} °C  "
          f"spread={res['spread']:.1f} °C  Cmax={Cmax:.2f}C")
    print(f"parasitics: pump={water_pump_power(d, g, loop)['P']:.2f} W  "
          f"stir@5cm/s={stirrer_power(d, g, fl, 0.05):.2f} W")
    adf = compare_architectures(d, g, fl, masses, d["T_amb"])
    print(adf[["Architecture", "T_cell", "Parasitic_W", "Thermal_mass_kg"]].round(1)
          .to_string(index=False))
    rw = runaway_screen(d, g, fl, masses)
    print(f"runaway: zone dT={rw['dT_zone']:.0f} °C  bulk dT={rw['dT_bulk']:.1f} °C  "
          f"headspace={rw['P_bar_g']:.1f} bar g")
    layout_figure(d, g)
    # calibration self-test: synthesise 'measured' data at cal=1.35, refit
    dd = dict(d); dd["cal_h"] = 1.35
    t_s = np.arange(0, 1201.0, 2.0); C_s = np.full_like(t_s, d["C1"])
    tr_s = solve_transient(dd, g, fl, masses, d["T_amb"], t_s, C_s)
    tm = t_s[::30]; Tm = np.interp(tm, tr_s["t"], tr_s["T_b"]) + np.random.default_rng(1).normal(0, 0.1, len(tm))
    c_fit, rmse = fit_calibration(d, g, fl, masses, d["T_amb"], t_s, C_s, tm, Tm)
    print(f"calibration self-test: true 1.35 -> fitted {c_fit:.2f} (RMSE {rmse:.2f} °C)")
    bm = benchmark_wang()
    print(f"benchmark: T_b(1800s)={bm['T_b'][-1]:.1f} (paper 32.3)")
    tp = dict(d); tp["tube_plane"] = "Below the cells"
    r_low = solve_steady(tp, g, fl, 1.0, d["T_amb"], C_rate=d["C1"])
    print(f"tubes-below check: u_ts={r_low['u_ts']*1000:.2f} mm/s  T_b={r_low['T_b']:.1f}")
    # ---- v3 checks ----
    print(f"weight: pack={masses['m_pack']:.0f} kg ({masses['whkg_pack']:.0f} Wh/kg, "
          f"{masses['whl_pack']:.0f} Wh/L)  enclosure={masses['m_struct']:.1f} kg "
          f"@{masses['enc']['t_mm']:.1f} mm  busbar={masses['m_bus']:.1f} kg  "
          f"holders={masses['m_holders']:.1f} kg")
    for nm in DRIVE_CYCLES:
        tc, vc, D = cycle_speed(nm)
        dist = np.trapezoid(vc, tc) / 1000
        err = 100 * (dist / D - 1)
        print(f"cycle {nm:<24s} {tc[-1]:5.0f} s  {dist:6.2f} km (official {D}, {err:+.2f}%)")
        assert abs(err) < 0.5, nm
    tc, vc, D = cycle_speed("WLTP Class 3b")
    P_b = vehicle_battery_power(tc, vc, d)
    tr_dc = simulate_pack(d, g, fl, masses, d["T_amb"], dict(kind="P", t=tc, P=P_b))
    print(f"WLTP: peak {tr_dc['C'].max():.2f}C / regen {-tr_dc['C'].min():.2f}C  "
          f"C_rms={tr_dc['C_rms']:.3f}  SoC {tr_dc['soc'][0]:.2f}->{tr_dc['soc'][-1]:.3f}  "
          f"T_b end {tr_dc['T_b'][-1]:.1f}")
    dch = dict(d); dch["soc0"] = 0.20; dch["T_start"] = 25.0
    tr_ch = simulate_pack(dch, g, fl, masses, d["T_amb"],
                          dict(kind="chg", t=np.arange(0.0, 7200.0, 2.0)))
    print(f"CC-CV 25C: SoC 0.20->{tr_ch['soc'][-1]:.3f}  peak chg {-tr_ch['C'].min():.2f}C  "
          f"end chg {-tr_ch['C'][-1]:.3f}C")
    assert tr_ch["soc"][-1] > 0.98 and -tr_ch["C"].min() <= d["C_chg"] + 1e-6
    dcold = dict(dch); dcold["T_start"] = 5.0
    tr_cd = simulate_pack(dcold, g, fl, masses, 5.0,
                          dict(kind="chg", t=np.arange(0.0, 3600.0, 2.0)))
    pk_cold = -tr_cd["C"].min()
    pk_early = -tr_cd["C"][: 40].min()          # first 80 s while still cold
    i08 = int(np.argmax(-tr_cd["C"] > 0.8))
    print(f"CC-CV 5°C start: first-80s chg {pk_early:.2f}C (plating cap "
          f"{plating_frac(5.0)*d['C_chg']:.2f}C at 5 °C); self-heats to 0.8C "
          f"by t={tr_cd['t'][i08]/60:.0f} min; run peak {pk_cold:.2f}C")
    assert pk_early < 0.35 * d["C_chg"] and pk_cold > pk_early
    xs_, Ts_ = goal_seek(d, fl, "flow_lpm", 0.5, 60.0, False, d["T_amb"], 2.0,
                         d["T_limit"], False)
    print(f"goal-seek flow for {d['T_limit']:.0f} °C at 2C: {xs_ if xs_ is None else round(xs_,1)} L/min -> {Ts_:.1f}")
    html = report_html(d, g, fl, res, masses, Cmax, [])
    assert "<html" in html and len(html) > 1500
    adf2 = compare_architectures(d, g, fl, masses, d["T_amb"], 2.0)
    assert len(adf2) == 4
    d46 = dict(d); d46.update(FORMATS["4680"], fmt="4680", Ns=96, Np=2,
                              pitch=0.054, cap_Ah=FORMATS["4680"]["cap"],
                              r_dc=FORMATS["4680"]["rdc"], m_cell=FORMATS["4680"]["mcell"])
    g46 = build_geometry(d46)
    r46 = solve_steady(d46, g46, fl, 1.0, d["T_amb"], C_rate=2.0)
    print(f"4680 sanity: {g46['N']} cells  T_b={r46['T_b']:.1f}  dT_core={r46['dT_core']:.1f} °C")
    assert 20 < res["T_b"] < 90 and 0.1 < Cmax < 12
    assert abs(bm["T_b"][-1] - 32.3) < 4.0, "benchmark drifted"
    assert abs(c_fit - 1.35) < 0.15, "calibration fit drifted"
    assert res["u_ts"] > 1e-4, "thermosiphon dead at defaults"
    # ---- v5 checks ----
    ch = chiller_model(res["Q_w"], d["T_water_in"], d["T_amb"])
    print(f"chiller: duty {res['Q_w']/1000:.2f} kW  COP {ch['COP']:.1f}  "
          f"P_el {ch['P_el']:.0f} W")
    assert 1.5 < ch["COP"] < 12
    Q_bus = (d['C1']*d['cap_Ah']*d['Np'])**2 * busbar_props(d, g)['R']
    fig_s = heat_sankey(res, Q, Q_bus, 0.1, 0.0, ch)
    assert len(fig_s.data) == 1
    sts, totdT = station_list(d, g, fl, res, masses, Q, Q_bus, 0.1, 0.0, ch)
    assert len(sts) == 7 and totdT > 5
    sc = scale_to_C(d, g, fl, masses, d["T_amb"], 4.0)
    print(f"scale-to-4C: T_b={sc['r0']['T_b']:.1f} °C  heat "
          f"{sc['r0']['Q_eff']/1000:.2f} kW  chiller {sc['chiller']['P_el']/1000:.2f} kW el")
    for nm, x, T in sc["fixes"]:
        print(f"   4C single lever {nm:<28s} "
              f"{'insufficient alone' if x is None else f'{x:.2f} -> {T:.1f} °C'}")
    c_sweep_fig(d, g, fl, d["T_amb"], 2.0); setpoint_trade(d, g, fl, 2.0, d["T_amb"])
    fx = thermal_xray_fig(d, g, fl, res, res["Q_eff"])
    assert len(fx.layout.shapes) >= 8
    adf6 = full_arch_study(d, g, fl, masses, 2.0, COST_DEFAULTS)
    assert len(adf6) == 6
    stat = adf6.iloc[0]; part = adf6[adf6["name"].str.startswith("Partial")].iloc[0]
    ind = adf6[adf6["name"].str.startswith("Indirect")].iloc[0]
    air = adf6[adf6["name"].str.startswith("Forced")].iloc[0]
    assert part["T_hot"] > stat["T_hot"] and ind["spread"] > 5.0
    assert air["maxC"] < stat["maxC"] and ind["mass"] < stat["mass"]
    print(f"arch study: static maxC {stat['maxC']:.2f} | indirect "
          f"{ind['T_hot']:.1f} °C hot, spread {ind['spread']:.1f} (fails 5) "
          f"| air maxC {air['maxC']:.2f} | partial hot {part['T_hot']:.1f}")
    figs = dict(sankey=fig_s, ladder=waterfall_chart(res, Q, d),
                resist=resistance_chart(res)[0], transient=fig_s,
                csweep=fig_s, setpoint=fig_s, pack3d=fig_s)
    spec_c = dict(kind="C", t=t, C=C)
    tr_s = simulate_pack(d, g, fl, masses, d["T_amb"], spec_c)
    figs["archs"] = fig_s
    secs = report_sections(d, g, fl, res, masses, tr_s, spec_c, Cmax, d["C1"],
                           Q, Q_bus, 0.1, 0.0, ch, figs, arch_df=adf6,
                           cost=COST_DEFAULTS)
    meta_s = dict(version=APP_VERSION, date="test", spec="smoke",
                  ok=True, T_gov=res["T_b"], T_limit=d["T_limit"],
                  kwh=masses["E_kwh"], mass=masses["m_pack"], Cmax=Cmax,
                  chil_el=ch["P_el"])
    html = export_report_html(secs, figs, meta_s)
    deep = sankey_deep_dive(d, g, fl, res, masses, Q, Q_bus, 0.1, 0.0, ch,
                            d["C1"])
    cp_pay_s = dict(
        design=dict(C1=2.0, T_amb=25.0, coolant=d["coolant"],
                    fmt=d["fmt"], Ns=d["Ns"], Np=d["Np"],
                    cap_Ah=d["cap_Ah"], r_dc=d["r_dc"],
                    k_dcir=d["k_dcir"], k_rad=d["k_rad"],
                    pitch=d["pitch"], arrangement=d["arrangement"],
                    tube_plane=d.get("tube_plane", "Top of pack"),
                    flow_lpm=10.0, T_water_in=20.0, n_tubes=16,
                    loop_fluid="Water", tube_od=d["tube_od"],
                    tube_wall=d["tube_wall"], tube_mat=d["tube_mat"],
                    fins_on=True, circ0="thermosiphon", u0=0.05,
                    plate_t=0.0015, plate_contact=0.8, T_limit=45.0,
                    limit_core=False, h_ext=d["h_ext"],
                    edge_margin=d["edge_margin"],
                    bottom_gap=d["bottom_gap"],
                    tube_zone=d["tube_zone"], gas_gap=d["gas_gap"],
                    manifold_margin=d["manifold_margin"],
                    passes=d["passes"], end_fraction=d["end_fraction"],
                    holder_block=d.get("holder_block", 0.25),
                    m_holder_g=d["m_holder_g"],
                    struct_mass=d["struct_mass"],
                    sigma_MPa=d["sigma_MPa"], stiff=d["stiff"],
                    p_des_bar=d["p_des_bar"], bus_J=d["bus_J"],
                    v_nom=d["v_nom"], fin_h=d["fin_h"],
                    fin_t=d["fin_t"], fin_p=d["fin_p"], k_fin=205.0,
                    fin_mat=d["fin_mat"]),
        formats={"18650": dict(d=0.0186, h=0.0652, cap=3.0, m=0.047,
                               r=35.0),
                 "21700": dict(d=0.0211, h=0.0703, cap=5.0, m=0.069,
                               r=25.0),
                 "4680": dict(d=0.046, h=0.080, cap=26.0, m=0.355,
                              r=6.0)},
        fluids=[dict(name=d["coolant"], rho=fl["rho"], cp=fl["cp"],
                     k=fl["k"], nu25=fl["nu25"], B=fl["B"],
                     beta=fl["beta"])],
        waters={k_: dict(rho=v_["rho"], cp=v_["cp"], k=v_["k"],
                         mu=v_["mu"]) for k_, v_ in WATER_LOOP.items()},
        base=dict(T_b=res["T_b"]),
        consts=dict(K_loop=5.0, cal=d.get("cal_h", 1.0), KT=K_TUBE,
                    RT=RHO_TUBE))
    if os.environ.get("CP_DUMP"):
        import json as _json
        open("/tmp/cp_payload.json", "w").write(_json.dumps(cp_pay_s))
    cp_html = cockpit_html(cp_pay_s)
    assert len(cp_html) > 25000 and "PACK COCKPIT" in cp_html
    assert "TEMPERATURE LADDER" in cp_html and "MAX-C" in cp_html
    print(f"cockpit: {len(cp_html)//1000} kB component")
    assert len(deep) > 4000 and "first law" in deep and "COP" in deep
    print(f"deep dive: {len(deep)} chars, {deep.count('**')//2} bold terms")
    assert len(secs) == 9 and "<html" in html
    assert "id='toc'" in html and "<table>" in html and APP_VERSION in html
    assert "callout" in html and "IntersectionObserver" in html
    print(f"export: {len(html)//1000} kB, TOC + tables + scrollspy present")
    print(f"report: {len(secs)} sections, HTML {len(html)//1000} kB")
    # ---- v7 checks: serpentine plates + circuit ----
    ds = dict(d, plate_on=True, u_oil=0.05, plate_t=0.0015,
              plate_contact=0.8, n_tubes=max(g["n_rows"] - 1, 1))
    gs = build_geometry(ds)          # geometry must match the design
    rs = solve_steady(ds, gs, fl, 1.0, d["T_amb"], C_rate=d["C1"])
    A_pl, eta_pl, m_pl = plate_fin_area(ds, gs, rs["h_tube"])
    # served-fraction coupling: halving the tubes must halve plate area
    A_half, _, _ = plate_fin_area(
        dict(ds, n_tubes=max((g["n_rows"] - 1) // 2, 1)), gs, rs["h_tube"])
    assert abs(A_half / A_pl - 0.5) < 0.03, "plate served-fraction"
    sp = serpentine_pump(ds, gs, fl, 0.05)
    dTp = sp["dT_est"](rs["Q_eff"])
    print(f"serpentine: T_b {res['T_b']:.1f} -> {rs['T_b']:.1f} °C  "
          f"plates +{A_pl:.1f} m² (eta {eta_pl:.2f}, {m_pl:.1f} kg)  "
          f"pump {sp['P']:.2f} W  channel dT {dTp:.1f} °C  "
          f"spread {res['spread']:.1f} -> {rs['spread']:.1f} °C")
    assert rs["T_b"] < res["T_b"] - 2.5 and sp["P"] < 5.0 and A_pl > 0.7
    assert dTp < 6.0, "parallel channels should carry the heat"
    # ---- v10 checks: tube shapes + external pump loop ----
    for shp, kw in [("Square", dict(tube_w=0.010)),
                    ("Rectangular", dict(tube_w=0.012, tube_h=0.008))]:
        dq = dict(d, tube_shape=shp, **kw)
        gq = build_geometry(dq)
        rq = solve_steady(dq, gq, fl, 1.0, d["T_amb"], C_rate=d["C1"])
        assert gq["tube_sec"]["shape"] == shp and rq["T_b"] < 60
        assert gq["tube_sec"]["fRe"] < 64.0 and             gq["tube_sec"]["lam_nu"] < 3.66
    dx = dict(ds, circ="External pump loop (closed)", pipe_id=0.019,
              pipe_len=2.5)
    xp1 = ext_loop_pump(dx, gs, fl, 0.05)
    xp2 = ext_loop_pump(dict(dx, pipe_id=0.032), gs, fl, 0.05)
    assert 0 < xp1["P"] < 200 and xp2["P"] < xp1["P"],         "bigger bore must cut the ext-pump power"
    assert xp1["dp_pipe"] > xp1["dp_pack"],         "at 19 mm the pipes should dominate the pack channels"
    print(f"tube shapes: square/rect solve; ext loop 19mm {xp1['P']:.1f} W"
          f" -> 32mm {xp2['P']:.1f} W ({xp1['Vdot_lpm']:.0f} L/min)")
    dp_ = dict(d, circ="Bottom propeller (axial, up)", u_oil=0.05)
    rp_ = solve_steady(dp_, g, fl, 1.0, d["T_amb"], C_rate=d["C1"])
    rs_ = solve_steady(dict(d, circ="Open stirring", u_oil=0.05), g, fl,
                       1.0, d["T_amb"], C_rate=d["C1"])
    pp_ = prop_power(dp_, g, fl, 0.05)
    assert rp_["T_b"] < res["T_b"] - 1.0, "propeller must beat passive"
    assert rp_["T_b"] > rs_["T_b"] + 0.3,         "axial film must be weaker than crossflow at the same u"
    assert 0 < pp_["P"] < 20
    print(f"propeller: axial T_b {rp_['T_b']:.1f} vs crossflow "
          f"{rs_['T_b']:.1f} °C at 5 cm/s; P {pp_['P']:.2f} W "
          f"({pp_['Vdot_lpm']:.0f} L/min swept)")
    # ---- v10.4 verification-ladder cases ----
    c1 = case_still_bath(fl, 3.0, 0.021, 0.070, 35.0)
    assert abs(c1["resid"]) < 0.01, "case1 energy balance"
    assert abs(c1["h"] * c1["A"] * c1["dT"] - 3.0) < 0.01
    ca = case_mode_temp(fl, 3.0, 0.021, 0.070, 35.0, 0.05, "axial")
    cc = case_mode_temp(fl, 3.0, 0.021, 0.070, 35.0, 0.05, "cross")
    assert cc["h"] > ca["h"] > c1["h"], "cross > axial > still"
    wm = case_bath_warmup(fl, 1080, 2.3, 63.0, 6.0, 5.0, 25.0, 120.0)
    assert abs(wm["T_ss"] - (25.0 + 1080 * 2.3 / 30.0)) < 1e-6
    td_ = dict(tube_shape="Round", tube_od=0.010, tube_wall=0.001,
               tube_mat="Copper", L_tube_case=0.85)
    ch_ = case_chain(fl, WATER_LOOP["Water"], 3.0, 0.021, 0.070,
                     20.0, 0.6, td_)
    lad = 3.0 * (ch_["R_ot"] + ch_["R_wall"] + ch_["R_in"])
    assert abs(lad - (ch_["T_oil"] - ch_["T_w_mean"])) < 0.05,         "case4 ladder must close"
    pd_ = dict(plate_on=True, plate_t=0.0015, plate_mat="Aluminium",
               plate_contact=0.9, n_tubes=1, manifold_margin=0.0,
               h_cell=0.07)
    ch_p = case_chain(fl, WATER_LOOP["Water"], 3.0, 0.021, 0.070,
                      20.0, 0.6, td_,
                      plate=dict(d=pd_, g=dict(n_rows=2, Lx=0.85)))
    assert ch_p["A_pl"] > 0 and ch_p["T_oil"] < ch_["T_oil"],         "a bonded plate must cool the oil"
    # ---- basic-module FEA export (WP3 report + sealed) ----
    import fea_export as _FX
    _w = fluid_dict(cool_df[cool_df["name"] ==
                            "Deionized water"].iloc[0])
    _pb = fea_p_basic(_w, 0.025, 0.030, 0.005, 0.0, 0.010, 0.300,
                      0.75, 202.4, 2719.0, 871.0, 25.0, 1800.0,
                      60.0, 0.0, 25.0)
    _qv = _pb["Q_cell"] / (_pb["a_cell"] ** 2 * _pb["L_z"])
    assert abs(_qv - 1e5) < 1e-6, "WP3 q_v must be exactly 1e5 W/m3"
    _pb.update(bar_name="aluminium", top_fixed=True, T_top=25.0,
               k_mult=1.0, steady=True, cls="ipl_wp3_report")
    _j1 = _FX.comsol_basic_2d(_pb)
    assert "TemperatureBoundary" in _j1 and "Stationary" in _j1 \
        and "intTop(ht.ntflux)" in _j1 and "Q_cell/L_z" in _j1
    assert "_upload_to_app.txt" in _j1 and "_summary.txt" in _j1
    _ps = dict(_pb, top_fixed=False, steady=False,
               cls="ipl_sealed_check")
    _j2 = _FX.comsol_basic_2d(_ps)
    assert "Transient" in _j2 and "dTdt_pred*t" in _j2 and \
        "TemperatureBoundary" not in _j2
    print(f"fea basic (WP3): q_v {_qv:.0f} W/m3, sealed slope "
          f"{_pb['dTdt_pred']*60:.4f} K/min, both variants OK")
    # ---- WP3 analytics roundtrip: series -> file -> parse -> check
    _x = np.linspace(0, 0.025, 61); _y = np.linspace(0, 0.030, 73)
    _T = wp3_series_field(_x, _y, 0.025, 0.030, 0.005, 0.005, 0.0125,
                          0.010, 1e5, 0.6, 25.0)
    _dy = _y[-1] - _y[-2]
    _Qp = np.trapezoid(-0.6 * (_T[-1, :] - _T[-2, :]) / _dy, _x)
    assert abs(_Qp - 2.5) < 0.12, "series must conserve energy"
    _rows = "\n".join(f"{xv:.9e} {yv:.9e} {_T[j, i]:.9e}"
                       for j, yv in enumerate(_y)
                       for i, xv in enumerate(_x))
    _xu, _yu, _Tp, _inf = parse_comsol_field("% x y T\n" + _rows)
    assert _inf["kind"] == "2d"
    assert _Tp.shape == _T.shape and \
        np.nanmax(np.abs(_Tp - _T)) < 1e-6
    _m = wp3_check_field(_xu, _yu, _Tp, 0.025, 0.030, 0.005, 0.005,
                         0.0125, 0.010, 1e5, 0.6, 25.0)
    assert _m["anchors_ok"] and _m["rms_out"] < 1e-6
    _X, _Y = np.meshgrid(_xu, _yu)
    _mb = wp3_check_field(_xu, _yu, _Tp + 0.3 * np.exp(
        -((_X - 0.02) ** 2 + (_Y - 0.02) ** 2) / 1e-5),
        0.025, 0.030, 0.005, 0.005, 0.0125, 0.010, 1e5, 0.6, 25.0)
    assert _mb["rms_out"] > 0.03, "checker must flag a bad field"
    # ---- 3D: generator markers + 4-column roundtrip with z-check
    _p3 = dict(_pb, cls="ipl_wp3_report_3d")
    _j3 = _FX.comsol_basic_3d(_p3)
    assert "Block" in _j3 and "regulargridz3" in _j3 and \
        "intTop(ht.ntflux) - Q_cell" in _j3 and \
        "throws IOException" in _j3
    _j3b = _FX.comsol_basic_3d(dict(_p3, build_only=True,
                                    cls="ipl_wp3_report_3d_build"))
    assert 'std1").run()' not in _j3b
    _yd = np.linspace(0, 0.30, 9)
    _rows3 = "\n".join(
        f"{xv:.9e} {yv:.9e} {zv:.9e} {_T[j, i]:.9e}"
        for yv in _yd for j, zv in enumerate(_y)
        for i, xv in enumerate(_x))
    _x3, _z3, _Ts3, _i3 = parse_comsol_field("% x y z T\n" + _rows3)
    assert _i3["kind"] == "3d" and _i3["depth_n"] == 9
    assert _i3["zvar"] < 1e-9, "tiled field must be z-invariant"
    assert np.nanmax(np.abs(_Ts3 - _T)) < 1e-6
    _m3 = wp3_check_field(_x3, _z3, _Ts3, 0.025, 0.030, 0.005, 0.005,
                          0.0125, 0.010, 1e5, 0.6, 25.0)
    assert _m3["anchors_ok"] and _m3["rms_out"] < 1e-6
    # perturb one depth plane and require zvar to flag it
    _pert = []
    for yv in _yd:
        for j, zv in enumerate(_y):
            for i, xv in enumerate(_x):
                t = _T[j, i] + (0.2 if abs(yv - _yd[4]) < 1e-12
                                else 0.0)
                _pert.append(f"{xv:.9e} {yv:.9e} {zv:.9e} {t:.9e}")
    _, _, _, _i3b = parse_comsol_field("% h\n" + "\n".join(_pert))
    assert _i3b["zvar"] > 0.19, "depth perturbation must be flagged"
    print(f"wp3 3d: markers OK, z-invariant roundtrip zvar "
          f"{_i3['zvar']:.1e}, perturbed plane flagged "
          f"{_i3b['zvar']:.2f} K")
    # ---- fan rung: closed form vs independent FD; energy identity;
    #      exporter markers in both dimensions
    def _fd(WW, HH, aa, gg, qq, kk, rr, cc, uu, Ti, N=4000):
        yy = np.linspace(0, HH, N); hh = yy[1] - yy[0]
        qb = np.where((yy >= gg) & (yy <= gg + aa),
                      qq * aa / WW, 0.0)
        lo = kk / hh ** 2 + rr * cc * uu / (2 * hh)
        di = -2 * kk / hh ** 2
        up = kk / hh ** 2 - rr * cc * uu / (2 * hh)
        A_ = np.zeros((N, N)); b_ = -qb.copy()
        idx = np.arange(1, N - 1)
        A_[idx, idx - 1] = lo; A_[idx, idx] = di; A_[idx, idx + 1] = up
        A_[0, 0] = 1; b_[0] = Ti
        A_[-1, -1] = 1; A_[-1, -2] = -1; b_[-1] = 0
        return yy, np.linalg.solve(A_, b_)
    for _u in (0.01, 0.0005):
        _yy, _Tfd = _fd(0.025, 0.030, 0.005, 0.010, 1e5, 0.6,
                        997.0, 4180.0, _u, 25.0)
        _Tan, _fi = fan_mean_exact(_yy, 0.025, 0.030, 0.005, 0.005,
                                   0.010, 1e5, 0.6, 997.0,
                                   4180.0, _u, 25.0)
        assert np.max(np.abs(_Tan - _Tfd)) < 1e-3
        _res = (_fi["adv_W_per_m"] + _fi["cond_bottom_W_per_m"]
                - 1e5 * 0.005 * 0.005)
        assert abs(_res) < 1e-9, "fan energy split must be exact"
    _pf = dict(_pb, top_fixed=False, steady=True, flow_mode="fan",
               u_fan=0.0005, T_in=25.0)
    _jf2 = _FX.comsol_basic_2d(dict(_pf, cls="ipl2d"))
    _jf3 = _FX.comsol_basic_3d(dict(_pf, cls="ipl3d"))
    assert "FluidHeatTransferModel" in _jf2 and \
        '{"0", "u_fan", "0"}' in _jf2 and "Tout_pred" in _jf2 \
        and '"ConvectiveOutflow"' in _jf2 \
        and 'feature("temp1")' not in _jf2
    assert '{"0", "0", "u_fan"}' in _jf3 and \
        "Tout_pred" in _jf3 and '"ConvectiveOutflow"' in _jf3
    # rectangular bar + water pipes: markers in both dims
    _pp = dict(_pb, b_w=0.021, b_h=0.070, n_pipes=3,
               d_pipe=0.008, pipe_drop=0.007, h_w=1800.0,
               T_w=20.0)
    _jp2 = _FX.comsol_basic_2d(dict(_pp, cls="iplX"))
    _jp3 = _FX.comsol_basic_3d(dict(_pp, cls="iplX"))
    assert _jp2.count('"Circle"') == 3 and '"Difference"' in _jp2
    assert "BALANCE - must be ~0" in _jp2 and \
        'set("h", "h_w")' in _jp2 and '"b_w", "b_h"' in _jp2
    assert _jp3.count('"Cylinder"') == 3 and \
        '"axistype", "y"' in _jp3 and "BALANCE" in _jp3
    # rectangular fan analytics: closed form vs FD, bw != bh
    def _fdr(WW, HH, bw_, bh_, gg, qq, kk, rr, cc, uu, Ti,
             N=3000):
        yy = np.linspace(0, HH, N); hh = yy[1] - yy[0]
        qb = np.where((yy >= gg) & (yy <= gg + bh_),
                      qq * bw_ / WW, 0.0)
        lo = kk / hh ** 2 + rr * cc * uu / (2 * hh)
        di = -2 * kk / hh ** 2
        up_ = kk / hh ** 2 - rr * cc * uu / (2 * hh)
        A_ = np.zeros((N, N)); b_ = -qb.copy()
        idx = np.arange(1, N - 1)
        A_[idx, idx - 1] = lo; A_[idx, idx] = di
        A_[idx, idx + 1] = up_
        A_[0, 0] = 1; b_[0] = Ti
        A_[-1, -1] = 1; A_[-1, -2] = -1; b_[-1] = 0
        return yy, np.linalg.solve(A_, b_)
    _yyr, _Tfr = _fdr(0.025, 0.030, 0.004, 0.012, 0.008, 1e5,
                      0.6, 997.0, 4180.0, 0.001, 25.0)
    _Tar, _fir = fan_mean_exact(_yyr, 0.025, 0.030, 0.004, 0.012,
                                0.008, 1e5, 0.6, 997.0, 4180.0,
                                0.001, 25.0)
    assert np.max(np.abs(_Tar - _Tfr)) < 2e-3
    _resr = (_fir["adv_W_per_m"] + _fir["cond_bottom_W_per_m"]
             - 1e5 * 0.004 * 0.012)
    assert abs(_resr) < 1e-9, "rect fan energy split must be exact"
    _jfb = _FX.comsol_basic_2d(dict(_pf, cls="ipl2d",
                                    build_only=True))
    assert 'std1").run()' not in _jfb
    _pbat = dict(_pb, heat_mode="battery", I_cell=15.0,
                 R0_cell=0.0055, R1_cell=0.0032)
    for _g in (_FX.comsol_basic_2d, _FX.comsol_basic_3d):
        _jb = _g(dict(_pbat, cls="iplX"))
        assert '"I_cell"' in _jb and \
            'I_cell^2*(R0_cell + R1_cell)' in _jb
    assert '"I_cell"' not in _FX.comsol_basic_2d(
        dict(_pb, cls="iplX"))
    print("fan rung: closed form vs FD < 1 mK at Pe 2084 and 104, "
          "energy split exact, exporter markers OK both dims; "
          "battery-current parameterisation verified")
    # ---- battery card: parser, round-trip, defaults
    import battery_card as _BC
    assert _BC.parse_condition("JP50_25deg_3C_2nd") == (25.0, 3.0, 2)
    assert _BC.parse_condition("X_35deg_1C") == (35.0, 1.0, 1)
    _cd = pd.DataFrame(dict(
        temp_C=[25.0, 25.0, 35.0], rate_C=[1.0, 1.0, 1.0],
        rep=[1, 1, 1], soc_pct=[30.0, 70.0, 50.0],
        dir=["dis", "dis", "dis"], i_pulse_A=[5, 5, 5],
        dur_s=[10, 10, 10], ocv_V=[3.5, 3.9, 3.7],
        r0_ohm=[0.006, 0.005, 0.0036], r1_ohm=[0.003, 0.003,
                                               0.002],
        tau1_s=[20, 22, 18], r2_ohm=[np.nan] * 3,
        tau2_s=[np.nan] * 3, q_pulse_W=[0.2, 0.2, 0.15]))
    _txt = _BC.card_write(_cd, dict(cell="T", capacity_Ah=5.0,
                                    validation=["x: rmse=1"]))
    _c2, _m2 = _BC.card_read(_txt)
    assert len(_c2) == 3 and _m2["capacity_Ah"] == 5.0
    _mdl = _BC.model_from_card(_c2, temp_C=25.0, rate_C=1.0)
    assert abs(float(_mdl["r0"](50.0)) - 0.0055) < 1e-6
    _mdl35 = _BC.model_from_card(_c2, temp_C=34.0, rate_C=1.0)
    assert abs(float(_mdl35["r0"](50.0)) - 0.0036) < 1e-9
    _dm = _BC.default_model()
    assert abs(_BC.q_steady(_dm, 5.0, 50.0) - 25 * 0.03) < 1e-9
    _tc = _BC.comsol_two_col(_c2, "r0_ohm", 25.0, 1.0)
    assert "soc_pct,r0_ohm" in _tc and _tc.count("\n") >= 4
    print("battery card: parser, round-trip, nearest-condition, "
          "default and 2-col export OK")
    # ---- HPPC battery model: ground truth must be recovered
    import battery_hppc as _BH
    _btr = _BH.synth_truth()
    _bout = _BH.run_pipeline(_BH.make_synth_hppc(noise_mv=0.3,
                                                 seed=1),
                             cap_Ah=4.5, order=1)
    _bsub = _bout["tab"].dropna(subset=["r1"])
    _r0e = float(np.nanmedian(np.abs(_bsub["r0"] -
                 _btr["r0"](_bsub["soc"])) /
                 _btr["r0"](_bsub["soc"])))
    _r1e = float(np.nanmedian(np.abs(_bsub["r1"] -
                 _btr["r1"](_bsub["soc"])) /
                 _btr["r1"](_bsub["soc"])))
    _te = float(np.nanmedian(np.abs(_bsub["tau1"] - 25.0) / 25.0))
    _gg = np.linspace(15, 95, 30)
    _oce = float(np.sqrt(np.mean((_bout["model"]["ocv"](_gg) -
                 _btr["ocv"](_gg)) ** 2)) * 1000)
    assert len(_bout["pulses"]) >= 18
    assert _r0e < 0.03 and _r1e < 0.10 and _te < 0.15
    assert _oce < 6 and _bout["rmse_mv"] < 4
    assert _bout["sim"]["ident_rel"] < 5e-3
    print(f"hppc: {len(_bout['pulses'])} pulses, R0 err "
          f"{100*_r0e:.1f}%, R1 {100*_r1e:.1f}%, tau "
          f"{100*_te:.1f}%, OCV {_oce:.1f} mV, sim RMSE "
          f"{_bout['rmse_mv']:.2f} mV, energy law "
          f"{100*_bout['sim']['ident_rel']:.2f}%")
    print(f"wp3 analytics: series peak {_T.max():.2f} degC, "
          f"top-flux {_Qp:.3f} W/m, roundtrip rms "
          f"{_m['rms_out']:.1e}, perturbation flagged "
          f"{_mb['rms_out']:.2f} K")
    print(f"cases: still {c1['T_s']:.1f} °C (h {c1['h']:.0f}) | "
          f"axial h {ca['h']:.0f} < cross h {cc['h']:.0f} | "
          f"bath ss {wm['T_ss']:.0f} °C tau {wm['tau_min']:.0f} min | "
          f"chain oil {ch_['T_oil']:.1f} -> {ch_p['T_oil']:.1f} °C "
          f"with plate (ladder closes {lad:.2f} °C)")
    fc = thermal_circuit_fig(d, g, fl, res, res["Q_eff"])
    assert len(fc.layout.shapes) >= 6
    nu_T_fig(fl)
    bmdf = load_benchmark()
    assert bmdf is not None and len(bmdf) >= 50
    assert bmdf["whkg"].notna().sum() >= 50 and bmdf["whl"].notna().sum() >= 50
    p_kg = 100 * (bmdf["whkg"].dropna() < masses["whkg_pack"]).mean()
    print(f"benchmark db: {len(bmdf)} packs  whkg {bmdf['whkg'].min():.0f}-"
          f"{bmdf['whkg'].max():.0f}  this design at {p_kg:.0f}th pct")
    # ---- zonal plate-channel model gates ----
    zbank = KernelBank(D=0.021, pitch=0.0215,
                       s_grid=(0.0018, 0.0024), n=60, nz=50)
    zd = dict(n_rows=8, n_cols=8, pitch=0.0215, d_cell=0.021,
              h_cell=0.070, cap_Ah=5.0, r_dc=25.0, k_dcir=0.012,
              C=2.0, plate_t=0.0015, plate_contact=0.8, s_nom=0.002,
              T_in=20.0, flow_lpm=10.0, tube_od=0.010,
              tube_wall=0.0008, k_tube=385.0, T_amb=25.0, h_ext=5.0,
              A_case=0.8, nu25=9e-6, B=3200.0, rho=920.0, cp=2000.0,
              k_oil=0.13, beta=7.5e-4, dp_extra=25.0, T_limit=45.0)
    zr = solve_zonal(zd, zbank, nz=8, iters=180)
    assert abs(zr["closure"]) < 0.01, f"zonal closure {zr['closure']}"
    assert 22 < zr["T_max"] < 90 and zr["spread"] < 3.0
    assert zr["T_core_max"] >= zr["T_max"], "core must exceed can"
    # buoyancy referenced to return -> near-zero net head at default
    _Tbar = zr["Tb"].mean(axis=(1, 2)).mean()
    assert _Tbar - zr["T_plen"] < 0.5, "buoyancy head not return-referenced"
    # optional oil->tube bypass is conservative-side (lowers T) and closes
    zr_b = solve_zonal(dict(zd, wetted_tube_frac=0.5), zbank, nz=8,
                       iters=180)
    assert zr_b["T_max"] < zr["T_max"] and zr_b["Q_bare"] > 0
    assert abs(zr_b["closure"]) < 0.01, "bypass broke closure"
    zl = derived_layout(zd, 150.0)
    assert zl["n_tubes"] >= 2 and 0.3 < zl["eta_bay"] < 0.95
    zmc = monte_carlo(zd, zbank, M=5, sigma_s=0.2e-3, nz=8, iters=25)
    assert 0.0 <= zmc["p_exceed"] <= 1.0 and "Tcore" in zmc
    assert zmc["n_exceed"] > 0 or zmc["p_ub95"] is not None
    assert zmc["closure_worst"] < 0.02
    print(f"zonal: Tmax {zr['T_max']:.1f} core {zr['T_core_max']:.1f} "
          f"closure {zr['closure']*100:.2f}% tubes {zr['lay']['n_tubes']} "
          f"eta {zr['lay']['eta_bay']:.2f} h_face "
          f"{zr['lay']['h_face']:.0f} | bypass {zr_b['UA_bare']:.0f}W/K "
          f"-> {zr_b['T_max']:.1f} | MC ok")
    print("SMOKE OK")


# ------------------------------------------------------------------ #
#  v3: Electrical model (OCV, entropic, plating derate, formats)      #
# ------------------------------------------------------------------ #
SOC_GRID = np.array([0, .05, .10, .20, .30, .40, .50, .60, .70, .80, .90, .95, 1.0])
OCV_GRID = np.array([3.00, 3.30, 3.45, 3.55, 3.61, 3.65, 3.69, 3.75, 3.83, 3.93,
                     4.03, 4.09, 4.20])          # generic graphite-NMC
DUDT_GRID = np.array([0.10, 0.05, 0.00, -0.05, -0.09, -0.11, -0.11, -0.10, -0.09,
                      -0.06, -0.03, -0.02, 0.00]) * 1e-3   # V/K, generic NMC shape
PLATE_T = np.array([-10, 0, 5, 10, 15, 20, 25, 45, 60])
PLATE_F = np.array([0.02, 0.08, 0.20, 0.45, 0.70, 0.90, 1.00, 1.00, 0.50])

def ocv(soc):  return float(np.interp(soc, SOC_GRID, OCV_GRID))
def dudt(soc): return float(np.interp(soc, SOC_GRID, DUDT_GRID))
def plating_frac(T_C): return float(np.interp(T_C, PLATE_T, PLATE_F))

# format presets: capacity/DCIR/mass from teardown literature
# (4680: About:Energy teardown 86.7 Wh -> ~24 Ah, 244 Wh/kg -> ~355 g;
#  DCIR ~5-6 mΩ mid-SoC is an estimate from published teardowns)
FORMATS = {
    "18650": dict(d_cell=0.018, h_cell=0.065, cap=3.4, rdc=30.0, mcell=0.048),
    "21700": dict(d_cell=0.021, h_cell=0.070, cap=5.0, rdc=25.0, mcell=0.070),
    "4680":  dict(d_cell=0.046, h_cell=0.080, cap=24.0, rdc=5.5, mcell=0.355),
}

# ------------------------------------------------------------------ #
#  v3: Drive-cycle library                                            #
#  Coarse (t [s], v [km/h]) breakpoints of the official profiles,     #
#  uniformly speed-scaled so the integrated distance matches the      #
#  official figure exactly. Adequate for thermal work (the pack       #
#  filters everything above ~0.01 Hz); swap in the official 1 Hz      #
#  trace via 'Custom CSV' for certification-grade inputs.             #
# ------------------------------------------------------------------ #
def _u(t0):   # one ECE-15 urban unit for NEDC
    return [(t0, 0), (t0+11, 0), (t0+15, 15), (t0+23, 15), (t0+28, 0), (t0+49, 0),
            (t0+54, 32), (t0+85, 32), (t0+96, 0), (t0+117, 0), (t0+122, 35),
            (t0+133, 50), (t0+155, 50), (t0+163, 35), (t0+176, 35), (t0+188, 0), (t0+195, 0)]

DRIVE_CYCLES = {
    "WLTP Class 3b": (23.27, [(0,0),(11,0),(26,25),(40,42),(48,48),(61,25),(96,0),(122,0),
        (135,30),(160,45),(175,56.5),(201,50),(240,25),(260,0),(300,0),(316,30),(345,49),
        (375,40),(390,20),(420,0),(465,0),(480,28),(511,46),(536,56.5),(556,45),(575,20),
        (589,0),(611,0),(640,40),(670,60),(700,76.6),(735,65),(770,50),(800,60),(840,70),
        (870,55),(900,40),(940,30),(980,45),(1005,25),(1022,0),(1060,0),(1090,40),(1130,65),
        (1180,85),(1230,97.4),(1290,90),(1350,80),(1400,60),(1440,30),(1477,0),(1500,20),
        (1530,60),(1570,90),(1620,110),(1660,125),(1700,131.3),(1740,120),(1765,90),
        (1785,40),(1800,0)]),
    "NEDC": (11.03, _u(0)+_u(195)+_u(390)+_u(585)+[(780,0),(790,35),(805,50),(830,70),
        (870,70),(880,50),(930,50),(940,70),(970,70),(985,100),(1035,100),(1050,120),
        (1090,120),(1120,80),(1160,30),(1180,0)]),
    "UDDS (FTP-75 city)": (12.07, [(0,0),(20,0),(35,40),(60,48),(90,25),(115,0),(125,0),
        (150,55),(185,75),(205,88),(230,91.2),(260,80),(300,60),(330,40),(345,0),(360,0),
        (380,45),(420,55),(450,40),(470,0),(500,30),(530,45),(560,30),(580,0),(610,40),
        (640,50),(670,35),(690,0),(720,40),(760,55),(800,45),(830,0),(860,35),(900,50),
        (930,30),(950,0),(980,40),(1020,55),(1060,45),(1090,25),(1110,0),(1140,35),
        (1180,48),(1220,40),(1260,55),(1300,45),(1340,20),(1369,0)]),
    "HWFET (highway)": (16.45, [(0,0),(30,40),(60,70),(100,80),(160,88),(220,78),(280,86),
        (340,92),(400,96.4),(460,88),(520,80),(580,86),(640,90),(700,75),(740,40),(765,0)]),
    "US06 (aggressive)": (12.89, [(0,0),(15,50),(30,90),(50,108),(70,112),(90,95),(110,60),
        (130,30),(145,0),(160,40),(180,80),(210,105),(250,120),(290,129.2),(330,125),
        (370,118),(410,125),(450,110),(480,80),(510,95),(540,60),(570,30),(600,0)]),
    "Artemis Motorway 130": (28.74, [(0,0),(30,60),(60,100),(100,118),(150,125),(200,131.8),
        (260,122),(320,128),(380,115),(430,125),(490,130),(550,118),(610,108),(660,90),
        (700,110),(760,125),(820,130),(880,120),(930,100),(980,70),(1030,30),(1068,0)]),
}

def cycle_speed(name):
    """1 s resampled speed [m/s], scaled to the official distance."""
    D_km, pts = DRIVE_CYCLES[name]
    tp = np.array([p[0] for p in pts], float)
    vp = np.array([p[1] for p in pts], float) / 3.6
    t = np.arange(0.0, tp[-1] + 1e-9, 1.0)
    v = np.interp(t, tp, vp)
    dist = np.trapezoid(v, t)
    v *= (D_km * 1000.0) / max(dist, 1.0)
    return t, v, D_km

def vehicle_battery_power(t, v, d):
    """Wheel power -> battery power [W] with drivetrain efficiency, capped
    regen, and constant accessory load. Positive = discharge."""
    a = np.gradient(v, t)
    P_wheel = d["veh_m"] * a * v + d["veh_m"] * G * d["Crr"] * v \
              + 0.5 * 1.20 * d["CdA"] * v ** 3
    P = np.where(P_wheel >= 0, P_wheel / max(d["eta_dt"], 0.05) + d["P_acc"],
                 P_wheel * d["eta_rg"] + d["P_acc"])
    return np.clip(P, -d["P_rg"] * 1000.0, None)

def busbar_props(d, g):
    """Total busbar resistance and mass. Auto: series run of length
    Ns x pitch at the design current density; override with R_bus > 0."""
    I_des = d["C1"] * d["cap_Ah"] * d["Np"]
    A_mm2 = max(I_des / d["bus_J"], 10.0)
    L = d["Ns"] * d["pitch"] * 1.15
    R_auto = 1.7e-8 * L / (A_mm2 * 1e-6)          # copper
    R = (d["R_bus"] * 1e-3) if d["R_bus"] > 0 else R_auto
    m = A_mm2 * 1e-6 * L * 8960.0 * 1.25          # + joints/terminals
    return dict(R=R, m=m, A_mm2=A_mm2)

def simulate_pack(d, g, fl, masses, T_amb, spec) -> dict:
    """v3 transient: electro-thermal. Tracks SoC, OCV, terminal voltage,
    CC-CV charge with plating-derated current, entropic heat, busbar heat.
    spec kinds: 'C' (array of C, thermal-only unless track_soc),
    'P' (battery power array, drive cycles), 'chg' (CC-CV), 'cyc' (cycling)."""
    loop = WATER_LOOP[d["loop_fluid"]]
    mdot = d["flow_lpm"] / 60.0 * loop["rho"] / 1000.0
    C_b, C_il = masses["C_batt"], masses["C_oil"]
    bus = busbar_props(d, g)
    t_arr = spec["t"]
    n = len(t_arr)
    T_b = np.full(n, d["T_start"]); T_il = np.full(n, d["T_start"])
    T_core = np.full(n, d["T_start"])
    soc = np.full(n, d["soc0"]); I_c = np.zeros(n); V_c = np.zeros(n)
    C_tr = np.zeros(n); Q_tr = np.zeros(n)
    rep = solve_steady(d, g, fl, 1.0, T_amb, C_rate=max(d["C1"], 0.5))
    R_wall, R_in, R_atm = rep["R_wall"], rep["R_in"], rep["R_atm"]
    ch = d.get("cal_h", 1.0)
    N, cap = g["N"], d["cap_Ah"]
    mode = "dis"; rest_t = 0.0; cyc_count = 0
    for i in range(1, n):
        dt = t_arr[i] - t_arr[i - 1]
        s = soc[i - 1]; Tb = T_b[i - 1]
        R_cell = r_of_T(d, Tb) * 1e-3
        U = ocv(s)
        # --- current demand per cell (discharge positive) ---
        if spec["kind"] == "C":
            I = spec["C"][i - 1] * cap * (1 if d["dirn"] == "Discharge" else -1)
        elif spec["kind"] == "P":
            P_cell = spec["P"][i - 1] / N
            Rq = d["chg_mult"] * R_cell if P_cell < 0 else R_cell
            disc = U * U - 4.0 * Rq * P_cell
            I = (U - math.sqrt(disc)) / (2.0 * Rq) if disc > 0 else U / (2.0 * Rq)
            if I < 0:   # regen: plating derate
                I = -min(-I, plating_frac(Tb) * d["C_chg"] * cap)
        elif spec["kind"] in ("chg", "cyc"):
            if spec["kind"] == "cyc":
                if mode == "dis" and s <= d["soc_min"]:
                    mode, rest_t = "rest1", 0.0
                elif mode == "rest1":
                    rest_t += dt
                    if rest_t >= d["cyc_rest"]: mode = "chg"
                elif mode == "chg" and s >= 0.999:
                    mode, rest_t, cyc_count = "rest2", 0.0, cyc_count + 1
                elif mode == "rest2":
                    rest_t += dt
                    if rest_t >= d["cyc_rest"]:
                        mode = "dis" if cyc_count < d["n_cyc"] else "done"
            else:
                mode = "chg" if s < 0.999 else "done"
            if mode == "dis":
                I = d["C1"] * cap
            elif mode == "chg":
                R_ch = d["chg_mult"] * R_cell
                I_cc = -plating_frac(Tb) * d["C_chg"] * cap
                V_at_cc = U - I_cc * R_ch
                if V_at_cc < d["v_max"]:
                    I = I_cc
                else:                       # CV phase
                    I = -(d["v_max"] - U) / R_ch
                    if -I < d["v_cut"] * cap:
                        I = 0.0
                        if spec["kind"] == "chg": mode = "done"
            else:
                I = 0.0
        # SoC bounds
        if spec["kind"] != "C" or d["track_soc"]:
            if (I > 0 and s <= 0.002) or (I < 0 and s >= 0.999 and spec["kind"] == "P"):
                I = 0.0
            soc[i] = min(max(s - I * dt / (3600.0 * cap), 0.0), 1.0)
        else:
            soc[i] = s
        Rq = d["chg_mult"] * R_cell if I < 0 else R_cell
        q_cell = I * I * Rq
        if d["entropic"]:
            q_cell -= I * (Tb + 273.15) * dudt(s)
        Q = N * q_cell + (I * d["Np"]) ** 2 * bus["R"]
        # --- thermal step (v2 core) ---
        u_ts, _ = thermosiphon_u(d, g, fl, max(Q, 1.0), 0.5 * (Tb + T_il[i-1]))
        u_eff = max(d["u_oil"], u_ts)
        cellf = h_cell_side(fl, Tb, T_il[i-1], d["h_cell"], d["d_cell"], g["gap_mm"], u_eff,
                            axial=d.get("circ", "").startswith("Bottom"))
        T_wall_est = T_il[i-1] - 0.6 * (T_il[i-1] - d["T_water_in"])
        tubf = h_tube_side(fl, T_il[i-1], T_wall_est, d["tube_od"], u_eff)
        R_b = 1.0 / max(ch * cellf["h"] * g["A_cells"], 1e-9)
        R_ot = 1.0 / max(ch * tubf["h"] * rep["A_oilside"], 1e-9)
        Q_bi = (Tb - T_il[i-1]) / R_b
        Q_w = (T_il[i-1] - (d["T_water_in"] + 0.5 * max(Q_bi, 0) / max(mdot * loop["cp"], 1e-9))) \
              / (R_ot + R_wall + R_in)
        Q_a = (T_il[i-1] - T_amb) / R_atm
        T_b[i] = Tb + dt * (Q - Q_bi) / C_b
        T_il[i] = T_il[i-1] + dt * (Q_bi - Q_w - Q_a) / C_il
        T_core[i] = T_b[i] + (Q / N) * r_core(d)
        I_c[i] = I; V_c[i] = U - I * Rq; C_tr[i] = I / cap; Q_tr[i] = Q
    return dict(t=t_arr, T_b=T_b, T_il=T_il, T_core=T_core, soc=soc, I=I_c, V=V_c,
                C=C_tr, Q=Q_tr, C_rms=float(np.sqrt(np.mean(C_tr ** 2))),
                bus=bus, cycles_done=cyc_count)

# ------------------------------------------------------------------ #
#  v3: goal-seek, report export, production benchmarks                #
# ------------------------------------------------------------------ #
GOAL_LEVERS = {"Total water flow [L/min]": ("flow_lpm", 0.5, 60.0, False),
               "Number of tubes": ("n_tubes", 1, 60, True),
               "Stirring [m/s]": ("u_oil", 0.0, 0.20, False),
               "Water inlet [°C]": ("T_water_in", 0.0, 40.0, False),
               "Cell pitch [mm]": ("pitch", 0.022, 0.035, False)}

def goal_seek(d, fl, lever_key, lo, hi, is_int, T_amb, C_duty, T_target, invert):
    """Bisect one lever so steady cell T meets T_target. invert=True for
    levers where increasing the value makes the pack hotter (T_water_in)."""
    def T_at(x):
        dd = dict(d); dd[lever_key] = int(round(x)) if is_int else float(x)
        gg = build_geometry(dd)
        return solve_steady(dd, gg, fl, 1.0, T_amb, C_rate=C_duty)["T_b"]
    T_lo, T_hi = T_at(lo), T_at(hi)
    cool_end_is_hi = T_hi < T_lo
    if min(T_lo, T_hi) > T_target:
        return None, max(T_lo, T_hi) if invert else min(T_lo, T_hi)
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if (T_at(mid) > T_target) == cool_end_is_hi:
            lo = mid
        else:
            hi = mid
    x = 0.5 * (lo + hi)
    return (int(math.ceil(x)) if is_int else x), T_at(x)

# Production packs for the benchmark table. Sources (retrieved July 2026):
# teardown/certification data, not OEM datasheets - Tesla do not publish
# pack masses. M3 2170: Rickard teardown via Teslarati/EVANNEX (478 kg,
# 4416 cells, ~75 kWh usable). MY 4680: EU certification 447 kg / ~79 kWh
# gross (Electrek May 2026); Munro-derived 445 kg (batterydesign.net).
# Plaid: Munro 181.5 Wh/kg at ~99 kWh -> ~545 kg (derived); 760 kW vehicle
# peak (EVKX). M3 LFP: batterydesign.net 438 kg, 55 kWh, 125 Wh/kg.
# HPB80: batterydesign.net. Peak power = vehicle rating (battery-side
# peaks unpublished).
PRODUCTION_PACKS = [
    dict(Pack="This design", kWh=None, kg=None, Whkg=None, WhL=None, kWpk=None,
         Cooling="static immersion + internal water HX"),
    dict(Pack="Mercedes AMG HPB80", kWh=6.1, kg=89, Whkg=68.5, WhL=None, kWpk=150,
         Cooling="pumped dielectric immersion + external HX, 45 °C set point"),
    dict(Pack="Tesla Model 3 LR (2170)", kWh=75.0, kg=478, Whkg=157, WhL=None, kWpk=377,
         Cooling="glycol ribbon/serpentine side cooling, 4416 cells, 96S46P"),
    dict(Pack="Tesla Model Y (4680 structural)", kWh=79.0, kg=447, Whkg=177, WhL=None,
         kWpk=331, Cooling="glycol side ribbons retained; cells glued in steel tub"),
    dict(Pack="Tesla Model S Plaid (18650)", kWh=99.0, kg=545, Whkg=181.5, WhL=None,
         kWpk=760, Cooling="micro-channel glycol ribbons, 7920 cells, 110S72P"),
    dict(Pack="Tesla Model 3 LFP (CATL prismatic)", kWh=55.0, kg=438, Whkg=125, WhL=None,
         kWpk=239, Cooling="bottom cold plate under prismatic cells"),
]

def report_html(d, g, fl, res, masses, Cmax, figs) -> str:
    rows = [
        ("Pack", f"{masses['E_kwh']:.1f} kWh, {d['Ns']}S{d['Np']}P = {g['N']} x "
                 f"{d['fmt']}, {d['Ns']*d['v_nom']:.0f} V"),
        ("Coolant", f"{fl['name']}, {masses['V_oil_L']:.0f} L / {masses['m_oil']:.0f} kg, "
                    f"stirring {d['u_oil']*100:.1f} cm/s"),
        ("Internal HX", f"{d['n_tubes']} x {d['tube_od']*1000:.0f} mm {d['tube_mat']} tubes, "
                        f"{d['tube_plane']}, fins {'on' if d['fins_on'] else 'off'}"),
        ("Water loop", f"{d['flow_lpm']:.0f} L/min at {d['T_water_in']:.0f} °C"),
        ("Steady at duty", f"can {res['T_b']:.1f} / core {res['T_core']:.1f} °C at "
                           f"C_rms, spread {res['spread']:.1f} °C, "
                           f"thermosiphon {res['u_ts']*1000:.1f} mm/s"),
        ("Capability", f"max continuous {Cmax:.2f}C to {d['T_limit']:.0f} °C "
                       f"({'core' if d['limit_core'] else 'can'})"),
        ("Mass", f"pack {masses['m_pack']:.0f} kg -> {masses['whkg_pack']:.0f} Wh/kg, "
                 f"{masses['whl_pack']:.0f} Wh/L (enclosure {masses['m_struct']:.0f} kg at "
                 f"{masses['enc']['t_mm']:.1f} mm eff.)"),
        ("Calibration", f"oil-film factor {d['cal_h']:.2f}"),
    ]
    tab = "".join(f"<tr><th style='text-align:left;padding:4px 12px 4px 0'>{k}</th>"
                  f"<td>{v}</td></tr>" for k, v in rows)
    parts = "".join(f.to_html(full_html=False, include_plotlyjs=False) for f in figs)
    return f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<title>Immersion Pack Lab report</title>
<style>body{{font-family:Georgia,serif;max-width:960px;margin:2em auto;color:#1F2933}}
h1{{border-bottom:3px solid #F59E0B}}table{{border-collapse:collapse;margin:1em 0}}</style>
</head><body><h1>Immersion Pack Lab - design report</h1>
<p>Generated by the app; model honesty: oil-side h +/-30% unless calibrated.
Wang et al. 2023 benchmark: 33.0 vs 32.3 °C measured.</p>
<table>{tab}</table>{parts}
<p style="font-size:0.85em;color:#666">Two-node network; Churchill-Chu / Churchill-Bernstein /
Hausen-Gnielinski correlations; DCIR(T), entropic heat, thermosiphon and busbar terms included.
Sources: Wang 2023 (J. Energy Storage 62, 106821); coolant_comparison_reviewed.xlsx.</p>
</body></html>"""


# ------------------------------------------------------------------ #
#  v4: 3D pack view                                                   #
# ------------------------------------------------------------------ #
def _add_cyl_z(V, F, I, xc, yc, z0, z1, r, val, n=10):
    """Append a vertical cylinder (side + top cap) to vertex/face lists."""
    b = len(V)
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.column_stack([xc + r * np.cos(th), yc + r * np.sin(th)])
    for x, y in ring: V.append((x, y, z0))
    for x, y in ring: V.append((x, y, z1))
    V.append((xc, yc, z1))
    I.extend([val] * (2 * n + 1))
    for k in range(n):
        k2 = (k + 1) % n
        F.append((b + k, b + k2, b + n + k))
        F.append((b + k2, b + n + k2, b + n + k))
        F.append((b + n + k, b + n + k2, b + 2 * n))

def _add_cyl_x(V, F, I, x0, x1, yc, zc, r, val, n=12):
    """Append a horizontal (along-x) open cylinder."""
    b = len(V)
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.column_stack([yc + r * np.cos(th), zc + r * np.sin(th)])
    for y, z in ring: V.append((x0, y, z))
    for y, z in ring: V.append((x1, y, z))
    I.extend([val] * (2 * n))
    for k in range(n):
        k2 = (k + 1) % n
        F.append((b + k, b + k2, b + n + k))
        F.append((b + k2, b + n + k2, b + n + k))

def _mesh(V, F, I, **kw):
    V = np.array(V); F = np.array(F)
    return go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2],
                     i=F[:, 0], j=F[:, 1], k=F[:, 2],
                     intensity=np.array(I), flatshading=True, **kw)

def pack_3d_figure(d, g, show_oil=True, show_tubes=True, show_box=True,
                   show_fins=True, height=480):
    fig = go.Figure()
    # --- cells, coloured by centre-vs-edge tendency ---
    V, F, I = [], [], []
    z0 = d["bottom_gap"]; z1 = z0 + d["h_cell"]
    cx, cy = g["Lx"] / 2, g["Ly"] / 2
    cnt = 0
    for r_ in range(g["n_rows"]):
        for c_ in range(g["n_cols"]):
            if cnt >= g["N"]: break
            x = d["edge_margin"] + (c_ + 0.5) * d["pitch"] \
                + (d["pitch"] / 2 if (d["arrangement"] == "Hexagonal" and r_ % 2) else 0)
            y = d["edge_margin"] + d["pitch"] / 2 + r_ * g["row_pitch"]
            val = 1.0 - math.hypot(x - cx, y - cy) / math.hypot(cx, cy)
            _add_cyl_z(V, F, I, x, y, z0, z1, d["d_cell"] / 2, val, n=10)
            cnt += 1
    fig.add_trace(_mesh(V, F, I, colorscale="RdYlBu_r", showscale=False,
                        name="cells", lighting=dict(ambient=0.72, diffuse=0.75, specular=0.15)))
    # --- tubes (+ translucent fin envelope) ---
    if show_tubes:
        Vt, Ft, It = [], [], []
        Vf, Ff, If_ = [], [], []
        inter = d.get("tube_plane") == "Interstitial (between rows)"
        x0, x1 = d["manifold_margin"], g["Lx"] - d["manifold_margin"]
        for j_ in range(d["n_tubes"]):
            if inter:
                gaps = max(g["n_rows"] - 1, 1)
                yj = d["edge_margin"] + d["pitch"] / 2 \
                     + (j_ % gaps + 0.5) * g["row_pitch"] * (g["n_rows"] - 1) / gaps
                zj = z0 + d["h_cell"] * (0.3 + 0.4 * ((j_ // gaps) % 2))
            else:
                yj = (j_ + 0.5) * g["Ly"] / d["n_tubes"]
                zj = z1 + d["tube_zone"] / 2
            _add_cyl_x(Vt, Ft, It, x0, x1, yj, zj, d["tube_od"] / 2, 1.0)
            if d["fins_on"] and show_fins:
                _add_cyl_x(Vf, Ff, If_, x0, x1, yj, zj,
                           d["tube_od"] / 2 + d["fin_h"], 1.0, n=10)
        fig.add_trace(_mesh(Vt, Ft, It, colorscale=[[0, "#D97706"], [1, "#D97706"]],
                            showscale=False, name="tubes"))
        if Vf:
            fig.add_trace(_mesh(Vf, Ff, If_, colorscale=[[0, "#A5B4CC"], [1, "#A5B4CC"]],
                                showscale=False, opacity=0.16, name="fin envelope"))
    # --- oil fill ---
    if show_oil:
        fz = g["fill_h"]
        xs = [0, g["Lx"], g["Lx"], 0, 0, g["Lx"], g["Lx"], 0]
        ys = [0, 0, g["Ly"], g["Ly"], 0, 0, g["Ly"], g["Ly"]]
        zs = [0, 0, 0, 0, fz, fz, fz, fz]
        fig.add_trace(go.Mesh3d(x=xs, y=ys, z=zs, alphahull=0, opacity=0.10,
                                color="#F59E0B", name="oil", hoverinfo="skip"))
    # --- enclosure wireframe ---
    if show_box:
        Lx, Ly, Lz = g["Lx"], g["Ly"], g["Lz"]
        E = [((0,0,0),(Lx,0,0)),((Lx,0,0),(Lx,Ly,0)),((Lx,Ly,0),(0,Ly,0)),((0,Ly,0),(0,0,0)),
             ((0,0,Lz),(Lx,0,Lz)),((Lx,0,Lz),(Lx,Ly,Lz)),((Lx,Ly,Lz),(0,Ly,Lz)),((0,Ly,Lz),(0,0,Lz)),
             ((0,0,0),(0,0,Lz)),((Lx,0,0),(Lx,0,Lz)),((Lx,Ly,0),(Lx,Ly,Lz)),((0,Ly,0),(0,Ly,Lz))]
        ex, ey, ez = [], [], []
        for (a, b) in E:
            ex += [a[0], b[0], None]; ey += [a[1], b[1], None]; ez += [a[2], b[2], None]
        fig.add_trace(go.Scatter3d(x=ex, y=ey, z=ez, mode="lines",
                                   line=dict(color=INK, width=3),
                                   name="enclosure", hoverinfo="skip"))
    fig.update_layout(
        height=height, margin=dict(l=0, r=0, t=30, b=0), showlegend=False,
        scene=dict(aspectmode="data",
                   xaxis=dict(visible=False), yaxis=dict(visible=False),
                   zaxis=dict(visible=False),
                   camera=dict(eye=dict(x=1.15, y=-1.35, z=0.95)),
                   bgcolor="rgba(0,0,0,0)"),
        paper_bgcolor="rgba(0,0,0,0)",
        title=f"{g['N']} x {d['fmt']} cells, {d['n_tubes']} tubes "
              f"({d['tube_plane'].lower()}), oil to {g['fill_h']*1000:.0f} mm")
    return fig


def thin_labels(xs, ys, names, logx=False, min_dx=0.10, min_dy=0.12):
    """Keep a label only if it does not crowd an already-kept one; everything
    remains identifiable on hover. Distances are fractions of axis span
    (log10 span when logx)."""
    xs = np.asarray(xs, float); ys = np.asarray(ys, float)
    xv = np.log10(np.maximum(xs, 1e-9)) if logx else xs
    xspan = max(xv.max() - xv.min(), 1e-9); yspan = max(ys.max() - ys.min(), 1e-9)
    order = np.argsort(-ys)                      # label the notable (high) first
    kept, out = [], [""] * len(xs)
    for idx in order:
        if all(abs(xv[idx] - xv[j]) / xspan > min_dx
               or abs(ys[idx] - ys[j]) / yspan > min_dy for j in kept):
            kept.append(idx); out[idx] = names[idx]
    return out


def pack_views_fig(d, g, masses, height=560):
    """Engineering triptych: plan (x-y) on top, side section (y-z, tubes as
    circles with fins and plates) and front section (x-z, tube runs) below.
    Shows oil level, headspace, bottom gap, wall thickness and key dims."""
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=2, cols=2, specs=[[{"colspan": 2}, None],
                                               [{}, {}]],
                        row_heights=[0.52, 0.48], horizontal_spacing=0.07,
                        vertical_spacing=0.14,
                        subplot_titles=("Plan  (x-y)",
                                        "Side section  (y-z)",
                                        "Front section  (x-z)"))
    Lx, Ly, Lz, fz = g["Lx"], g["Ly"], g["Lz"], g["fill_h"]
    tw = masses["enc"]["t_mm"] / 1000
    OIL, CELL, TUBE, FINC, PLATE, WALL = ("rgba(245,158,11,0.20)", "#9AA7B8",
                                          "#D97706", "rgba(165,180,204,0.55)",
                                          "#6366F1", "#0F172A")
    inter = d.get("tube_plane") == "Interstitial (between rows)"
    z_tube = (d["bottom_gap"] + d["h_cell"] * 0.5) if inter else              (d["bottom_gap"] + d["h_cell"] + d["tube_zone"] / 2)
    r_t, r_f = d["tube_od"] / 2, d["tube_od"] / 2 + (d["fin_h"] if d["fins_on"] else 0)

    def wallrect(row, col, W, H):
        fig.add_shape(type="rect", x0=-tw, y0=-tw, x1=W + tw, y1=H + tw,
                      line=dict(color=WALL, width=1.5),
                      fillcolor="rgba(15,23,42,0.05)", row=row, col=col)
        fig.add_shape(type="rect", x0=0, y0=0, x1=W, y1=H,
                      line=dict(color=WALL, width=1), fillcolor="white",
                      row=row, col=col)

    # ---------- plan ----------
    wallrect(1, 1, Lx, Ly)
    fig.add_shape(type="rect", x0=0, y0=0, x1=Lx, y1=Ly, fillcolor=OIL,
                  line=dict(width=0), row=1, col=1)
    step = max(1, g["N"] // 400)          # cap plan markers for speed
    xs, ys = [], []
    cnt = 0
    for r_ in range(g["n_rows"]):
        for c_ in range(g["n_cols"]):
            if cnt >= g["N"]: break
            if cnt % step == 0:
                xs.append(d["edge_margin"] + (c_ + 0.5) * d["pitch"]
                          + (d["pitch"] / 2 if (d["arrangement"] == "Hexagonal"
                                                and r_ % 2) else 0))
                ys.append(d["edge_margin"] + d["pitch"] / 2 + r_ * g["row_pitch"])
            cnt += 1
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", hoverinfo="skip",
                             marker=dict(size=max(3, d["d_cell"] / Lx * 430),
                                         color=CELL, line=dict(width=0)),
                             showlegend=False), row=1, col=1)
    for j in range(d["n_tubes"]):
        yj = ((j % max(g["n_rows"] - 1, 1) + 0.5) * g["row_pitch"]
              * (g["n_rows"] - 1) / max(g["n_rows"] - 1, 1)
              + d["edge_margin"] + d["pitch"] / 2) if inter else              (j + 0.5) * Ly / d["n_tubes"]
        fig.add_shape(type="line", x0=d["manifold_margin"], y0=yj,
                      x1=Lx - d["manifold_margin"], y1=yj,
                      line=dict(color=TUBE, width=2.5), row=1, col=1)

    # ---------- side section (y-z): tubes as circles ----------
    wallrect(2, 1, Ly, Lz)
    fig.add_shape(type="rect", x0=0, y0=0, x1=Ly, y1=fz, fillcolor=OIL,
                  line=dict(width=0), row=2, col=1)
    for r_ in range(g["n_rows"]):
        yc = d["edge_margin"] + d["pitch"] / 2 + r_ * g["row_pitch"]
        fig.add_shape(type="rect", x0=yc - d["d_cell"] / 2, y0=d["bottom_gap"],
                      x1=yc + d["d_cell"] / 2, y1=d["bottom_gap"] + d["h_cell"],
                      line=dict(color="#7C8797", width=0.8), fillcolor=CELL,
                      row=2, col=1)
    if d.get("plate_on"):
        for r_ in range(g["n_rows"] - 1):
            yp = d["edge_margin"] + d["pitch"] / 2 + (r_ + 0.5) * g["row_pitch"]
            fig.add_shape(type="line", x0=yp, y0=d["bottom_gap"],
                          x1=yp, y1=d["bottom_gap"] + d["h_cell"],
                          line=dict(color=PLATE, width=2.5), row=2, col=1)
    for j in range(d["n_tubes"]):
        yj = ((j % max(g["n_rows"] - 1, 1) + 0.5) * g["row_pitch"]
              + d["edge_margin"] + d["pitch"] / 2) if inter else              (j + 0.5) * Ly / d["n_tubes"]
        zj = z_tube if not inter else              d["bottom_gap"] + d["h_cell"] * (0.3 + 0.4 * ((j // max(g["n_rows"]-1,1)) % 2))
        if d["fins_on"]:
            fig.add_shape(type="circle", x0=yj - r_f, y0=zj - r_f,
                          x1=yj + r_f, y1=zj + r_f, fillcolor=FINC,
                          line=dict(width=0), row=2, col=1)
        fig.add_shape(type="circle", x0=yj - r_t, y0=zj - r_t,
                      x1=yj + r_t, y1=zj + r_t, fillcolor=TUBE,
                      line=dict(color="#B45309", width=1), row=2, col=1)

    # ---------- front section (x-z): tube runs ----------
    wallrect(2, 2, Lx, Lz)
    fig.add_shape(type="rect", x0=0, y0=0, x1=Lx, y1=fz, fillcolor=OIL,
                  line=dict(width=0), row=2, col=2)
    fig.add_shape(type="rect", x0=d["edge_margin"], y0=d["bottom_gap"],
                  x1=Lx - d["edge_margin"], y1=d["bottom_gap"] + d["h_cell"],
                  fillcolor=CELL, line=dict(color="#7C8797", width=0.8),
                  opacity=0.85, row=2, col=2)
    band = r_f if d["fins_on"] else r_t
    fig.add_shape(type="rect", x0=d["manifold_margin"], y0=z_tube - band,
                  x1=Lx - d["manifold_margin"], y1=z_tube + band,
                  fillcolor=FINC, line=dict(width=0), row=2, col=2)
    fig.add_shape(type="rect", x0=d["manifold_margin"], y0=z_tube - r_t,
                  x1=Lx - d["manifold_margin"], y1=z_tube + r_t,
                  fillcolor=TUBE, line=dict(width=0), row=2, col=2)

    # dims + oil level annotations
    for rc, W, H in [((1, 1), Lx, Ly), ((2, 1), Ly, Lz), ((2, 2), Lx, Lz)]:
        fig.add_annotation(x=W / 2, y=-tw - 0.012, text=f"{W*1000:.0f} mm",
                           showarrow=False, font=dict(size=10, color="#64748B"),
                           row=rc[0], col=rc[1])
    for col in (1, 2):
        fig.add_annotation(x=0.012, y=fz, text=f"oil {fz*1000:.0f}",
                           showarrow=False, xanchor="left", yanchor="bottom",
                           font=dict(size=10, color="#B45309"),
                           row=2, col=col)
    for rc in [(1, 1), (2, 1), (2, 2)]:
        fig.update_xaxes(visible=False, row=rc[0], col=rc[1])
        fig.update_yaxes(visible=False, scaleanchor="x" if rc == (1, 1) else None,
                         row=rc[0], col=rc[1])
    fig.update_yaxes(scaleanchor="x", scaleratio=1, row=1, col=1)
    fig.update_yaxes(scaleanchor="x3", scaleratio=1, row=2, col=1)
    fig.update_yaxes(scaleanchor="x4", scaleratio=1, row=2, col=2)
    fig.update_layout(height=height, margin=dict(l=8, r=8, t=34, b=8),
                      showlegend=False)
    for a in fig.layout.annotations[:3]:
        a.font = dict(size=12, color="#334155")
    return fig

# ------------------------------------------------------------------ #
#  v4: input panels (workbench layout)                                #
# ------------------------------------------------------------------ #
def design_inputs(cool_df) -> dict:
    d = {}
    c1 = c2 = c3 = contextlib.nullcontext()
    ss = st.session_state
    t_cells = (f"Cells - {ss.get('w_Ns',108)}S{ss.get('w_Np',10)}P x "
               f"{ss.get('w_fmt','21700')}, {ss.get('w_rdc',25.0):.0f} mΩ")
    with c1, st.expander(t_cells, expanded=False):
        d["fmt"] = _w(st.selectbox, "Format", "fmt", "21700",
                      options=list(FORMATS) + ["Custom"],
                      help="Presets from teardown data; 4680 DCIR is an estimate.")
        if d["fmt"] != "Custom" and st.button(f"Apply {d['fmt']} preset"):
            f_ = FORMATS[d["fmt"]]
            st.session_state["_pending"] = dict(
                w_dcell=f_["d_cell"]*1000, w_hcell=f_["h_cell"]*1000,
                w_cap=f_["cap"], w_rdc=f_["rdc"], w_mcell=f_["mcell"],
                w_pitch=max(st.session_state.get("w_pitch", 27.0),
                            f_["d_cell"]*1000 + 6.0))
            st.rerun()
        d["d_cell"] = _w(st.slider, "Diameter [mm]", "dcell", 21.0, min_value=10.0, max_value=60.0, step=0.5) / 1000
        d["h_cell"] = _w(st.slider, "Height [mm]", "hcell", 70.0, min_value=40.0, max_value=130.0, step=1.0) / 1000
        d["cap_Ah"] = _w(st.number_input, "Capacity [Ah]", "cap", 5.0, min_value=1.0, max_value=30.0, step=0.1)
        d["v_nom"] = _w(st.number_input, "Nominal V", "vnom", 3.7, min_value=3.0, max_value=4.0, step=0.05)
        d["r_dc"] = _w(st.number_input, "DCIR at 25 °C [mΩ]", "rdc", 25.0, min_value=2.0, max_value=80.0, step=0.5)
        d["k_dcir"] = _w(st.slider, "DCIR fall [%/K]", "kdcir", 1.2, min_value=0.0, max_value=3.0, step=0.1) / 100.0
        d["m_cell"] = _w(st.number_input, "Mass [kg]", "mcell", 0.070, min_value=0.03, max_value=0.6, step=0.001, format="%.3f")
        d["cp_cell"] = _w(st.number_input, "cp [J/kgK]", "cpcell", 950.0, min_value=700.0, max_value=1200.0, step=10.0)
        d["k_rad"] = _w(st.slider, "Radial k_r [W/mK]", "krad", 0.9, min_value=0.3, max_value=2.0, step=0.1)
        with st.expander("Electrical detail"):
            d["v_max"] = _w(st.slider, "Charge V limit", "vmax", 4.2, min_value=4.0, max_value=4.4, step=0.05)
            d["chg_mult"] = _w(st.slider, "Charge DCIR x", "chgm", 1.10, min_value=1.0, max_value=1.5, step=0.05)
            d["v_cut"] = _w(st.slider, "CV cut-off [C]", "vcut", 0.05, min_value=0.02, max_value=0.2, step=0.01)
            d["entropic"] = _w(st.checkbox, "Entropic heat", "entro", True)
        d["Ns"] = _w(st.number_input, "Series (S)", "Ns", 108, min_value=1, max_value=300)
        d["Np"] = _w(st.number_input, "Parallel (P)", "Np", 10, min_value=1, max_value=60)
        with st.popover("Suggest S x P"):
            tE = st.number_input("Target kWh", 1.0, 200.0, 20.0, 0.5)
            tV = st.number_input("Target V", 48.0, 900.0, 400.0, 10.0)
            if st.button("Apply suggestion"):
                vn = st.session_state.get("w_vnom", 3.7); cp_ = st.session_state.get("w_cap", 5.0)
                ns = max(int(round(tV / vn)), 1)
                st.session_state["_pending"] = dict(
                    w_Ns=ns, w_Np=max(int(round(tE * 1000 / (ns * vn * cp_))), 1))
                st.rerun()
    t_lay = (f"Layout and coolant - {ss.get('w_pitch',27.0):.0f} mm pitch, "
             f"{ss.get('w_fluid','MIVOLT DF7')}, "
             f"{ss.get('w_circ','Thermosiphon only').split(' (')[0].lower()}")
    with c2, st.expander(t_lay, expanded=False):
        d["arrangement"] = _w(st.radio, "Arrangement", "arr", "Square",
                              options=["Square", "Hexagonal"], horizontal=True)
        d["pitch"] = _w(st.slider, "Cell pitch [mm]", "pitch", 27.0, min_value=22.0, max_value=60.0, step=0.5) / 1000
        d["edge_margin"] = _w(st.slider, "Edge margin [mm]", "edge", 10.0, min_value=5.0, max_value=40.0, step=1.0) / 1000
        d["tube_zone"] = _w(st.slider, "Tube zone height [mm]", "tz", 35.0, min_value=15.0, max_value=80.0, step=1.0) / 1000
        d["gas_gap"] = _w(st.slider, "Headspace [mm]", "gas", 10.0, min_value=0.0, max_value=40.0, step=1.0) / 1000
        d["end_fraction"] = _w(st.slider, "End-caps wetted", "ends", 0.0, min_value=0.0, max_value=1.0, step=0.1)
        names = list(cool_df["name"])
        d["coolant"] = _w(st.selectbox, "Fluid", "fluid",
                          "MIVOLT DF7" if "MIVOLT DF7" in names else names[0], options=names)
        with st.expander("ν(T) - viscosity falls with temperature (in the solver)"):
            st.plotly_chart(nu_T_fig(fluid_dict(
                cool_df[cool_df["name"] == d["coolant"]].iloc[0])),
                use_container_width=True, key="nuT")
        d["circ"] = _w(st.selectbox, "Circulation method", "circ",
                       "Thermosiphon only",
                       options=["Thermosiphon only", "Open stirring",
                                "Serpentine plates (guided)",
                                "External pump loop (closed)",
                                "Bottom propeller (axial, up)"],
                       help="Realistic options are compared in Learn and "
                            "Ideas. The external pump keeps the dielectric "
                            "in a closed loop through an outside pump - no "
                            "external heat exchanger; the internal water "
                            "tubes still remove the heat, so oil and water "
                            "never meet.")
        d["u_oil"] = 0.0
        d["plate_on"] = False
        _ext = d["circ"].startswith("External")
        if d["circ"] == "Open stirring":
            d["u_oil"] = _w(st.slider, "Stirring velocity [m/s]", "uoil", 0.0,
                            min_value=0.0, max_value=0.20, step=0.005)
        if d["circ"].startswith("Bottom"):
            d["u_oil"] = _w(st.slider, "Upward velocity through the "
                            "array [m/s]", "upr", 0.05, min_value=0.0,
                            max_value=0.15, step=0.005,
                            help="A shrouded axial impeller under the "
                                 "array pushes the oil straight up, "
                                 "aligned with buoyancy. The film uses "
                                 "the AXIAL flat-plate correlation "
                                 "along the can, not crossflow - at "
                                 "equal velocity that is a weaker film "
                                 "than a horizontal sweep, which the "
                                 "solver shows honestly.")
        if _ext:
            d["u_loop"] = _w(st.slider, "Loop velocity through the pack "
                             "[m/s]", "ulp", 0.05, min_value=0.0,
                             max_value=0.20, step=0.005,
                             help="Set by the external pump's operating "
                                  "point; drives both oil films exactly "
                                  "like stirring/guided flow.")
            d["u_oil"] = d["u_loop"]
            d["pipe_id"] = _w(st.slider, "External pipe bore [mm]", "pid",
                              19.0, min_value=6.0, max_value=50.0,
                              step=1.0) / 1000
            d["pipe_len"] = _w(st.slider, "External pipe length, out and "
                               "back [m]", "plen", 2.5, min_value=0.5,
                               max_value=8.0, step=0.25)
            d["plate_on"] = _w(st.checkbox, "Guide the flow with serpentine "
                               "plates", "extpl", True,
                               help="Plates form the parallel channels the "
                                    "pump feeds AND act as fins bonded to "
                                    "the water tubes.")
        if d["circ"] == "Serpentine plates (guided)":
            d["plate_on"] = True
            d["u_guided"] = _w(st.slider, "Guided channel velocity [m/s]", "ugd",
                               0.05, min_value=0.0, max_value=0.15, step=0.005)
            d["u_oil"] = d["u_guided"]
        if d["plate_on"]:
            d["plate_t"] = _w(st.select_slider, "Plate thickness [mm]", "plt", 1.5,
                              options=[1.0, 1.5, 2.0]) / 1000
            d["plate_mat"] = _w(st.selectbox, "Plate material", "plm", "Aluminium",
                                options=["Aluminium", "Copper"])
            d["plate_contact"] = _w(st.slider, "Plate-to-tube contact factor",
                                    "plc", 0.8, min_value=0.4, max_value=1.0,
                                    step=0.05,
                                    help="Brazed round tube ~0.9, clamped "
                                         "~0.6; a flat-faced square or "
                                         "rectangular tube bonds better "
                                         "(~0.95).")
        with st.expander("Holders"):
            d["m_holder_g"] = _w(st.slider, "Holder mass [g/cell]", "mhold", 8.0, min_value=0.0, max_value=25.0, step=1.0)
            d["holder_block"] = _w(st.slider, "Gap-flow blockage", "hblk", 0.20, min_value=0.0, max_value=0.6, step=0.05)
    t_hx = (f"Heat exchanger, water, structure - {ss.get('w_ntub',16)} tubes, "
            f"{ss.get('w_flow',10.0):.0f} L/min at {ss.get('w_twin',20.0):.0f}°C")
    with c3, st.expander(t_hx, expanded=False):
        d["n_tubes"] = _w(st.slider, "Tubes", "ntub", 16, min_value=1, max_value=60)
        _npl_hint = max(ss.get("w_Ns", 108) * ss.get("w_Np", 10), 1)
        d["tube_shape"] = _w(st.selectbox, "Tube cross-section", "tshape",
                             "Round",
                             options=["Round", "Square", "Rectangular"],
                             help="Square/rectangular tubes present a "
                                  "flat face to the serpentine plates - a "
                                  "far better bond than a round tube's "
                                  "line contact - and more wetted "
                                  "perimeter per tube. The solver uses "
                                  "the proper duct physics (Shah-London "
                                  "laminar Nu and f-Re by aspect ratio, "
                                  "Gnielinski on the hydraulic diameter).")
        if d["tube_shape"] == "Round":
            d["tube_od"] = _w(st.slider, "Tube OD [mm]", "tod", 10.0, min_value=4.0, max_value=25.0, step=0.5) / 1000
        elif d["tube_shape"] == "Square":
            d["tube_w"] = _w(st.slider, "Tube side [mm]", "tw", 10.0,
                             min_value=4.0, max_value=25.0, step=0.5) / 1000
            d["tube_od"] = d["tube_w"]
        else:
            d["tube_w"] = _w(st.slider, "Tube width (into plate) [mm]",
                             "tw", 12.0, min_value=4.0, max_value=30.0,
                             step=0.5) / 1000
            d["tube_h"] = _w(st.slider, "Tube height (along plate) [mm]",
                             "th", 8.0, min_value=4.0, max_value=30.0,
                             step=0.5) / 1000
            d["tube_od"] = max(d["tube_w"], d["tube_h"])
        d["tube_wall"] = _w(st.slider, "Wall [mm]", "twall", 1.0, min_value=0.5, max_value=3.0, step=0.25) / 1000
        d["tube_mat"] = _w(st.selectbox, "Tube material", "tmat", "Copper", options=list(K_TUBE))
        if ss.get("w_circ", "").startswith(("Serpentine", "External")) or \
           ss.get("w_extpl", False):
            _ncol = math.ceil(math.sqrt(_npl_hint))
            _npl = max(math.ceil(_npl_hint / _ncol) - 1, 1)
            _sug = sorted({_npl} | {_npl // k for k in (2, 3, 4)
                                    if _npl // k >= 1}, reverse=True)
            _nt_now = ss.get("w_ntub", 16)
            if _nt_now < _npl:
                st.caption(f"With ~{_npl} plate gaps, {_nt_now} tubes "
                           f"leaves {_npl - _nt_now} plates without a "
                           f"bonded tube - those act only as passive "
                           f"spreaders and the model derates the plate "
                           f"area to {100*min(1,_nt_now/_npl):.0f}%. "
                           f"Symmetric counts: "
                           f"{', '.join(str(s) for s in _sug)}.")
            elif _nt_now > _npl:
                st.caption(f"More tubes ({_nt_now}) than plate gaps "
                           f"(~{_npl}); the extras act as bare tubes in "
                           f"the tube zone. One per plate ({_npl}) is "
                           f"the natural count.")
        d["passes"] = _w(st.slider, "Passes", "pass", 1, min_value=1, max_value=4)
        d["tube_plane"] = _w(st.selectbox, "Tube plane", "tplane", "Top of pack",
                             options=["Top of pack", "Interstitial (between rows)",
                                      "Mid-height", "Below the cells"])
        if d["tube_shape"] == "Round":
            d["fins_on"] = _w(st.checkbox, "Annular fins", "fins", True)
        else:
            d["fins_on"] = False
            st.caption("Annular fins fit round tubes only; flat-sided "
                       "tubes rely on the plates for extended area.")
        if d["fins_on"]:
            d["fin_h"] = _w(st.slider, "Fin height [mm]", "finh", 8.0, min_value=2.0, max_value=20.0, step=0.5) / 1000
            d["fin_t"] = _w(st.slider, "Fin thickness [mm]", "fint", 0.5, min_value=0.2, max_value=1.5, step=0.1) / 1000
            d["fin_p"] = _w(st.slider, "Fin pitch [mm]", "finp", 4.0, min_value=2.0, max_value=12.0, step=0.5) / 1000
            d["fin_mat"] = _w(st.selectbox, "Fin material", "finm", "Aluminium", options=["Aluminium", "Copper"])
        d["loop_fluid"] = _w(st.selectbox, "Water loop fluid", "loopf", "Water", options=list(WATER_LOOP))
        d["flow_lpm"] = _w(st.slider, "Flow [L/min]", "flow", 10.0, min_value=0.5, max_value=60.0, step=0.5)
        d["T_water_in"] = _w(st.slider, "Inlet T [°C]", "twin", 20.0, min_value=0.0, max_value=40.0, step=1.0)
        with st.expander("Structure and busbars"):
            d["R_bus"] = _w(st.number_input, "Busbar R [mΩ] (0=auto)", "rbus", 0.0, min_value=0.0, max_value=20.0, step=0.1)
            d["bus_J"] = _w(st.slider, "Busbar J [A/mm²]", "busj", 5.0, min_value=2.0, max_value=10.0, step=0.5)
            d["sigma_MPa"] = _w(st.slider, "Allowable stress [MPa]", "sigma", 80.0, min_value=30.0, max_value=200.0, step=5.0)
            d["stiff"] = _w(st.slider, "Stiffening knock-down", "stiff", 0.45, min_value=0.2, max_value=1.0, step=0.05)
            d["p_des_bar"] = _w(st.slider, "Design pressure [bar g]", "pdes", 0.5, min_value=0.1, max_value=2.0, step=0.1)
    return d

def duty_inputs() -> dict:
    d = {}
    d["duty"] = _w(st.selectbox, "Duty profile", "duty", "Constant C",
                   options=["Constant C", "Drive cycle", "Charge (CC-CV)",
                            "Cycling (dis/chg x N)", "Fast charge then rest",
                            "Pulse train", "CSV upload"])
    d["C1"] = _w(st.slider, "C-rate (primary / discharge)", "c1", 2.0, min_value=0.2, max_value=8.0, step=0.1)
    if d["duty"] == "Constant C":
        d["dirn"] = _w(st.radio, "Direction", "dirn", "Discharge",
                       options=["Discharge", "Charge"], horizontal=True)
        d["track_soc"] = _w(st.checkbox, "Track SoC (off = thermal-only)", "tsoc", False)
    if d["duty"] == "Drive cycle":
        d["cycle"] = _w(st.selectbox, "Cycle", "cycle", "WLTP Class 3b", options=list(DRIVE_CYCLES),
                        help="Breakpoint profiles scaled to official distance; CSV upload for 1 Hz traces.")
        d["repeat_cyc"] = _w(st.checkbox, "Repeat to fill simulation", "repcyc", True)
        with st.expander("Vehicle"):
            d["veh_m"] = _w(st.number_input, "Mass [kg]", "vehm", 1900.0, min_value=600.0, max_value=4000.0, step=50.0)
            d["CdA"] = _w(st.number_input, "Cd x A [m2]", "cda", 0.62, min_value=0.3, max_value=1.5, step=0.01)
            d["Crr"] = _w(st.number_input, "Crr", "crr", 0.009, min_value=0.005, max_value=0.02, step=0.001, format="%.3f")
            d["eta_dt"] = _w(st.slider, "Drivetrain eff.", "etad", 0.92, min_value=0.7, max_value=0.98, step=0.01)
            d["eta_rg"] = _w(st.slider, "Regen recovery", "etar", 0.65, min_value=0.0, max_value=0.95, step=0.05)
            d["P_rg"] = _w(st.slider, "Regen cap [kW]", "prg", 60.0, min_value=0.0, max_value=300.0, step=5.0)
            d["P_acc"] = _w(st.number_input, "Accessories [W]", "pacc", 500.0, min_value=0.0, max_value=5000.0, step=50.0)
    if d["duty"] == "Cycling (dis/chg x N)":
        d["n_cyc"] = _w(st.slider, "Cycles", "ncyc", 3, min_value=1, max_value=10)
        d["cyc_rest"] = _w(st.slider, "Rest [s]", "crest", 600.0, min_value=0.0, max_value=3600.0, step=60.0)
        d["soc_min"] = _w(st.slider, "Discharge to SoC", "socmin", 0.10, min_value=0.0, max_value=0.5, step=0.05)
    if d["duty"] in ("Fast charge then rest", "Pulse train"):
        d["t1"] = _w(st.slider, "Primary phase [s]", "t1", 900.0, min_value=60.0, max_value=3600.0, step=30.0)
    if d["duty"] == "Pulse train":
        d["C2"] = _w(st.slider, "C-rate (secondary)", "c2", 0.5, min_value=0.0, max_value=4.0, step=0.1)
        d["t2"] = _w(st.slider, "Secondary phase [s]", "t2", 600.0, min_value=30.0, max_value=3600.0, step=30.0)
    if d["duty"] == "CSV upload":
        dup = st.file_uploader("Duty CSV: t_s and C (or P_kW)", type="csv", key="duty_up")
        if dup is not None:
            try:
                probe = dict(DEFAULTS)
                probe.update({k: st.session_state.get(f"w_{k2}", probe[k]) for k, k2 in
                              [("Ns", "Ns"), ("Np", "Np"), ("v_nom", "vnom"), ("cap_Ah", "cap")]})
                st.session_state["duty_csv"] = duty_from_csv(dup, probe)
                tt, cc = st.session_state["duty_csv"]
                st.caption(f"{len(tt)} points, {tt[-1]:.0f} s, peak {cc.max():.1f}C.")
            except Exception as e:
                st.error(f"Could not parse: {e}")
    if d["duty"] in ("Charge (CC-CV)", "Cycling (dis/chg x N)", "Drive cycle"):
        d["C_chg"] = _w(st.slider, "Charge rating (CC) [C]", "cchg", 1.0, min_value=0.2, max_value=4.0, step=0.1,
                        help="Plating map derates below 25 °C; also caps regen.")
        d["soc0"] = _w(st.slider, "Start SoC", "soc0", 0.90 if d["duty"] != "Charge (CC-CV)" else 0.20,
                       min_value=0.0, max_value=1.0, step=0.05)
    st.markdown("---")
    d["duration"] = _w(st.slider, "Simulation length [s]", "dur", 3600.0, min_value=300.0, max_value=14400.0, step=300.0)
    d["T_start"] = _w(st.slider, "Start temperature [°C]", "tstart", 25.0, min_value=-10.0, max_value=45.0, step=1.0)
    d["T_limit"] = _w(st.slider, "Cell limit [°C]", "tlim", 45.0, min_value=35.0, max_value=60.0, step=1.0)
    d["limit_core"] = _w(st.checkbox, "Apply limit to core", "limcore", False)
    d["T_amb"] = _w(st.slider, "Ambient [°C]", "tamb", 25.0, min_value=-10.0, max_value=45.0, step=1.0)
    return d

# ------------------------------------------------------------------ #
#  v4 main: workbench                                                 #
# ------------------------------------------------------------------ #
def main():
    st.set_page_config(page_title="Immersion Pack Lab", layout="wide", page_icon=None)
    st.markdown(CSS, unsafe_allow_html=True)
    for _k, _v in st.session_state.pop("_pending", {}).items():
        st.session_state[_k] = _v
    hero_box = st.container()
    kpi_box = st.container()

    if "cool_df" not in st.session_state:
        st.session_state.cool_df = _read_coolants()
    cool_df = st.session_state.cool_df

    tabs = st.tabs(["Design", "Duty", "Results", "Cockpit", "Zones",
                    "Improve", "Ideas", "Safety", "Compare",
                    "Learn", "Cases", "FEA", "Validate", "System",
                    "Report", "Battery"])

    with tabs[0]:
        colL, colR = st.columns([1.15, 1], gap="large")
        with colL:
            st.caption("Define the pack. The live view on the right follows "
                       "every change and stays with you as you scroll.")
            d = dict(DEFAULTS)
            d.update(design_inputs(cool_df))
    with tabs[1]:
        cduty, cplot = st.columns([1, 1.4], gap="large")
        with cduty:
            st.caption("What you ask of the pack. Steady state uses the duty's RMS C-rate.")
            d.update(duty_inputs())
    # model-tuning keys (widgets live in tab 8, read here pre-compute)
    for key, wkey in [("cal_h", "cal_h"), ("K_loop", "kloop"), ("h_ext", "hext"),
                      ("bottom_gap", "bgap"), ("manifold_margin", "mman"),
                      ("T_service_max", "tserv")]:
        default = DEFAULTS[key]
        v = st.session_state.get(f"w_{wkey}", default * (1000 if key in
                                 ("bottom_gap", "manifold_margin") else 1))
        d[key] = v / 1000 if key in ("bottom_gap", "manifold_margin") else v

    # ---------------- compute ---------------- #
    fl = fluid_dict(cool_df[cool_df["name"] == d["coolant"]].iloc[0])
    g = build_geometry(d)
    res0 = solve_steady(d, g, fl, 1.0, d["T_amb"], C_rate=d["C1"])
    masses = build_masses(d, g, fl, res0["fin"])
    Cmax = max_continuous_C(d, g, fl, d["T_amb"], d["T_limit"])
    if d["duty"] == "Drive cycle":
        tc, vc, D_km = cycle_speed(d["cycle"])
        if d["repeat_cyc"] and tc[-1] < d["duration"]:
            reps = int(math.ceil(d["duration"] / tc[-1]))
            vc = np.tile(vc, reps)[: int(d["duration"]) + 1]
            tc = np.arange(0.0, len(vc), 1.0)
        spec = dict(kind="P", t=tc, P=vehicle_battery_power(tc, vc, d), v=vc, D_km=D_km)
    elif d["duty"] == "Charge (CC-CV)":
        spec = dict(kind="chg", t=np.arange(0.0, d["duration"] + 1e-9, 2.0))
    elif d["duty"] == "Cycling (dis/chg x N)":
        spec = dict(kind="cyc", t=np.arange(0.0, d["duration"] + 1e-9, 2.0))
    else:
        csv_tc = st.session_state.get("duty_csv") if d["duty"] == "CSV upload" else None
        t0, C0 = duty_profile(d["duty"], d["duration"], d["C1"], d["t1"], d["C2"],
                              d["t2"], csv_tc)
        spec = dict(kind="C", t=t0, C=C0)
    tr = simulate_pack(d, g, fl, masses, d["T_amb"], spec)
    t_arr, C_arr = tr["t"], tr["C"]
    C_steady = max(tr["C_rms"], 0.05)
    res = solve_steady(d, g, fl, 1.0, d["T_amb"], C_rate=C_steady)
    Q_duty = res["Q_eff"]
    loop = WATER_LOOP[d["loop_fluid"]]
    T_gov = res["T_core"] if d["limit_core"] else res["T_b"]
    ok = T_gov <= d["T_limit"]
    P_pump = water_pump_power(d, g, loop)["P"]
    if d["circ"].startswith("External"):
        _xp = ext_loop_pump(d, g, fl, d["u_oil"])
        P_stir = _xp["P"]
        res["ext"] = _xp
    elif d["circ"].startswith("Bottom"):
        _pp = prop_power(d, g, fl, d["u_oil"])
        P_stir = _pp["P"]
        res["prop"] = _pp
    elif d.get("plate_on"):
        _sp = serpentine_pump(d, g, fl, d["u_oil"])
        P_stir = _sp["P"]
    else:
        P_stir = stirrer_power(d, g, fl, d["u_oil"])
    Q_bus = (C_steady * d["cap_Ah"] * d["Np"]) ** 2 * tr["bus"]["R"]
    cost_now = dict(COST_DEFAULTS)
    cost_now.update(st.session_state.get("cost_edit", {}))
    arch_df = _arch_study_cached(json.dumps(d, sort_keys=True, default=float),
                                 float(round(C_steady, 3)),
                                 json.dumps(cost_now, sort_keys=True))
    chil = chiller_model(res["Q_w"] + P_pump, d["T_water_in"], d["T_amb"])

    # ---------------- hero + KPI cards ---------------- #
    T_gov = res["T_core"] if d["limit_core"] else res["T_b"]
    ok = T_gov <= d["T_limit"]
    with hero_box:
        st.markdown(
            f"<div class='hero'><div class='hero-top'>"
            f"<h1>Immersion Pack Lab <span style='font-size:.6em;"
            f"color:#94A3B8'>{APP_VERSION}</span></h1>"
            f"<span class='sub'>{masses['E_kwh']:.1f} kWh / "
            f"{d['Ns']*d['v_nom']:.0f} V - {g['N']} x {d['fmt']} in "
            f"{fl['name'].split('(')[0].strip()} - {d['duty']}"
            f"{(' / ' + d['cycle']) if d['duty'] == 'Drive cycle' else ''}"
            f"</span>"
            f"<span class='chip'>{'WITHIN LIMIT' if ok else 'OVER LIMIT'} - "
            f"{T_gov:.1f} / {d['T_limit']:.0f} °C at {C_steady:.2f}C rms"
            f"</span></div>"
            f"<div class='hero-stats'>"
            f"<div class='hstat'>can / core<b class='{'ok' if ok else 'bad'}'>"
            f"{res['T_b']:.1f} / {res['T_core']:.1f}°C</b></div>"
            f"<div class='hstat'>max continuous<b>{Cmax:.2f} C</b></div>"
            f"<div class='hstat'>heat at duty<b>{Q_duty/1000:.2f} kW</b></div>"
            f"<div class='hstat'>chiller<b>{chil['P_el']/1000:.2f} kW el</b></div>"
            f"<div class='hstat'>parasitics<b>{P_pump+P_stir:.0f} W</b></div>"
            f"<div class='hstat'>coolant<b>{masses['V_oil_L']:.0f} L / "
            f"{masses['m_oil']:.0f} kg</b></div>"
            f"<div class='hstat'>pack<b>{masses['m_pack']:.0f} kg - "
            f"{masses['whkg_pack']:.0f} Wh/kg</b></div>"
            f"</div></div>", unsafe_allow_html=True)

    # ---------------- sidebar: status, save/load, report ---------------- #
    sb = st.sidebar
    sb.markdown(
f"<div class='kpi'><div class='l'>Design status - {APP_VERSION}</div>"
        f"<div class='v {'ok' if ok else 'bad'}'>"
        f"{'Within limit' if ok else 'Over limit'}</div>"
        f"<div class='s'>{res['T_b']:.1f} °C at {C_steady:.2f}C rms - "
        f"thermosiphon {res['u_ts']*1000:.1f} mm/s - spread "
        f"{res['spread']:.1f} °C</div></div>", unsafe_allow_html=True)
    if not ok:
        sb.error("Over limit - see Improve tab.")
    sb.markdown("---")
    with sb.expander("Save / load design"):
        state = {k: v for k, v in st.session_state.items()
                 if k.startswith("w_") and isinstance(v, (int, float, str, bool))}
        st.download_button("Download design (.json)", data=pd.Series(state).to_json(),
                           file_name="pack_design.json", use_container_width=True)
        up = st.file_uploader("Load design", type="json", key="design_up")
        if up is not None:
            sig = up.name + str(up.size)
            if st.session_state.get("applied_design") != sig:
                for k, v in pd.read_json(up, typ="series").to_dict().items():
                    if k.startswith("w_"):
                        st.session_state[k] = v
                st.session_state["applied_design"] = sig
                st.rerun()
    # (full narrative report lives in tab 9; sidebar keeps a shortcut)
    sb.caption("Full narrative report: Report tab.")
    sb.caption("Wang et al. 2023 benchmark: 33.0 vs 32.3 °C. Oil-side h honest to "
               "+/-30% until calibrated (tab 8).")

    # ---------------- Design tab: sticky live view ---------------- #
    with colR:
        with st.container(key="liveview"):
            with st.container(border=True):
                v1, v2, v3 = st.columns([2.2, 0.8, 0.8])
                view = v1.radio("View", ["Views", "3D", "Plan"],
                                horizontal=True,
                                label_visibility="collapsed")
                s_oil = v2.toggle("Oil", True)
                s_box = v3.toggle("Box", True)
                s_tub, s_fin = True, True
                if view == "Views":
                    st.plotly_chart(pack_views_fig(d, g, masses, height=520),
                                    use_container_width=True, key="liveviews",
                                    config=PLOTCFG)
                elif view == "3D":
                    st.plotly_chart(pack_3d_figure(d, g, s_oil, s_tub, s_box,
                                                   s_fin, height=470),
                                    use_container_width=True, key="live3d",
                                    config=PLOTCFG)
                else:
                    st.plotly_chart(layout_figure(d, g, height=470),
                                    use_container_width=True, key="liveplan")
                T_gov0 = res["T_core"] if d["limit_core"] else res["T_b"]
                ok0 = T_gov0 <= d["T_limit"]
                pin = st.session_state.get("pinned")
                def _dlt(cur, key, fmt="{:+.1f}"):
                    if not pin or key not in pin: return ""
                    return f" ({fmt.format(cur - pin[key])} vs A)"
                kpi_cards([
                    ("Can / core", f"{res['T_b']:.1f} / {res['T_core']:.1f}°C",
                     f"{T_gov0-d['T_limit']:+.1f} °C margin"
                     + _dlt(res['T_b'], 'T_can'), "ok" if ok0 else "bad"),
                    ("Cell spread (max-min)", f"{res['spread']:.1f} °C",
                     "hottest vs coldest cell; keep under 5 °C" + _dlt(res['spread'], 'Spread_K'), ""),
                    ("Mass", f"{masses['m_pack']:.0f} kg",
                     f"{masses['whkg_pack']:.0f} Wh/kg"
                     + _dlt(masses['m_pack'], 'Pack_kg', '{:+.0f} kg'), ""),
                ])
                st.caption(f"Box {g['Lx']*1000:.0f}×{g['Ly']*1000:.0f}×"
                           f"{g['Lz']*1000:.0f} mm - gap {g['gap_mm']:.1f} mm - "
                           f"circulation: {d['circ'].lower()}"
                           + (f" - plates +{res.get('A_plate',0):.1f} m² eff."
                              if d.get('plate_on') else ""))
        with st.expander("Mass and packaging audit"):
            mdf2 = pd.DataFrame([
                ["Cells", masses["m_cells"]], ["Coolant", masses["m_oil"]],
                ["Tubes", masses["m_tubes"]], ["Fins", masses["m_fins"]],
                ["Busbars", masses["m_bus"]], ["Cell holders", masses["m_holders"]],
                [f"Enclosure ({masses['enc']['t_mm']:.1f} mm eff.)",
                 masses["m_struct"]],
            ], columns=["Item", "kg"])
            mdf2["% of pack"] = 100 * mdf2["kg"] / masses["m_pack"]
            st.dataframe(mdf2.round(1), hide_index=True, use_container_width=True)
            st.caption(f"Worst cell ~{res['T_worst']:.1f} °C, best "
                       f"~{res['T_best']:.1f} (spread {res['spread']:.1f} °C vs "
                       f"5 K). Buffer "
                       f"{(masses['C_oil']+masses['C_batt'])/1e3:.0f} kJ/K.")

    # ---------------- Duty tab: response ---------------- #
    with tabs[1]:
        with cplot:
            if spec["kind"] == "P":
                figV = go.Figure(go.Scatter(x=spec["t"]/60, y=spec["v"]*3.6,
                                            line=dict(color="#6366F1", width=2)))
                figV.update_layout(height=200, margin=dict(l=10, r=10, t=30, b=10),
                                   title=f"{d['cycle']} speed trace",
                                   xaxis_title="min", yaxis_title="km/h",
                                   plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)")
                st.plotly_chart(figV, use_container_width=True)
            sc1, sc2 = st.columns([2.2, 1])
            t_end = float(tr["t"][-1])
            t_pick = sc1.slider("Inspect time [s]", 0.0, t_end,
                                min(t_end, t_end * 0.5), step=max(t_end/200, 2.0),
                                key="t_scrub")
            band = sc2.checkbox("±30% oil-film band", False)
            i_p = int(np.argmin(np.abs(tr["t"] - t_pick)))
            kpi_cards([
                ("At t = {:.0f} s".format(tr["t"][i_p]),
                 f"{tr['T_b'][i_p]:.1f} / {tr['T_core'][i_p]:.1f}°C",
                 "can / core", "bad" if tr["T_b"][i_p] > d["T_limit"] else "ok"),
                ("C-rate", f"{tr['C'][i_p]:+.2f} C",
                 f"SoC {tr['soc'][i_p]*100:.0f}%", ""),
                ("Oil", f"{tr['T_il'][i_p]:.1f}°C",
                 f"heat {tr['Q'][i_p]/1000:.2f} kW", ""),
            ])
            figT = go.Figure()
            if band:
                trs = []
                for cf in (0.7, 1.3):
                    dd = dict(d); dd["cal_h"] = d["cal_h"] * cf
                    trs.append(simulate_pack(dd, g, fl, masses, d["T_amb"], spec))
                figT.add_trace(go.Scatter(
                    x=np.concatenate([tr["t"], tr["t"][::-1]])/60,
                    y=np.concatenate([trs[0]["T_b"], trs[1]["T_b"][::-1]]),
                    fill="toself", fillcolor="rgba(192,57,43,0.12)",
                    line=dict(width=0), name="h +/-30%", hoverinfo="skip"))
            figT.add_trace(go.Scatter(x=tr["t"]/60, y=tr["T_core"], name="Core",
                                      line=dict(color="#B91C1C", width=2, dash="dot")))
            figT.add_trace(go.Scatter(x=tr["t"]/60, y=tr["T_b"], name="Cell surface",
                                      line=dict(color="#EF4444", width=3)))
            figT.add_trace(go.Scatter(x=tr["t"]/60, y=tr["T_il"], name="Bulk oil",
                                      line=dict(color=ACCENT, width=3)))
            figT.add_hline(y=d["T_limit"], line_dash="dash", line_color="#B91C1C")
            figT.add_vline(x=t_pick / 60, line_color="#6366F1", line_width=2,
                           opacity=0.7)
            figT.add_trace(go.Scatter(x=t_arr/60, y=C_arr, name="C-rate", yaxis="y2",
                                      line=dict(color="#94A3B8", dash="dot")))
            if d["duty"] != "Constant C" or d.get("track_soc"):
                figT.add_trace(go.Scatter(x=tr["t"]/60, y=tr["soc"]*4, yaxis="y2",
                                          name="SoC (x4)", line=dict(color="#10B981", width=2)))
            figT.update_layout(height=380, xaxis_title="Time [min]",
                               yaxis_title="Temperature [°C]",
                               yaxis2=dict(title="C / SoCx4", overlaying="y", side="right",
                                           showgrid=False,
                                           range=[min(0, float(C_arr.min())*1.2),
                                                  max(float(np.abs(C_arr).max())*1.6, 4.2)]),
                               plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
                               hovermode="x unified",
                               margin=dict(l=10, r=10, t=10, b=10))
            st.plotly_chart(figT, use_container_width=True, config=PLOTCFG)
            if spec["kind"] == "P":
                st.caption(f"Official distance {spec['D_km']:.2f} km/cycle; peak "
                           f"{tr['C'].max():.2f}C / regen {-tr['C'].min():.2f}C; "
                           f"C_rms {tr['C_rms']:.2f}; SoC {tr['soc'][0]*100:.0f} -> "
                           f"{tr['soc'][-1]*100:.0f}%.")
            elif spec["kind"] in ("chg", "cyc"):
                st.caption(f"SoC {tr['soc'][0]*100:.0f} -> {tr['soc'][-1]*100:.0f}%; peak "
                           f"charge {-tr['C'].min():.2f}C (plating-derated); cycles "
                           f"{tr['cycles_done']}.")

    # ---------------- Results ---------------- #
    with tabs[2]:
        stations, totdT = station_list(d, g, fl, res, masses, Q_duty, Q_bus,
                                       P_pump, P_stir, chil)
        films = [("Can to oil", res["R_b"]), ("Oil to tube+fins", res["R_ot"]),
                 ("Water film", res["R_in"])]
        weak_name = max(films, key=lambda x: x[1])[0]
        lp_state = dict(
            T_b=res["T_b"], T_core=res["T_core"], T_w_in=d["T_water_in"],
            dT_water=res["dT_water"], T_limit=d["T_limit"],
            T_amb=d["T_amb"], q_kw=Q_duty / 1000, c_rms=C_steady,
            spread=res["spread"],
            u_mm_s=max(d["u_oil"], res["u_ts"]) * 1000,
            flow_lpm=d["flow_lpm"], flow_norm=min(d["flow_lpm"] / 20, 2.0),
            mode=("serpentine" if d.get("plate_on") else
                  ("stirred" if d["u_oil"] > 0 else "thermosiphon")),
            interstitial=d.get("tube_plane") == "Interstitial (between rows)",
            n_rows_draw=min(g["n_rows"], 14),
            n_tubes_draw=min(d["n_tubes"], 12),
            d_over_p=d["d_cell"] / d["pitch"],
            fill_frac=g["fill_h"] / g["Lz"],
            cell_top_frac=(d["bottom_gap"] + d["h_cell"]) / g["Lz"],
            cell_bot_frac=d["bottom_gap"] / g["Lz"],
            headspace_mm=d["gas_gap"] * 1000,
            weak=weak_name.lower(),
            weak_region=("tubes" if weak_name in
                         ("Water film", "Oil to tube+fins") else "cells"),
            chil_duty_kw=(res["Q_w"] + P_pump) / 1000,
            chil_cop=chil["COP"], chil_el_w=chil["P_el"])
        _q_cell = Q_duty / g["N"]
        lp_stations = [
            dict(id="cells",
                 title=f"Cells - {g['N']} x {d['fmt']}",
                 hint=(f"can {res['T_b']:.1f} °C, core {res['T_core']:.1f} "
                       f"°C, {_q_cell:.2f} W each"),
                 body=(f"<p>Each cell makes <b>{_q_cell:.2f} W</b> at "
                       f"{C_steady:.2f}C rms; heat scales with C². The can "
                       f"sits at <b>{res['T_b']:.1f} °C</b>, the core "
                       f"<b>{res['T_core']:.1f} °C</b> (radial conduction, "
                       f"k_r = {d['k_rad']:.1f} W/m·K).</p>"
                       f"<p>The first oil film costs "
                       f"<b>{Q_duty*res['R_b']:.1f} °C</b>: h = "
                       f"{res['h_cell']:.0f} W/m²·K over "
                       f"{g['A_cells']:.1f} m². Move the oil to thin it.</p>"
                       f"<p>DCIR now {r_of_T(d, res['T_b']):.1f} mΩ vs "
                       f"{d['r_dc']:.0f} at 25 °C - the hot pack makes "
                       f"{100*(1-r_of_T(d,res['T_b'])/d['r_dc']):.0f}% "
                       "less heat.</p>")),
            dict(id="oil",
                 title="Bulk oil - mixer and flywheel",
                 hint=(f"{res['T_il']:.1f} °C, circulation "
                       f"{max(d['u_oil'], res['u_ts'])*1000:.1f} mm/s"),
                 body=(f"<p>Bulk oil at <b>{res['T_il']:.1f} °C</b>. "
                       f"Circulation {max(d['u_oil'],res['u_ts'])*1000:.1f} "
                       f"mm/s ({lp_state['mode']}) keeps the spread at "
                       f"<b>{res['spread']:.1f} °C</b> (criterion 5).</p>"
                       f"<p>Thermal buffer "
                       f"{(masses['C_oil']+masses['C_batt'])/1e3:.0f} kJ/K "
                       "- short peaks never reach the steady picture.</p>"
                       f"<p>Casing sheds {res['Q_atm']:.0f} W free.</p>")),
            dict(id="tubes",
                 title=f"Tubes and fins - {d['n_tubes']} x "
                       f"{d['tube_mat'].lower()}",
                 hint=(f"h_oil {res['h_tube']:.0f}, water Re "
                       f"{res['Re_water']:.0f} ({res['water_regime']})"),
                 body=(f"<p>Second oil film: h = {res['h_tube']:.0f} "
                       f"W/m²·K, A_eff = {res['A_oilside']:.1f} m² -> "
                       f"<b>{Q_duty*res['R_ot']:.1f} °C</b>. Wall: "
                       f"negligible.</p>"
                       f"<p>Water film: {res['water_regime']}, Re = "
                       f"{res['Re_water']:.0f}, h = {res['h_water']:.0f} "
                       f"-> <b>{Q_duty*res['R_in']:.1f} °C</b>. In laminar "
                       "flow more speed does nothing - cross Re 3000.</p>"
                       f"<p>Stream warms {res['dT_water']:.1f} °C along "
                       "the tubes; counterflow plumbing removes that "
                       "spread for free.</p>")),
            dict(id="head",
                 title="Headspace and enclosure",
                 hint=f"{d['gas_gap']*1000:.0f} mm gas, burst "
                      f"{d['p_des_bar']:.1f} bar g",
                 body=(f"<p>{d['gas_gap']*1000:.0f} mm nitrogen blanket "
                       "absorbs oil expansion and vent gas; burst disc at "
                       f"{d['p_des_bar']:.1f} bar g sizes the "
                       f"{masses['enc']['t_mm']:.1f} mm effective wall "
                       f"({masses['m_struct']:.0f} kg - the price of "
                       "abuse tolerance).</p>")),
            dict(id="chiller",
                 title="Chiller",
                 hint=(f"{(res['Q_w']+P_pump)/1000:.2f} kW duty, COP "
                       f"{chil['COP']:.1f}, {chil['P_el']:.0f} W el"),
                 body=(f"<p>Re-cools the loop to {d['T_water_in']:.0f} °C. "
                       f"Duty <b>{(res['Q_w']+P_pump)/1000:.2f} kW</b>, "
                       f"COP <b>{chil['COP']:.1f}</b> (45% of Carnot on a "
                       f"{chil['lift']:.0f} °C lift) -> "
                       f"<b>{chil['P_el']:.0f} W</b> electricity. Room "
                       f"receives {(res['Q_w']+P_pump+chil['P_el'])/1000:.2f}"
                       " kW - heat plus the compressor's wage.</p>"
                       f"<p>Warmer set point: double win (COP up, DCIR "
                       f"down) to the {d['T_limit']:.0f} °C limit; above "
                       f"~{d['T_amb']+5:.0f} °C a dry cooler replaces the "
                       "compressor entirely.</p>")),
        ]
        components.html(live_pack_html(lp_state, lp_stations),
                        height=436)
        st.caption("Live Pack: oil particles move at the solved circulation "
                   "speed, water beads travel and warm along the tubes, and "
                   "colours are the real solved temperatures. Its controls "
                   "run in the browser - instant, no recompute.")

        st.markdown("#### The heat journey")
        ev = st.plotly_chart(thermal_circuit_fig(d, g, fl, res, Q_duty),
                             use_container_width=True, key="res_circuit",
                             config=PLOTCFG, on_select="rerun",
                             selection_mode="points")
        names_s = [s_[0] for s_ in stations]
        if ev and ev.selection and ev.selection.points:
            idx = ev.selection.points[0].get("point_index", 0)
            cmap = {0: 0, 1: 1, 2: 3, 3: 4, 4: 5, 5: 6}
            st.session_state["seg_station"] = names_s[cmap.get(idx, 0)]
        sel = st.segmented_control("Station", names_s, key="seg_station",
                                   default=names_s[0],
                                   label_visibility="collapsed")             if hasattr(st, "segmented_control") else             st.radio("Station", names_s, horizontal=True, key="seg_station")
        for name, dT, nums, imp in stations:
            if name != (sel or names_s[0]):
                continue
            with st.container(border=True):
                h1, h2 = st.columns([3, 1])
                h1.markdown(f"**{name}**")
                if dT > 0:
                    h2.progress(min(dT / max(totdT, 1e-9), 1.0),
                                text=f"{dT:.1f} °C ({100*dT/max(totdT,1e-9):.0f}%)")
                st.markdown(nums)
                st.markdown(f"*Improve:* {imp}")
        with st.expander("Full flow map, ladder and set-point trade"):
            st.plotly_chart(heat_sankey(res, Q_duty, Q_bus, P_pump, P_stir,
                                        chil), use_container_width=True,
                            key="res_sankey", config=PLOTCFG)
            st.markdown(sankey_deep_dive(d, g, fl, res, masses, Q_duty,
                                         Q_bus, P_pump, P_stir, chil,
                                         C_steady))
            st.plotly_chart(thermal_xray_fig(d, g, fl, res, Q_duty),
                            use_container_width=True, key="res_xray",
                            config=PLOTCFG)
            st.plotly_chart(setpoint_trade(d, g, fl, C_steady, d["T_amb"]),
                            use_container_width=True, key="res_setpoint")

    # ---------------- Cockpit ---------------- #
    with tabs[3]:
        st.caption("The whole design space on the rails: pack "
                   "architecture, geometry, fluids, tubes, circulation and "
                   "limits, every readout linked, instant. The physics AND "
                   "the geometry/mass build now run in your browser as "
                   "faithful ports of the app's solver, audited against it "
                   "on load (chip, top left). Instruments: temperature "
                   "ladder, flight recorder, SNAP reference deltas, MAX-C "
                   "autothrottle, turbulence trim, dry-cooler set, and "
                   "scenario presets. Nothing changes the design until you "
                   "copy and apply below.")
        _H_loop = (d["h_cell"] / 2 + d["tube_zone"] / 2
                   if d.get("tube_plane", "Top of pack") == "Top of pack"
                   else 0.008 if d.get("tube_plane") ==
                   "Interstitial (between rows)" else d["h_cell"] / 2)
        cp_fluids = []
        for _, rw in cool_df.iterrows():
            fdd = fluid_dict(rw)
            cp_fluids.append(dict(name=fdd["name"], rho=fdd["rho"],
                                  cp=fdd["cp"], k=fdd["k"],
                                  nu25=fdd["nu25"], B=fdd["B"],
                                  beta=fdd["beta"]))
        cp_pay = dict(
            design=dict(C1=float(C_steady), T_amb=d["T_amb"],
                        coolant=d["coolant"], fmt=d["fmt"],
                        Ns=d["Ns"], Np=d["Np"], cap_Ah=d["cap_Ah"],
                        r_dc=d["r_dc"], k_dcir=d["k_dcir"],
                        k_rad=d["k_rad"], pitch=d["pitch"],
                        arrangement=d["arrangement"],
                        tube_plane=d.get("tube_plane", "Top of pack"),
                        flow_lpm=d["flow_lpm"],
                        T_water_in=d["T_water_in"],
                        n_tubes=d["n_tubes"],
                        loop_fluid=d["loop_fluid"],
                        ver=APP_VERSION,
                        tube_od=d["tube_od"], tube_wall=d["tube_wall"],
                        tube_mat=d["tube_mat"], fins_on=d["fins_on"],
                        tshape=d.get("tube_shape", "Round").lower()[:5]
                        .replace("recta", "rect"),
                        tw=d.get("tube_w", 0.012),
                        th=d.get("tube_h", 0.008),
                        pipe_id=d.get("pipe_id", 0.019),
                        pipe_len=d.get("pipe_len", 2.5),
                        circ0=("extpump" if d["circ"].startswith("External")
                               else "prop" if d["circ"].startswith("Bottom")
                               else "serpentine" if d.get("plate_on") else
                               "stirred" if d["u_oil"] > 0 else
                               "thermosiphon"),
                        u0=(d.get("u_guided", 0.05) if d.get("plate_on")
                            else max(d["u_oil"], 0.05)),
                        plate_t=d.get("plate_t", 0.0015),
                        plate_contact=d.get("plate_contact", 0.8),
                        T_limit=d["T_limit"],
                        limit_core=bool(d["limit_core"]),
                        h_ext=d["h_ext"],
                        edge_margin=d["edge_margin"],
                        bottom_gap=d["bottom_gap"],
                        tube_zone=d["tube_zone"], gas_gap=d["gas_gap"],
                        manifold_margin=d["manifold_margin"],
                        passes=d["passes"],
                        end_fraction=d["end_fraction"],
                        holder_block=d.get("holder_block", 0.25),
                        m_holder_g=d["m_holder_g"],
                        struct_mass=d["struct_mass"],
                        sigma_MPa=d["sigma_MPa"], stiff=d["stiff"],
                        p_des_bar=d["p_des_bar"], bus_J=d["bus_J"],
                        v_nom=d["v_nom"], fin_h=d["fin_h"],
                        fin_t=d["fin_t"], fin_p=d["fin_p"],
                        k_fin=205.0 if d["fin_mat"] == "Aluminium"
                        else 385.0, fin_mat=d["fin_mat"]),
            formats={"18650": dict(d=0.0186, h=0.0652, cap=3.0,
                                   m=0.047, r=35.0),
                     "21700": dict(d=0.0211, h=0.0703, cap=5.0,
                                   m=0.069, r=25.0),
                     "4680": dict(d=0.046, h=0.080, cap=26.0,
                                  m=0.355, r=6.0)},
            fluids=cp_fluids,
            waters={k_: dict(rho=v_["rho"], cp=v_["cp"], k=v_["k"],
                             mu=v_["mu"])
                    for k_, v_ in WATER_LOOP.items()},
            base=dict(T_b=res["T_b"]),
            consts=dict(K_loop=d.get("K_loop", 5.0),
                        cal=d.get("cal_h", 1.0),
                        KT=K_TUBE, RT=RHO_TUBE))
        components.html(cockpit_html(cp_pay), height=820)

        st.markdown("##### 3D view of the serpentine architecture")
        st.caption("The geometry the flat instruments cannot show: cells "
                   "in a grid, a metal cooling plate in each gap, water "
                   "tubes bonded along the plates and running "
                   "**through-plane** (into the page), and coolant "
                   "flowing front to back as it collects heat. Drag to "
                   "rotate; the flow animates. This is the actual heat "
                   "path - cell → oil film → plate → tube wall → water.")
        _dTw3 = res.get("dT_water", 2.0)
        p3 = dict(D=round(d["d_cell"] * 1000, 1),
                  Ht=round(d["h_cell"] * 1000, 1),
                  pitch=round(d["pitch"] * 1000, 1),
                  gap=round(g["gap_mm"], 2),
                  plate_t=round(d.get("plate_t", 0.0015) * 1000, 2),
                  tube_od=round(d["tube_od"] * 1000, 1),
                  nx=4, ny=3, tubes_z=2,
                  Tcell=round(res["T_b"], 1),
                  Tw_in=round(d["T_water_in"], 1),
                  Tw_out=round(d["T_water_in"] + _dTw3, 1),
                  limit=round(d["T_limit"], 1),
                  ns_np=f"{d['Ns']}S{d['Np']}P")
        components.html(_P3.pack3d_html(p3), height=560)
        if not d.get("plate_on"):
            _cm = "stirred" if d["u_oil"] > 0 else "static thermosiphon"
            st.caption(
                f"Note: your Design uses a {_cm} circulation method, not "
                "the guided serpentine plates - the 3D view above shows "
                "the serpentine concept for reference. Switch Circulation "
                "method to *Serpentine plates* in Design (or the cockpit) "
                "to make it the live architecture.")

        with st.expander("What the cockpit computes (the ported "
                         "equations)"):
            st.markdown("The cockpit is a faithful browser port of the "
                        "app's steady two-node solver, audited against "
                        "the Python result on load (the chip, top "
                        "left). Every watt crosses the same resistance "
                        "chain from cell core to coolant:")
            st.latex(r"T_\mathrm{core}=T_w+Q\big(R_\mathrm{core}"
                     r"+R_\mathrm{cell\,film}+R_\mathrm{tube\,film}"
                     r"+R_\mathrm{wall}+R_\mathrm{water}\big)")
            st.latex(r"R=\dfrac{1}{hA},\qquad "
                     r"h_\mathrm{film}=\dfrac{Nu\,k_\mathrm{oil}}{L},"
                     r"\qquad Q=(C\,\mathrm{Ah})^{2}R_{dc}(T)\,N")
            st.markdown("The autothrottle (MAX-C) inverts this chain "
                        "for the C-rate that lands the governing "
                        "temperature exactly on the limit, by bisection "
                        "on the same equations; the turbulence trim "
                        "reports the flow that trips Re = 2300 in the "
                        "tubes. Full derivations with the film "
                        "correlations are in the Learn tab; the "
                        "duct-level version is in Zones.")
        with st.expander("Apply cockpit settings to the design"):
            st.caption("Press 'Copy settings for Design' in the cockpit, "
                       "paste here, and apply. Values land on the real "
                       "Design widgets.")
            cp_txt = st.text_area("Paste settings JSON", "", height=68,
                                  key="cp_paste",
                                  label_visibility="collapsed")
            if st.button("Apply to design", key="cp_apply") and cp_txt:
                try:
                    j = json.loads(cp_txt)
                    pen = st.session_state.setdefault("_pending", {})
                    mapping = {"flow": "w_flow", "twin": "w_twin",
                               "ntub": "w_ntub", "pitch": "w_pitch",
                               "tamb": "w_tamb", "fluid": "w_fluid",
                               "fins": "w_fins", "plt": "w_plt",
                               "plc": "w_plc", "c1": "w_c1",
                               "ns": "w_Ns", "np": "w_Np",
                               "cap": "w_cap", "rdc": "w_rdc",
                               "kdcir": "w_kdcir", "arr": "w_arr",
                               "tplane": "w_tplane", "tod": "w_tod",
                               "twall": "w_twall", "tmat": "w_tmat",
                               "loop": "w_loopf", "tlim": "w_tlim",
                               "limc": "w_limcore", "hext": "w_hext",
                               "tw": "w_tw", "th": "w_th",
                               "pipeD": "w_pid", "pipeL": "w_plen"}
                    ints = {"ntub", "ns", "np"}
                    for k_, wk in mapping.items():
                        if k_ in j:
                            pen[wk] = (int(j[k_]) if k_ in ints
                                       else j[k_])
                    circ_map = {"thermosiphon": "Thermosiphon only",
                                "stirred": "Open stirring",
                                "serpentine":
                                "Serpentine plates (guided)",
                                "extpump":
                                "External pump loop (closed)",
                                "prop":
                                "Bottom propeller (axial, up)"}
                    if "tshape" in j:
                        pen["w_tshape"] = {"round": "Round",
                                           "squar": "Square",
                                           "rect": "Rectangular"}.get(
                            j["tshape"], "Round")
                    if "circ" in j:
                        pen["w_circ"] = circ_map.get(j["circ"],
                                                     "Thermosiphon only")
                        if j["circ"] == "extpump" and "u" in j:
                            pen["w_ulp"] = float(j["u"])
                        elif j["circ"] == "prop" and "u" in j:
                            pen["w_upr"] = float(j["u"])
                        elif j["circ"] == "serpentine" and "u" in j:
                            pen["w_ugd"] = float(j["u"])
                        elif j["circ"] == "stirred" and "u" in j:
                            pen["w_uoil"] = float(j["u"])
                    st.success("Applied - the whole app now reflects the "
                               "cockpit settings.")
                    st.rerun()
                except Exception as e_:
                    st.error(f"Could not parse that: {e_}")

    # ---------------- Zones ---------------- #
    with tabs[4]:
        st.markdown("#### Zonal plate-channel model - every cell, "
                    "every channel, every root")
        st.caption("The plate-channel architecture solved as a full "
                   "thermal-hydraulic network: exact developing-flow "
                   "kernels (fea4, on the true lens cross-section) for "
                   "every oil slot, per-cell DCIR-coupled balances, "
                   "plate bays rooted to water tubes at the derived "
                   "pitch, water marched tube by tube, and the "
                   "recirculation plenum solved in closed form. Kernel "
                   "gates: fRe 95.7 vs 96.0, one-wall Nu 5.387 vs "
                   "5.385, two-wall Nu 7.541 vs 7.541 exact.")
        zc1, zc2, zc3, zc4, zc5, zc6 = st.columns(6)
        z_s = zc1.slider("Slot s [mm]", 1.5, 3.0, 2.0, 0.1,
                         key="zn_s")
        z_gx = zc2.slider("In-row crevice [mm]", 0.4, 1.4, 0.5, 0.1,
                          key="zn_gx")
        z_dp = zc3.slider("Pump head [Pa]", 2.0, 100.0, 25.0, 1.0,
                          key="zn_dp")
        z_C = zc4.slider("C-rate", 0.5, 4.0, float(C_steady), 0.1,
                         key="zn_C")
        z_hm = zc5.selectbox("Heat map", ["Uniform", "Busbar end",
                                          "Centre hot"], key="zn_hm")
        z_amp = zc6.slider("Map amplitude [%]", 0, 30, 10,
                           key="zn_amp")
        zk1, zk2, zk3, zk4 = st.columns(4)
        z_hc = zk1.slider("Contact h_c [W/m²·K]", 1000, 20000, 8000,
                          500, key="zn_hc",
                          help="Tube-to-plate collar contact "
                               "conductance. Unmeasured pending braze "
                               "data - the headline margin is "
                               "conditional on this (see sweep below).")
        z_col = zk2.slider("Collar factor (x plate_t)", 2.0, 10.0, 6.0,
                           0.5, key="zn_col",
                           help="Effective collar engagement length as "
                                "a multiple of plate thickness; sets the "
                                "contact area per crossing.")
        z_byp = zk3.checkbox("Direct oil->tube path", value=False,
                             key="zn_byp",
                             help="Adds the bare-tube oil-to-water path "
                                  "the root-only model omits. Magnitude "
                                  "depends on header geometry, so it is "
                                  "off by default and conservative when "
                                  "off.")
        z_wet = zk4.slider("Wetted tube fraction", 0.0, 1.0, 0.5, 0.1,
                           key="zn_wet", disabled=not z_byp)

        @st.cache_resource(show_spinner=False)
        def _zonal_bank(Dc, px):
            return KernelBank(D=Dc, pitch=px,
                              s_grid=(0.0014, 0.0018, 0.0022,
                                      0.0027, 0.0032),
                              k_oil=0.13, n=100, nz=100)

        px_row = d["d_cell"] + z_gx / 1000.0
        with st.spinner("Computing duct kernels for this geometry "
                        "(one-off, cached)..."):
            zbank = _zonal_bank(round(d["d_cell"], 4),
                                round(px_row, 4))

        zd = dict(n_rows=int(g["n_rows"]), n_cols=int(g["n_cols"]),
                  pitch=px_row, d_cell=d["d_cell"],
                  h_cell=d["h_cell"], cap_Ah=d["cap_Ah"],
                  r_dc=d["r_dc"], k_dcir=d["k_dcir"], C=z_C,
                  plate_t=d.get("plate_t", 0.0015),
                  plate_contact=d.get("plate_contact", 0.8),
                  s_nom=z_s / 1000.0, T_in=d["T_water_in"],
                  flow_lpm=d["flow_lpm"], tube_od=d["tube_od"],
                  tube_wall=d["tube_wall"],
                  k_tube=K_TUBE[d["tube_mat"]], T_amb=d["T_amb"],
                  h_ext=d["h_ext"], A_case=float(g["A_box_ext"]),
                  nu25=fl["nu25"], B=fl["B"], rho=fl["rho"],
                  cp=fl["cp"], k_oil=fl["k"], beta=fl["beta"],
                  dp_extra=z_dp, T_limit=d["T_limit"],
                  h_contact=float(z_hc), collar_factor=float(z_col),
                  k_rad=d.get("k_rad", 0.9),
                  wetted_tube_frac=(float(z_wet) if z_byp else 0.0),
                  h_oil_tube=150.0)
        nrz, ncz = zd["n_rows"], zd["n_cols"]
        if z_hm == "Busbar end":
            hmap = np.tile(np.linspace(0, z_amp / 100.0, ncz),
                           (nrz, 1))
        elif z_hm == "Centre hot":
            yy, xx = np.mgrid[0:nrz, 0:ncz]
            rr = np.hypot((yy - nrz / 2) / (nrz / 2),
                          (xx - ncz / 2) / (ncz / 2))
            hmap = z_amp / 100.0 * np.exp(-2.5 * rr ** 2)
        else:
            hmap = None

        zkey = f"zres_{json.dumps([round(v, 6) if isinstance(v, float) else v for v in zd.values()])}_{z_hm}_{z_amp}"
        if st.session_state.get("zn_key") != zkey:
            with st.spinner("Solving the zonal network..."):
                zr = solve_zonal(zd, zbank, heat_map=hmap, nz=10,
                                 iters=260)
                zr2 = solve_zonal(zd, zbank, heat_map=hmap, nz=16,
                                  iters=40, init=zr)
                # F3/V8: matched cross-check - lumped serpentine solver
                # at the SAME architecture and the zonal's own velocity.
                u_z = float(zr["ubar"].mean())
                d_match = dict(d, plate_on=True, u_oil=u_z,
                               plate_t=zd["plate_t"],
                               plate_contact=zd["plate_contact"])
                try:
                    res_match = solve_steady(d_match, g, fl, 1.0,
                                             d["T_amb"], C_rate=z_C)
                except Exception:
                    res_match = None
            st.session_state["zn_key"] = zkey
            st.session_state["zn_res"] = (zr, zr2)
            st.session_state["zn_match"] = (res_match, u_z)
        zr, zr2 = st.session_state["zn_res"]
        res_match, u_z = st.session_state.get("zn_match", (None, 0.0))
        lay = zr["lay"]

        zm = st.columns(6)
        zm[0].metric("Hottest can", f"{zr['T_max']:.1f} °C",
                     f"{zr['T_max'] - d['T_limit']:+.1f} vs limit",
                     delta_color="inverse")
        zm[1].metric("Hottest core", f"{zr['T_core_max']:.1f} °C",
                     f"{zr['T_core_max'] - d['T_limit']:+.1f} vs limit",
                     delta_color="inverse",
                     help="Can temperature plus the jellyroll core-to-"
                          "can rise R_core = 1/(4 pi k_r H). If the "
                          "45 °C limit is a plating/core limit, this is "
                          "the number that must clear it.")
        zm[2].metric("Spread", f"{zr['spread']:.2f} °C")
        zm[3].metric("Energy closure",
                     f"{abs(zr['closure']) * 100:.2f} %")
        zm[4].metric("Pump (slots)", f"{zr['P_pump']:.1f} W")
        zm[5].metric("z-resolution check",
                     f"{abs(zr['T_max'] - zr2['T_max']):.02f} °C",
                     "nz 10 vs 16")
        st.caption(f"Derived tube layout from the fin rule P = 2/m, fed "
                   f"the **kernel-implied** plate film "
                   f"h_face = {lay['h_face']:.0f} W/m²·K (a11/pitch at "
                   f"the design point, so the layout is self-consistent "
                   f"with the duct it feeds): m = {lay['m']:.1f} /m, "
                   f"rule P = {lay['P_rule'] * 1000:.0f} mm, snapped to "
                   f"{lay['P_snap'] * 1000:.1f} mm = every "
                   f"{lay['cells_per_tube']} cells, so {lay['n_tubes']} "
                   f"tubes per plate; bay fin efficiency eta = "
                   f"{lay['eta_bay']:.2f} (product rule, conservative by "
                   f"~3-16% vs a 2D plate solve, so T is an upper bound "
                   f"in that respect).")

        if res_match is not None:
            gap = zr["T_mean"] - res_match["T_b"]
            st.caption(
                f"**Independent cross-check (matched):** the lumped "
                f"two-node solver, configured for this same serpentine "
                f"architecture at the zonal's own velocity "
                f"(ū = {u_z * 1000:.0f} mm/s), gives "
                f"T_cell = {res_match['T_b']:.1f} °C against the zonal "
                f"mean {zr['T_mean']:.1f} °C ({gap:+.1f} °C). They use "
                f"genuinely different machinery - the lumped uses "
                f"crossflow cell films and an area-credit plate "
                f"coupling; the zonal routes every watt through "
                f"per-crossing collars and (with the bypass off) omits "
                f"the direct oil-to-tube path - so a residual gap of a "
                f"few °C is expected and brackets the modelling "
                f"uncertainty. Completing each model with the other's "
                f"missing physics (the zonal with the oil-to-tube "
                f"bypass fully on, the lumped with the collar chain "
                f"added) collapses the gap: both land near **38-40 °C** "
                f"at central parameters, so the shipped default with the "
                f"bypass off is most likely ~2-4 °C conservative, and "
                f"the rig point buys back exactly that margin.")

        with st.expander("Contact-conductance sensitivity (the headline "
                         "margin is conditional on h_c)"):
            st.caption("Because the contact toll scales as 1/h_c, the "
                       "sensitivity is hyperbolic: the margin degrades "
                       "fast below the nominal 8000 W/m²·K. Run the "
                       "sweep to see where this design crosses its "
                       "limit.")
            if st.button("Run h_c sweep", key="zn_hcsweep"):
                hc_list = [20000, 12000, 8000, 6000, 4000, 3000, 2000]
                with st.spinner("Sweeping contact conductance..."):
                    rows = []
                    for hc in hc_list:
                        rr = solve_zonal(dict(zd, h_contact=float(hc)),
                                         zbank, heat_map=hmap, nz=10,
                                         iters=200, init=zr)
                        rows.append((hc, rr["T_max"], rr["T_core_max"]))
                st.session_state["zn_hcsw"] = rows
            if "zn_hcsw" in st.session_state:
                rows = st.session_state["zn_hcsw"]
                figh = go.Figure()
                figh.add_trace(go.Scatter(
                    x=[r[0] for r in rows], y=[r[1] for r in rows],
                    name="Hottest can", mode="lines+markers",
                    line=dict(color="#6366F1", width=3)))
                figh.add_trace(go.Scatter(
                    x=[r[0] for r in rows], y=[r[2] for r in rows],
                    name="Hottest core", mode="lines+markers",
                    line=dict(color="#B91C1C", width=3, dash="dot")))
                figh.add_hline(y=d["T_limit"], line_dash="dash",
                               line_color="#B91C1C",
                               annotation_text=f"{d['T_limit']:.0f} °C "
                                               "limit")
                figh.add_vline(x=z_hc, line_dash="dot",
                               line_color="#10B981",
                               annotation_text="current")
                figh.update_layout(
                    height=300, xaxis_title="contact h_c [W/m²·K]",
                    yaxis_title="°C", legend=dict(orientation="h",
                                                  y=1.15),
                    margin=dict(l=10, r=10, t=40, b=10))
                st.plotly_chart(figh, width='stretch', key="zn_hcfig")

                def _cross(vals):
                    xs = np.array([r[0] for r in rows], float)
                    ys = np.array(vals, float)
                    inv = 1.0 / xs
                    for i in range(len(ys) - 1):
                        if (ys[i] - d["T_limit"]) * \
                           (ys[i + 1] - d["T_limit"]) < 0:
                            f = (d["T_limit"] - ys[i]) / (ys[i + 1]
                                                          - ys[i])
                            return 1.0 / (inv[i] + f * (inv[i + 1]
                                                        - inv[i]))
                    return None
                xc = _cross([r[1] for r in rows])
                xco = _cross([r[2] for r in rows])
                msg = []
                if xc:
                    msg.append(f"the **can** clears the limit down to "
                               f"h_c ≈ {xc:.0f} W/m²·K")
                else:
                    msg.append("the **can** clears the limit across "
                               "this whole range")
                if xco:
                    msg.append(f"but the **core** only clears above "
                               f"h_c ≈ {xco:.0f} W/m²·K - within "
                               f"~{100 * (1 - xco / z_hc):.0f}% of the "
                               f"current {z_hc:.0f}")
                st.caption(
                    ". ".join(msg) + ". So on a can criterion the braze "
                    "can be quite imperfect, but on a core/plating "
                    "criterion the design passes only if the joint is "
                    "essentially as good as assumed - which is exactly "
                    "what the one-tube rig point would measure.")

        if z_byp and zr.get("UA_bare", 0) > 0:
            st.caption(f"Direct oil-to-tube path ON: bare-tube "
                       f"UA ≈ {zr['UA_bare']:.0f} W/K carrying "
                       f"{zr['Q_bare']:.0f} W "
                       f"({100 * zr['Q_bare'] / max(zr['Q_gen'], 1):.0f}% "
                       f"of the duty) straight from the oil to the "
                       f"water, bypassing the collar chain. This is the "
                       f"largest known modelling omission when off, and "
                       f"it is conservative-side, so switching it on "
                       f"*lowers* the predicted temperature. Its exact "
                       f"magnitude depends on how much tube length the "
                       f"header leaves wetted by oil.")

        zg1, zg2 = st.columns([3, 2])
        with zg1:
            figz = go.Figure(go.Heatmap(
                z=zr["T_cell"], colorscale="RdYlBu_r",
                colorbar=dict(title="°C", thickness=12)))
            for xk in lay["xk"]:
                figz.add_vline(x=xk / zd["pitch"] - 0.5,
                               line=dict(color="rgba(80,90,120,.55)",
                                         width=1, dash="dot"))
            figz.update_layout(height=380,
                               title="Per-cell temperature map "
                                     "(dotted: derived tube lines)",
                               xaxis_title="column",
                               yaxis_title="row",
                               margin=dict(l=10, r=10, t=50, b=10))
            st.plotly_chart(figz, width='stretch', key="zn_map")
        with zg2:
            share = zr["mdot_col"] / zr["mdot_col"].mean() * 100
            figf = go.Figure(go.Bar(
                y=share, marker_color="#6366F1",
                hovertemplate="channel %{x}: %{y:.0f}%"
                              "<extra></extra>"))
            figf.update_layout(height=185,
                               title="Channel flow share [% of mean]",
                               margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(figf, width='stretch', key="zn_flow")
            figw = go.Figure(go.Bar(
                y=zr["T_wat"][:, -1], marker_color="#38BDF8",
                hovertemplate="tube %{x}: %{y:.2f} °C"
                              "<extra></extra>"))
            figw.update_layout(height=185,
                               title="Water outlet per tube [°C]",
                               margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(figw, width='stretch', key="zn_wat")

        st.markdown("##### Tolerance Monte Carlo")
        zt = st.columns(4)
        z_ss = zt[0].slider("sigma slot [mm]", 0.05, 0.50, 0.25,
                            0.05, key="zn_ss")
        z_sc = zt[1].slider("sigma contact", 0.0, 0.25, 0.10, 0.01,
                            key="zn_sc")
        z_M = zt[2].slider("Samples", 30, 300, 90, 30, key="zn_M")
        run_mc = zt[3].button("Run Monte Carlo", key="zn_run",
                              type="primary")
        if run_mc:
            pb = st.progress(0.0)
            mc = monte_carlo(zd, zbank, M=int(z_M),
                             sigma_s=z_ss / 1000.0, sigma_c=z_sc,
                             heat_map=hmap,
                             progress=lambda f: pb.progress(f))
            pb.empty()
            st.session_state["zn_mc"] = mc
        if "zn_mc" in st.session_state:
            mc = st.session_state["zn_mc"]
            mm = st.columns(4)
            mm[0].metric("Tmax, mean of samples",
                         f"{mc['Tmax'].mean():.2f} °C",
                         f"sd {mc['Tmax'].std():.3f}")
            core_txt = (f"{mc['Tcore'].max():.1f} °C"
                        if "Tcore" in mc else "-")
            mm[1].metric("Worst can / core",
                         f"{mc['Tmax'].max():.2f} °C",
                         f"core {core_txt}")
            if mc.get("n_exceed", 0) == 0 and mc.get("p_ub95"):
                mm[2].metric(f"Exceed {zd['T_limit']:.0f} °C",
                             f"0 / {mc.get('M', len(mc['Tmax']))}",
                             f"≤{mc['p_ub95'] * 100:.1f}% (95% UB)",
                             delta_color="off")
            else:
                mm[2].metric(f"P(exceed {zd['T_limit']:.0f} °C)",
                             f"{mc['p_exceed']:.3f}")
            mm[3].metric("Worst closure",
                         f"{mc['closure_worst'] * 100:.2f} %")
            # verification obs 3: the core is the razor-thin node - show
            # its exceedance count beside the can-based one.
            if "n_exceed_core" in mc:
                n_ec = mc["n_exceed_core"]
                if n_ec == 0 and mc.get("p_ub95_core"):
                    st.success(
                        f"Core node: 0 / {mc['M']} samples over "
                        f"{zd['T_limit']:.0f} °C (≤"
                        f"{mc['p_ub95_core'] * 100:.1f}% 95% UB). Both "
                        f"the can and the core clear the limit under "
                        f"tolerance.")
                else:
                    can_clause = (
                        f"The can still clears comfortably "
                        f"(0 / {mc['M']}), but "
                        if mc.get("n_exceed", 0) == 0 else
                        f"The can exceeds too "
                        f"({mc.get('n_exceed', 0)} / {mc['M']}), and ")
                    st.warning(
                        f"**Core node: {n_ec} / {mc['M']} samples put "
                        f"the core over {zd['T_limit']:.0f} °C** (worst "
                        f"core {mc['Tcore'].max():.2f} °C, point "
                        f"estimate {mc['p_exceed_core'] * 100:.0f}%). "
                        f"{can_clause}on a core/plating criterion the "
                        f"rule of three no longer applies and this "
                        f"design sits at the limit as drawn - see the "
                        f"h_c sweep for why the core margin is thin.")
            zh1, zh2 = st.columns(2)
            fh1 = go.Figure(go.Histogram(x=mc["Tmax"], nbinsx=24,
                                         marker_color="#6366F1"))
            fh1.update_layout(height=220,
                              title="Hottest cell across samples",
                              margin=dict(l=10, r=10, t=40, b=10))
            zh1.plotly_chart(fh1, width='stretch', key="zn_h1")
            fh2 = go.Figure(go.Histogram(x=mc["spread"], nbinsx=24,
                                         marker_color="#F59E0B"))
            fh2.update_layout(height=220,
                              title="Pack spread across samples",
                              margin=dict(l=10, r=10, t=40, b=10))
            zh2.plotly_chart(fh2, width='stretch', key="zn_h2")
            sd = mc["Tmax"].std()
            n_ex = mc.get("n_exceed", int(np.sum(
                mc["Tmax"] > zd["T_limit"])))
            ub_txt = (f"With {n_ex}/{mc.get('M', len(mc['Tmax']))} "
                      f"exceedances the 95% upper bound on the exceed "
                      f"probability is {mc['p_ub95'] * 100:.1f}% (rule "
                      f"of three)" if n_ex == 0 and mc.get("p_ub95")
                      else f"{n_ex} samples exceeded the limit")
            st.caption(
                f"Reading, stated honestly: the sample spread "
                f"(sd {sd:.3f} °C) is at or below the solver's own "
                f"convergence residue at this tolerance, so it should "
                f"be read as *scatter indistinguishable from "
                f"convergence noise, < 0.1 °C*, not a tight physical "
                f"prediction. {ub_txt}. The physics underneath is real "
                f"and is the point: with the lens-shaped slots the flow "
                f"sensitivity to width is nearly linear (the flat-slot "
                f"cubic law does not apply here), and the plates plus "
                f"the mixed plenum homogenise the pack, so slot-width "
                f"and contact *scatter* barely move the hottest cell. "
                f"What the Monte Carlo does **not** sample is the "
                f"epistemic band on the contact conductance itself "
                f"(see the h_c sweep) - that, not tolerance, is where "
                f"the margin risk lives.")
        st.markdown("**Model notes, stated plainly:** cell cans are "
                    "treated as isothermal (high-k casing) and the "
                    "core-to-can rise is superposed on top (reported as "
                    "*Hottest core*); kernels assume laminar developing "
                    "flow (checked: Re stays far below transition); the "
                    "buoyancy head is referenced to the loop return "
                    "(plenum), so at the default the pump dominates and "
                    "the flow is conservative; the water film uses the "
                    "same transition-bridged correlation as the lumped "
                    "solver. Three numbers are engineering choices, not "
                    "derived quantities, and are exposed above: the "
                    "contact conductance h_c (unmeasured pending braze "
                    "data - see the sweep), the collar engagement "
                    "factor, and (off by default) the wetted fraction "
                    "of the direct oil-to-tube path. Radiation is "
                    "excluded (order 1 W). The staircase boundary in "
                    "fea4 was validated against exact limits before "
                    "use.")

        with st.expander("Methods, physics, and how every number is "
                         "reached (equations beside the results)",
                         expanded=True):
            st.markdown("Almost everything below is closed-form or "
                        "provably contractive, and every correlation "
                        "traces to a validated FEA gate. Three inputs "
                        "are honest engineering assumptions rather than "
                        "derived quantities, and the model is upfront "
                        "about them: the contact conductance h_c "
                        "(unmeasured pending braze data), the collar "
                        "engagement factor that sets the contact area, "
                        "and the plate fin efficiency (a conservative "
                        "product rule). Read the chain top to bottom: "
                        "it is exactly what the solver walks for each "
                        "watt.")

            sc1, sc2 = st.columns([3, 2])
            with sc1:
                st.plotly_chart(_S.schematic_unit(zd, g, lay),
                                width='stretch', key="zn_sch_unit")
            with sc2:
                st.plotly_chart(_S.schematic_lens(zd, g),
                                width='stretch', key="zn_sch_lens")

            st.markdown("##### A. Layout from the fin rule")
            st.markdown("The plate is a fin cooled at its rooted "
                        "line and loaded over its face. The fin "
                        "parameter and the optimal root pitch are")
            st.latex(r"m=\sqrt{\dfrac{2\,h_\mathrm{face}}"
                     r"{k_p\,t_p}},\qquad P=\dfrac{2}{m},\qquad "
                     r"\eta_\mathrm{bay}="
                     r"\dfrac{\tanh(mP/2)}{mP/2}\cdot"
                     r"\dfrac{\tanh(mH/2)}{mH/2}")
            st.caption(f"Live: m = {lay['m']:.1f} /m, rule pitch "
                       f"{lay['P_rule']*1000:.0f} mm snapped to the "
                       f"lattice at {lay['P_snap']*1000:.1f} mm "
                       f"(= {lay['cells_per_tube']} cells), giving "
                       f"{lay['n_tubes']} tubes per plate at mid-height "
                       f"and bay efficiency eta = {lay['eta_bay']:.2f}. "
                       "Tubes land on lattice points, never chosen by "
                       "hand.")

            st.markdown("##### B. Hydraulics: laminar slot against "
                        "buoyancy")
            st.markdown("Each slot column balances the available head "
                        "(pump plus thermal buoyancy) against fully "
                        "developed laminar friction; viscosity is "
                        "evaluated at the local bulk temperature "
                        "(Andrade law).")
            st.latex(r"\Delta p=\Delta p_\mathrm{pump}"
                     r"+\rho\,\beta\,g\,H\,(\bar T-T_\mathrm{in}),"
                     r"\qquad \bar u=\Delta p\,"
                     r"\dfrac{2D_h^{2}}{f\!Re\;\mu\,H},\qquad "
                     r"\nu(T)=\nu_{25}\,e^{\,B\left(1/T-1/298\right)}")
            _Dh, _A, _fRe = zbank.props(zd["s_nom"])
            _u = float(zr["ubar"].mean())
            st.caption(f"Live: fRe = {_fRe:.1f} (lens duct, not 96), "
                       f"mean slot velocity {_u*1000:.1f} mm/s, so the "
                       "flow sensitivity to slot width is nearly linear "
                       "here, not cubic (see the choke test in the "
                       "results above).")

            st.markdown("##### C. Duct kernels (fea4): the exact "
                        "developing-flow coupling")
            st.markdown("On the true lens cross-section the velocity "
                        "field is a Poisson solve; the developing "
                        "energy problem is solved twice (one wall hot "
                        "at a time) and superposed into an exact 2x2 "
                        "transport matrix linking each wall's line flux "
                        "to the moving bulk:")
            st.latex(r"\nabla^2 w=-1\ \ (\text{no-slip}),\qquad "
                     r"f\!Re=\dfrac{2D_h^{2}}{\bar w};\qquad "
                     r"\begin{bmatrix}q'_c\\ q'_p\end{bmatrix}"
                     r"=\mathbf{a}(z^{*},s)"
                     r"\begin{bmatrix}T_c-T_b\\ T_p-T_b\end{bmatrix},"
                     r"\quad z^{*}=\dfrac{z}{D_h\,Re\,Pr}")
            _alpha = zd["k_oil"] / (zd["rho"] * zd["cp"])
            _RePr = _u * _Dh / _alpha
            _zop = (zd["h_cell"] * 0.5) / (_Dh * max(_RePr, 1e-9))
            _a = zbank.a_of(zd["s_nom"], np.array([_zop]))[0]
            st.latex(r"\mathbf{a}(z^{*}\!=\!%.2g)=\begin{bmatrix}"
                     r"%.1f & %.2f\\ %.2f & %.1f\end{bmatrix}\ "
                     r"\mathrm{W\,m^{-1}K^{-1}}" %
                     (_zop, _a[0, 0], _a[0, 1], _a[1, 0], _a[1, 1]))
            st.caption("Validation gates (parallel-plate limits with "
                       "known answers): fRe 95.7 vs 96.0, one-wall "
                       "Nu 5.387 vs 5.385, two-wall Nu 7.541 vs 7.541. "
                       "The whole 70 mm slot runs in the thermal "
                       "entrance (Pr about 127), so the cell film is "
                       "roughly 3x the stagnant-thermosiphon value.")

            st.markdown("##### D. Channel energy march (exact "
                        "exponential, energy-consistent)")
            st.markdown("Marching bulk enthalpy up each segment, the "
                        "walls act through the conductance-weighted mix "
                        "T_eff. The update is the exact solution of the "
                        "linear segment (unconditionally stable), and "
                        "the segment-mean driver uses the exact factor "
                        "f-bar so the flux and the enthalpy rise agree "
                        "to machine precision at any flow, including "
                        "starved channels.")
            st.latex(r"\dot m c_p\dfrac{dT_b}{dz}=G\,(T_\mathrm{eff}"
                     r"-T_b),\quad G=\textstyle\sum a,\ "
                     r"T_\mathrm{eff}=\dfrac{G_cT_c+G_pT_p}{G}")
            st.latex(r"T_b(z{+}dz)=T_\mathrm{eff}+(T_b-T_\mathrm{eff})"
                     r"e^{-x},\quad \bar f=\dfrac{1-e^{-x}}{x},\quad "
                     r"x=\dfrac{G\,dz}{\dot m c_p}")

            st.markdown("##### E. Cells: heat that fights back "
                        "(DCIR feedback)")
            st.markdown("Each cell is an isothermal balance; generation "
                        "falls as it warms. It is solved by Newton "
                        "using the march's own Jacobian J, so the cell "
                        "and its channels are consistent, not lagged.")
            st.latex(r"q_\mathrm{gen}(T)=(C\,\mathrm{Ah})^{2}R_{dc}\,"
                     r"e^{-k_{dc}(T-25)}=\hat q+J\,(T-T_\mathrm{old})")

            st.markdown("##### F. Plate face and root network "
                        "(closed form)")
            sc3, sc4 = st.columns([2, 3])
            with sc3:
                st.plotly_chart(_S.schematic_network(zd, lay),
                                width='stretch', key="zn_sch_net")
            with sc4:
                st.markdown("The plate face is not iterated. Its heat "
                            "is linear in its own temperature, so it is "
                            "solved directly against the series root "
                            "resistance (contact, then tube wall, then "
                            "the water film), and the water is marched "
                            "tube by tube afterwards.")
                st.latex(r"R_\mathrm{root}=\dfrac{1}{h_cA_\mathrm{ct}}"
                         r"+\dfrac{\ln(d_o/d_i)}{2\pi k_t L}"
                         r"+\dfrac{1}{h_w\pi d_i L}")
                st.latex(r"T_f=\dfrac{T_w/R_\mathrm{root}-S_0}"
                         r"{W+1/R_\mathrm{root}},\qquad "
                         r"q_\mathrm{root}=\dfrac{T_f-T_w}"
                         r"{R_\mathrm{root}}")
                st.caption("W and S0 are accumulated from the same "
                           "march (per bay), so the plate closure is "
                           "exact given the bulk field. Water film h_w "
                           "from Gnielinski (turbulent) or Hausen "
                           "(laminar developing).")

            st.markdown("##### G. Recirculation plenum (closed-form "
                        "fixed point)")
            st.markdown("The mixed exit re-enters as the common inlet. "
                        "Because each column exit is affine in the "
                        "inlet, the recirculation fixed point is solved "
                        "in one line instead of ratcheted, which is "
                        "what removed the slow mode that used to stall "
                        "convergence.")
            st.latex(r"T_\mathrm{exit}=A+P\,T_\mathrm{plen},\quad "
                     r"P=\prod_z e^{-x}\ \Rightarrow\ "
                     r"T_\mathrm{plen}^{*}=\dfrac{\langle A\rangle}"
                     r"{1-\langle P\rangle}")

            st.markdown("##### H. Tolerance Monte Carlo")
            st.latex(r"s_i\sim\mathcal N(s_\mathrm{nom},\sigma_s)\ "
                     r"[\text{trunc.}],\quad c_k\sim"
                     r"\mathcal N(c,\sigma_c),\quad "
                     r"P_\mathrm{exceed}=\Pr\!\big(T_\mathrm{max}"
                     r">T_\mathrm{lim}\big)")
            st.caption("Every sample is a full network solve, "
                       "warm-started from the converged nominal case "
                       "(30 s for 60-90 samples). The result is the "
                       "physics of homogenisation: the plates and the "
                       "mixed plenum flatten manufacturing scatter, so "
                       "the water film at the roots, not tolerance, "
                       "governs the temperature level.")

            st.markdown("##### Validation battery (this run)")
            _match_txt = (
                f"the lumped serpentine solver **matched** to this "
                f"velocity **{res_match['T_b']:.1f} °C** "
                f"({zr['T_mean'] - res_match['T_b']:+.1f} °C; different "
                f"machinery, so a few °C is expected)"
                if res_match is not None else
                "an independent lumped solver in the matched mode")
            st.markdown(
                f"- Energy closure in vs out: "
                f"**{abs(zr['closure'])*100:.2f} %** "
                f"(generation {zr['Q_gen']:.0f} W = water "
                f"{zr['Q_water']:.0f} W + casing {zr['Q_case']:.0f} W).\n"
                f"- z-resolution independence: nz 10 vs 16 differ by "
                f"**{abs(zr['T_max']-zr2['T_max']):.02f} °C**.\n"
                f"- Independent cross-check: zonal mean "
                f"**{zr['T_mean']:.1f} °C** against {_match_txt}.\n"
                f"- Collapse test: with uniform inputs the pack spread "
                f"goes to zero, confirming the per-channel bookkeeping "
                f"introduces no spurious asymmetry.")

    # ---------------- Improve ---------------- #
    with tabs[5]:
        with st.container(border=True):
            st.markdown("#### Predictor - what happens if...")
            st.caption("Play with the levers; nothing is saved to the design. "
                       "Every slider move re-solves instantly.")
            p1, p2, p3, p4, p5 = st.columns(5)
            pf_flow = p1.slider("Flow [L/min]", 0.5, 60.0, float(d["flow_lpm"]), 0.5)
            pf_tubes = p2.slider("Tubes", 1, 60, int(d["n_tubes"]))
            pf_u = p3.slider("Oil velocity [m/s]", 0.0, 0.15, float(d["u_oil"]), 0.005)
            pf_tin = p4.slider("Inlet [°C]", 0.0, 40.0, float(d["T_water_in"]), 1.0)
            pf_C = p5.slider("C-rate", 0.2, 6.0, float(round(C_steady, 1)), 0.1)
            dp_ = dict(d, flow_lpm=pf_flow, n_tubes=int(pf_tubes), u_oil=pf_u,
                       T_water_in=pf_tin)
            gp_ = build_geometry(dp_)
            rp_ = solve_steady(dp_, gp_, fl, 1.0, d["T_amb"], C_rate=pf_C)
            okp = rp_["T_b"] <= d["T_limit"]
            kpi_cards([
                ("Predicted can / core",
                 f"{rp_['T_b']:.1f} / {rp_['T_core']:.1f}°C",
                 f"{rp_['T_b']-res['T_b']:+.1f} °C vs current",
                 "ok" if okp else "bad"),
                ("Margin to limit", f"{d['T_limit']-rp_['T_b']:+.1f} °C",
                 f"at {pf_C:.1f}C continuous", "ok" if okp else "bad"),
                ("Cell spread (max-min)", f"{rp_['spread']:.1f} °C",
                 f"{rp_['spread']-res['spread']:+.1f} vs current", ""),
                ("Heat", f"{rp_['Q_eff']/1000:.2f} kW",
                 f"chiller ~{chiller_model(rp_['Q_w'], pf_tin, d['T_amb'])['P_el']:.0f} W el",
                 ""),
            ])
        cpin1, cpin2, _sp = st.columns([1, 1, 4])
        summary = dict(Fluid=fl["name"], T_can=round(res["T_b"], 1),
                       T_core=round(res["T_core"], 1), Cmax=round(Cmax, 2),
                       Spread_K=round(res["spread"], 1), Pack_kg=round(masses["m_pack"]),
                       Whkg=round(masses["whkg_pack"]), WhL=round(masses["whl_pack"]),
                       Parasitic_W=round(P_pump + P_stir, 1), Tubes=d["n_tubes"],
                       Flow_lpm=d["flow_lpm"], Stir=d["u_oil"])
        if cpin1.button("Pin design as A"):
            st.session_state["pinned"] = summary
        if "pinned" in st.session_state and cpin2.button("Clear pin"):
            del st.session_state["pinned"]
        if "pinned" in st.session_state:
            st.dataframe(pd.DataFrame([st.session_state["pinned"], summary],
                                      index=["A (pinned)", "B (current)"]),
                         use_container_width=True)
        with st.container(border=True):
            st.markdown("#### Scale the system to a target C-rate")
            st.caption("Pick a continuous C and see what each subsystem "
                       "needs - chiller, flow, tubes, stirring, busbars, "
                       "core, charging.")
            C_t = st.slider("Target continuous C", 0.5, 6.0,
                            float(min(max(d["C1"], 0.5), 6.0)), 0.5)
            if st.button("Analyse target C"):
                with st.spinner("Solving current design and single-lever "
                                "fixes at target C..."):
                    sc = scale_to_C(d, g, fl, masses, d["T_amb"], C_t)
                st.dataframe(pd.DataFrame(sc["rows"],
                             columns=["Subsystem", "At this design", "Change needed / note"]),
                             hide_index=True, use_container_width=True)
            st.plotly_chart(c_sweep_fig(d, g, fl, d["T_amb"], C_steady),
                            use_container_width=True, key="imp_csweep")
        improve_core(d, g, fl, masses, res, Q_duty, C_steady, cool_df)

    # ---------------- Ideas ---------------- #
    with tabs[6]:
        st.caption("Concepts tried against the live design. Baseline = this "
                   "design with no forced circulation. Adopt a winner via "
                   "Design - Circulation method.")
        with st.container(border=True):
            st.markdown("#### Idea A - serpentine plate circulation (yours)")
            st.markdown(
                "Thin plates hang from the tubes between cell rows: they act "
                "as **conduction fins** (both faces wetted, no water inside) "
                "and guide a **small pump** to push oil along every row "
                "channel - forced flow exactly where the gaps are tightest, "
                "immune to the sub-6 mm gap penalty. Channels are **parallel "
                "and manifolded**; a single long serpentine would heat the "
                "oil tens of degrees end to end.")
            i1, i2, i3 = st.columns(3)
            iu = i1.slider("Channel velocity [m/s]", 0.01, 0.15, 0.05, 0.005,
                           key="id_u")
            it = i2.select_slider("Plate thickness [mm]", [1.0, 1.5, 2.0], 1.5,
                                  key="id_t")
            ic = i3.slider("Plate-tube contact", 0.4, 1.0, 0.8, 0.05, key="id_c")
            base = dict(d, plate_on=False, u_oil=0.0)
            stir = dict(d, plate_on=False, u_oil=iu)
            serp = dict(d, plate_on=True, u_oil=iu, plate_t=it / 1000,
                        plate_contact=ic)
            rows_i, bars = [], []
            for nm, dd_ in [("Baseline (thermosiphon)", base),
                            ("Open stirring", stir),
                            ("Serpentine plates", serp)]:
                r_ = solve_steady(dd_, g, fl, 1.0, d["T_amb"], C_rate=C_steady)
                if dd_.get("plate_on"):
                    P_ = serpentine_pump(dd_, g, fl, iu)["P"]
                    dm_ = plate_fin_area(dd_, g, r_["h_tube"])[2]
                else:
                    P_ = stirrer_power(dd_, g, fl, dd_["u_oil"])
                    dm_ = 0.0
                rows_i.append(dict(Concept=nm, T_can=round(r_["T_b"], 1),
                                   T_core=round(r_["T_core"], 1),
                                   Spread_K=round(r_["spread"], 1),
                                   Pump_W=round(P_, 1),
                                   Added_kg=round(dm_, 1),
                                   A_sink_m2=round(r_["A_oilside"], 1)))
                bars.append((nm, r_["T_b"], P_, dm_))
            idf = pd.DataFrame(rows_i)
            st.dataframe(idf, hide_index=True, use_container_width=True)
            figI = go.Figure()
            names_i = [b[0] for b in bars]
            figI.add_trace(go.Bar(name="Cell T [°C]", x=names_i,
                                  y=[b[1] for b in bars],
                                  marker_color="#F0655F",
                                  offsetgroup=1,
                                  text=[f"{b[1]:.1f}" for b in bars],
                                  textposition="outside"))
            figI.add_trace(go.Bar(name="Pump [W] (right axis)", x=names_i,
                                  y=[b[2] for b in bars],
                                  marker_color="#7C88F8", yaxis="y2",
                                  offsetgroup=2,
                                  text=[f"{b[2]:.1f}" for b in bars],
                                  textposition="outside"))
            pmax = max(max(b[2] for b in bars), 0.1)
            figI.update_layout(barmode="group", height=340,
                               yaxis=dict(title="Cell T [°C]"),
                               yaxis2=dict(title="Pump [W]",
                                           overlaying="y", side="right",
                                           range=[0, pmax * 1.5],
                                           showgrid=False),
                               title=f"Concepts at {C_steady:.2f}C rms "
                                     f"(limit {d['T_limit']:.0f} °C)",)
            st.plotly_chart(figI, use_container_width=True, key="ideas_bar")
            dT_serp = rows_i[0]["T_can"] - rows_i[2]["T_can"]
            dT_stir = rows_i[0]["T_can"] - rows_i[1]["T_can"]
            st.markdown(
                f"**Verdict at this operating point:** serpentine buys "
                f"**{dT_serp:.1f} °C** over baseline for "
                f"{rows_i[2]['Pump_W']:.1f} W and {rows_i[2]['Added_kg']:.0f} kg "
                f"of plates (open stirring buys {dT_stir:.1f} °C for "
                f"{rows_i[1]['Pump_W']:.1f} W). The plates also cut the spread "
                f"({rows_i[2]['Spread_K']:.1f} vs {rows_i[0]['Spread_K']:.1f} °C) "
                "because the sink is distributed. Risks: plate-to-tube joint "
                "quality (the contact slider), channel blockage, manifold "
                "design, and serviceability.")

    # ---------------- Safety ---------------- #
    with tabs[7]:
        st.markdown("Order-of-magnitude screening plus the engineering checklist. "
                    "Nothing here replaces abuse testing.")
        runaway_ui(d, g, fl, masses, res)
        exp_L = fl["beta"] * masses["V_oil_L"] * (d["T_service_max"] - 0.0)
        st.markdown(f"""
**Expansion and leaks.**
* Oil expansion over a 0 to {d['T_service_max']:.0f} °C service band:
  **~{exp_L:.1f} L** on {masses['V_oil_L']:.0f} L - bellows or bladder, never free air.
* Hold **oil pressure above water pressure** so any tube leak goes oil-to-water
  (transformer practice); or double-walled tubes with interstitial leak detection.
* Burst disc set at {d['p_des_bar']:.1f} bar g (the enclosure is sized for this in
  Design); route the vent away from occupants and expect oil ejection.
* Ester fluids: monitor moisture; copper: use inhibited oil or plated tubes.""")

    # ---------------- Compare ---------------- #
    with tabs[8]:
        st.subheader("Architectures")
        arch_tab(d, g, fl, masses, C_steady)
        st.markdown("---")
        st.subheader("Coolant shoot-out")
        coolant_tab(d, g, cool_df)
        st.markdown("---")
        st.subheader("Cooling architecture trade study")
        st.caption("The same cells and duty under six cooling approaches - "
                   "full immersion (three circulation variants), partial "
                   "immersion, indirect cold plate, and forced air. Immersion "
                   "rows use the full solver; the dry architectures use "
                   "closed-form conduction chains (bottom-cooled axial "
                   "gradient, pad and plate resistances). Costs are the "
                   "thermal system only, from the editable assumptions below "
                   "(cells themselves: BNEF Dec 2025 average $79/kWh cell, "
                   "$99/kWh BEV pack; Europe typically +56%).")
        with st.expander("Cost assumptions (edit and the table updates)"):
            ce = {}
            cc1, cc2, cc3, cc4 = st.columns(4)
            ce["oil_gbp_L"] = cc1.number_input("Dielectric [£/L]", 1.0, 40.0,
                float(cost_now["oil_gbp_L"]), 0.5, key="c_oil")
            ce["chiller_gbp_kWel"] = cc2.number_input("Chiller [£/kW el]",
                40.0, 400.0, float(cost_now["chiller_gbp_kWel"]), 10.0,
                key="c_chil")
            ce["coldplate_gbp"] = cc3.number_input("Cold plate [£]", 80.0,
                800.0, float(cost_now["coldplate_gbp"]), 20.0, key="c_cp")
            ce["fab_factor"] = cc4.number_input("Fabrication factor", 1.0,
                4.0, float(cost_now["fab_factor"]), 0.1, key="c_fab")
            st.session_state["cost_edit"] = ce
        adf = arch_df.copy()
        adf["verdict"] = [
            ("PASS" if (h <= d["T_limit"] and s <= 5) else
             "hotspot" if h > d["T_limit"] and s <= 5 else
             "uniformity" if h <= d["T_limit"] else "both")
            for h, s in zip(adf["T_hot"], adf["spread"])]
        show_a = adf[["name", "T_can", "T_hot", "spread", "maxC", "mass",
                      "parasitic", "chiller_el", "cost_th", "verdict",
                      "note"]].round(1)
        show_a.columns = ["Architecture", "T can [°C]", "Hotspot [°C]",
                          "Spread [°C]", "Max cont. C", "Pack [kg]",
                          "Parasitic [W]", "Chiller [W el]",
                          "Thermal cost [£]", f"At {C_steady:.1f}C", "Note"]
        st.dataframe(show_a, hide_index=True, use_container_width=True)
        ca1, ca2 = st.columns(2)
        with ca1:
            st.plotly_chart(arch_chart(arch_df, d["T_limit"]),
                            use_container_width=True, key="arch_hot")
        with ca2:
            figM = go.Figure()
            figM.add_trace(go.Bar(name="Pack mass [kg]",
                x=[n.split(" - ")[-1].split(" (")[0] for n in arch_df["name"]],
                y=arch_df["mass"], marker_color="#94A3B8", offsetgroup=1,
                text=arch_df["mass"].round(0), textposition="outside"))
            figM.add_trace(go.Bar(name="Thermal cost [£]",
                x=[n.split(" - ")[-1].split(" (")[0] for n in arch_df["name"]],
                y=arch_df["cost_th"], marker_color="#5BC8E8", yaxis="y2",
                offsetgroup=2, text=arch_df["cost_th"].round(0),
                textposition="outside"))
            figM.update_layout(barmode="group", height=400,
                               yaxis=dict(title="kg"),
                               yaxis2=dict(title="£", overlaying="y",
                                           side="right", showgrid=False),
                               title="What each approach weighs and costs "
                                     "(thermal system)")
            st.plotly_chart(figM, use_container_width=True, key="arch_mass")
        best_u = arch_df.loc[arch_df["spread"].idxmin(), "name"]
        light = arch_df.loc[arch_df["mass"].idxmin(), "name"]
        st.markdown(
            f"**Reading it honestly:** the indirect cold plate is the mass "
            f"and cost winner but its axial-conduction spread fails the 5 °C "
            f"criterion; forced air is a sub-2C architecture; partial "
            f"immersion cooks the dry cell tops. Full immersion pays "
            f"{masses['m_oil']:.0f} kg of oil and the burst-rated enclosure "
            f"to buy uniformity ({best_u.split(' - ')[-1]} is best at "
            f"{arch_df['spread'].min():.1f} °C), abuse tolerance and the "
            f"thermal flywheel; {light} is lightest at "
            f"{arch_df['mass'].min():.0f} kg. Dry-architecture rows are "
            f"closed-form estimates (±20-30%); immersion rows are the full "
            f"validated solver.")
        st.markdown("---")
        st.subheader("BEV road-car pack database")
        fig_bm = benchmark_db_tab(d, g, masses, Cmax, res)
        st.markdown("---")
        st.subheader("Teardown cooling references")
        bench_prod_tab(masses, Cmax)

    # ---------------- Learn ---------------- #
    with tabs[9]:
        learn_tab(d, g, fl, res, masses, cool_df, loop,
                  Q_duty, chil)

    # ---------------- Report ---------------- #
    figs_r = dict(
        sankey=heat_sankey(res, Q_duty, Q_bus, P_pump, P_stir, chil),
        ladder=thermal_xray_fig(d, g, fl, res, Q_duty),
        resist=resistance_chart(res)[0],
        transient=None,   # filled below from the Duty chart data
        csweep=c_sweep_fig(d, g, fl, d["T_amb"], C_steady),
        setpoint=setpoint_trade(d, g, fl, C_steady, d["T_amb"]),
        pack3d=pack_3d_figure(d, g),
        archs=arch_chart(arch_df, d["T_limit"]),
    )
    figT_r = go.Figure()
    figT_r.add_trace(go.Scatter(x=tr["t"]/60, y=tr["T_core"], name="Core",
                                line=dict(color="#B91C1C", width=2, dash="dot")))
    figT_r.add_trace(go.Scatter(x=tr["t"]/60, y=tr["T_b"], name="Cell",
                                line=dict(color="#EF4444", width=3)))
    figT_r.add_trace(go.Scatter(x=tr["t"]/60, y=tr["T_il"], name="Oil",
                                line=dict(color=ACCENT, width=3)))
    figT_r.add_hline(y=d["T_limit"], line_dash="dash", line_color="#B91C1C")
    figT_r.update_layout(height=320, xaxis_title="min",
                         yaxis_title="°C", plot_bgcolor="rgba(255,255,255,0)",
                         paper_bgcolor="rgba(0,0,0,0)",)
    figs_r["transient"] = figT_r
    secs = report_sections(d, g, fl, res, masses, tr, spec, Cmax, C_steady,
                           Q_duty, Q_bus, P_pump, P_stir, chil, figs_r,
                           arch_df=arch_df, cost=cost_now)
    import datetime as _dt
    meta_r = dict(version=APP_VERSION,
                  date=_dt.date.today().strftime("%d %b %Y"),
                  spec=(f"{masses['E_kwh']:.1f} kWh / "
                        f"{d['Ns']*d['v_nom']:.0f} V - {g['N']} x "
                        f"{d['fmt']} in "
                        f"{fl['name'].split('(')[0].strip()} - {d['duty']}"),
                  ok=ok, T_gov=T_gov, T_limit=d["T_limit"],
                  kwh=masses["E_kwh"], mass=masses["m_pack"], Cmax=Cmax,
                  chil_el=chil["P_el"])
    with tabs[13]:
        system_tab(d, g, fl, res, masses, loop, chil,
                   Q_duty, C_steady)

    with tabs[14]:
        st.caption("Take this report with you: the HTML file is "
                   "instant; Build renders every section, table and "
                   "figure into real Word, PowerPoint and PDF files "
                   "of the design exactly as it stands now.")
        cbt, cfg, cbo = st.columns([1.1, 0.9, 1.1])
        cbt.download_button("Download this report (.html)",
                            data=export_report_html(secs, figs_r, meta_r),
                            file_name="pack_design_report.html",
                            mime="text/html", use_container_width=True)
        inc_figs = cfg.checkbox("Include figures", True, key="rx_figs",
                                help="Figures render via the kaleido "
                                     "package (in requirements.txt).")
        if cbo.button("Build Word / PowerPoint / PDF", key="rx_build",
                      type="primary", use_container_width=True):
            import report_export as _RX
            with st.spinner("Rendering figures and laying out "
                            "documents..."):
                _pngs = _RX.figs_to_png(figs_r) if inc_figs else {}
                if inc_figs and not _pngs:
                    st.info("kaleido not available - exporting "
                            "text-only.")
                st.session_state["rx_files"] = dict(
                    docx=_RX.build_docx(secs, meta_r, _pngs),
                    pptx=_RX.build_pptx(secs, meta_r, _pngs),
                    pdf=_RX.build_pdf(secs, meta_r, _pngs))
        if "rx_files" in st.session_state:
            _fx = st.session_state["rx_files"]
            d1, d2, d3 = st.columns(3)
            d1.download_button(
                "Word report (.docx)", data=_fx["docx"],
                file_name="pack_design_report.docx",
                mime="application/vnd.openxmlformats-officedocument"
                     ".wordprocessingml.document",
                use_container_width=True, key="rx_dl_docx")
            d2.download_button(
                "PowerPoint deck (.pptx)", data=_fx["pptx"],
                file_name="pack_design_report.pptx",
                mime="application/vnd.openxmlformats-officedocument"
                     ".presentationml.presentation",
                use_container_width=True, key="rx_dl_pptx")
            d3.download_button(
                "PDF report (.pdf)", data=_fx["pdf"],
                file_name="pack_design_report.pdf",
                mime="application/pdf",
                use_container_width=True, key="rx_dl_pdf")
        render_report_tab(secs, figs_r, meta_r)

    with tabs[15]:
        battery_tab()

    # ---------------- Validate and tune ---------------- #
    with tabs[10]:
        cases_tab(d, g, fl, res, cool_df, loop)

    with tabs[11]:
        fea_tab(d, g, fl, cool_df, loop)

    with tabs[12]:
        st.subheader("Benchmark: Wang et al. 2023")
        bench_wang_tab()
        st.markdown("---")
        st.subheader("Calibrate to your rig")
        mup = st.file_uploader("Measured CSV: columns t_s, T_cell", type="csv", key="meas_up")
        if mup is not None:
            try:
                mdf = pd.read_csv(mup)
                cols = {c.lower().strip(): c for c in mdf.columns}
                st.session_state["meas_data"] = (
                    mdf[cols["t_s"]].to_numpy(float), mdf[cols["t_cell"]].to_numpy(float))
                st.success(f"Loaded {len(mdf)} points.")
            except Exception as e:
                st.error(f"Could not parse: {e}")
        meas = st.session_state.get("meas_data")
        if meas is not None:
            figM = go.Figure()
            figM.add_trace(go.Scatter(x=tr["t"]/60, y=tr["T_b"], name="Model",
                                      line=dict(color="#EF4444", width=3)))
            figM.add_trace(go.Scatter(x=meas[0]/60, y=meas[1], mode="markers",
                                      name="Measured", marker=dict(color=INK, size=5)))
            figM.update_layout(height=300, xaxis_title="min", yaxis_title="°C",
                               plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)")
            st.plotly_chart(figM, use_container_width=True)
            if st.button("Fit oil-film calibration factor"):
                with st.spinner("Sweeping 0.4-2.2..."):
                    cb, rm = fit_calibration(d, g, fl, masses, d["T_amb"],
                                             tr["t"], np.abs(tr["C"]), meas[0], meas[1])
                st.session_state["_pending"] = dict(w_cal_h=cb)
                st.success(f"cal_h = {cb:.2f} (RMSE {rm:.2f} °C) - applied below.")
                st.rerun()
        st.markdown("---")
        st.subheader("FEA studies (pack_fea_v1)")
        f1, f2, f3 = st.columns(3)
        with f1:
            figF1 = go.Figure(go.Bar(
                x=["Flat lid (FE)", "Best grillage (FE)", "App assumption"],
                y=[25.9, 14.6, 11.7], marker_color=["#94A3B8", "#6366F1", "#F59E0B"],
                text=["25.9", "14.6", "11.7"], textposition="outside"))
            figF1.update_layout(height=260, title="Lid mass [kg] at 0.5 bar g",
                                margin=dict(l=10, r=10, t=44, b=10))
            st.plotly_chart(figF1, use_container_width=True, key="fea1")
            st.caption("Roark validated to +0.0%. Honest knock-down 0.56 vs "
                       "the app's 0.45 (composite rib action would close the "
                       "gap) - consider 0.50.")
        with f2:
            figF2 = go.Figure()
            figF2.add_trace(go.Bar(name="App formula", x=["21700 2C", "4680 2C"],
                                   y=[2.85, 14.0], marker_color="#F59E0B"))
            figF2.add_trace(go.Bar(name="FE (mandrel + wetted ends)",
                                   x=["21700 2C", "4680 2C"], y=[2.28, 9.30],
                                   marker_color="#6366F1"))
            figF2.update_layout(barmode="group", height=260,
                                title="Core-to-can ΔT [K]",
                                margin=dict(l=10, r=10, t=44, b=10))
            st.plotly_chart(figF2, use_container_width=True, key="fea2")
            st.caption("Exact-solution validation -0.1%. The app is 20-34% "
                       "conservative on real cells.")
        with f3:
            figF3 = go.Figure(go.Bar(
                x=["App estimate", "Section model, co-flow", "Counterflow"],
                y=[3.25, 1.08, 0.06],
                marker_color=["#F59E0B", "#6366F1", "#10B981"],
                text=["3.25", "1.08", "0.06"], textposition="outside"))
            figF3.update_layout(height=260,
                                title="Water-rise cell spread [K]",
                                margin=dict(l=10, r=10, t=44, b=10))
            st.plotly_chart(figF3, use_container_width=True, key="fea3")
            st.caption("Energy closure 0.00%. Counterflow plumbing of "
                       "alternate tubes removes the term for free.")
        st.caption("Full self-validating scripts and figures: pack_fea_v1.zip "
                   "(fea1_lid.py, fea2_cell.py, fea3_pack.py, REPORT.md).")
        st.markdown("---")
        st.subheader("FEA4: plate-channel duct kernels (fea4_channel)")
        st.caption("The engine behind the Zones tab. Fully developed "
                   "laminar analysis on the true lens cross-section: a "
                   "velocity Poisson solve for the friction factor and a "
                   "two-case developing-flow (Graetz) march for the exact "
                   "2x2 wall-to-bulk transport matrix. Verified against "
                   "three limits with known closed-form or tabulated "
                   "answers before it was trusted anywhere.")

        @st.cache_resource(show_spinner=False)
        def _fea4_gates():
            import fea4_channel as F4
            v = F4.duct_coefficients(rect=True, s=0.002,
                                     pitch_x=0.02, n=180)
            Nu1 = (2 * 0.002) / (v["Rpp"] * 0.13 * v["Pp"])
            gg = F4.graetz_kernel(rect=True, s=0.002, pitch_x=0.02,
                                  n=120, nz=90)
            aFD = gg["a"][-1]
            Nu2 = ((aFD[0, 0] + aFD[0, 1]) * (2 * 0.002)
                   / (0.13 * gg["Pc"]))
            return v["fRe"], Nu1, Nu2

        with st.spinner("Running duct-kernel gates (cached)..."):
            _fRe, _Nu1, _Nu2 = _fea4_gates()
        g1, g2 = st.columns([3, 2])
        with g1:
            figF4 = go.Figure()
            figF4.add_trace(go.Bar(
                name="fea4 (this model)",
                x=["fRe (channel)", "Nu, one wall heated",
                   "Nu, both walls fixed-T"],
                y=[_fRe, _Nu1, _Nu2], marker_color="#6366F1",
                text=[f"{_fRe:.1f}", f"{_Nu1:.3f}", f"{_Nu2:.3f}"],
                textposition="outside"))
            figF4.add_trace(go.Bar(
                name="analytic target",
                x=["fRe (channel)", "Nu, one wall heated",
                   "Nu, both walls fixed-T"],
                y=[96.0, 5.385, 7.541], marker_color="#94A3B8",
                text=["96.0", "5.385", "7.541"],
                textposition="outside"))
            figF4.update_layout(
                barmode="group", height=300,
                title="Parallel-plate validation limits: model vs "
                      "analytic",
                legend=dict(orientation="h", y=1.12),
                margin=dict(l=10, r=10, t=54, b=10))
            st.plotly_chart(figF4, use_container_width=True, key="fea4")
            st.caption(f"All three land on target: fRe "
                       f"{_fRe:.1f} vs 96.0, one-wall Nu {_Nu1:.3f} vs "
                       f"5.385, two-wall Nu {_Nu2:.3f} vs 7.541. Two "
                       "operator bugs (a masked-cell energy term and a "
                       "wall-flux extrapolation) were caught precisely "
                       "because these gates failed first.")
        with g2:
            st.plotly_chart(_S.schematic_lens(
                dict(d_cell=0.021, pitch=0.0215, s_nom=0.002),
                dict(gap_mm=0.5)), use_container_width=True,
                key="fea4_lens")
        st.caption("Self-contained module with the gates as a runnable "
                   "__main__: fea4_channel.py. The earlier stagnant-slot "
                   "study returned an arc effectiveness near 0.13, which "
                   "is what redirected the concept to advective primary "
                   "cooling (oil carries, plates sink).")
        st.markdown("---")
        st.subheader("Model tuning")
        t1_, t2_, t3_ = st.columns(3)
        with t1_:
            _w(st.slider, "Oil-film calibration factor", "cal_h", 1.0, min_value=0.4, max_value=2.2, step=0.01)
            _w(st.slider, "Thermosiphon minor-loss K", "kloop", 5.0, min_value=0.0, max_value=20.0, step=0.5)
        with t2_:
            _w(st.slider, "Casing-to-ambient h [W/m²·K]", "hext", 5.0, min_value=0.0, max_value=25.0, step=0.5)
            _w(st.slider, "Gap under cells [mm]", "bgap", 5.0, min_value=0.0, max_value=20.0, step=1.0)
        with t3_:
            _w(st.slider, "Manifold margin [mm]", "mman", 20.0, min_value=0.0, max_value=60.0, step=5.0)
            _w(st.slider, "Max service T (expansion) [°C]", "tserv", 60.0, min_value=40.0, max_value=90.0, step=5.0)
        if st.checkbox("Edit fluid property table", False):
            st.session_state.cool_df = st.data_editor(cool_df, num_rows="dynamic", height=260)
        st.caption("These knobs change the physics for every tab. cal_h multiplies both "
                   "oil films; fit it above once rig data exists (spray-app pattern).")

# ------------------------------------------------------------------ #
#  v5: chiller model, heat Sankey, station cards, scale-to-C          #
# ------------------------------------------------------------------ #
def chiller_model(Q_w, T_w_in, T_amb, eta=0.45, approach=5.0, cond_dT=10.0):
    """Vapour-compression estimate: COP = eta_carnot x T_evap / lift.
    Evaporator = inlet - approach; condenser = ambient + cond_dT."""
    T_c = T_w_in - approach + 273.15
    T_h = T_amb + cond_dT + 273.15
    lift = max(T_h - T_c, 3.0)
    COP = max(eta * T_c / lift, 0.4)
    return dict(COP=COP, P_el=Q_w / COP, lift=lift)

def heat_sankey(res, Q_duty, Q_bus, P_pump, P_stir, chil):
    """Left = sources, right = sinks. Ribbon width = watts."""
    kw = lambda v: f"{v/1000:.2f} kW" if v >= 100 else f"{v:.0f} W"
    Qc = max(Q_duty - Q_bus, 1e-3)
    lab = [f"Cells  {kw(Qc)}", f"Busbars  {kw(Q_bus)}",
           f"Bulk oil  {kw(Q_duty)}", f"Water loop  {kw(res['Q_w'])}",
           f"Casing loss  {kw(res['Q_atm'])}",
           f"Chiller  {kw(res['Q_w']+P_pump)}",
           f"To ambient  {kw(res['Q_w']+P_pump+chil['P_el'])}",
           f"Chiller electricity  {kw(chil['P_el'])}"]
    node_c = ["#EF4444", "#D97706", "#F59E0B", "#0EA5E9", "#94A3B8",
              "#6366F1", "#64748B", "#10B981"]
    x = [0.01, 0.01, 0.30, 0.55, 0.72, 0.78, 0.99, 0.55]
    y = [0.36, 0.86, 0.45, 0.38, 0.93, 0.44, 0.40, 0.78]
    srcs = [0, 1, 2, 2, 3, 7, 5]
    dsts = [2, 2, 3, 4, 5, 5, 6]
    vals = [Qc, max(Q_bus, 1e-3), max(res["Q_w"], 1e-3),
            max(res["Q_atm"], 1e-3), max(res["Q_w"] + P_pump, 1e-3),
            max(chil["P_el"], 1e-3),
            max(res["Q_w"] + P_pump + chil["P_el"], 1e-3)]
    def rgba(hx, a=0.32):
        return f"rgba({int(hx[1:3],16)},{int(hx[3:5],16)},{int(hx[5:7],16)},{a})"
    fig = go.Figure(go.Sankey(
        arrangement="fixed", valueformat=".0f", valuesuffix=" W",
        node=dict(label=lab, x=x, y=y, pad=26, thickness=24, color=node_c,
                  line=dict(color="rgba(255,255,255,.9)", width=1),
                  hovertemplate="%{label}<extra></extra>"),
        link=dict(source=srcs, target=dsts, value=vals,
                  color=[rgba(node_c[s]) for s in srcs],
                  hovertemplate="%{source.label}  →  %{target.label}"
                                "<br>%{value:.0f} W<extra></extra>"),
        textfont=dict(family="Inter, -apple-system, sans-serif", size=14,
                      color="#0F172A")))
    fig.update_layout(height=340, margin=dict(l=10, r=10, t=52, b=16),
                      title="Heat flow map - every watt from the cells to "
                            "the air outside")
    return fig


def thermal_circuit_fig(d, g, fl, res, Q):
    """The cell-to-water chain as an interactive circuit. Box height = share
    of total resistance; colour = transport mode; hover = the governing law
    with live numbers; the weakest link is flagged."""
    el = [
        ("Cell interior", "conduction", res["dT_core"] / max(Q, 1e-9),
         f"Radial conduction through the jellyroll<br>"
         f"R = 1/(4*pi*k_r*H) = {r_core(d)*1000:.1f} mK/W per cell<br>"
         f"k_r = {d['k_rad']:.1f} W/m*K -> DT = {res['dT_core']:.1f} °C"),
        ("Can to oil film", "convection", res["R_b"],
         f"Natural/forced convection film<br>h = {res['h_cell']:.0f} W/m²·K "
         f"over A = {g['A_cells']:.1f} m²<br>R = 1/(hA) = "
         f"{res['R_b']*1000:.2f} mK/W -> DT = {Q*res['R_b']:.1f} °C"),
        ("Oil to tube+fins" + ("+plates" if d.get("plate_on") else ""),
         "convection", res["R_ot"],
         f"Convection onto the finned sink<br>h = {res['h_tube']:.0f} W/m²·K, "
         f"A_eff = {res['A_oilside']:.1f} m²"
         + (f" (plates add {res.get('A_plate',0):.1f} m2, eta "
            f"{res.get('eta_plate',0):.2f})" if d.get("plate_on") else "")
         + f"<br>R = {res['R_ot']*1000:.2f} mK/W -> DT = {Q*res['R_ot']:.1f} °C"),
        ("Tube wall", "conduction", res["R_wall"],
         f"Conduction through {d['tube_mat'].lower()}<br>"
         f"R = ln(r_o/r_i)/(2*pi*k*L*n) = {res['R_wall']*1000:.3f} mK/W "
         f"-> DT = {Q*res['R_wall']*1000:.0f} mK"),
        ("Water film", "convection", res["R_in"],
         f"{res['water_regime']} internal flow, Re = {res['Re_water']:.0f}<br>"
         f"h = {res['h_water']:.0f} W/m²·K -> R = {res['R_in']*1000:.2f} mK/W "
         f"-> DT = {Q*res['R_in']:.1f} °C"),
        ("Water stream", "advection", 0.5 * res["dT_water"] / max(Q, 1e-9),
         f"Heat carried away: DT = Q/(m_dot*c_p)<br>rise "
         f"{res['dT_water']:.1f} °C inlet to outlet (+{0.5*res['dT_water']:.1f} "
         "K at the mean cell)"),
    ]
    Rtot = sum(e[2] for e in el)
    films = [e for e in el if e[1] == "convection"]
    weakest = max(films, key=lambda e: e[2])[0]
    mode_c = dict(conduction="#94A3B8", convection="#6366F1", advection="#0EA5E9")
    fig = go.Figure()
    x = 0.0
    xw = 1.0 / len(el)
    hover_x, hover_y, hover_t = [], [], []
    for name, mode, R, tip in el:
        share = R / max(Rtot, 1e-12)
        h = 0.16 + 0.74 * share
        is_w = name == weakest
        fig.add_shape(type="rect", x0=x + 0.06 * xw, x1=x + 0.94 * xw,
                      y0=0.5 - h / 2, y1=0.5 + h / 2,
                      line=dict(color="#EF4444" if is_w else "#E2E8F0",
                                width=3 if is_w else 1),
                      fillcolor=mode_c[mode], opacity=0.9 if is_w else 0.75)
        fig.add_annotation(x=x + 0.5 * xw, y=0.5 + h / 2 + 0.07,
                           text=f"<b>{Q*R:.1f} °C</b>" if Q * R >= 0.05
                                else f"{Q*R*1000:.0f} mK",
                           showarrow=False, font=dict(size=12, color="#0F172A"))
        fig.add_annotation(x=x + 0.5 * xw, y=0.5 - h / 2 - 0.08,
                           text=name.replace(" ", "<br>", 1), showarrow=False,
                           font=dict(size=10.5, color="#475569"))
        if is_w:
            fig.add_annotation(x=x + 0.5 * xw, y=0.97, text="WEAKEST LINK",
                               showarrow=False,
                               font=dict(size=10, color="#B91C1C"))
        hover_x.append(x + 0.5 * xw); hover_y.append(0.5)
        hover_t.append(f"<b>{name}</b> - {mode}<br>{tip}<br>"
                       f"share of chain: {100*share:.0f}%")
        if x > 0:
            pass
        x += xw
    for i in range(len(el) - 1):
        fig.add_annotation(x=(i + 1) * xw, y=0.5, text="", showarrow=True,
                           ax=(i + 0.94) * xw, ay=0.5, axref="x", ayref="y",
                           arrowhead=2, arrowwidth=2, arrowcolor="#CBD5E1")
    fig.add_trace(go.Scatter(x=hover_x, y=hover_y, mode="markers",
                             marker=dict(size=42, opacity=0.0),
                             customdata=list(range(len(hover_x))),
                             hovertext=hover_t, hoverinfo="text",
                             showlegend=False))
    for mode, c in mode_c.items():
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                 marker=dict(size=12, color=c, symbol="square"),
                                 name=mode))
    fig.update_layout(height=330, xaxis=dict(visible=False, range=[0, 1]),
                      yaxis=dict(visible=False, range=[0, 1.05]),
                      legend=dict(orientation="h", y=-0.05),
                      margin=dict(l=10, r=10, t=54, b=10),
                      title=f"The thermal circuit at {Q/1000:.2f} kW: box "
                            "height = share of resistance, hover for the law")
    return fig

def nu_T_fig(fl):
    Ts = np.linspace(0, 80, 60)
    nus = [film_props(fl, t)["nu"] * 1e6 for t in Ts]
    fig = go.Figure(go.Scatter(x=Ts, y=nus, line=dict(color="#F59E0B", width=3)))
    fig.update_layout(height=210, xaxis_title="Oil temperature [°C]",
                      yaxis_title="ν [cSt]",
                      title=f"{fl['name']}: viscosity vs temperature "
                            "(Andrade law - the solver uses this at film T)",
                      margin=dict(l=10, r=10, t=44, b=10))
    return fig


def thermal_xray_fig(d, g, fl, res, Q):
    """The flagship: a schematic slice from cell core to water, drawn as
    physical strata (hardware), with the live temperature profile T(x)
    crossing them (results) and every drop annotated with its mechanism and
    law (physics). Widths are schematic; temperatures are real."""
    T_wm = d["T_water_in"] + 0.5 * res["dT_water"]
    T_si = T_wm + Q * res["R_in"]              # wall inner (water side)
    T_so = T_si + Q * res["R_wall"]            # wall outer (oil side)
    T_il = T_so + Q * res["R_ot"]              # bulk oil
    T_b, T_core = res["T_b"], res["T_core"]
    films = [("can-oil", res["R_b"]), ("oil-tube", res["R_ot"]),
             ("water", res["R_in"])]
    weak = max(films, key=lambda x: x[1])[0]

    Z = [  # (key, width, fill, label, glyph)
        ("cell",  .17, "rgba(51,65,85,0.92)",  "Cell", "cyl"),
        ("film1", .07, "rgba(245,158,11,0.38)", "First oil film", None),
        ("oil",   .20, "rgba(245,158,11,0.16)", "Bulk oil", "wave"),
        ("film2", .07, "rgba(245,158,11,0.38)", "Second oil film", None),
        ("fin",   .06, "rgba(165,180,204,0.55)", "Fins + tube", "comb"),
        ("wall",  .025, "rgba(148,163,184,0.9)", "", None),
        ("film3", .06, "rgba(14,165,233,0.30)", "Water film", None),
        ("water", .16, "rgba(14,165,233,0.14)", "Water stream", "drop"),
    ]
    xs, x0 = {}, 0.0
    for k, w, *_ in Z:
        xs[k] = (x0, x0 + w); x0 += w
    T_lo = min(d["T_water_in"], T_wm) - 2
    T_hi = max(T_core, T_b) + 3

    fig = go.Figure()
    for k, w, c, lab, gly in Z:
        a, b = xs[k]
        fig.add_shape(type="rect", x0=a, x1=b, y0=T_lo, y1=T_hi,
                      fillcolor=c, line=dict(width=0), layer="below")
        if lab:
            fig.add_annotation(x=(a + b) / 2, y=T_lo + 0.5, text=lab,
                               showarrow=False, yanchor="bottom",
                               font=dict(size=10.5, color="#475569"))
    # water inlet-outlet band + mean
    a, b = xs["water"]
    fig.add_shape(type="rect", x0=a, x1=b, y0=d["T_water_in"],
                  y1=d["T_water_in"] + res["dT_water"],
                  fillcolor="rgba(14,165,233,0.25)", line=dict(width=0))
    fig.add_annotation(x=(a + b) / 2, y=d["T_water_in"], yanchor="top",
                       text=f"in {d['T_water_in']:.0f} °C", showarrow=False,
                       font=dict(size=10, color="#0369A1"))
    fig.add_annotation(x=(a + b) / 2,
                       y=d["T_water_in"] + res["dT_water"], yanchor="bottom",
                       text=f"out {d['T_water_in']+res['dT_water']:.1f} °C",
                       showarrow=False, font=dict(size=10, color="#0369A1"))

    # temperature profile
    px_, py_, tips = [], [], []
    a, b = xs["cell"]
    for i in range(13):                      # parabola core -> can
        xi = a + (b - a) * i / 12
        py_.append(T_b + (T_core - T_b) * (1 - (i / 12) ** 2))
        px_.append(xi)
        tips.append(f"Inside the cell (conduction, k_r={d['k_rad']:.1f} "
                    f"W/m·K)<br>T = {py_[-1]:.1f} °C")
    def seg(k, Ta, Tb_, tip, n=2):
        a_, b_ = xs[k]
        for i in range(n + 1):
            px_.append(a_ + (b_ - a_) * i / n)
            py_.append(Ta + (Tb_ - Ta) * i / n)
            tips.append(tip + f"<br>T = {py_[-1]:.1f} °C")
    seg("film1", T_b, T_il,
        f"First oil film - convection, h = {res['h_cell']:.0f} W/m²·K over "
        f"{g['A_cells']:.1f} m² -> {Q*res['R_b']:.1f} °C")
    seg("oil", T_il, T_il, "Bulk oil - the mixer; nearly no drop")
    seg("film2", T_il, T_so,
        f"Second oil film - convection onto the finned sink, "
        f"h = {res['h_tube']:.0f}, A_eff = {res['A_oilside']:.1f} m² -> "
        f"{Q*res['R_ot']:.1f} °C")
    seg("fin", T_so, T_so, "Fin and tube metal - conduction, negligible")
    seg("wall", T_so, T_si,
        f"Tube wall - conduction, {Q*res['R_wall']*1000:.0f} m°C")
    seg("film3", T_si, T_wm,
        f"Water film - {res['water_regime']} flow, Re = "
        f"{res['Re_water']:.0f}, h = {res['h_water']:.0f} -> "
        f"{Q*res['R_in']:.1f} °C")
    seg("water", T_wm, T_wm,
        f"Water stream carries {res['Q_w']/1000:.2f} kW away "
        f"(rise {res['dT_water']:.1f} °C)")
    fig.add_trace(go.Scatter(x=px_, y=py_, mode="lines",
                             line=dict(color="#1F2937", width=3.5),
                             hovertext=tips, hoverinfo="text",
                             name="temperature"))
    # drop labels + weakest flag
    for k, Ta, Tb_, nm in [("film1", T_b, T_il, "can-oil"),
                           ("film2", T_il, T_so, "oil-tube"),
                           ("film3", T_si, T_wm, "water")]:
        a_, b_ = xs[k]
        red = nm == weak
        fig.add_annotation(x=(a_ + b_) / 2, y=(Ta + Tb_) / 2,
                           text=(f"<b>{abs(Ta-Tb_):.1f} °C</b>"
                                 + ("<br>WEAKEST" if red else "")),
                           showarrow=False,
                           font=dict(size=11,
                                     color="#B91C1C" if red else "#1F2937"),
                           bgcolor="rgba(255,255,255,.75)")
    fig.add_annotation(x=xs["cell"][0] + 0.012, y=T_core,
                       text=f"core {T_core:.1f}", showarrow=False,
                       xanchor="left", font=dict(size=10, color="#FCA5A5"))
    fig.add_annotation(x=xs["cell"][1], y=T_b, text=f"can {T_b:.1f}",
                       showarrow=False, xanchor="left",
                       font=dict(size=10, color="#7C2D12"))
    fig.add_hline(y=d["T_limit"], line_dash="dash", line_color="#B91C1C",
                  annotation_text=f"limit {d['T_limit']:.0f} °C")
    fig.update_layout(height=470,
                      xaxis=dict(visible=False, range=[0, 1]),
                      yaxis=dict(title="Temperature [°C]",
                                 range=[T_lo, T_hi + 2]),
                      showlegend=False,
                      title="Thermal X-ray - the hardware slice, the live "
                            "temperature profile, and every toll (hover the "
                            "line)")
    return fig


def sankey_deep_dive(d, g, fl, res, masses, Q_duty, Q_bus, P_pump, P_stir,
                     chil, C_steady):
    """Extensive live commentary on the heat-flow map. Every number is
    computed from the current solve; two extra solves (1C, 4C) feed the
    scaling discussion."""
    Q = Q_duty
    Qc = Q - Q_bus
    T_wm = d["T_water_in"] + 0.5 * res["dT_water"]
    R_chain = res["R_ot"] + res["R_wall"] + res["R_in"]
    UA_w = 1.0 / max(R_chain, 1e-9)
    mdot = d["flow_lpm"] / 60 * 1.0
    mcp = mdot * 4180.0
    NTU = UA_w / max(mcp, 1e-9)
    eff = 1.0 - math.exp(-NTU)
    C_th = masses["C_oil"] + masses["C_batt"]
    tau_min = C_th * (res["R_b"] + R_chain) / 60.0
    dT_case = max(res["T_il"] - d["T_amb"], 0.1)
    UA_case = res["Q_atm"] / dT_case
    resid = Q - (res["Q_w"] + res["Q_atm"])
    P_par = P_pump + P_stir
    sysCOP = res["Q_w"] / max(chil["P_el"] + P_par, 1e-9)
    q_cell = Q / g["N"]
    q_flux = Q / g["A_cells"] / 10.0  # mW/cm2
    r_hot = r_of_T(d, res["T_b"]); r_cold = d["r_dc"]
    ch10 = chiller_model(res["Q_w"], 10.0, d["T_amb"])
    ch30 = chiller_model(res["Q_w"], 30.0, d["T_amb"])
    r1 = solve_steady(d, g, fl, 1.0, d["T_amb"], C_rate=1.0)
    r4 = solve_steady(d, g, fl, 1.0, d["T_amb"], C_rate=4.0)
    weak = max([("the water film", res["R_in"]),
                ("the oil-to-tube film", res["R_ot"]),
                ("the can-to-oil film", res["R_b"])], key=lambda x: x[1])[0]
    dry_ok = d["T_water_in"] >= d["T_amb"] + 5.0
    md = f"""
#### Reading the heat-flow map in depth

**The encoding first.** Ribbon width is power in watts, nothing else: no
temperatures, no efficiencies, just the first law drawn to scale. Every
node must balance - what enters equals what leaves, because at steady
state nothing inside the pack is warming up or cooling down. Colour
follows the physical stream: red for heat born in the cells, amber while
the oil carries it, blue once the water has it, indigo through the
chiller, green for electricity entering, grey for anything leaving to the
room. If you ever see a node where the widths visibly do not balance, the
model is broken; here the ledger closes to {abs(resid):.1f} W
({100*abs(resid)/Q:.2f}% of the total), which is numerical rounding, not
missing physics.

**Cells - {Qc/1000:.2f} kW.** This is Joule heat: pack current squared
times internal resistance. At the duty's {C_steady:.2f}C RMS each of the
{g['N']} cells dissipates {q_cell:.2f} W, a surface intensity of about
{q_flux:.0f} mW/cm² - gentle by electronics standards, which is exactly
why the thin-film convection problem, not boiling or flux limits,
dominates immersion design at this C-rate. Two scaling laws hide in this
one ribbon. First, heat grows with the **square** of C-rate: halve the
duty to 1C and this solver says the ribbon shrinks to
{r1['Q_eff']/1000:.2f} kW; push to 4C and it becomes
{r4['Q_eff']/1000:.2f} kW at an (over-limit) equilibrium of
{r4['T_b']:.1f} °C. Second, the DCIR feedback softens everything: at the
current {res['T_b']:.1f} °C the cells' resistance has fallen from
{r_cold:.1f} to {r_hot:.1f} mΩ, so the pack makes
{100*(1-r_hot/r_cold):.0f}% less heat than the same current would generate
at 25 °C. A hot pack partially protects itself; a cold one is doubly
punished. This is why the 4C ribbon is {r4['Q_eff']/max(r1['Q_eff'],1):.1f}x
the 1C ribbon rather than the naive 16x.

**Busbars - {Q_bus:.0f} W.** Sized automatically at {d['bus_J']:.0f}
A/mm², they contribute {100*Q_bus/Q:.1f}% of the total. The ribbon is
almost invisible and that is the design intent: copper is cheap insurance
here. But note it carries the same C² law as the cells - at 4C this
ribbon alone becomes ~{Q_bus*(4/max(C_steady,1e-6))**2:.0f} W, which is
why the scale-to-C tool resizes the section.

**Bulk oil - the junction, {Q/1000:.2f} kW through.** Every watt passes
through this node because every heat-making surface is submerged; the oil
is not a bypass, it is the only road. Two properties matter and neither
is visible as width. It is the **mixer**: circulation at
{max(d['u_oil'], res['u_ts'])*1000:.1f} mm/s smears the heat sideways so
the best-to-worst cell spread stays at {res['spread']:.1f} °C. And it is
the **flywheel**: {C_th/1000:.0f} kJ/K of thermal mass gives the pack a
first-order time constant of roughly {tau_min:.0f} minutes, which is why
launch-and-recover transients in the Duty tab barely dent the steady
picture - the Sankey is the time-average this buffer relaxes to. One
honest omission: the {P_stir:.1f} W of circulation work also ends up in
the oil as heat, but at {100*P_stir/Q:.2f}% of the load it is not drawn.

**Water loop - {res['Q_w']/1000:.2f} kW.** The oil hands
{100*res['Q_w']/Q:.0f}% of the heat to the internal tube bank. The chain
resistance (oil film + wall + water film) gives UA =
{UA_w:.0f} W/K; against the water's capacity rate of {mcp:.0f} W/K that
is NTU = {NTU:.2f}, an exchanger effectiveness of {100*eff:.0f}% - the
water leaves having closed {100*eff:.0f}% of the gap to the oil
temperature. The stream warms only {res['dT_water']:.1f} °C from inlet to
outlet, which tells you transport is easy and the fight is in the films:
right now the biggest single toll in the chain is **{weak}**. A crucial
misreading to avoid: doubling the water flow does **not** widen this
ribbon. At steady state the ribbon width is set by how much heat the
cells make, full stop; better cooling lowers the *temperatures* needed to
drive the same watts across. The Sankey shows the traffic, the thermal
X-ray shows the price paid to move it.

**Casing loss - {res['Q_atm']/1000:.2f} kW.** The enclosure sheds
{100*res['Q_atm']/Q:.0f}% of the load straight to the room, driven by the
{dT_case:.1f} °C between bulk oil and ambient through an effective
UA of {UA_case:.1f} W/K (natural convection plus a little radiation,
lumped). It is free cooling, so resist the instinct to insulate the box;
this design wants the leak. The share grows in a cold room and vanishes -
or reverses - in a hot one: park at 45 °C ambient and this ribbon flips
direction, heat soaking *in* through the walls while the chiller pays for
it.

**Chiller - {(res['Q_w']+P_pump)/1000:.2f} kW duty in,
{chil['P_el']/1000:.2f} kW electricity in, everything out.** The loop
water returns warm and must be re-cooled to {d['T_water_in']:.0f} °C, and
because that set point sits {'below' if not dry_ok else 'above'} ambient
this is genuine refrigeration: a vapour-compression circuit whose
evaporator chills your water, whose compressor does the work you see
entering as the green ribbon, and whose condenser dumps the sum to the
air. First law at the node: {res['Q_w']/1000:.2f} kW of heat plus
{P_pump:.0f} W of pump work plus {chil['P_el']/1000:.2f} kW of
electricity equals the {(res['Q_w']+P_pump+chil['P_el'])/1000:.2f} kW
leaving to ambient - the room receives *more* heat than the battery made,
the surcharge being the compressor's wage. The COP of {chil['COP']:.1f}
(estimated at 45% of Carnot on the actual {chil['lift']:.0f} °C lift)
means each electrical watt moves {chil['COP']:.1f} thermal watts. The set
point moves this violently: at a 10 °C inlet the same duty would cost
{ch10['P_el']:.0f} W of electricity (COP {ch10['COP']:.1f}); at 30 °C it
falls to {ch30['P_el']:.0f} W (COP {ch30['COP']:.1f}) - and the warmer
inlet *also* shrinks the cell ribbon through DCIR, the rare double win,
bounded only by the {d['T_limit']:.0f} °C cell limit.
{"Since your set point is already at or above ambient plus a sensible approach, a plain dry cooler (radiator and fan) could replace the refrigeration circuit entirely." if dry_ok else f"Raise the set point to about {d['T_amb']+5:.0f} °C or beyond and the compressor becomes unnecessary - a dry cooler (radiator and fan) rejects the heat for fan power alone, deleting the {chil['P_el']:.0f} W, the refrigerant circuit, and a slice of the cost table."}

**The efficiency ledger.** Total cooling overhead - chiller electricity
plus pump plus circulation - is {(chil['P_el']+P_par)/1000:.2f} kW against
{res['Q_w']/1000:.2f} kW removed: a system cooling COP of
**{sysCOP:.1f}**, or {100*(chil['P_el']+P_par)/max(Q,1e-9):.0f}% of the
battery's own heat spent on moving it. For perspective, the
{(res['Q_w']+P_pump+chil['P_el'])/1000:.2f} kW arriving at ambient is a
domestic fan heater running flat out - worth remembering for lab
ventilation when the rig runs at duty for hours.

**What this map deliberately cannot tell you.** No temperatures (the
X-ray's job), no uniformity (the section model's job), no transients (the
Duty tab's job - the flywheel means peak ribbons never appear here), and
the chiller is a physics estimate, not a datasheet. Read it for *where
the watts go and what they cost to move*; read its neighbours for
*how hot anything gets*.
"""
    return md

def station_list(d, g, fl, res, masses, Q_duty, Q_bus, P_pump, P_stir, chil,
                 cf_tip=True):
    """The heat journey as (name, dT_or_metric, numbers_md, improve_md)."""
    Q = Q_duty
    dTs = dict(core=res["dT_core"], f1=Q * res["R_b"], f2=Q * res["R_ot"],
               wall=Q * res["R_wall"], f3=Q * res["R_in"],
               rise=0.5 * res["dT_water"])
    tot = sum(dTs.values())
    S = []
    S.append(("1. Heat source - inside the cell",
        dTs["core"],
        f"**{Q/1000:.2f} kW** at {res['T_b']:.1f} °C "
        f"(DCIR {r_of_T(d, res['T_b']):.1f} mΩ vs {d['r_dc']:.0f} at 25; the "
        f"hot pack makes {100*(1-r_of_T(d,res['T_b'])/d['r_dc']):.0f}% less heat). "
        f"Busbars add {Q_bus:.0f} W. Core runs **{res['dT_core']:.1f} °C** above "
        f"the can (R = 1/(4 pi k_r H)).",
        f"Lower-DCIR cells cut heat at the source; a warmer set point does too "
        f"(the AMG 45 °C logic). Busbars: thicker section (now sized at "
        f"{d['bus_J']:.0f} A/mm²). The core term ignores the mandrel and wetted "
        f"ends - FEA shows the real figure is ~20% lower (21700) to ~34% (4680); "
        f"only lower current or axial extraction move it further."))
    S.append(("2. Can to oil - the first film",
        dTs["f1"],
        f"h_cell = **{res['h_cell']:.0f} W/m²·K** over {g['A_cells']:.1f} m² -> "
        f"**{dTs['f1']:.1f} °C**. Gap factor {res['gapf']:.2f} at "
        f"{g['gap_mm']:.1f} mm; flow component from "
        f"{'stirring' if d['u_oil']>res['u_ts'] else 'thermosiphon'} "
        f"({max(d['u_oil'],res['u_ts'])*1000:.1f} mm/s).",
        "Levers in order: keep the gap at or above 6 mm; stir a few cm/s "
        f"({stirrer_power(d,g,fl,0.05):.1f} W buys the biggest single film gain); "
        "a lower-viscosity fluid adds ~x1.3 at most (Ra ~ 1/nu, h ~ Ra^0.25)."))
    S.append(("3. Bulk oil - the mixer and flywheel",
        0.0,
        f"Self-circulation **{res['u_ts']*1000:.1f} mm/s**, stratification "
        f"~{res['dT_loop']:.1f} °C, best-to-worst cell spread "
        f"**{res['spread']:.1f} °C** (criterion 5 °C). Thermal buffer "
        f"{(masses['C_oil']+masses['C_batt'])/1e3:.0f} kJ/K absorbs peaks.",
        "Keep the cold plane high (tube placement drives the loop head). "
        + ("**Plumb alternate tubes in opposite directions**: the section model "
           "shows counterflow collapses the water-rise spread from ~1.1 K to "
           "~0.06 K, free. " if cf_tip else "")
        + "Holders currently block "
        f"{d['holder_block']*100:.0f}% of the riser - open them up."))
    S.append(("4. Oil to tube - the finned second film",
        dTs["f2"],
        f"h_tube = **{res['h_tube']:.0f} W/m²·K**; fins eta "
        f"{res['fin']['eta']:.2f}, area x{res['fin']['area_gain']:.1f} -> "
        f"A_eff {res['A_oilside']:.1f} m² -> **{dTs['f2']:.1f} °C**.",
        "Fin harder (oil's low h keeps long fins ~90% efficient), add tubes, "
        "or go interstitial to distribute the sink. This film shares the "
        "same stirring lever as station 2."))
    S.append(("5. Tube wall",
        dTs["wall"],
        f"{d['tube_mat']}, {d['tube_wall']*1000:.1f} mm -> "
        f"**{dTs['wall']*1000:.0f} mK**. Negligible.",
        f"Material is free here: aluminium instead of copper adds only "
        f"~{dTs['wall']*1000*(K_TUBE['Copper']/K_TUBE['Aluminium']-1):.0f} mK. "
        "Choose tubes for corrosion and joining, not conductivity."))
    S.append(("6. Water film and flow",
        dTs["f3"] + dTs["rise"],
        f"**{res['water_regime']}**, Re = {res['Re_water']:.0f}, h = "
        f"{res['h_water']:.0f} -> film **{dTs['f3']:.1f} °C**; inlet-to-outlet "
        f"rise {res['dT_water']:.1f} °C (mean +{dTs['rise']:.1f} °C). Pump "
        f"{P_pump:.1f} W.",
        "In laminar flow more velocity does nothing - trip turbulence "
        "(Re > 3000) with more flow or fewer/narrower parallel paths, or "
        "split the water into counterflowing halves to also fix uniformity."))
    S.append(("7. Chiller - closing the loop",
        0.0,
        f"Duty **{(res['Q_w']+P_pump)/1000:.2f} kW** at {d['T_water_in']:.0f} °C "
        f"inlet; est. COP **{chil['COP']:.1f}** (lift {chil['lift']:.0f} °C) -> "
        f"electrical **{chil['P_el']:.0f} W**. Casing sheds {res['Q_atm']:.0f} W "
        "for free.",
        "Warmer inlet is a double win up to the cell limit: COP rises and "
        "DCIR(T) cuts the heat itself - see the set-point trade curve below. "
        "Size the chiller for continuous duty and let the oil flywheel eat "
        "the peaks."))
    return S, tot

def setpoint_trade(d, g, fl, C_duty, T_amb):
    rows = []
    for Tin in (5, 10, 15, 20, 25, 30):
        dd = dict(d); dd["T_water_in"] = float(Tin)
        r = solve_steady(dd, g, fl, 1.0, T_amb, C_rate=C_duty)
        ch = chiller_model(r["Q_w"], Tin, T_amb)
        rows.append((Tin, r["T_b"], r["Q_eff"], ch["P_el"], ch["COP"]))
    a = np.array(rows)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=a[:, 0], y=a[:, 3], name="Chiller electrical [W]",
                             line=dict(color="#6366F1", width=3)))
    fig.add_trace(go.Scatter(x=a[:, 0], y=a[:, 1], name="Cell T [°C]",
                             yaxis="y2", line=dict(color="#EF4444", width=3)))
    fig.add_trace(go.Scatter(x=[a[0, 0], a[-1, 0]], y=[d["T_limit"]] * 2,
                             yaxis="y2", mode="lines", name="limit",
                             line=dict(color="#B91C1C", dash="dash")))
    fig.update_layout(height=300, xaxis_title="Water inlet set point [°C]",
                      yaxis_title="Chiller electrical [W]",
                      yaxis2=dict(title="Cell T [°C]", overlaying="y",
                                  side="right", showgrid=False),
                      title=f"Set-point trade at {C_duty:.2f}C: COP vs DCIR",
                      plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
                      margin=dict(l=10, r=10, t=50, b=10))
    return fig

def c_sweep_fig(d, g, fl, T_amb, C_now):
    Cs = np.linspace(0.5, 6.0, 12)
    Tb, Qs, Pel = [], [], []
    for C in Cs:
        r = solve_steady(d, g, fl, 1.0, T_amb, C_rate=float(C))
        Tb.append(r["T_b"]); Qs.append(r["Q_eff"] / 1000)
        Pel.append(chiller_model(r["Q_w"], d["T_water_in"], T_amb)["P_el"] / 1000)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Cs, y=Qs, name="Heat [kW]",
                             line=dict(color="#EF4444", width=3)))
    fig.add_trace(go.Scatter(x=Cs, y=Pel, name="Chiller electrical [kW]",
                             line=dict(color="#6366F1", width=3)))
    fig.add_trace(go.Scatter(x=Cs, y=Tb, name="Cell T [°C]", yaxis="y2",
                             line=dict(color="#F59E0B", width=3)))
    fig.add_trace(go.Scatter(x=[Cs[0], Cs[-1]], y=[d["T_limit"]] * 2,
                             yaxis="y2", mode="lines", name="limit",
                             line=dict(color="#B91C1C", dash="dash")))
    fig.add_vline(x=C_now, line_dash="dot", annotation_text=f"duty {C_now:.2f}C")
    fig.update_layout(height=340, xaxis_title="Continuous C-rate",
                      yaxis_title="kW",
                      yaxis2=dict(title="Cell T [°C]", overlaying="y",
                                  side="right", showgrid=False),
                      title="Scaling the same hardware with C "
                            "(heat ~ C^2, softened by DCIR(T))",
                      plot_bgcolor="rgba(255,255,255,0)", paper_bgcolor="rgba(0,0,0,0)",
                      margin=dict(l=10, r=10, t=60, b=10))
    return fig

def scale_to_C(d, g, fl, masses, T_amb, C_t):
    """What must change, subsystem by subsystem, to run C_t continuously."""
    r0 = solve_steady(d, g, fl, 1.0, T_amb, C_rate=C_t)
    ok0 = r0["T_b"] <= d["T_limit"]
    ch = chiller_model(r0["Q_w"], d["T_water_in"], T_amb)
    I_pack = C_t * d["cap_Ah"] * d["Np"]
    bus = busbar_props(dict(d, C1=C_t), g)
    rows = [("Heat to remove", f"{r0['Q_eff']/1000:.2f} kW at {r0['T_b']:.1f} °C "
             f"({'OK' if ok0 else 'OVER LIMIT'})",
             "size everything below for this, continuously")]
    rows.append(("Chiller", f"duty {(r0['Q_w'])/1000:.2f} kW, COP {ch['COP']:.1f} "
                 f"-> {ch['P_el']/1000:.2f} kW electrical",
                 "or raise the set point and re-check (COP + DCIR both improve)"))
    fixes = []
    for name, (key, lo, hi, isint) in GOAL_LEVERS.items():
        inv = key == "T_water_in"
        x, T = goal_seek(d, fl, key, lo, hi, isint, T_amb, C_t, d["T_limit"], inv)
        fixes.append((name, x, T))
        cur = d[key] * (1000 if key == "pitch" else 1)
        shown = "insufficient alone" if x is None else \
                f"{x*1000:.1f}" if key == "pitch" else f"{x:.1f}"
        rows.append((f"Single lever: {name}",
                     f"now {cur:.1f}", f"needs {shown}"
                     + ("" if x is None else f" -> {T:.1f} °C")))
    # stirred combination if static fails everywhere useful
    dd = dict(d, u_oil=max(d["u_oil"], 0.05))
    r_st = solve_steady(dd, g, fl, 1.0, T_amb, C_rate=C_t)
    rows.append(("Combination: stir 5 cm/s",
                 f"{r_st['T_b']:.1f} °C ({stirrer_power(d,g,fl,0.05):.1f} W)",
                 "then re-run the single levers from there"))
    rows.append(("Busbars", f"{I_pack:.0f} A pack current",
                 f"section {bus['A_mm2']:.0f} mm² at {d['bus_J']:.0f} A/mm2, "
                 f"{bus['m']:.1f} kg, adds {(I_pack**2*bus['R']):.0f} W"))
    rows.append(("Water side", f"Re = {r0['Re_water']:.0f} ({r0['water_regime']})",
                 "target Re > 3000; laminar extra flow is wasted"))
    rows.append(("Cell core", f"+{r0['dT_core']:.1f} °C core-to-can at {C_t:.1f}C "
                 "(formula; FEA says ~20% less for a real 21700)",
                 "no coolant touches this - format and current only"))
    chg_ok = plating_frac(25.0) * d["C_chg"]
    rows.append(("Charging at this C", f"cell rating {d['C_chg']:.1f}C CC; at "
                 f"{C_t:.1f}C the cell must be rated for it and warm "
                 f"(map allows {plating_frac(15)*C_t:.1f}C at 15 °C)",
                 "preheat to 25 °C+ or accept CC derating; CV tail unchanged"))
    return dict(rows=rows, r0=r0, chiller=ch, fixes=fixes)

# ------------------------------------------------------------------ #
#  v5: full narrative report (in-app tab + HTML export)               #
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
#  v7.1: BEV road-car pack benchmark database (user-supplied)         #
# ------------------------------------------------------------------ #
_BM_COLS = {
 "Model/Pack": "model", "Energy": "E_use", "Unnamed: 2": "E_tot",
 "Nominal Capacity": "cap_Ah", "Power": "P10s_kW", "Pack Mass": "m_pack",
 "Cell Mass": "m_cells", "Pack-Cells": "m_noncell",
 "Cell to Pack Ratio": "ctp_m", "Pack Dimensions": "dx", "Unnamed: 10": "dy",
 "Unnamed: 11": "dz", "Pack Volume": "V_pack", "Cell Total Volume": "V_cells",
 "Cell to Pack Ratio.1": "ctp_v", "Pack-Cells.1": "V_noncell",
 "Gravimetric Energy Density": "whkg", "Volumetric Energy Density": "whl",
 "Gravimetric Power Density": "w10s_kg", "Unnamed: 19": "wcont_kg_dis",
 "Unnamed: 20": "wcont_kg_chg", "Volumetric Power Density": "w10s_l",
 "Unnamed: 22": "wcont_l"}

def load_benchmark(path=None, uploaded=None):
    path = path or str(_APPDIR / "pack_benchmark.xlsx")
    try:
        import io
        srcx = io.BytesIO(uploaded) if uploaded is not None else path
        df = pd.read_excel(srcx, sheet_name="Sheet1").rename(columns=_BM_COLS)
    except Exception:
        return None
    df = df.iloc[1:].reset_index(drop=True)
    for c in df.columns:
        if c != "model":
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["whkg"] = df["whkg"].fillna(df["E_tot"] * 1000 / df["m_pack"])
    df["whl"] = df["whl"].fillna(df["E_tot"] * 1000 / df["V_pack"])
    df["kW_kg_10s"] = df["P10s_kW"] / df["m_pack"] * 1000
    df["C10s"] = df["P10s_kW"] / df["E_tot"]
    return df

def _pct(series, val):
    s = series.dropna()
    return 100.0 * (s < val).mean() if len(s) else float("nan")

def benchmark_db_tab(d, g, masses, Cmax, res):
    up_bytes = st.session_state.get("bm_upload_bytes")
    bm = load_benchmark(uploaded=up_bytes)
    if bm is None:
        st.error("pack_benchmark.xlsx not found next to app.py. Drop the "
                 "file in the app folder, or load it here for this session:")
        upf = st.file_uploader("Load benchmark xlsx", type="xlsx",
                               key="bm_upl")
        if upf is not None:
            st.session_state["bm_upload_bytes"] = upf.getvalue()
            st.rerun()
        return None
    srcnote = ("uploaded this session" if up_bytes is not None
               else "pack_benchmark.xlsx next to app.py")
    st.markdown(f"<span class='chip ok'>DATA LOADED - {len(bm)} BEV "
                f"road-car packs ({srcnote})</span>", unsafe_allow_html=True)

    f1, f2, f3 = st.columns([1.2, 1, 1])
    e_rng = f1.slider("Pack energy filter [kWh]", 0.0,
                      float(np.ceil(bm["E_tot"].max() / 10) * 10),
                      (0.0, float(np.ceil(bm["E_tot"].max() / 10) * 10)),
                      key="bm_e")
    q = f2.text_input("Search model", "", key="bm_q")
    only_p = f3.toggle("Only packs with 10 s power", False, key="bm_p")
    v = bm[(bm["E_tot"] >= e_rng[0]) & (bm["E_tot"] <= e_rng[1])]
    if q:
        v = v[v["model"].str.contains(q, case=False, na=False)]
    if only_p:
        v = v[v["P10s_kW"].notna()]

    def eclass(e):
        return ("City < 50 kWh" if e < 50 else "Mid 50-80 kWh" if e < 80
                else "Large 80-110 kWh" if e < 110 else "Flagship > 110 kWh")
    v = v.copy(); v["grp"] = v["E_tot"].apply(eclass)
    GRP_C = {"City < 50 kWh": "#10B981", "Mid 50-80 kWh": "#0EA5E9",
             "Large 80-110 kWh": "#6366F1", "Flagship > 110 kWh": "#8B5CF6"}
    my_whkg, my_whl = masses["whkg_pack"], masses["whl_pack"]
    m_rib = masses["m_pack"] - masses["m_struct"] \
            + masses["m_struct"] * 0.56 / 0.45
    whkg_rib = masses["E_kwh"] * 1000 / m_rib

    kpi_cards([
        ("Packs shown", f"{len(v)}", f"of {len(bm)} in the file", ""),
        ("Wh/kg percentile", f"{_pct(bm['whkg'], my_whkg):.0f}%",
         f"this design: {my_whkg:.0f} Wh/kg", "brand"),
        ("Wh/L percentile", f"{_pct(bm['whl'], my_whl):.0f}%",
         f"this design: {my_whl:.0f} Wh/L", "brand"),
        ("Cell-to-pack mass", f"{100*masses['m_cells']/masses['m_pack']:.0f}%",
         f"database median {100*bm['ctp_m'].median():.0f}%", ""),
    ])
    VIEWS = ["Density map", "Power", "Rankings", "Explorer",
             "Head-to-head", "Full table"]
    sub = st.segmented_control("Benchmark view", VIEWS, default=VIEWS[0],
                               key="bm_view", label_visibility="collapsed") \
        if hasattr(st, "segmented_control") else \
        st.radio("Benchmark view", VIEWS, horizontal=True, key="bm_view")
    sub = sub or VIEWS[0]

    if sub == "Density map":
        z1, z2 = st.columns([1, 3])
        zones = z1.toggle("Median zones", True, key="bm_zones")
        figA = go.Figure()
        if zones:
            q1x, mx, q3x = (bm["whkg"].quantile(.25), bm["whkg"].median(),
                            bm["whkg"].quantile(.75))
            q1y, my_, q3y = (bm["whl"].quantile(.25), bm["whl"].median(),
                             bm["whl"].quantile(.75))
            figA.add_vrect(x0=q1x, x1=q3x, fillcolor="rgba(99,102,241,0.05)",
                           line_width=0)
            figA.add_vline(x=mx, line_dash="dot", line_color="#94A3B8",
                           annotation_text=f"median {mx:.0f}")
            figA.add_hrect(y0=q1y, y1=q3y, fillcolor="rgba(99,102,241,0.05)",
                           line_width=0)
            figA.add_hline(y=my_, line_dash="dot", line_color="#94A3B8",
                           annotation_text=f"median {my_:.0f}")
        for gname, gg_ in v.groupby("grp"):
            figA.add_trace(go.Scatter(
                x=gg_["whkg"], y=gg_["whl"], mode="markers", name=gname,
                marker=dict(size=7 + 20 * gg_["E_tot"] / bm["E_tot"].max(),
                            color=GRP_C[gname], opacity=0.85,
                            line=dict(color="white", width=1)),
                hovertext=[f"<b>{mn}</b><br>{e:.0f} kWh, {mm:.0f} kg<br>"
                           f"{wk:.0f} Wh/kg, {wl:.0f} Wh/L"
                           for mn, e, mm, wk, wl in zip(
                               gg_["model"], gg_["E_tot"], gg_["m_pack"],
                               gg_["whkg"], gg_["whl"])],
                hoverinfo="text"))
        figA.add_trace(go.Scatter(x=[my_whkg], y=[my_whl],
            mode="markers+text", text=["This design"],
            textposition="top center", name="This design",
            marker=dict(symbol="star", size=22, color="#EF4444",
                        line=dict(color="white", width=2))))
        figA.add_trace(go.Scatter(x=[whkg_rib], y=[my_whl], mode="markers",
            name="FEA-honest enclosure",
            marker=dict(symbol="star-open", size=18, color="#EF4444",
                        line=dict(width=2))))
        figA.update_layout(height=470, xaxis_title="Gravimetric [Wh/kg]",
                           yaxis_title="Volumetric [Wh/L]",
                           title="Energy-density map - bubble = pack kWh; "
                                 "click legend classes on/off",)
        st.plotly_chart(figA, use_container_width=True, key="bm_A",
                        config=PLOTCFG)
        st.session_state["bm_figA"] = figA
        z2.caption("Double-click a legend entry to isolate one class.")

    elif sub == "Power":
        vp = v[v["kW_kg_10s"].notna()]
        figB = go.Figure()
        for gname, gg_ in vp.groupby("grp"):
            figB.add_trace(go.Scatter(x=gg_["whkg"], y=gg_["kW_kg_10s"],
                mode="markers", name=gname,
                marker=dict(size=10, color=GRP_C[gname], opacity=0.85),
                hovertext=gg_["model"], hoverinfo="text+x+y"))
        vc = v[v["wcont_kg_dis"].notna()]
        figB.add_trace(go.Scatter(x=vc["whkg"], y=vc["wcont_kg_dis"],
            mode="markers", name="continuous (database)",
            marker=dict(size=11, color="#0F172A", symbol="diamond"),
            hovertext=vc["model"], hoverinfo="text+x+y"))
        my_cont = Cmax * masses["E_kwh"] * 1000 / masses["m_pack"]
        figB.add_trace(go.Scatter(x=[my_whkg], y=[my_cont],
            mode="markers+text", text=["this (continuous)"],
            textposition="top center", name="This design",
            marker=dict(symbol="star", size=20, color="#EF4444")))
        figB.add_trace(go.Scatter(x=[68.5], y=[150000 / 89],
            mode="markers+text", text=["AMG HPB80 (peak)"],
            textposition="bottom center", name="HPB80",
            marker=dict(symbol="x", size=12, color="#D97706")))
        figB.update_layout(height=460, xaxis_title="Wh/kg",
                           yaxis_title="W/kg",
                           title="Power vs energy density - database 10 s "
                                 "ratings vs this design's continuous",)
        st.plotly_chart(figB, use_container_width=True, key="bm_B",
                        config=PLOTCFG)

    elif sub == "Rankings":
        cC, cD = st.columns(2)
        with cC:
            vr = v.dropna(subset=["whkg"]).sort_values("whkg").tail(18)
            figC = go.Figure(go.Bar(y=vr["model"], x=vr["whkg"],
                                    orientation="h", marker_color="#94A3B8"))
            figC.add_vline(x=my_whkg, line_color="#EF4444", line_width=3,
                           annotation_text=f"this design {my_whkg:.0f}")
            figC.add_vline(x=68.5, line_color="#D97706", line_dash="dot",
                           annotation_text="HPB80")
            figC.update_layout(height=460, xaxis_title="Wh/kg",
                               title="Top of the field vs this design")
            st.plotly_chart(figC, use_container_width=True, key="bm_C")
        with cD:
            vm = bm.dropna(subset=["ctp_m"]).sort_values("ctp_m")
            figD = go.Figure(go.Bar(x=vm["model"], y=100 * vm["ctp_m"],
                                    marker_color="#0EA5E9"))
            figD.add_hline(y=100 * masses["m_cells"] / masses["m_pack"],
                           line_color="#EF4444", line_width=3,
                           annotation_text=f"this design "
                           f"{100*masses['m_cells']/masses['m_pack']:.0f}%")
            figD.update_layout(height=460,
                               yaxis_title="cell mass / pack mass [%]",
                               title="Cell-to-pack mass ratio",
                               xaxis_tickangle=-40)
            st.plotly_chart(figD, use_container_width=True, key="bm_D")

    elif sub == "Explorer":
        NUMCOLS = {"Wh/kg": "whkg", "Wh/L": "whl",
                   "Pack energy [kWh]": "E_tot", "Pack mass [kg]": "m_pack",
                   "10 s power [kW]": "P10s_kW", "10 s C-rate": "C10s",
                   "10 s W/kg": "kW_kg_10s", "Cell-to-pack mass": "ctp_m",
                   "Pack volume [L]": "V_pack", "Capacity [Ah]": "cap_Ah"}
        e1, e2, e3, e4 = st.columns([1.2, 1.2, 0.7, 0.9])
        xk = e1.selectbox("X", list(NUMCOLS), 2, key="bm_x")
        yk = e2.selectbox("Y", list(NUMCOLS), 0, key="bm_y")
        logx = e3.toggle("log X", False, key="bm_lx")
        fitl = e4.toggle("Trend line", True, key="bm_fit")
        ve = v.dropna(subset=[NUMCOLS[xk], NUMCOLS[yk]])
        figE = go.Figure()
        for gname, gg_ in ve.groupby("grp"):
            figE.add_trace(go.Scatter(x=gg_[NUMCOLS[xk]],
                                      y=gg_[NUMCOLS[yk]], mode="markers",
                                      name=gname,
                                      marker=dict(size=9,
                                                  color=GRP_C[gname],
                                                  opacity=0.85),
                                      hovertext=gg_["model"],
                                      hoverinfo="text+x+y"))
        if fitl and len(ve) > 3:
            xx = np.log10(ve[NUMCOLS[xk]]) if logx else ve[NUMCOLS[xk]]
            a1, a0 = np.polyfit(xx, ve[NUMCOLS[yk]], 1)
            xs_ = np.linspace(xx.min(), xx.max(), 40)
            figE.add_trace(go.Scatter(x=10 ** xs_ if logx else xs_,
                                      y=a0 + a1 * xs_, mode="lines",
                                      name="trend",
                                      line=dict(color="#64748B",
                                                dash="dash")))
        mine = {"whkg": my_whkg, "whl": my_whl, "E_tot": masses["E_kwh"],
                "m_pack": masses["m_pack"],
                "ctp_m": masses["m_cells"] / masses["m_pack"],
                "V_pack": masses["V_outer_L"]}
        if NUMCOLS[xk] in mine and NUMCOLS[yk] in mine:
            figE.add_trace(go.Scatter(x=[mine[NUMCOLS[xk]]],
                                      y=[mine[NUMCOLS[yk]]],
                                      mode="markers+text",
                                      text=["This design"],
                                      textposition="top center",
                                      name="This design",
                                      marker=dict(symbol="star", size=18,
                                                  color="#EF4444")))
        figE.update_layout(height=460, xaxis_title=xk, yaxis_title=yk,
                           xaxis_type="log" if logx else "linear",)
        st.plotly_chart(figE, use_container_width=True, key="bm_E",
                        config=PLOTCFG)

    elif sub == "Head-to-head":
        opts = bm.dropna(subset=["whkg", "whl"])["model"].tolist()
        pick = st.multiselect("Pick up to four packs to duel this design",
                              opts,
                              default=[o for o in opts if "Model 3" in o][:1]
                              or opts[:1], max_selections=4, key="bm_duel")
        METRICS = [("Wh/kg", "whkg", my_whkg), ("Wh/L", "whl", my_whl),
                   ("kWh", "E_tot", masses["E_kwh"]),
                   ("10 s W/kg", "kW_kg_10s",
                    Cmax * masses["E_kwh"] * 1000 / masses["m_pack"]),
                   ("Cell-to-pack %", "ctp_m",
                    masses["m_cells"] / masses["m_pack"])]
        cats = [m_[0] for m_ in METRICS]
        def _norm(key, val):
            s = bm[key].dropna()
            rng = max(s.max() - s.min(), 1e-9)
            return max(min((val - s.min()) / rng, 1.05), 0.0)
        figH = go.Figure()
        figH.add_trace(go.Scatterpolar(
            r=[_norm(k, mv) for _, k, mv in METRICS]
              + [_norm("whkg", my_whkg)],
            theta=cats + [cats[0]], name="This design",
            line=dict(color="#EF4444", width=3), fill="toself",
            fillcolor="rgba(239,68,68,0.10)"))
        for pm in pick:
            row = bm[bm["model"] == pm].iloc[0]
            rr = [_norm(k, row[k]) if row[k] == row[k] else 0
                  for _, k, _ in METRICS]
            figH.add_trace(go.Scatterpolar(r=rr + [rr[0]],
                                           theta=cats + [cats[0]],
                                           name=pm, fill="toself",
                                           opacity=0.75))
        figH.update_layout(height=470, polar=dict(radialaxis=dict(
            visible=True, range=[0, 1.05], showticklabels=False)),
            legend=dict(orientation="h", y=-0.08),
            title="Normalised to the database range - outer edge = best in "
                  "file; this design's power axis is continuous, theirs 10 s")
        st.plotly_chart(figH, use_container_width=True, key="bm_H")

    else:
        show = v[["model", "E_use", "E_tot", "m_pack", "whkg", "whl",
                  "P10s_kW", "C10s", "ctp_m", "V_pack"]].round(1)
        st.dataframe(show, hide_index=True, use_container_width=True,
                     height=460)
        st.download_button("Download filtered set (.csv)",
                           show.to_csv(index=False),
                           "pack_benchmark_view.csv", key="bm_dl")

    near = bm.dropna(subset=["whkg", "whl"]).copy()
    near["dist"] = np.hypot((near["whkg"] - my_whkg) / bm["whkg"].std(),
                            (near["whl"] - my_whl) / bm["whl"].std())
    nn = near.nsmallest(3, "dist")["model"].tolist()
    st.markdown(f"**Where this design sits:** {my_whkg:.0f} Wh/kg is the "
                f"**{_pct(bm['whkg'], my_whkg):.0f}th percentile** of these "
                f"packs and {my_whl:.0f} Wh/L the "
                f"**{_pct(bm['whl'], my_whl):.0f}th**; nearest neighbours: "
                f"{', '.join(nn)}. The gap is the immersion tax - oil "
                f"({masses['m_oil']:.0f} kg) and the pressure-rated "
                f"enclosure ({masses['m_struct']:.0f} kg). What this file "
                "cannot show is the thermal case: abuse tolerance, "
                "uniformity, the flywheel.")
    return st.session_state.get("bm_figA")


# ------------------------------------------------------------------ #
#  v8.6: full cooling-architecture trade study with costs             #
# ------------------------------------------------------------------ #
COST_DEFAULTS = dict(cell_usd_kwh=79.0, usd_gbp=0.79, eu_premium=1.56,
                     oil_gbp_L=9.0, alu_gbp_kg=3.2, cu_gbp_kg=8.5,
                     fab_factor=1.8, pump_gbp=28.0, chiller_gbp_kWel=140.0,
                     coldplate_gbp=280.0, pads_gbp=55.0, blower_gbp=45.0,
                     glycol_gbp=15.0, manifold_gbp=40.0)

_WATER = dict(rho=1000.0, cp=4180.0, k=0.60, mu=8.9e-4)

def _cell_axials(d, Q_cell, f_dry):
    A_cs = math.pi * (d["d_cell"] / 2) ** 2
    kz = 30.0
    dT_top_dry = Q_cell * f_dry ** 2 * d["h_cell"] / (2 * kz * A_cs)
    dT_ax_mean = Q_cell * d["h_cell"] / (3 * kz * A_cs)
    dT_ax_top = Q_cell * d["h_cell"] / (2 * kz * A_cs)
    return dT_top_dry, dT_ax_mean, dT_ax_top

def arch_solve(name, d, g, fl, masses, C, cost):
    """Return one architecture's honest numbers at continuous C."""
    Q_cell_25 = (C * d["cap_Ah"]) ** 2 * d["r_dc"] / 1000
    mat = cost["cu_gbp_kg"] if d["tube_mat"] == "Copper" else cost["alu_gbp_kg"]
    base_therm_cost = ((masses["m_tubes"] + masses["m_fins"]) * mat
                       * cost["fab_factor"] + cost["pump_gbp"]
                       + cost["manifold_gbp"])
    r = dict(name=name, note="")
    if name.startswith(("Full immersion", "Partial")):
        dd = dict(d)
        if "serpentine" in name:
            dd.update(plate_on=True, u_oil=0.05)
        elif "stirred" in name:
            dd.update(plate_on=False, u_oil=0.05)
        else:
            dd.update(plate_on=False, u_oil=0.0)
        s = solve_steady(dd, g, fl, 1.0, d["T_amb"], C_rate=C)
        f_w = 0.6 if name.startswith("Partial") else 1.0
        extra = s["Q_eff"] * s["R_b"] * (1 / f_w - 1)
        Q_cell = s["Q_eff"] / g["N"]
        dT_dry, _, _ = _cell_axials(d, Q_cell, 1 - f_w)
        r["T_can"] = s["T_b"] + extra
        r["T_hot"] = max(s["T_core"] + extra, r["T_can"] + dT_dry)
        r["spread"] = s["spread"] * (1.4 if f_w < 1 else 1.0)
        P_c = (serpentine_pump(dd, g, fl, 0.05)["P"] if "serpentine" in name
               else stirrer_power(dd, g, fl, dd["u_oil"]))
        r["parasitic"] = water_pump_power(d, g, _WATER)["P"] + P_c
        r["Q"] = s["Q_eff"]; r["Q_w"] = s["Q_w"]
        m_oil = masses["m_oil"] * f_w
        m_pl = plate_fin_area(dd, g, s["h_tube"])[2] if "serpentine" in name else 0
        r["mass"] = masses["m_pack"] - masses["m_oil"] + m_oil + m_pl
        r["cost_th"] = (base_therm_cost + m_oil / fl["rho"] * 1000
                        * cost["oil_gbp_L"]
                        + m_pl * cost["alu_gbp_kg"] * cost["fab_factor"])
        r["note"] = ("distributed sink, best uniformity" if "serpentine"
                     in name else
                     "top of cells runs dry - hotspot above the oil line"
                     if f_w < 1 else "the baseline architecture")
    elif name.startswith("Indirect"):
        Q = Q_cell_25 * g["N"] * math.exp(-d["k_dcir"] * 10)  # ~35 degC est
        Q_cell = Q / g["N"]
        _, ax_m, ax_t = _cell_axials(d, Q_cell, 0.0)
        A_cs = math.pi * (d["d_cell"] / 2) ** 2
        R_pad = 0.0015 / (3.0 * A_cs)
        R_pl = 0.30
        rise = Q / ((d["flow_lpm"] / 60) * 1000 * 4180)
        r["T_can"] = d["T_water_in"] + 0.5 * rise \
            + Q_cell * (R_pad + R_pl) + ax_m
        r["T_hot"] = d["T_water_in"] + rise + Q_cell * (R_pad + R_pl) \
            + ax_t + res_core_extra(d, Q_cell)
        r["spread"] = ax_t + rise
        r["parasitic"] = water_pump_power(d, g, _WATER)["P"] + 2
        r["Q"] = Q; r["Q_w"] = Q
        m_plate = g["Lx"] * g["Ly"] * 0.003 * 2700 + 4.0
        enc_light = masses["m_struct"] * 0.45
        r["mass"] = (masses["m_pack"] - masses["m_oil"]
                     - masses["m_struct"] + enc_light
                     - masses["m_tubes"] - masses["m_fins"] + m_plate + 1.5)
        r["cost_th"] = (cost["coldplate_gbp"] + cost["pads_gbp"]
                        + cost["pump_gbp"] + cost["glycol_gbp"]
                        + cost["manifold_gbp"])
        r["note"] = ("bottom cold plate + pads; heat must conduct down the "
                     "jellyroll (k_z) - axial gradient dominates")
    else:  # Forced air
        h_a, util = 45.0, 0.7
        Q = Q_cell_25 * g["N"] * math.exp(-d["k_dcir"] * 14)
        dT_air = 8.0
        R = 1 / (h_a * g["A_cells"] * util)
        r["T_can"] = d["T_amb"] + 0.5 * dT_air + Q * R
        r["T_hot"] = d["T_amb"] + dT_air + Q * R \
            + res_core_extra(d, Q / g["N"])
        r["spread"] = dT_air + 1.5
        mdot = Q / (1005 * dT_air)
        r["parasitic"] = mdot / 1.2 * 250 / 0.4
        r["Q"] = Q; r["Q_w"] = 0.0
        r["mass"] = (masses["m_pack"] - masses["m_oil"]
                     - masses["m_struct"] * 0.5 - masses["m_tubes"]
                     - masses["m_fins"] + 2.5)
        r["cost_th"] = cost["blower_gbp"] * 2
        r["note"] = ("no chiller possible - tied to ambient; uniformity "
                     "fails the 5 °C criterion")
    ch = chiller_model(max(r["Q_w"], 1e-3), d["T_water_in"], d["T_amb"]) \
        if r["Q_w"] > 0 else dict(P_el=0.0, COP=float("nan"))
    r["chiller_el"] = ch["P_el"]
    r["cost_th"] += cost["chiller_gbp_kWel"] * ch["P_el"] / 1000
    return r

def res_core_extra(d, Q_cell):
    return Q_cell / (4 * math.pi * d["k_rad"] * d["h_cell"])

def arch_maxC(name, d, g, fl, masses, cost, T_lim):
    lo, hi = 0.3, 6.0
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        r = arch_solve(name, d, g, fl, masses, mid, cost)
        if r["T_hot"] > T_lim: hi = mid
        else: lo = mid
    return lo

ARCH_LIST = ["Full immersion - static (this project)",
             "Full immersion - serpentine plates",
             "Full immersion - stirred",
             "Partial immersion (60% fill)",
             "Indirect cold plate (dry)",
             "Forced air (dry)"]

def full_arch_study(d, g, fl, masses, C, cost):
    rows = []
    for nm in ARCH_LIST:
        r = arch_solve(nm, d, g, fl, masses, C, cost)
        r["maxC"] = arch_maxC(nm, d, g, fl, masses, cost, d["T_limit"])
        rows.append(r)
    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False)
def _arch_study_cached(d_json, C, cost_json):
    dd = json.loads(d_json); cc = json.loads(cost_json)
    gg = build_geometry(dd)
    cool = _read_coolants()
    fl_ = fluid_dict(cool[cool["name"] == dd["coolant"]].iloc[0])
    rr = solve_steady(dd, gg, fl_, 1.0, dd["T_amb"], C_rate=C)
    mm = build_masses(dd, gg, fl_, rr["fin"])
    return full_arch_study(dd, gg, fl_, mm, C, cc)


@st.cache_data(show_spinner=False)
def _sweep_cached(d_json, kx, x0, x1, ky, y0, y1, C):
    dd0 = json.loads(d_json)
    cool = _read_coolants()
    fl_ = fluid_dict(cool[cool["name"] == dd0["coolant"]].iloc[0])
    xs = np.linspace(x0, x1, 7); ys = np.linspace(y0, y1, 7)
    Z = np.zeros((7, 7))
    for i, yv in enumerate(ys):
        for j, xv in enumerate(xs):
            dd = dict(dd0)
            dd[kx] = int(round(xv)) if kx == "n_tubes" else float(xv)
            dd[ky] = int(round(yv)) if ky == "n_tubes" else float(yv)
            gg = build_geometry(dd)
            Z[i, j] = solve_steady(dd, gg, fl_, 1.0, dd0["T_amb"],
                                   C_rate=C)["T_b"]
    return xs, ys, Z

def arch_chart(df, T_lim):
    short = [n.replace("Full immersion - ", "FI ").replace(" (this project)",
             "").replace(" (dry)", "").replace(" (60% fill)", " 60%")
             for n in df["name"]]
    fig = go.Figure()
    fig.add_trace(go.Bar(name="Hotspot [°C]", x=short, y=df["T_hot"],
                         marker_color="#F0655F", offsetgroup=1,
                         text=df["T_hot"].round(1), textposition="outside"))
    fig.add_trace(go.Bar(name="Spread [°C]", x=short, y=df["spread"],
                         marker_color="#7C88F8", offsetgroup=2, yaxis="y2",
                         text=df["spread"].round(1), textposition="outside"))
    fig.add_hline(y=T_lim, line_dash="dash", line_color="#B91C1C",
                  annotation_text=f"limit {T_lim:.0f}")
    fig.add_hline(y=5, line_dash="dot", line_color="#4A54D8", yref="y2",
                  annotation_text="5 °C uniformity")
    fig.update_layout(barmode="group", height=400,
                      yaxis=dict(title="Hotspot [°C]"),
                      yaxis2=dict(title="Spread [°C]", overlaying="y",
                                  side="right", showgrid=False,
                                  range=[0, max(df["spread"].max(), 6) * 1.4]),
                      title="Six ways to cool the same cells - hotspot vs "
                            "uniformity at this duty")
    return fig

def report_sections(d, g, fl, res, masses, tr, spec, Cmax, C_steady, Q_duty,
                    Q_bus, P_pump, P_stir, chil, figs, arch_df=None,
                    cost=None):
    ok = (res["T_core"] if d["limit_core"] else res["T_b"]) <= d["T_limit"]
    stations, totdT = station_list(d, g, fl, res, masses, Q_duty, Q_bus,
                                   P_pump, P_stir, chil)
    st_md = ""
    for name, dT, nums, imp in stations:
        share = f" - **{dT:.1f} °C** ({100*dT/max(totdT,1e-9):.0f}% of the ladder)" \
                if dT > 0 else ""
        st_md += f"\n**{name}**{share}\n\n{nums}\n\n*Improve:* {imp}\n"
    S = []
    S.append(("Executive summary", f"""
{masses['E_kwh']:.1f} kWh / {d['Ns']*d['v_nom']:.0f} V pack of {g['N']} x
{d['fmt']} cells in static {fl['name']}, cooled by {d['n_tubes']} internal
water tubes ({d['tube_plane'].lower()}). At the duty's {C_steady:.2f}C RMS the
cells sit at **{res['T_b']:.1f} °C can / {res['T_core']:.1f} core** against a
{d['T_limit']:.0f} °C limit: **{'WITHIN LIMIT' if ok else 'OVER LIMIT'}**.
Max continuous **{Cmax:.2f}C**. Pack **{masses['m_pack']:.0f} kg**
({masses['whkg_pack']:.0f} Wh/kg, {masses['whl_pack']:.0f} Wh/L). Chiller duty
{res['Q_w']/1000:.2f} kW -> ~{chil['P_el']/1000:.2f} kW electrical at COP
{chil['COP']:.1f}. Parasitics {P_pump+P_stir:.0f} W.""", ["sankey"]))
    S.append(("The design", f"""
Box {g['Lx']*1000:.0f} x {g['Ly']*1000:.0f} x {g['Lz']*1000:.0f} mm,
{g['n_cols']} x {g['n_rows']} {d['arrangement'].lower()} grid at
{d['pitch']*1000:.1f} mm pitch (gap {g['gap_mm']:.1f} mm), oil to
{g['fill_h']*1000:.0f} mm with {d['gas_gap']*1000:.0f} mm headspace.
Enclosure sized for {d['p_des_bar']:.1f} bar g: {masses['enc']['t_mm']:.1f} mm
effective wall, {masses['m_struct']:.0f} kg. Mass ledger: cells
{masses['m_cells']:.0f}, oil {masses['m_oil']:.0f}, enclosure
{masses['m_struct']:.0f}, holders {masses['m_holders']:.1f}, tubes+fins
{masses['m_tubes']+masses['m_fins']:.1f}, busbars {masses['m_bus']:.1f} kg.
*Improve:* the enclosure is the second-heaviest non-cell item; FEA supports a
stiffened knock-down of ~0.50-0.56 (ribbed), or drop the burst set pressure
and re-check runaway venting.""", ["pack3d"]))
    S.append(("The heat journey, station by station", st_md, ["ladder", "resist"]))
    S.append(("Duty and transient response", f"""
Duty: **{d['duty']}**{(' - ' + d['cycle']) if d['duty']=='Drive cycle' else ''},
{spec['t'][-1]:.0f} s simulated, RMS {C_steady:.2f}C, peak
{np.abs(tr['C']).max():.2f}C. SoC {tr['soc'][0]*100:.0f} ->
{tr['soc'][-1]*100:.0f}%. The oil flywheel
({(masses['C_oil']+masses['C_batt'])/1e3:.0f} kJ/K) filters peaks; size steady
hardware for the RMS, not the spike.
*Improve:* if transients ever graze the limit, buffer harder (more oil) before
buying chiller.""", ["transient"]))
    S.append(("Scaling with C-rate", f"""
Heat scales ~C^2, softened by DCIR(T). This hardware holds
{d['T_limit']:.0f} °C up to **{Cmax:.2f}C** continuous. The chart shows
heat, chiller electricity and cell temperature as C rises; the Improve tab's
scale-to-C tool lists the cheapest change per subsystem for any target.
*Improve at 4C-class duties:* stirring becomes mandatory (static is a 1-2C
architecture), water must be turbulent, busbars resize with I, and the core
term grows to ~{4*4/(C_steady**2+1e-9)*res['dT_core']:.0f} °C-class - check the
core limit, and note charging additionally needs the plating map satisfied
(warm cells).""", ["csweep", "setpoint"]))
    if arch_df is not None:
        cost = cost or COST_DEFAULTS
        Q_45_4C = (Q_duty * (4.0 / max(C_steady, 1e-6)) ** 2
                   * math.exp(-d["k_dcir"] * (45.0 - res["T_b"])))
        m_therm = (masses["m_oil"] + masses["m_tubes"] + masses["m_fins"]
                   + masses.get("m_plates", 0.0)
                   + masses["m_struct"] * 0.55)
        cell_cost_usd = masses["E_kwh"] * 1000 * cost["cell_usd_kwh"] / 1000
        tbl = "| Architecture | Hotspot [°C] | Spread [°C] | Max cont. C | " \
              "Pack [kg] | Parasitic [W] | Chiller [W el] | Thermal £ |\n" \
              "|---|---|---|---|---|---|---|---|\n"
        for _, rr in arch_df.iterrows():
            tbl += (f"| {rr['name']} | {rr['T_hot']:.1f} | "
                    f"{rr['spread']:.1f} | {rr['maxC']:.2f} | "
                    f"{rr['mass']:.0f} | {rr['parasitic']:.0f} | "
                    f"{rr['chiller_el']:.0f} | {rr['cost_th']:.0f} |\n")
        i_best = arch_df["spread"].idxmin(); i_light = arch_df["mass"].idxmin()
        S.append(("Cooling requirement, impacts, and the architecture trade",
                  f"""
**What the cooling system must do.** Remove **{Q_duty/1000:.2f} kW
continuously** at the duty's {C_steady:.2f}C RMS while holding the cells at
or below {d['T_limit']:.0f} °C, with best-to-worst cell spread under 5 °C.
Water side: {d['flow_lpm']:.0f} L/min at {d['T_water_in']:.0f} °C
({res['water_regime']}, Re = {res['Re_water']:.0f}; target Re > 3000).
Chiller: **{(res['Q_w']+P_pump)/1000:.2f} kW duty ->
{chil['P_el']/1000:.2f} kW electrical** at COP {chil['COP']:.1f}. For a
4C-class duty held at the 45 °C limit the load rises to roughly
**{Q_45_4C/1000:.1f} kW**, and no single lever suffices - stirring or the
serpentine plates plus turbulent water become mandatory.

**Weight impact.** The thermal system costs
**{m_therm:.0f} kg ({100*m_therm/masses['m_pack']:.0f}% of the pack)**: oil
{masses['m_oil']:.0f} kg, tubes and fins
{masses['m_tubes']+masses['m_fins']:.1f} kg, plates
{masses.get('m_plates',0):.1f} kg, and the burst-rating premium on the
enclosure ~{masses['m_struct']*0.55:.0f} kg. That is the direct cause of
the 4th-percentile Wh/kg position against the BEV database.

**Cost impact.** Thermal-system costs per architecture are tabled below
from editable assumptions (dielectric {cost['oil_gbp_L']:.0f} £/L, chiller
{cost['chiller_gbp_kWel']:.0f} £/kW el, fabrication x{cost['fab_factor']:.1f}).
For scale: the cells themselves are ~**{cell_cost_usd:,.0f} USD** at the BNEF
December-2025 average of 79 USD/kWh (cell level; BEV packs averaged 99 USD/kWh,
Europe typically +56%). Against that, the flooded coolant volume is the
largest thermal line item for the immersion rows, so the thermal system is a
meaningful fraction of cell cost, not a rounding error - the full bill of
materials, item by item, is in the System tab.

**Performance and the honest comparison.** The same cells and duty under
six approaches - immersion rows from the full validated solver, dry rows
from closed-form conduction chains (bottom-cooled axial gradient, pad and
plate resistances; treat as +/-20-30%):

{tbl}
{arch_df.loc[i_best,'name']} gives the best uniformity
({arch_df['spread'].min():.1f} °C); {arch_df.loc[i_light,'name']} is
lightest ({arch_df['mass'].min():.0f} kg). The indirect cold plate wins
mass, cost and mean temperature but its axial-conduction spread fails the
5 °C criterion; forced air is a sub-2C architecture; partial immersion
overheats the dry cell tops. Full immersion buys uniformity, abuse
tolerance and the thermal flywheel with oil mass and the pressure-rated
box - and the serpentine plates are its best variant.
*Improve:* fit the rig data to the two oil films (cal factor), then revisit
this table; the dry-architecture constants (pad k, plate R) deserve a
vendor quote before any down-select.""", ["archs"]))

    S.append(("Reading the heat-flow map in depth",
              sankey_deep_dive(d, g, fl, res, masses, Q_duty, Q_bus, P_pump,
                               P_stir, chil, C_steady), ["sankey"]))

    S.append(("Uniformity and safety", f"""
Best-to-worst spread estimate **{res['spread']:.1f} °C** (criterion 5 °C); the
section model shows the water-rise part is ~3x conservative and that
**counterflow plumbing removes it almost entirely**. Runaway screening: one
cell's {d['E_tr']:.0f} kJ heats its local zone modestly and the whole pack by
under a kelvin, but venting pressurises the headspace - burst disc at
{d['p_des_bar']:.1f} bar g, oil-above-water pressure rule, expansion bellows
sized for ~{fl['beta']*masses['V_oil_L']*(d['T_service_max']):.1f} L.""", []))
    S.append(("Validation and honesty", f"""
Wang et al. 2023 rebuild: 33.0 vs 32.3 °C measured, resistances within
30-50%, zero tuning. FD/FE studies: lid coefficients match Roark (+0.0%),
cell core matches the exact solution (-0.1%) and shows the app conservative
by 20-34% on real geometry; the pack-section model closes energy to 0.00%.
Oil-side h honest to +/-30% until the calibration factor
(currently {d['cal_h']:.2f}) is fitted to rig data. Constants open to
challenge are listed in the README.""", []))
    return S

def render_report_tab(secs, figs, meta):
    words = sum(len(md.split()) for _, md, _ in secs)
    st.markdown(
        f"<div class='hero'><div class='hero-top'>"
        f"<h1>Design report</h1>"
        f"<span class='sub'>{meta['spec']}</span>"
        f"<span class='chip'>{meta['version']} - {meta['date']}</span>"
        f"<span class='chip {'ok' if meta['ok'] else 'bad'}'>"
        f"{'WITHIN LIMIT' if meta['ok'] else 'OVER LIMIT'} - "
        f"{meta['T_gov']:.1f} / {meta['T_limit']:.0f} °C</span>"
        f"</div><div class='hero-stats'>"
        f"<div class='hstat'>sections<b>{len(secs)}</b></div>"
        f"<div class='hstat'>reading time<b>~{max(words//200,1)} min</b></div>"
        f"<div class='hstat'>energy<b>{meta['kwh']:.1f} kWh</b></div>"
        f"<div class='hstat'>mass<b>{meta['mass']:.0f} kg</b></div>"
        f"<div class='hstat'>max continuous<b>{meta['Cmax']:.2f} C</b></div>"
        f"<div class='hstat'>chiller<b>{meta['chil_el']/1000:.2f} kW el</b>"
        f"</div></div></div>", unsafe_allow_html=True)
    names = [f"{i+1}. {t}" for i, (t, _, _) in enumerate(secs)]
    opts = ["Full report"] + names
    pick = (st.pills("Contents", opts, default="Full report",
                     key="rep_nav", label_visibility="collapsed")
            if hasattr(st, "pills") else
            st.radio("Contents", opts, horizontal=True, key="rep_nav"))
    pick = pick or "Full report"
    show = (list(range(len(secs))) if pick == "Full report"
            else [opts.index(pick) - 1])
    for i in show:
        title, md, fkeys = secs[i]
        with st.container(border=True):
            st.markdown(f"### {i+1}. {title}")
            st.markdown(md)
            for k in fkeys:
                if k in figs:
                    st.plotly_chart(figs[k], use_container_width=True,
                                    key=f"rep_{i}_{k}")
    if pick != "Full report":
        i = show[0]
        cprev, cnext, _ = st.columns([1, 1, 5])
        if i > 0 and cprev.button(f"<- {i}. {secs[i-1][0][:22]}",
                                  key="rep_prev"):
            st.session_state.setdefault("_pending", {})["rep_nav"] = opts[i]
            st.rerun()
        if i < len(secs) - 1 and cnext.button(
                f"{i+2}. {secs[i+1][0][:22]} ->", key="rep_next"):
            st.session_state.setdefault("_pending", {})["rep_nav"] = \
                opts[i + 2]
            st.rerun()

def export_report_html(secs, figs, meta) -> str:
    try:
        import markdown as _md
        conv = lambda t: _md.markdown(t, extensions=["tables"])
    except Exception:
        conv = lambda t: "<p>" + t.replace("\n\n", "</p><p>") + "</p>"
    toc, body = "", ""
    for i, (title, md, fkeys) in enumerate(secs):
        sid = f"sec{i+1}"
        toc += f"<a href='#{sid}' data-sec='{sid}'>{i+1}. {title}</a>"
        html_md = conv(md)
        html_md = html_md.replace("<p><em>Improve:</em>",
                                  "<p class='callout'><em>Improve:</em>")
        body += f"<section id='{sid}'><h2>{i+1}. {title}</h2>{html_md}"
        for k in fkeys:
            if k in figs:
                body += ("<div class='fig'>"
                         + figs[k].to_html(full_html=False,
                                           include_plotlyjs=False)
                         + "</div>")
        body += "</section>"
    chip = ("<span class='chip ok'>WITHIN LIMIT</span>" if meta["ok"]
            else "<span class='chip bad'>OVER LIMIT</span>")
    stats = "".join(
        f"<div class='stat'><span>{k}</span><b>{v}</b></div>" for k, v in [
            ("energy", f"{meta['kwh']:.1f} kWh"),
            ("mass", f"{meta['mass']:.0f} kg"),
            ("governing T", f"{meta['T_gov']:.1f} / {meta['T_limit']:.0f} °C"),
            ("max continuous", f"{meta['Cmax']:.2f} C"),
            ("chiller", f"{meta['chil_el']/1000:.2f} kW el"),
        ])
    css = ("body{font-family:Inter,-apple-system,'Segoe UI',sans-serif;"
           "color:#1F2937;margin:0;background:#fff;line-height:1.55}"
           "*{box-sizing:border-box}"
           ".wrap{display:grid;grid-template-columns:230px 1fr;gap:34px;"
           "max-width:1180px;margin:0 auto;padding:26px}"
           "nav{position:sticky;top:20px;align-self:start;"
           "border:1px solid #E7EAF0;border-radius:14px;padding:12px;"
           "font-size:.85rem}"
           "nav a{display:block;color:#64748B;text-decoration:none;"
           "padding:5px 8px;border-radius:8px;"
           "border-left:3px solid transparent}"
           "nav a:hover{background:#F1F5F9}"
           "nav a.on{color:#6E77F0;border-left-color:#6E77F0;"
           "background:#EEF2FF;font-weight:600}"
           "header.cover{grid-column:1/-1;background:"
           "linear-gradient(120deg,#EEF2FF,#E6F7FD);"
           "border:1px solid #E3E8F4;border-radius:18px;padding:22px 26px}"
           "header.cover h1{margin:0;font-size:1.6rem;"
           "letter-spacing:-.02em}"
           "header .sub{color:#64748B;margin-top:4px}"
           ".chips{margin-top:10px}"
           ".chip{display:inline-block;padding:3px 12px;"
           "border-radius:999px;font-size:.76rem;font-weight:600;"
           "margin-right:8px;background:#fff;border:1px solid #E3E8F4;"
           "color:#4A54D8}"
           ".chip.ok{background:#DCFCE7;color:#15803D;"
           "border-color:#86EFAC}"
           ".chip.bad{background:#FEE2E2;color:#B91C1C;"
           "border-color:#FCA5A5}"
           ".stats{display:flex;gap:10px;flex-wrap:wrap;margin-top:14px}"
           ".stat{background:rgba(255,255,255,.8);"
           "border:1px solid #E3E8F4;border-radius:12px;padding:8px 14px;"
           "font-size:.72rem;color:#64748B}"
           ".stat b{display:block;font-size:1rem;color:#1F2937}"
           "section{border:1px solid #E7EAF0;border-radius:16px;"
           "padding:6px 22px;margin:0 0 18px 0;"
           "box-shadow:0 1px 2px rgba(16,24,40,.04)}"
           "h2{color:#1F2937;letter-spacing:-.02em}"
           ".callout{background:#EEF2FF;border-left:4px solid #6E77F0;"
           "border-radius:8px;padding:10px 14px}"
           "table{border-collapse:collapse;width:100%;font-size:.88rem;"
           "margin:12px 0}"
           "th,td{border-bottom:1px solid #E7EAF0;padding:7px 9px;"
           "text-align:left}"
           "th{background:#F8FAFC;font-weight:600}"
           "tr:hover td{background:#F8FAFC}"
           ".fig{margin:10px 0}"
           "#top{position:fixed;right:22px;bottom:22px;background:#6E77F0;"
           "color:#fff;border:none;border-radius:999px;width:42px;"
           "height:42px;font-size:20px;cursor:pointer;opacity:.85;"
           "display:none}"
           "footer{grid-column:1/-1;color:#64748B;font-size:.8rem;"
           "border-top:1px solid #E7EAF0;padding-top:14px;margin-top:6px}"
           "@media print{nav,#top{display:none!important}"
           ".wrap{display:block}section{break-inside:avoid;border:none;"
           "box-shadow:none;padding:0}header.cover{border:none}}"
           "@media (max-width:900px){.wrap{display:block}"
           "nav{position:static;margin-bottom:16px}}")
    js = ("const links=[...document.querySelectorAll('nav a')];"
          "const secs=[...document.querySelectorAll('main section')];"
          "const io=new IntersectionObserver(es=>{es.forEach(e=>{"
          "if(e.isIntersecting){links.forEach(l=>l.classList.toggle('on',"
          "l.dataset.sec===e.target.id));}});},"
          "{rootMargin:'-20% 0px -70% 0px'});"
          "secs.forEach(s=>io.observe(s));"
          "links.forEach(l=>l.addEventListener('click',ev=>{"
          "ev.preventDefault();document.getElementById(l.dataset.sec)"
          ".scrollIntoView({behavior:'smooth'});}));"
          "const topb=document.getElementById('top');"
          "addEventListener('scroll',()=>{topb.style.display="
          "scrollY>600?'block':'none';});")
    return ("<!DOCTYPE html><html><head><meta charset='utf-8'>"
            "<meta name='viewport' "
            "content='width=device-width,initial-scale=1'>"
            "<script src='https://cdn.plot.ly/plotly-2.32.0.min.js'>"
            "</script>"
            f"<title>Immersion Pack Lab report {meta['version']}</title>"
            f"<style>{css}</style></head><body><div class='wrap'>"
            "<header class='cover'>"
            "<h1>Immersion Pack Lab - design report</h1>"
            f"<div class='sub'>{meta['spec']}</div>"
            f"<div class='chips'><span class='chip'>{meta['version']}"
            f"</span><span class='chip'>{meta['date']}</span>{chip}</div>"
            f"<div class='stats'>{stats}</div></header>"
            "<nav id='toc'><b style='display:block;margin:2px 8px 8px'>"
            "Contents</b>" + toc + "</nav>"
            f"<main>{body}</main>"
            "<button id='top' "
            "onclick='scrollTo({top:0,behavior:\"smooth\"})'>^</button>"
            f"<footer>Generated by Immersion Pack Lab {meta['version']} "
            f"on {meta['date']}. Physics anchored to Wang et al. 2023 "
            "(J. Energy Storage 62, 106821) and the pack_fea_v1 studies; "
            "benchmark data: pack_benchmark.xlsx (58 BEV packs); cell "
            "costs: BNEF December 2025 survey. Interactive figures "
            "require the plotly CDN; everything else is self-contained."
            "</footer></div>"
            f"<script>{js}</script></body></html>")


if __name__ == "__main__":
    if os.environ.get("SMOKE"):
        smoke()
    else:
        main()
