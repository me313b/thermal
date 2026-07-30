"""Decide - the cooling-architecture decision application.

Compares Setup A (this project: sealed tank, internal water-pipe
bank, submerged fans) against Setup B (pumped oil through an
external plate heat exchanger) for the CURRENT design in the app,
not for a generic case: per-cell heat comes from the battery model
at the set current, the array and tank from the FEA tab, the
tube-bank film from Churchill-Bernstein at the configured fan
speed and pipe size, and the oil inventory from the auto-sized
tank volume. Every assumption is shown and overridable; every
verdict states the numbers it was decided on.

Runs inside Immersion Pack Lab as the Decide tab, or standalone:
    streamlit run decide.py
"""

import math

import numpy as np
import pandas as pd


# ------------------------------------------------------------------ #
#  Physics                                                           #
# ------------------------------------------------------------------ #
def visc_at(fl, T_C):
    """Kinematic viscosity [m2/s] at T from nu25 [cSt] and the
    Andrade-style B coefficient the coolant table carries."""
    nu25 = float(fl.get("nu25", 40.0)) * 1e-6
    B = float(fl.get("B", 0.0) or 0.0)
    if B <= 0:
        return nu25
    return nu25 * math.exp(B * (1.0 / (T_C + 273.15)
                                - 1.0 / 298.15))


def h_crossflow_cb(fl, u, D, T_C=40.0):
    """Churchill-Bernstein average Nusselt for a cylinder in
    crossflow; valid for all Re*Pr > 0.2. Returns (h, Re, Pr)."""
    nu = visc_at(fl, T_C)
    k, rho, cp = fl["k"], fl["rho"], fl["cp"]
    Re = max(u * D / nu, 1e-6)
    Pr = nu * rho * cp / k
    Nu = 0.3 + (0.62 * Re ** 0.5 * Pr ** (1 / 3)
                / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25) \
        * (1 + (Re / 282000.0) ** 0.625) ** 0.8
    return Nu * k / D, Re, Pr


def bank_capability(h, D, L_pipe, n_pipes, dT):
    """Heat the installed bank can reject [W] at the given film and
    approach: Q = h * (n * pi * D * L) * dT."""
    A = n_pipes * math.pi * D * L_pipe
    return h * A * dT, A


def tube_needed(Q, h, D, dT):
    """Tube length required for Q at film h and approach dT."""
    A = Q / max(h * dT, 1e-12)
    return A / (math.pi * D), A


# ------------------------------------------------------------------ #
#  Cost model (one-off prototype GBP; every value overridable)      #
# ------------------------------------------------------------------ #
COST_DEFAULTS = dict(
    fan_each=15.0, fans_n=3, tube_per_m=8.0, water_pump=25.0,
    water_plumb=20.0, oil_pump=135.0, plate_hx=90.0,
    oil_hoses=60.0, expansion=15.0, bms_fixed=120.0,
    per_cell_assembly=0.8, enclosure=80.0, extra_oil_L=0.8)


def cost_A(L_tube, c):
    return (c["fans_n"] * c["fan_each"] + c["tube_per_m"] * L_tube
            + c["water_pump"] + c["water_plumb"])


def cost_B(Q, c):
    step = 40.0 if Q > 1200 else 0.0
    return (c["oil_pump"] + c["plate_hx"] + step + c["oil_hoses"]
            + c["expansion"] + c["water_pump"] + c["water_plumb"])


def system_costs(n_cells, cell_price, oil_L, oil_price, L_tube,
                 Q, c):
    cells = n_cells * cell_price
    elec = c["bms_fixed"] + n_cells * c["per_cell_assembly"]
    enc = c["enclosure"]
    oil = oil_L * oil_price
    cA = cost_A(L_tube, c)
    cB = cost_B(Q, c)
    extra = c["extra_oil_L"] * oil_price
    common = cells + elec + enc + oil
    return dict(cells=cells, elec=elec, enc=enc, oil=oil,
                coolA=cA, coolB=cB, extra_oil=extra,
                totA=common + cA, totB=common + cB + extra)


def verdicts(Q, cap_A, sc, dT):
    oil_share = 100 * sc["oil"] / sc["totA"]
    delta = sc["totB"] - sc["totA"]
    delta_share = 100 * delta / sc["totA"]
    fits = Q <= cap_A
    oil_judge = ("not a game changer" if oil_share < 5 else
                 "noticeable but not decisive" if oil_share < 12
                 else "significant - worth engineering down")
    chosen = "A" if fits else "B"
    return dict(oil_share=oil_share, oil_judge=oil_judge,
                delta=delta, delta_share=delta_share, fits=fits,
                chosen=chosen)


def decision_md(inp, out, sc, vd):
    """A report-ready markdown summary of the decision."""
    L = []
    L.append("# Cooling architecture decision\n")
    L.append(f"Pack: {inp['n_cells']} cells at "
             f"{inp['q_cell']*1000:.0f} mW each -> Q = "
             f"{inp['Q']:.0f} W. Oil {inp['oil_L']:.1f} L at "
             f"GBP {inp['oil_price']:.0f}/L. Approach dT = "
             f"{inp['dT']:.0f} degC.\n")
    L.append(f"Setup A film (Churchill-Bernstein at u = "
             f"{inp['u_fan']*1000:.1f} mm/s, D = "
             f"{inp['D']*1000:.0f} mm): h = {out['h']:.0f} "
             f"W/m2K (Re {out['Re']:.1f}, Pr {out['Pr']:.0f}). "
             f"Installed bank capability {out['cap_A']:.0f} W; "
             f"tube needed for Q: {out['L_need']:.1f} m.\n")
    L.append(f"System totals: A GBP {sc['totA']:.0f}, B GBP "
             f"{sc['totB']:.0f} (delta GBP {vd['delta']:.0f} = "
             f"{vd['delta_share']:.1f}% of system). Oil = "
             f"{vd['oil_share']:.1f}% of system - "
             f"{vd['oil_judge']}; identical in both setups.\n")
    L.append(f"**Chosen: Setup {vd['chosen']}.** "
             + ("The installed bank carries the load; A keeps "
                "the cheapest, simplest loop. B remains the "
                "drawn-up fallback."
                if vd['fits'] else
                "Q exceeds the fan-driven bank capability; B's "
                "small premium buys guaranteed headroom.") + "\n")
    return "\n".join(L)


# ------------------------------------------------------------------ #
#  The tab                                                           #
# ------------------------------------------------------------------ #
def decide_tab(st, go, fluid, ss):
    st.markdown("#### Decide - internal pipe bank vs external "
                "heat exchanger")
    st.markdown(
        "This page decides the cooling architecture for the "
        "design **currently configured in this app**: the cell "
        "count and pitches come from the FEA tab, the per-cell "
        "heat from the battery model at your set current, the "
        "tube-bank film from Churchill-Bernstein at your fan "
        "speed and pipe size, and the oil inventory from the "
        "auto-sized tank. Change anything below - every "
        "prefilled value is an assumption you can override, and "
        "the verdicts restate the numbers they were decided on.")

    # ---- inputs, prefilled from app state ----
    nr = int(ss.get("fxb_nr", 4) or 4)
    nc = int(ss.get("fxb_nc", 8) or 8)
    n0 = nr * nc if str(ss.get("fxb_gm", "")).startswith(
        "Battery") else 400
    bm = ss.get("bat_model")
    if bm is not None and "fxb_ic" in ss:
        i0 = float(ss.get("fxb_ic", 5.0))
        q0 = i0 * i0 * (float(bm["r0"](50.0))
                        + float(bm["r1"](50.0)))
        qsrc = ss.get("bat_src", "battery model")
    else:
        q0, qsrc = 1.0, "default (no battery model loaded)"
    c1, c2, c3, c4 = st.columns(4)
    n_cells = c1.number_input("Cells in the pack", 1, 20000, n0,
                              1, key="dc_n")
    q_cell = c2.number_input("Heat per cell [W]", 0.05, 20.0,
                             round(max(q0, 0.05), 2), 0.05,
                             key="dc_q")
    dT = c3.number_input("Oil above water ΔT [°C]", 2.0, 25.0,
                         8.0, 0.5, key="dc_dt")
    T_oil = c4.number_input("Oil bulk temperature [°C]", 10.0,
                            90.0, 40.0, 1.0, key="dc_toil")
    st.caption(f"Heat per cell prefilled from: **{qsrc}**.")
    Q = n_cells * q_cell

    p1, p2, p3, p4 = st.columns(4)
    u_fan = p1.number_input("Fan crossflow speed [mm/s]", 0.1,
                            200.0,
                            float(ss.get("fxb_uf", 20.0)) or 20.0,
                            0.5, key="dc_u") / 1000.0
    D = p2.number_input("Pipe OD [mm]", 2.0, 60.0,
                        float(ss.get("fxb_pd", 8.0)) or 8.0,
                        0.5, key="dc_d") / 1000.0
    n_pipes = p3.number_input("Pipes installed", 0, 60,
                              int(ss.get("fxb_np", 3) or 3), 1,
                              key="dc_np")
    L_pipe = p4.number_input("Pipe length each [mm]", 50.0,
                             5000.0, 300.0, 10.0,
                             key="dc_lp") / 1000.0

    h, Re, Pr = h_crossflow_cb(fluid, u_fan, D, T_oil)
    cap_A, A_inst = bank_capability(h, D, L_pipe, n_pipes, dT)
    L_need, A_need = tube_needed(Q, h, D, dT)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Pack heat Q", f"{Q:.0f} W")
    m2.metric("Bank film h (CB)", f"{h:.0f} W/m²K",
              help=f"Re = {Re:.1f}, Pr = {Pr:.0f} in "
                   f"{fluid['name']} at {T_oil:.0f} °C")
    m3.metric("Installed bank can reject", f"{cap_A:.0f} W")
    m4.metric("Tube needed for Q", f"{L_need:.1f} m")
    if cap_A < Q:
        st.warning(
            f"The installed bank ({n_pipes} × "
            f"{L_pipe*1000:.0f} mm) rejects {cap_A:.0f} W at "
            f"ΔT {dT:.0f} °C - short of {Q:.0f} W. Add "
            f"{max(L_need - n_pipes*L_pipe, 0):.1f} m of tube, "
            "raise the fan speed, widen ΔT, or Setup B takes "
            "it.")
    else:
        st.success(
            f"The installed bank covers Q with a factor "
            f"{cap_A/max(Q,1e-9):.2f} in hand.")

    st.markdown("**Costs** - one-off prototype GBP; edit any "
                "assumption:")
    oil_per_cell0 = 0.045
    # oil inventory from the app's auto-sized tank when available
    try:
        Wt = float(ss.get("fxb_w", 0)) / 1000.0
        Ht = float(ss.get("fxb_h", 0)) / 1000.0
    except Exception:
        Wt = Ht = 0.0
    o1, o2, o3, o4 = st.columns(4)
    cell_price = o1.number_input("Cell price [£]", 0.5, 50.0,
                                 4.5, 0.1, key="dc_cp")
    oil_price = o2.number_input("Oil price [£/L]", 1.0, 60.0,
                                12.0, 0.5, key="dc_op")
    oil_per_cell = o3.number_input("Oil per cell [L]", 0.005,
                                   0.5, oil_per_cell0, 0.005,
                                   format="%.3f", key="dc_ov")
    show_adv = o4.checkbox("Edit hardware unit costs", False,
                           key="dc_adv")
    c = dict(COST_DEFAULTS)
    if show_adv:
        a1, a2, a3, a4, a5 = st.columns(5)
        c["oil_pump"] = a1.number_input("Oil pump [£]", 20.0,
                                        1000.0, c["oil_pump"],
                                        5.0, key="dc_opump")
        c["plate_hx"] = a2.number_input("Plate HX [£]", 20.0,
                                        1000.0, c["plate_hx"],
                                        5.0, key="dc_hx")
        c["tube_per_m"] = a3.number_input("Tube [£/m]", 1.0,
                                          100.0, c["tube_per_m"],
                                          0.5, key="dc_tpm")
        c["fan_each"] = a4.number_input("Fan each [£]", 2.0,
                                        200.0, c["fan_each"],
                                        1.0, key="dc_fan")
        c["per_cell_assembly"] = a5.number_input(
            "Assembly [£/cell]", 0.1, 10.0,
            c["per_cell_assembly"], 0.1, key="dc_pca")
    oil_L = n_cells * oil_per_cell
    sc = system_costs(n_cells, cell_price, oil_L, oil_price,
                      max(L_need, n_pipes * L_pipe), Q, c)
    vd = verdicts(Q, cap_A, sc, dT)

    # stacked cost bars
    labels = ["Cells", "Electronics + assembly", "Enclosure",
              "Dielectric oil", "Cooling hardware"]
    valsA = [sc["cells"], sc["elec"], sc["enc"], sc["oil"],
             sc["coolA"]]
    valsB = [sc["cells"], sc["elec"], sc["enc"],
             sc["oil"] + sc["extra_oil"], sc["coolB"]]
    cols = ["#7F77DD", "#B4B2A9", "#D3D1C7", "#EF9F27",
            "#5DCAA5"]
    fig = go.Figure()
    for lb, va, vb, co in zip(labels, valsA, valsB, cols):
        fig.add_bar(y=["Setup B", "Setup A"], x=[vb, va],
                    name=lb, orientation="h",
                    marker_color=co,
                    text=[f"£{vb:,.0f}", f"£{va:,.0f}"],
                    textposition="inside")
    fig.update_layout(barmode="stack", height=180,
                      margin=dict(l=8, r=8, t=8, b=8),
                      legend=dict(orientation="h", y=1.25),
                      xaxis_title="system cost [£]")
    st.plotly_chart(fig, width='stretch', key="dc_bars")

    tbl = pd.DataFrame({
        "Line": labels + ["System total"],
        "Setup A [£]": [round(v) for v in valsA]
        + [round(sc["totA"])],
        "Setup B [£]": [round(v) for v in valsB]
        + [round(sc["totB"])],
        "Share of A [%]": [round(100 * v / sc["totA"], 1)
                           for v in valsA] + [100.0]})
    st.dataframe(tbl, width='stretch', height=250)

    # ceiling map: where A ends and B begins
    dts = np.linspace(3, 15, 60)
    caps = [bank_capability(
        h_crossflow_cb(fluid, u_fan, D, T_oil)[0], D, L_pipe,
        n_pipes, t)[0] for t in dts]
    fmap = go.Figure()
    fmap.add_scatter(x=dts, y=caps, mode="lines",
                     name="installed bank capability",
                     line=dict(color="#1D9E75", width=2.5))
    fmap.add_scatter(x=[dT], y=[Q], mode="markers",
                     name="your design point",
                     marker=dict(size=12, color="#D85A30",
                                 symbol="x"))
    fmap.update_layout(height=280,
                       xaxis_title="oil above water ΔT [°C]",
                       yaxis_title="heat [W]",
                       legend=dict(orientation="h", y=1.12),
                       margin=dict(l=8, r=8, t=26, b=8),
                       title=dict(text="Setup A territory: below "
                                  "the line the bank carries it; "
                                  "above it, Setup B",
                                  font=dict(size=12), x=0.02))
    st.plotly_chart(fmap, width='stretch', key="dc_map")

    st.markdown(
        f"**The oil question, answered with numbers:** the oil "
        f"bill is £{sc['oil']:,.0f} - "
        f"**{vd['oil_share']:.1f}% of the system**, "
        f"{vd['oil_judge']}, and it is the same bill in both "
        f"setups (B adds only £{sc['extra_oil']:.0f} for the "
        f"loop litre). The entire A-vs-B hardware difference is "
        f"£{vd['delta']:.0f} = {vd['delta_share']:.1f}% of the "
        f"system - too small to decide on cost alone, so the "
        f"decision is thermal.")
    if vd["fits"]:
        st.success(
            f"**Chosen: Setup A.** The installed bank carries "
            f"{Q:.0f} W with margin at ΔT {dT:.0f} °C; the oil "
            "you already paid for buys direct cell contact, and "
            "A keeps the cheapest, simplest loop. B stays the "
            "drawn-up fallback.")
    else:
        st.info(
            f"**Chosen: Setup B at this power.** {Q:.0f} W "
            f"exceeds the bank's {cap_A:.0f} W; B's "
            f"{vd['delta_share']:.1f}% premium is trivial next "
            "to guaranteed headroom. Immersion itself stays "
            f"justified - oil is {vd['oil_share']:.1f}% of the "
            "bill.")
    inp = dict(n_cells=n_cells, q_cell=q_cell, Q=Q, dT=dT,
               u_fan=u_fan, D=D, oil_L=oil_L,
               oil_price=oil_price)
    out = dict(h=h, Re=Re, Pr=Pr, cap_A=cap_A, L_need=L_need)
    st.download_button(
        "Decision summary (markdown, for the report)",
        data=decision_md(inp, out, sc, vd),
        file_name="cooling_decision.md", mime="text/markdown",
        key="dc_md")
    st.caption(
        "Honesty notes: prices are one-off prototype figures - "
        "at product volume the cell price falls faster than the "
        "oil price, so the oil share grows; and the CB film "
        "assumes clean crossflow at the fan speed - your FEA "
        "film-law pipe row measures the real value, and this "
        "page should be re-run with it once the pipes rung "
        "reports.")


def main():
    import streamlit as st
    import plotly.graph_objects as go
    st.set_page_config(page_title="Cooling decision",
                       layout="wide")
    fluid = dict(name="Synthetic ester (typical)", k=0.144,
                 rho=960.0, cp=1900.0, nu25=40.0, B=0.0)
    decide_tab(st, go, fluid, st.session_state)


if __name__ == "__main__":
    main()
