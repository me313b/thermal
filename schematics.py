"""schematics.py - labelled, interactive schematics for the plate-channel
architecture and the zonal solve topology. Plotly figures (hover/zoom,
house style) so the geometry symbols and the network structure have a
picture to sit beside the equations.
"""
import numpy as np
import plotly.graph_objects as go

INK = "#0F172A"
BRAND = "#6366F1"
OIL = "#38BDF8"
PLATE = "#94A3B8"
CELL = "#FDE68A"
CELL_L = "#F59E0B"
WATER = "#0EA5E9"
GRIDL = "rgba(80,90,120,.35)"


def _base(fig, h, title):
    fig.update_layout(
        height=h, title=dict(text=title, x=0.01, xanchor="left"),
        plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=8, r=8, t=44, b=8), showlegend=False,
        font=dict(family="Inter, sans-serif", size=12, color="#334155"))
    fig.update_xaxes(visible=False, scaleanchor="y", scaleratio=1)
    fig.update_yaxes(visible=False)
    return fig


def schematic_unit(d, g, lay):
    """Plan-view slice through the pack: cells, plates in the gaps, the
    oil slot s, the in-row crevice g_x, tube cross-sections at the
    derived pitch, with a zoom callout of the lens subchannel fea4
    solves."""
    D = d["d_cell"]; R = D / 2
    p = d["pitch"]; s = d["s_nom"]; tp = d["plate_t"]
    gx = p - D
    ncol, nrow = 4, 3
    row_gap = D + s * 2 + tp                 # cell + two slots + plate
    fig = go.Figure()

    # plates between rows (horizontal bands)
    for r in range(nrow - 1):
        yc = (r + 0.5) * row_gap + R
        fig.add_shape(type="rect", x0=-0.3 * p, x1=(ncol - 0.7) * p,
                      y0=yc - tp / 2, y1=yc + tp / 2,
                      fillcolor=PLATE, line=dict(width=0), layer="below")
        # tubes piercing this plate at derived xk
        for xk in lay["xk"]:
            if xk <= (ncol - 0.5) * p:
                fig.add_shape(type="circle", xref="x", yref="y",
                              x0=xk - d["tube_od"] / 2,
                              x1=xk + d["tube_od"] / 2,
                              y0=yc - d["tube_od"] / 2,
                              y1=yc + d["tube_od"] / 2,
                              fillcolor=WATER,
                              line=dict(color="#0369A1", width=1))

    # cells
    for r in range(nrow):
        yc = r * row_gap + R
        for c in range(ncol):
            xc = c * p
            fig.add_shape(type="circle", x0=xc - R, x1=xc + R,
                          y0=yc - R, y1=yc + R, fillcolor=CELL,
                          line=dict(color=CELL_L, width=1.5))

    # oil-up arrows in the slots
    for c in (0, ncol - 1):
        xc = c * p
        fig.add_annotation(x=xc, y=row_gap * 0.5 + R, ax=xc,
                           ay=row_gap * 0.5 + R + 0.5 * p,
                           xref="x", yref="y", axref="x", ayref="y",
                           showarrow=True, arrowhead=2, arrowwidth=2,
                           arrowcolor=OIL)

    # dimension labels
    y0 = R
    fig.add_annotation(x=0, y=-R - 0.28 * p, ax=p, ay=-R - 0.28 * p,
                       xref="x", yref="y", axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowside="end+start",
                       arrowwidth=1.2, arrowcolor=INK)
    fig.add_annotation(x=0.5 * p, y=-R - 0.42 * p,
                       text=f"pitch p = {p*1000:.1f} mm",
                       showarrow=False, font=dict(size=11, color=INK))
    # crevice gx between two cells in a row
    fig.add_annotation(x=0.5 * p, y=y0, text=f"g_x={gx*1000:.1f}",
                       showarrow=False, font=dict(size=10, color="#B91C1C"))
    # slot s between cell top and plate
    yplate = 0.5 * row_gap + R
    fig.add_annotation(x=(ncol - 1) * p + R + 0.02 * p,
                       y=(R + (yplate - tp / 2)) / 2 + 0.5 * R,
                       text=f"slot s={s*1000:.1f}", showarrow=False,
                       font=dict(size=10, color=OIL), textangle=0,
                       xanchor="left")
    fig.add_annotation(x=lay["xk"][0], y=yplate + tp / 2 + 0.12 * p,
                       text=f"water tube (OD {d['tube_od']*1000:.0f})",
                       showarrow=False, font=dict(size=10, color="#0369A1"))
    fig.add_annotation(x=(ncol - 1.5) * p, y=(nrow - 1) * row_gap + R,
                       text="cells (isothermal cans)", showarrow=False,
                       font=dict(size=10, color=CELL_L), xanchor="center")

    _base(fig, 430, "Plate-channel unit (plan view): oil rises through "
                    "the slots, water crosses through tubes in the plates")
    fig.update_xaxes(range=[-0.6 * p, (ncol - 0.2) * p])
    fig.update_yaxes(range=[-R - 0.6 * p, (nrow - 1) * row_gap + R
                            + 0.3 * p])
    return fig


def schematic_lens(d, g):
    """Zoom of the true lens subchannel fea4 meshes: one pitch wide,
    cell crown vs flat plate, showing the masked oil domain and the two
    heated walls."""
    D = d["d_cell"]; R = D / 2; p = d["pitch"]; s = d["s_nom"]
    fig = go.Figure()
    # plate (top wall)
    fig.add_shape(type="rect", x0=-0.5 * p, x1=0.5 * p,
                  y0=R + s, y1=R + s + 0.18 * p, fillcolor=PLATE,
                  line=dict(width=0))
    # two neighbour cell crowns
    for dxc in (-p, 0, p):
        fig.add_shape(type="circle", x0=dxc - R, x1=dxc + R,
                      y0=-0.4 * p, y1=R, fillcolor=CELL,
                      line=dict(color=CELL_L, width=1.5))
    # oil lens (shaded region between crowns and plate, one pitch)
    xs = np.linspace(-0.5 * p, 0.5 * p, 60)
    crown = np.sqrt(np.maximum(R**2 - xs**2, 0.0))
    fig.add_trace(go.Scatter(
        x=np.concatenate([xs, xs[::-1]]),
        y=np.concatenate([crown, np.full_like(xs, R + s)[::-1]]),
        fill="toself", fillcolor="rgba(56,189,248,.28)",
        line=dict(color=OIL, width=1), hoverinfo="skip"))
    fig.add_annotation(x=0, y=R + s + 0.09 * p, text="plate wall  T_p",
                       showarrow=False, font=dict(size=11, color=INK))
    fig.add_annotation(x=0, y=0.25 * R, text="cell wall  T_c",
                       showarrow=False, font=dict(size=11, color="#92400E"))
    fig.add_annotation(x=0.28 * p, y=R + 0.35 * s, text="oil (lens)",
                       showarrow=False, font=dict(size=11, color="#0369A1"))
    fig.add_annotation(x=-0.5 * p, y=-0.3 * p, ax=0.5 * p, ay=-0.3 * p,
                       xref="x", yref="y", axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowside="end+start",
                       arrowwidth=1.2, arrowcolor=INK)
    fig.add_annotation(x=0, y=-0.44 * p, text="one pitch (periodic)",
                       showarrow=False, font=dict(size=10, color=INK))
    _base(fig, 300, "The lens subchannel fea4 solves (cross-section)")
    fig.update_xaxes(range=[-0.62 * p, 0.62 * p])
    fig.update_yaxes(range=[-0.5 * p, R + s + 0.24 * p])
    return fig


def schematic_network(d, lay):
    """Solve topology: the chain each watt follows and the closed-form /
    contractive relation on every edge."""
    fig = go.Figure()
    nodes = [
        (0.5, 5.2, "CELL i,j", CELL_L,
         "q_gen(T)=I²R·e^(−k(T−25))"),
        (0.5, 4.0, "OIL BULK  T_b(z)", OIL,
         "exact exp. march"),
        (0.5, 2.8, "PLATE FACE  T_f", PLATE,
         "closed-form vs root"),
        (0.5, 1.6, "ROOT (tube)", "#38BDF8",
         "R_ct+R_wall+R_film"),
        (0.5, 0.4, "WATER  T_w(tube)", WATER,
         "marched tube-by-tube"),
        (3.4, 4.0, "PLENUM", "#A78BFA",
         "T*=⟨A⟩/(1−⟨P⟩)"),
    ]
    for x, y, lab, col, sub in nodes:
        fig.add_shape(type="rect", x0=x - 0.42, x1=x + 0.42,
                      y0=y - 0.28, y1=y + 0.28, fillcolor=col,
                      opacity=0.22, line=dict(color=col, width=2))
        fig.add_annotation(x=x, y=y + 0.08, text=f"<b>{lab}</b>",
                           showarrow=False, font=dict(size=11, color=INK))
        fig.add_annotation(x=x, y=y - 0.13, text=sub, showarrow=False,
                           font=dict(size=9.5, color="#475569"))

    def edge(y0, y1, lab):
        fig.add_annotation(x=0.5, y=y1 + 0.28, ax=0.5, ay=y0 - 0.28,
                           xref="x", yref="y", axref="x", ayref="y",
                           showarrow=True, arrowhead=2, arrowwidth=1.8,
                           arrowcolor=INK)
        fig.add_annotation(x=0.98, y=(y0 + y1) / 2, text=lab,
                           showarrow=False, xanchor="left",
                           font=dict(size=9.5, color="#334155"))
    edge(5.2, 4.0, "a₀₀(z*)·(T_c−T_b)")
    edge(4.0, 2.8, "a₁₁·η_bay·(T_f−T_b)")
    edge(2.8, 1.6, "q=(T_f−T_w)/R_root")
    edge(1.6, 0.4, "ṁ_w c_p dT_w")
    # recirculation loop
    fig.add_annotation(x=2.98, y=3.7, ax=0.92, ay=0.4,
                       xref="x", yref="y", axref="x", ayref="y",
                       showarrow=True, arrowhead=2, arrowwidth=1.5,
                       arrowcolor="#A78BFA")
    fig.add_annotation(x=0.92, y=4.0, ax=2.98, ay=4.2,
                       xref="x", yref="y", axref="x", ayref="y",
                       showarrow=True, arrowhead=2, arrowwidth=1.5,
                       arrowcolor="#A78BFA")
    fig.add_annotation(x=2.0, y=4.45, text="oil recirculates",
                       showarrow=False, font=dict(size=9.5,
                                                  color="#7C3AED"))
    fig.add_annotation(x=2.0, y=0.05, text="mixed exit feeds plenum",
                       showarrow=False, font=dict(size=9.5,
                                                  color="#7C3AED"))
    _base(fig, 430, "Zonal solve topology: every edge is closed-form or "
                    "contractive (no ratcheting)")
    fig.update_xaxes(range=[-0.1, 4.1])
    fig.update_yaxes(range=[-0.1, 5.7])
    return fig


def schematic_oil_loop(d, g, xp=None):
    """Closed dielectric loop with an EXTERNAL pump and NO external heat
    exchanger: oil leaves the pack, passes the pump, and returns; the
    internal water tubes (crossing the pack) do the heat rejection, so
    the oil and water circuits never meet."""
    fig = go.Figure()
    _base(fig, 320, "Closed oil loop - external pump, internal water HX")
    # pack box
    fig.add_shape(type="rect", x0=0.30, y0=0.12, x1=0.94, y1=0.88,
                  line=dict(color="#94A3B8", width=2),
                  fillcolor="rgba(245,158,11,.10)")
    fig.add_annotation(x=0.62, y=0.93, text="pack (flooded)",
                       showarrow=False, font=dict(size=11, color="#64748B"))
    # cells + plates + embedded tubes (three units)
    for i in range(3):
        x0 = 0.38 + i * 0.17
        fig.add_shape(type="rect", x0=x0, y0=0.22, x1=x0 + 0.075, y1=0.78,
                      line=dict(color="#B45309", width=1),
                      fillcolor="rgba(241,82,82,.55)")
        if i < 2:
            px = x0 + 0.095
            fig.add_shape(type="rect", x0=px, y0=0.22, x1=px + 0.035,
                          y1=0.78, line=dict(color="#CBD5E1", width=1),
                          fillcolor="rgba(203,213,225,.75)")
            fig.add_shape(type="rect", x0=px + 0.006, y0=0.46,
                          x1=px + 0.029, y1=0.54,
                          line=dict(color="#0369A1", width=1),
                          fillcolor="#38BDF8")
    fig.add_annotation(x=0.62, y=0.50, ax=0.62, ay=0.50, text="",
                       showarrow=False)
    fig.add_annotation(x=0.655, y=0.36, text="water tubes<br>(into page)",
                       showarrow=False, font=dict(size=9, color="#0369A1"))
    # oil channel arrows through the pack (left to right)
    for yy in (0.30, 0.50, 0.70):
        fig.add_annotation(x=0.90, y=yy, ax=0.34, ay=yy,
                           xref="x", yref="y", axref="x", ayref="y",
                           showarrow=True, arrowhead=3, arrowwidth=1.6,
                           arrowcolor="#F59E0B")
    # external loop: supply (bottom) and return (top) to the pump
    fig.add_shape(type="line", x0=0.30, y0=0.20, x1=0.10, y1=0.20,
                  line=dict(color="#F59E0B", width=4))
    fig.add_shape(type="line", x0=0.10, y0=0.20, x1=0.10, y1=0.80,
                  line=dict(color="#F59E0B", width=4))
    fig.add_shape(type="line", x0=0.10, y0=0.80, x1=0.30, y1=0.80,
                  line=dict(color="#F59E0B", width=4))
    fig.add_annotation(x=0.20, y=0.84, ax=0.13, ay=0.84, showarrow=True,
                       arrowhead=3, arrowcolor="#F59E0B", text="")
    fig.add_annotation(x=0.13, y=0.16, ax=0.20, ay=0.16, showarrow=True,
                       arrowhead=3, arrowcolor="#F59E0B", text="")
    # pump symbol
    fig.add_shape(type="circle", x0=0.055, y0=0.42, x1=0.145, y1=0.58,
                  line=dict(color="#F59E0B", width=3),
                  fillcolor="rgba(245,158,11,.15)")
    fig.add_annotation(x=0.10, y=0.50, text="P", showarrow=False,
                       font=dict(size=15, color="#B45309"))
    lab = "external oil pump"
    if xp:
        lab += (f"<br>{xp['Vdot_lpm']:.0f} L/min - "
                f"{xp['dp']/1000:.1f} kPa - {xp['P']:.0f} W")
    fig.add_annotation(x=0.10, y=0.32, text=lab, showarrow=False,
                       font=dict(size=10, color="#B45309"))
    fig.add_annotation(x=0.10, y=0.66, text="no external HX",
                       showarrow=False, font=dict(size=9, color="#64748B"))
    fig.add_annotation(
        x=0.62, y=0.05,
        text="oil loop (amber) and water loop (blue) never meet - the "
             "tubes reject the heat inside the pack",
        showarrow=False, font=dict(size=10, color="#475569"))
    return fig
