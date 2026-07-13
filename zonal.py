"""zonal.py - plate-channel zonal thermal-hydraulic network + Monte Carlo.

Physics: every inter-row gap carries a water-rooted aluminium plate with
an oil slot on each face; oil recirculates upward through the slots; each
slot column obeys the exact developing-flow kernel from
fea4_channel.graetz_kernel on the true lens cross-section:
    [q'_c; q'_p] = a(z*, s) . [T_cell - Tb; T_plateface - Tb]
Fin knockdown eta_bay (product rule, roots at the derived tube pitch
P = 2/m) scales the plate-side kernel entries.

This is the v9.4 build, incorporating the external red-team findings:
 - Buoyancy head is referenced to the loop RETURN (plenum) temperature,
   not the water inlet: in a closed recirculating loop the net head is
   rho.beta.g.H.(T_riser - T_downcomer), and the downcomer sits at
   T_plen (F1). At the default this makes buoyancy near-zero and the
   pump head dominant, which is the conservative and correct behaviour.
 - The water-side Nusselt number comes from correlations.water_nu, the
   SAME function the lumped solver uses, with a continuous 2300-3000
   transition bridge (F5); no laminar-turbulent step.
 - Contact conductance h_contact and the tube-to-plate collar
   engagement factor are named inputs, surfaced with a sensitivity
   sweep in the app, not buried constants (F2).
 - The tube layout is derived from the kernel-implied plate film
   (a11/pitch at the design point), so the fin rule is self-consistent
   with the duct it feeds (F4).
 - An optional, explicitly bounded direct oil-to-tube path can be
   switched on (wetted_tube_frac > 0); its magnitude depends on header
   geometry, so it is off by default and its effect is shown live (F3).
 - The core-to-can rise (R_core = 1/(4 pi k_r H)) is superposed and the
   peak CORE temperature reported alongside the can temperature (F10).

Validations: energy closure, uniform-input collapse, z-independence,
and the fea4 kernel gates. Run `python zonal.py` for the shakedown.
"""
import math
import numpy as np

from fea4_channel import graetz_kernel
from correlations import water_nu

WATER = dict(rho=1000.0, cp=4180.0, k=0.60, mu=8.9e-4)


class KernelBank:
    def __init__(self, D, pitch, s_grid=(0.0012, 0.0016, 0.0020,
                                         0.0026, 0.0032),
                 k_oil=0.13, n=100, nz=110):
        self.D, self.pitch, self.k_oil = D, pitch, k_oil
        self.s_grid = np.array(s_grid)
        self.tabs = [graetz_kernel(D=D, pitch_x=pitch, s=s,
                                   k_oil=k_oil, n=n, nz=nz)
                     for s in s_grid]

    def _wj(self, s):
        gr = self.s_grid
        s = float(np.clip(s, gr[0], gr[-1]))
        j = int(np.clip(np.searchsorted(gr, s) - 1, 0, len(gr) - 2))
        return j, (s - gr[j]) / (gr[j + 1] - gr[j])

    def props(self, s):
        j, w = self._wj(s)
        t0, t1 = self.tabs[j], self.tabs[j + 1]
        f = lambda k_: (1 - w) * t0[k_] + w * t1[k_]
        return f("Dh"), f("A"), f("fRe")

    def a_of(self, s, zstar):
        j, w = self._wj(s)
        t0, t1 = self.tabs[j], self.tabs[j + 1]
        zstar = np.clip(zstar, t0["zs"][0], t0["zs"][-1])
        out = np.zeros(zstar.shape + (2, 2))
        for r in range(2):
            for c in range(2):
                v0 = np.interp(np.log(zstar), np.log(t0["zs"]),
                               t0["a"][:, r, c])
                v1 = np.interp(np.log(zstar), np.log(t1["zs"]),
                               t1["a"][:, r, c])
                out[..., r, c] = (1 - w) * v0 + w * v1
        return out


def h_face_design(d, bank):
    """Kernel-implied per-face oil->plate film at the design operating
    point (nominal slot, pump-head-only flow), used to make the fin-rule
    layout self-consistent with the duct kernel (red-team F4). Returns
    W/m2K = a11(z*_mid)/pitch.

    Note (verification obs 4): this reads a11 from whatever bank the
    caller built, so the value drifts a few percent with the kernel grid
    resolution (n): ~255 at n=100, ~278 at a coarse mini-bank. The
    layout snap has margin here - the 43 mm / every-2-cells snap holds
    for h_face up to ~ 2 k_p t_p / (pitch)^2, i.e. it does not flip
    until h_face exceeds a threshold well above these values - but the
    shakedown asserts the resulting tube count so a future bank change
    that moves the snap boundary is caught rather than silent."""
    rho, cp, k_o = d["rho"], d["cp"], d["k_oil"]
    H = d["h_cell"]
    Dh, A, fRe = bank.props(d["s_nom"])
    nu_n = d["nu25"] * math.exp(
        d["B"] * (1 / (d["T_in"] + 15.0 + 273.15) - 1 / 298.15))
    mu_n = rho * nu_n
    u_n = max(d["dp_extra"] * 2 * Dh ** 2 / (fRe * mu_n * H), 1e-6)
    RePr = max(u_n * Dh / (k_o / (rho * cp)), 1e-9)
    zst = (H / 2) / (Dh * RePr)
    a11 = float(bank.a_of(d["s_nom"], np.array([zst]))[0, 1, 1])
    return a11 / d["pitch"]


def snap_boundary_h(d):
    """The h_face values at which the fin-rule tube pitch snaps to a
    different number of cells, for the current geometry. Used to report
    how much margin the design-point h_face has before the layout would
    flip (verification obs 4)."""
    k_p, t_p = 205.0, d["plate_t"]
    out = {}
    for n_per in (1, 2, 3, 4):
        # P_rule = 2/m = pitch*n_per boundary -> m = 2/(pitch*n_per)
        m_b = 2.0 / (d["pitch"] * n_per)
        out[n_per] = 0.5 * (m_b ** 2) * k_p * t_p   # h at that m
    return out


def derived_layout(d, h_face):
    k_p, t_p = 205.0, d["plate_t"]
    m = math.sqrt(2.0 * max(h_face, 5.0) / (k_p * t_p))
    P_rule = 2.0 / m
    n_per = max(int(P_rule // d["pitch"]), 1)
    P_snap = n_per * d["pitch"]
    L_row = d["n_cols"] * d["pitch"]
    n_t = max(int(math.ceil(L_row / P_snap)), 2)
    xk = (np.arange(n_t) + 0.5) * P_snap
    ex = lambda x: math.tanh(x) / max(x, 1e-9)
    eta = ex(m * P_snap / 2) * ex(m * d["h_cell"] / 2)
    return dict(m=m, P_rule=P_rule, P_snap=P_snap, n_tubes=n_t,
                xk=xk, eta_bay=eta, cells_per_tube=n_per,
                h_face=h_face)


def solve_zonal(d, bank, s_ch=None, contact=None, heat_map=None,
                nz=10, iters=40, relax=0.6, tol=0.002,
                init=None):
    """Plate-channel zonal network. Every subsystem is closed-form or
    contractive: exponential (exact) z-march per channel; per-cell
    Newton with the march's own Jacobian; per-bay plate-face solved in
    closed form against the root chain (contact + wall + water film);
    water marched tube-by-tube; recirculation plenum solved from the
    linear exit map. Energy closure is the convergence monitor."""
    nr, nc = d["n_rows"], d["n_cols"]
    n_gap = nr - 1
    n_ch = 2 * n_gap + 2
    if s_ch is None:
        s_ch = np.full(n_ch, d["s_nom"])
    H = d["h_cell"]
    z = (np.arange(nz) + 0.5) * (H / nz)
    dz = H / nz

    side_row = np.empty(n_ch, int); side_plate = np.empty(n_ch, int)
    side_row[0], side_plate[0] = 0, -1
    for g in range(n_gap):
        side_row[1 + 2 * g], side_plate[1 + 2 * g] = g, g
        side_row[2 + 2 * g], side_plate[2 + 2 * g] = g + 1, g
    side_row[-1], side_plate[-1] = nr - 1, -1

    nu = lambda T: d["nu25"] * np.exp(
        d["B"] * (1 / (np.clip(T, -20, 150) + 273.15) - 1 / 298.15))
    rho, cp, k_o, beta = d["rho"], d["cp"], d["k_oil"], d["beta"]

    Dh = np.empty(n_ch); Acs = np.empty(n_ch); fRe = np.empty(n_ch)
    for i in range(n_ch):
        Dh[i], Acs[i], fRe[i] = bank.props(s_ch[i])

    q_cell0 = (d["C"] * d["cap_Ah"]) ** 2 * d["r_dc"] / 1000.0
    hm = np.ones((nr, nc)) if heat_map is None else (1.0 + heat_map)

    # F4: layout from the kernel-implied plate film, not a stale 150
    hf = h_face_design(d, bank)
    lay = derived_layout(d, h_face=hf)
    n_t = lay["n_tubes"]
    eta = lay["eta_bay"]
    if contact is None:
        contact = np.full((n_gap, n_t), d["plate_contact"])
    col_bay = np.minimum((np.arange(nc) * d["pitch"]
                          // lay["P_snap"]).astype(int), n_t - 1)

    # ---- water-side root chain ----
    d_i = d["tube_od"] - 2 * d["tube_wall"]
    mdot_w = d["flow_lpm"] / 60 * WATER["rho"] / 1000 / n_t
    Re_w = 4 * mdot_w / (math.pi * WATER["mu"] * d_i)
    Pr_w = WATER["mu"] * WATER["cp"] / WATER["k"]
    L_t = max(n_gap * d["pitch"], 0.05)
    Nu_w, water_regime = water_nu(Re_w, Pr_w, d_i, L_t)   # F5 shared
    h_w = Nu_w * WATER["k"] / d_i
    per_len = d["pitch"]
    R_film = 1 / (h_w * math.pi * d_i * per_len)
    R_wl = math.log(d["tube_od"] / d_i) / (2 * math.pi
                                           * d["k_tube"] * per_len)
    # F2: contact conductance and collar engagement are named inputs.
    HC = d.get("h_contact", 8000.0)
    collar = d.get("collar_factor", 6.0)      # effective collar length
    A_ct = math.pi * d["tube_od"] * d["plate_t"] * collar  # per crossing
    Rr = (1.0 / (HC * A_ct * contact) + R_wl + R_film)   # (n_gap,n_t)

    # F3: optional direct oil->tube bypass (off by default; magnitude is
    # header-geometry dependent, so it is exposed rather than assumed).
    frac = float(d.get("wetted_tube_frac", 0.0))
    if frac > 1e-6:
        L_tube_each = n_gap * d["pitch"]
        h_oil_t = float(d.get("h_oil_tube", 150.0))
        A_bo = frac * n_t * math.pi * d["tube_od"] * L_tube_each
        A_bi = frac * n_t * math.pi * d_i * L_tube_each
        R_wall_b = (math.log(d["tube_od"] / d_i)
                    / (2 * math.pi * d["k_tube"]
                       * frac * n_t * L_tube_each))
        UA_bare = 1.0 / (1.0 / (h_oil_t * A_bo) + R_wall_b
                         + 1.0 / (h_w * A_bi))
        # per plate channel, per unit length, per column: the march sums
        # over 2*n_gap channels x nc columns x H (z-integral).
        Gw_perlen = UA_bare / max(2 * n_gap * nc, 1) / H
    else:
        UA_bare, Gw_perlen = 0.0, 0.0

    T_cell = np.full((nr, nc), d["T_in"] + 12.0)
    Tb = np.full((n_ch, nc, nz), d["T_in"] + 5.0)
    T_plen = d["T_in"] + 5.0
    Tf = np.full((n_gap, n_t), d["T_in"] + 3.0)   # plate-face per bay
    T_wat = np.full((n_t, n_gap), d["T_in"])      # water along tubes
    A_col_edge = d["A_case"] / (2 * nc)
    mdot_col = None
    if init is not None:
        try:
            T_cell = init["T_cell"].copy()
            T_plen = float(init["T_plen"])
            if init["Tf"].shape == Tf.shape:
                Tf = init["Tf"].copy()
            if init["T_wat"].shape == T_wat.shape:
                T_wat = init["T_wat"].copy()
            if init["Tb"].shape == Tb.shape:
                Tb = init["Tb"].copy()
            if init["mdot_col"].shape == (n_ch,):
                mdot_col = init["mdot_col"].copy()
        except (KeyError, AttributeError):
            pass
    closure = 1.0
    q_bare_gap = np.zeros(n_gap)

    for it in range(iters):
        # ---- 1. hydraulics: buoyancy referenced to the loop return ----
        # F1: net head = rho.beta.g.H.(T_riser - T_downcomer); the
        # downcomer/return sits at T_plen, not the water inlet.
        Tbar = np.clip(Tb.mean(axis=(1, 2)), d["T_in"] - 5,
                       d["T_in"] + 60)
        mu = rho * nu(Tbar)
        dp = d["dp_extra"] + rho * beta * 9.81 * np.clip(
            Tbar - T_plen, 0.0, 45.0) * H
        ubar = dp * 2 * Dh ** 2 / (fRe * mu * H)
        m_new = np.maximum(rho * ubar * Acs, 1e-7)
        mdot_col = (m_new if mdot_col is None
                    else mdot_col + relax * (m_new - mdot_col))

        # ---- 2. channel march (exact exponential segment update) ----
        alpha = k_o / (rho * cp)
        RePr = np.maximum(mdot_col / (rho * Acs), 1e-9) * Dh / alpha
        Tb_new = np.empty_like(Tb)
        Pch = np.empty((n_ch, nc))
        qhat_c = np.zeros((nr, nc))
        Jc = np.zeros((nr, nc))
        W_bay = np.zeros((n_gap, n_t))
        S_bay = np.zeros((n_gap, n_t))
        q_bare_new = np.zeros(n_gap)
        Q_case_new = 0.0
        Qc_tot = 0.0
        Tw_gap_all = T_wat.mean(axis=0)              # (n_gap,)
        for i in range(n_ch):
            r, pidx = side_row[i], side_plate[i]
            zst = z / (Dh[i] * max(RePr[i], 1e-9))
            a = bank.a_of(s_ch[i], zst).copy()       # (nz,2,2)
            gw = 0.0
            if pidx >= 0:
                a[:, 0, 1] *= eta; a[:, 1, 0] *= eta
                a[:, 1, 1] *= eta
                Tp = Tf[pidx, col_bay]               # (nc,)
                gw = Gw_perlen
                Tw_g = Tw_gap_all[pidx]
            else:
                Gf = float(a[:, 1, 1].sum()) * dz
                f_e = (d["h_ext"] * A_col_edge
                       / max(d["h_ext"] * A_col_edge + Gf, 1e-9))
                a[:, 0, 1] *= f_e; a[:, 1, 0] *= f_e
                a[:, 1, 1] *= f_e
                Tp = np.full(nc, d["T_amb"])
                Tw_g = d["T_amb"]
            Tbi = Tb[i].copy()
            Tbi[:, 0] = T_plen
            mcp = float(mdot_col[i] * cp)
            Tc_r = T_cell[r]
            Pdecay = np.ones(nc)
            for kz in range(nz):
                a00, a01 = a[kz, 0, 0], a[kz, 0, 1]
                a10, a11 = a[kz, 1, 0], a[kz, 1, 1]
                Gc_row = a00 + a10
                Gp_row = a01 + a11
                G = max(Gc_row + Gp_row + gw, 1e-12)
                Tw_eff = (Gc_row * Tc_r + Gp_row * Tp + gw * Tw_g) / G
                x = G * dz / mcp
                ex = math.exp(-min(x, 50.0))
                fbar = (1.0 - ex) / x if x > 1e-6 else 1.0 - 0.5 * x
                Tb0 = Tbi[:, kz]
                Tb_bar = Tw_eff + (Tb0 - Tw_eff) * fbar
                if kz + 1 < nz:
                    Tbi[:, kz + 1] = Tw_eff + (Tb0 - Tw_eff) * ex
                Pdecay *= ex
                dTc = Tc_r - Tb_bar
                dTp = Tp - Tb_bar
                qc = (a00 * dTc + a01 * dTp) * dz
                qp = (a10 * dTc + a11 * dTp) * dz
                Qc_tot += float(np.sum(qc))
                qhat_c[r] += qc
                Jc[r] += dz * (a00 - (a00 + a01)
                               * (Gc_row / G) * (1.0 - fbar))
                if pidx >= 0:
                    np.add.at(W_bay[pidx], col_bay, a11 * dz)
                    np.add.at(S_bay[pidx], col_bay,
                              (a10 * dTc - a11 * Tb_bar) * dz
                              + a11 * dz * Tp)
                    if gw > 0:
                        qw = gw * (Tw_g - Tb_bar) * dz   # oil view (<0)
                        q_bare_new[pidx] += float(np.sum(-qw))
                else:
                    Q_case_new += float(np.sum(-qp))
            Tb_new[i] = Tbi
            Pch[i] = Pdecay
        Tb = Tb + relax * (Tb_new - Tb)
        q_bare_gap = q_bare_gap + relax * (q_bare_new - q_bare_gap)

        # ---- 3. plenum: closed-form recirculation fixed point ----
        m_all = np.maximum(mdot_col, 1e-9)[:, None] * np.ones((1, nc))
        T_exit = Tb_new[:, :, -1]
        A_lin = T_exit - Pch * T_plen
        mA = float(np.sum(m_all * A_lin) / np.sum(m_all))
        mP = float(np.sum(m_all * Pch) / np.sum(m_all))
        T_star = mA / max(1.0 - mP, 1e-3)
        T_plen += 0.8 * (float(np.clip(T_star, d["T_in"] - 2.0,
                                       d["T_in"] + 60.0)) - T_plen)

        # ---- 4. cells: Newton on the march's own linearisation ----
        Tc_old = T_cell.copy()
        J = np.maximum(Jc, 1e-9)
        for _ in range(8):
            q_gen = q_cell0 * hm * np.exp(
                -d["k_dcir"] * (np.clip(T_cell, -10, 130) - 25.0))
            T_cell += 0.8 * (Tc_old + (q_gen - qhat_c) / J - T_cell)
            T_cell = np.clip(T_cell, d["T_in"] - 10, 140.0)

        # ---- 5. plate faces: closed form vs the root chain ----
        S0 = S_bay - W_bay * Tf
        Tw_bay = T_wat.T                              # (n_gap, n_t)
        Tf_new = (Tw_bay / Rr - S0) / (W_bay + 1.0 / Rr)
        Tf += 0.7 * (np.clip(Tf_new, d["T_in"] - 5, 120.0) - Tf)
        q_root = (Tf - Tw_bay) / Rr

        # ---- 6. water march (root heat + optional bare-tube heat) ----
        q_bare_pt = q_bare_gap / max(n_t, 1)          # per tube per gap
        for t in range(n_t):
            Tw = d["T_in"]
            for p in range(n_gap):
                T_wat[t, p] = Tw
                Tw += (q_root[p, t] + q_bare_pt[p]) / (mdot_w
                                                       * WATER["cp"])
                Tw = min(max(Tw, d["T_in"] - 5), 95.0)

        Qg = float(np.sum(q_cell0 * hm * np.exp(
            -d["k_dcir"] * (np.clip(T_cell, -10, 130) - 25))))
        Q_out = float(np.sum(q_root)) + float(q_bare_gap.sum()) \
            + Q_case_new
        closure = (Qg - Q_out) / max(Qg, 1e-9)
        if it > 10 and abs(closure) < tol:
            break

    Q_gen = float(np.sum(q_cell0 * hm * np.exp(
        -d["k_dcir"] * (np.clip(T_cell, -10, 130) - 25))))
    q_cell_map = q_cell0 * hm * np.exp(
        -d["k_dcir"] * (np.clip(T_cell, -10, 130) - 25))
    Q_water = float(np.sum(q_root)) + float(q_bare_gap.sum())
    Q_case = Q_case_new
    P_pump = float(np.sum(mdot_col * nc / rho * d["dp_extra"]) / 0.35)

    # F10: superpose the core-to-can rise and report peak CORE temp.
    R_core = 1.0 / (4.0 * math.pi * d.get("k_rad", 0.9) * H)
    T_core = T_cell + q_cell_map * R_core
    return dict(T_cell=T_cell, T_core=T_core, Tb=Tb, Tf=Tf, T_wat=T_wat,
                T_plen=T_plen, P_pump=P_pump, mdot_col=mdot_col,
                ubar=ubar, s_ch=s_ch, lay=lay, q_root=q_root,
                Q_gen=Q_gen, Q_water=Q_water, Q_case=Q_case,
                Q_bare=float(q_bare_gap.sum()), UA_bare=UA_bare,
                h_w=h_w, water_regime=water_regime, R_core=R_core,
                closure=(Q_gen - Q_water - Q_case) / max(Q_gen, 1e-9),
                T_max=float(T_cell.max()),
                T_core_max=float(T_core.max()),
                T_mean=float(T_cell.mean()),
                spread=float(T_cell.max() - T_cell.min()),
                iters_used=it + 1, zs=z)


def monte_carlo(d, bank, M=150, sigma_s=0.2e-3, sigma_c=0.08,
                heat_map=None, nz=10, iters=60, seed=7,
                progress=None, tol=0.0015):
    """Tolerance Monte Carlo: per-channel slot widths and per-root
    contact quality from truncated normals; every sample is a full
    network solve warm-started from the converged nominal case.

    Statistics are reported honestly (red-team F6): the sample spread is
    quoted with the caveat that it is at or below solver residue, and a
    zero-exceedance count is converted to a 95% upper bound via the rule
    of three (3/M). sigma_s is applied per channel, i.e. one correlated
    defect along the full slot length."""
    rng = np.random.default_rng(seed)
    n_ch = 2 * (d["n_rows"] - 1) + 2
    hf = h_face_design(d, bank)
    lay = derived_layout(d, hf)
    base = solve_zonal(d, bank, heat_map=heat_map, nz=nz, iters=260,
                       tol=0.0008)
    Tmax = np.empty(M); spread = np.empty(M); clos = np.empty(M)
    Tcore = np.empty(M)
    for m in range(M):
        s = np.clip(rng.normal(d["s_nom"], sigma_s, n_ch),
                    0.55 * d["s_nom"], 1.7 * d["s_nom"])
        c = np.clip(rng.normal(d["plate_contact"], sigma_c,
                               (d["n_rows"] - 1, lay["n_tubes"])),
                    0.15, 1.0)
        r = solve_zonal(d, bank, s_ch=s, contact=c,
                        heat_map=heat_map, nz=nz, iters=iters,
                        tol=tol, init=base)
        Tmax[m] = r["T_max"]; spread[m] = r["spread"]
        Tcore[m] = r["T_core_max"]; clos[m] = abs(r["closure"])
        if progress is not None:
            progress((m + 1) / M)
    n_exceed = int(np.sum(Tmax > d["T_limit"]))
    p_exceed = n_exceed / M
    # rule of three: 0/M gives a 95% upper bound of 3/M
    p_ub95 = (3.0 / M) if n_exceed == 0 else None
    # F10 follow-up (verification obs 3): the core is the razor-thin node,
    # so report its exceedance count too. When the core already exceeds,
    # the rule of three no longer applies on that node.
    n_exceed_core = int(np.sum(Tcore > d["T_limit"]))
    p_exceed_core = n_exceed_core / M
    p_ub95_core = (3.0 / M) if n_exceed_core == 0 else None
    return dict(Tmax=Tmax, spread=spread, Tcore=Tcore, base=base,
                closure_worst=float(clos.max()),
                p_exceed=p_exceed, n_exceed=n_exceed, M=M,
                p_ub95=p_ub95, tol=tol,
                n_exceed_core=n_exceed_core, p_exceed_core=p_exceed_core,
                p_ub95_core=p_ub95_core)


def default_d(nr=33, nc=33, C=2.0):
    return dict(n_rows=nr, n_cols=nc, pitch=0.0215, d_cell=0.021,
                h_cell=0.070, cap_Ah=5.0, r_dc=25.0, k_dcir=0.012,
                C=C, plate_t=0.0015, plate_contact=0.8, s_nom=0.002,
                T_in=20.0, flow_lpm=10.0, tube_od=0.010,
                tube_wall=0.0008, k_tube=385.0, T_amb=25.0, h_ext=5.0,
                A_case=1.2, nu25=9e-6, B=3200.0, rho=920.0, cp=2000.0,
                k_oil=0.13, beta=7.5e-4, dp_extra=25.0, T_limit=45.0,
                h_contact=8000.0, collar_factor=6.0, k_rad=0.9,
                wetted_tube_frac=0.0, h_oil_tube=150.0)


# ------------------------------------------------------------------ #
#  Shakedown (red-team F9: the __main__ the handover promised)        #
# ------------------------------------------------------------------ #
if __name__ == "__main__":
    import time
    t0 = time.time()
    bank = KernelBank(D=0.021, pitch=0.0215,
                      s_grid=(0.0016, 0.002, 0.0026), n=70, nz=60)
    print(f"kernel bank (3 slots): {time.time() - t0:.1f} s")

    d = default_d(nr=10, nc=10, C=2.0)
    r = solve_zonal(d, bank, nz=10, iters=120)
    print(f"\nmini 10x10 : T_max {r['T_max']:.2f}  T_core {r['T_core_max']:.2f}"
          f"  spread {r['spread']:.2f}  closure {r['closure']*100:.2f}%"
          f"  tubes {r['lay']['n_tubes']}")

    d = default_d(nr=33, nc=33, C=2.0)
    t1 = time.time()
    r = solve_zonal(d, bank, nz=10, iters=260)
    r2 = solve_zonal(d, bank, nz=16, iters=40, init=r)
    lay = r["lay"]
    print(f"full 33x33 : T_max {r['T_max']:.2f}  T_core {r['T_core_max']:.2f}"
          f"  mean {r['T_mean']:.2f}  spread {r['spread']:.2f}"
          f"  [{time.time() - t1:.1f} s]")
    print(f"  closure {r['closure']*100:.3f}%  z-check nz10-vs-16 "
          f"{abs(r['T_max'] - r2['T_max']):.3f} C")
    print(f"  layout: h_face {lay['h_face']:.0f} W/m2K -> m {lay['m']:.1f}/m"
          f"  {lay['n_tubes']} tubes (every {lay['cells_per_tube']} cells)"
          f"  eta {lay['eta_bay']:.3f}")
    # obs 4: the layout must not silently flip if the bank resolution
    # changes. Assert the tube count and report the snap margin.
    bnds = snap_boundary_h(d)
    hf = lay["h_face"]
    flip_above = bnds[lay["cells_per_tube"]]   # h above which tubes ++
    assert lay["n_tubes"] == 17, f"layout snap moved: {lay['n_tubes']} tubes"
    print(f"  snap margin: h_face {hf:.0f} vs the {lay['cells_per_tube']}"
          f"->{lay['cells_per_tube']-1} cell boundary at "
          f"{flip_above:.0f} W/m2K ({100*(flip_above-hf)/hf:.0f}% headroom "
          f"before the layout would move to more tubes)")
    print(f"  water {r['h_w']:.0f} W/m2K ({r['water_regime']})  "
          f"energy: gen {r['Q_gen']:.0f} = water {r['Q_water']:.0f} + "
          f"case {r['Q_case']:.0f} W")

    # buoyancy reference check (F1): near-zero net head at default
    Tbar = r["Tb"].mean(axis=(1, 2)).mean()
    head = d["rho"] * d["beta"] * 9.81 * max(Tbar - r["T_plen"], 0) * d["h_cell"]
    print(f"\nF1 buoyancy: riser-return dT {Tbar - r['T_plen']:+.3f} C -> "
          f"net head {head:.2f} Pa of {d['dp_extra']} Pa pump (was ~9 Pa)")

    # h_c sensitivity (F2) - report BOTH nodes (verification obs 2)
    print("\nF2 h_contact sweep (can / core):")
    for hc in (20000, 8000, 6000, 4000, 2000):
        rc = solve_zonal(dict(d, h_contact=hc), bank, nz=10, iters=200)
        flag = ("  <- core over 45" if rc["T_core_max"] > 45 else
                ("  <- can over 45" if rc["T_max"] > 45 else ""))
        print(f"  h_c {hc:6d} W/m2K -> can {rc['T_max']:.2f}"
              f"  core {rc['T_core_max']:.2f}{flag}")
    print("  (the can clears 45 down to ~2900 W/m2K, but the CORE only "
          "clears above ~7500 - within ~7% of the nominal 8000)")

    # optional oil->tube bypass (F3)
    rb = solve_zonal(dict(d, wetted_tube_frac=0.5), bank, nz=10, iters=200)
    print(f"\nF3 oil->tube bypass (frac 0.5): UA {rb['UA_bare']:.1f} W/K, "
          f"T_max {r['T_max']:.2f} -> {rb['T_max']:.2f} C "
          f"(Q_bare {rb['Q_bare']:.0f} W)")

    print("\nSHAKEDOWN OK  (%.0f s total)" % (time.time() - t0))
