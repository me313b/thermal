"""zonal.py - plate-channel zonal thermal-hydraulic network + Monte Carlo.

Physics: every inter-row gap carries a water-rooted aluminium plate with
an oil slot on each face; oil flows upward; each slot column obeys the
exact developing-flow kernel from fea4_channel.graetz_kernel on the true
lens cross-section:
    [q'_c; q'_p] = a(z*, s) . [T_cell - Tb; T_plateface - Tb]
Fin knockdown eta_bay (product rule, roots at the derived tube pitch
P = 2/m, tubes at mid-height) scales the plate-side kernel entries.
Hydraulics: laminar fRe(s) at mu(T) with buoyancy assist against a
shared manifold head -> the cubic-law tolerance sensitivity. Water is
marched tube-by-tube across the plates through contact + wall + film.
Validations: energy closure, uniform-input collapse, z-independence.
"""
import math
import numpy as np

from fea4_channel import graetz_kernel

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
                xk=xk, eta_bay=eta, cells_per_tube=n_per)


def solve_zonal(d, bank, s_ch=None, contact=None, heat_map=None,
                nz=10, iters=40, relax=0.6, tol=0.002,
                init=None):
    """Plate-channel zonal network. Every subsystem is closed-form or
    contractive: exponential (exact) z-march per channel; per-cell
    Newton with the march's own Jacobian; per-bay plate-face solved in
    closed form against the root chain (contact + wall + water film);
    water marched tube-by-tube; recirculation plenum solved from the
    linear exit map T_exit = A + P.T_plen. Energy closure is the
    convergence monitor and is returned."""
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

    lay = derived_layout(d, h_face=150.0)
    n_t = lay["n_tubes"]
    eta = lay["eta_bay"]
    if contact is None:
        contact = np.full((n_gap, n_t), d["plate_contact"])
    col_bay = np.minimum((np.arange(nc) * d["pitch"]
                          // lay["P_snap"]).astype(int), n_t - 1)
    ncol_bay = np.bincount(col_bay, minlength=n_t).astype(float)

    d_i = d["tube_od"] - 2 * d["tube_wall"]
    mdot_w = d["flow_lpm"] / 60 * WATER["rho"] / 1000 / n_t
    Re_w = 4 * mdot_w / (math.pi * WATER["mu"] * d_i)
    Pr_w = WATER["mu"] * WATER["cp"] / WATER["k"]
    L_t = max(n_gap * d["pitch"], 0.05)
    if Re_w < 2300:
        gz = (d_i / L_t) * Re_w * Pr_w
        Nu_w = 3.66 + 0.0668 * gz / (1 + 0.04 * gz ** (2 / 3))
    else:
        f = (0.790 * math.log(max(Re_w, 3000)) - 1.64) ** -2
        Nu_w = ((f / 8) * (Re_w - 1000) * Pr_w
                / (1 + 12.7 * math.sqrt(f / 8)
                   * (Pr_w ** (2 / 3) - 1)))
    h_w = Nu_w * WATER["k"] / d_i
    per_len = d["pitch"]
    R_film = 1 / (h_w * math.pi * d_i * per_len)
    R_wl = math.log(d["tube_od"] / d_i) / (2 * math.pi
                                           * d["k_tube"] * per_len)
    A_ct = math.pi * d["tube_od"] * d["plate_t"] * 6.0
    HC = 8000.0
    Rr = (1.0 / (HC * A_ct * contact) + R_wl + R_film)   # (n_gap,n_t)

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
            Tf = init["Tf"].copy()
            T_wat = init["T_wat"].copy()
            if init["Tb"].shape == Tb.shape:
                Tb = init["Tb"].copy()
            mdot_col = init["mdot_col"].copy()
        except (KeyError, AttributeError):
            pass
    closure = 1.0

    for it in range(iters):
        # ---- 1. hydraulics (clamped, under-relaxed) ----
        Tbar = np.clip(Tb.mean(axis=(1, 2)), d["T_in"] - 5,
                       d["T_in"] + 60)
        mu = rho * nu(Tbar)
        dp = d["dp_extra"] + rho * beta * 9.81 * np.clip(
            Tbar - d["T_in"], 0.0, 45.0) * H
        ubar = dp * 2 * Dh ** 2 / (fRe * mu * H)
        m_new = np.maximum(rho * ubar * Acs, 1e-7)
        mdot_col = (m_new if mdot_col is None
                    else mdot_col + relax * (m_new - mdot_col))

        # ---- 2. channel march (exact segment update) ----
        alpha = k_o / (rho * cp)
        RePr = np.maximum(mdot_col / (rho * Acs), 1e-9) * Dh / alpha
        Tb_new = np.empty_like(Tb)
        Pch = np.empty((n_ch, nc))
        qhat_c = np.zeros((nr, nc))
        Jc = np.zeros((nr, nc))
        W_bay = np.zeros((n_gap, n_t))
        S_bay = np.zeros((n_gap, n_t))
        Q_case_new = 0.0
        Qc_tot = 0.0
        for i in range(n_ch):
            r, pidx = side_row[i], side_plate[i]
            zst = z / (Dh[i] * max(RePr[i], 1e-9))
            a = bank.a_of(s_ch[i], zst).copy()       # (nz,2,2)
            if pidx >= 0:
                a[:, 0, 1] *= eta; a[:, 1, 0] *= eta
                a[:, 1, 1] *= eta
                Tp = Tf[pidx, col_bay]               # (nc,)
            else:
                Gf = float(a[:, 1, 1].sum()) * dz
                f_e = (d["h_ext"] * A_col_edge
                       / max(d["h_ext"] * A_col_edge + Gf, 1e-9))
                a[:, 0, 1] *= f_e; a[:, 1, 0] *= f_e
                a[:, 1, 1] *= f_e
                Tp = np.full(nc, d["T_amb"])
            Tbi = Tb[i].copy()
            Tbi[:, 0] = T_plen
            mcp = float(mdot_col[i] * cp)
            Tc_r = T_cell[r]
            Pdecay = np.ones(nc)
            for kz in range(nz):
                a00, a01 = a[kz, 0, 0], a[kz, 0, 1]
                a10, a11 = a[kz, 1, 0], a[kz, 1, 1]
                Gc_row, Gp_row = a00 + a10, a01 + a11
                G = max(Gc_row + Gp_row, 1e-12)
                Tw_eff = (Gc_row * Tc_r + Gp_row * Tp) / G
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
                    # S collects q_p at the CURRENT Tf; convert to the
                    # Tf-independent part: q_p = S0 + W.Tf with
                    # S0 = S - W.Tf_current, handled after the sweep.
                else:
                    Q_case_new += float(np.sum(-qp))
            Tb_new[i] = Tbi
            Pch[i] = Pdecay
        Tb = Tb + relax * (Tb_new - Tb)

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
        # q_root(Tf) = -(S0 + W Tf) = (Tf - Tw)/Rr
        S0 = S_bay - W_bay * Tf
        Tw_bay = T_wat.T                              # (n_gap, n_t)
        Tf_new = (Tw_bay / Rr - S0) / (W_bay + 1.0 / Rr)
        Tf += 0.7 * (np.clip(Tf_new, d["T_in"] - 5, 120.0) - Tf)
        q_root = (Tf - Tw_bay) / Rr

        # ---- 6. water march ----
        for t in range(n_t):
            Tw = d["T_in"]
            for p in range(n_gap):
                T_wat[t, p] = Tw
                Tw += q_root[p, t] / (mdot_w * WATER["cp"])
                Tw = min(max(Tw, d["T_in"] - 5), 95.0)

        Qg = float(np.sum(q_cell0 * hm * np.exp(
            -d["k_dcir"] * (np.clip(T_cell, -10, 130) - 25))))
        Q_out = float(np.sum(q_root)) + Q_case_new
        closure = (Qg - Q_out) / max(Qg, 1e-9)
        if it > 10 and abs(closure) < tol:
            break

    Q_gen = float(np.sum(q_cell0 * hm * np.exp(
        -d["k_dcir"] * (np.clip(T_cell, -10, 130) - 25))))
    Q_water = float(np.sum(q_root))
    Q_case = Q_case_new
    P_pump = float(np.sum(mdot_col * nc / rho * d["dp_extra"]) / 0.35)
    return dict(T_cell=T_cell, Tb=Tb, Tf=Tf, T_wat=T_wat,
                T_plen=T_plen, P_pump=P_pump, mdot_col=mdot_col,
                ubar=ubar, s_ch=s_ch, lay=lay, q_root=q_root,
                Q_gen=Q_gen, Q_water=Q_water, Q_case=Q_case,
                closure=(Q_gen - Q_water - Q_case) / max(Q_gen, 1e-9),
                T_max=float(T_cell.max()),
                T_mean=float(T_cell.mean()),
                spread=float(T_cell.max() - T_cell.min()),
                iters_used=it + 1, zs=z)


def monte_carlo(d, bank, M=150, sigma_s=0.2e-3, sigma_c=0.08,
                heat_map=None, nz=10, iters=35, seed=7,
                progress=None):
    """Tolerance Monte Carlo: per-channel slot widths and per-root
    contact quality sampled from truncated normals; every sample is a
    full network solve warm-started from the converged nominal case."""
    rng = np.random.default_rng(seed)
    n_ch = 2 * (d["n_rows"] - 1) + 2
    lay = derived_layout(d, 150.0)
    base = solve_zonal(d, bank, heat_map=heat_map, nz=nz, iters=240,
                       tol=0.002)
    Tmax = np.empty(M); spread = np.empty(M); clos = np.empty(M)
    for m in range(M):
        s = np.clip(rng.normal(d["s_nom"], sigma_s, n_ch),
                    0.55 * d["s_nom"], 1.7 * d["s_nom"])
        c = np.clip(rng.normal(d["plate_contact"], sigma_c,
                               (d["n_rows"] - 1, lay["n_tubes"])),
                    0.15, 1.0)
        r = solve_zonal(d, bank, s_ch=s, contact=c,
                        heat_map=heat_map, nz=nz, iters=iters,
                        tol=0.004, init=base)
        Tmax[m] = r["T_max"]; spread[m] = r["spread"]
        clos[m] = abs(r["closure"])
        if progress is not None:
            progress((m + 1) / M)
    return dict(Tmax=Tmax, spread=spread, base=base,
                closure_worst=float(clos.max()),
                p_exceed=float(np.mean(Tmax > d["T_limit"])))


def default_d(nr=33, nc=33, C=2.0):
    return dict(n_rows=nr, n_cols=nc, pitch=0.0215, d_cell=0.021,
                h_cell=0.070, cap_Ah=5.0, r_dc=25.0, k_dcir=0.012,
                C=C, plate_t=0.0015, plate_contact=0.8, s_nom=0.002,
                T_in=20.0, flow_lpm=10.0, tube_od=0.010,
                tube_wall=0.0008, k_tube=385.0, T_amb=25.0, h_ext=5.0,
                A_case=1.2, nu25=9e-6, B=3200.0, rho=920.0, cp=2000.0,
                k_oil=0.13, beta=7.5e-4, dp_extra=25.0, T_limit=45.0)
