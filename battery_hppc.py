"""HPPC battery model extraction for the Immersion Pack Lab.

Energy law used as the built-in check: the integral of
I*(OCV - V) equals the integral of I^2*R0 + V1^2/R1 (+V2^2/R2)
over any profile that starts and ends at rest (to the quadrature
accuracy of the sampling - trapezoid across instantaneous current
steps leaves a few-per-mille residual at 1 Hz); pointwise the two
differ by the RC capacitor's storage rate. The instantaneous
dissipation series q is the thermally correct one.

From a Hybrid Pulse Power Characterization test (rest - pulse -
rest at stepped SOC) this module extracts an equivalent-circuit
model: OCV(SOC), ohmic R0(SOC) from the instantaneous voltage step,
and one or two RC polarisation branches (R1, tau1 [, R2, tau2])
from the relaxation tails - then PROVES the result by re-simulating
the entire measured profile and reporting the millivolt RMS error,
plus the heat identity I^2*R0 + V1^2/R1 (+V2^2/R2) == I*(OCV - V),
which is what the thermal app actually consumes.

Sign convention inside this module: discharge current POSITIVE.
Pure numpy/pandas; the RC fits use a deterministic tau-grid with
linear least squares for the amplitudes (no fragile optimiser).
"""

import numpy as np
import pandas as pd


# ------------------------------------------------------------------ #
#  Parsing                                                           #
# ------------------------------------------------------------------ #
def map_columns(df):
    """Find time / current / voltage (and optional soc, temp)
    columns by fuzzy name match. Returns dict of column names."""
    cols = {c.lower().strip(): c for c in df.columns}
    out = {}

    def find(keys):
        for k in keys:
            for lc, orig in cols.items():
                if k in lc:
                    return orig
        return None

    out["t"] = find(["time", "t_s", "t (", "sec"])
    out["i"] = find(["current", "i_a", "i (", "amp"])
    out["v"] = find(["voltage", "v_", "v (", "volt"])
    out["soc"] = find(["soc"])
    out["temp"] = find(["temp"])
    if out["t"] is None and "t" in cols:
        out["t"] = cols["t"]
    if out["i"] is None and "i" in cols:
        out["i"] = cols["i"]
    if out["v"] is None and "v" in cols:
        out["v"] = cols["v"]
    missing = [k for k in ("t", "i", "v") if out[k] is None]
    if missing:
        raise ValueError(
            "could not identify column(s) %s - expected names "
            "containing time / current / voltage" % missing)
    return out


def detect_sign(t, i, v):
    """+1 if discharge is already positive in the data, -1 if it must
    be flipped. Uses the largest current step: real cells drop in
    voltage when discharge current is applied."""
    di = np.diff(i)
    k = int(np.argmax(np.abs(di)))
    dv = v[k + 1] - v[k]
    if di[k] == 0:
        return 1
    return 1 if dv * di[k] < 0 else -1


def prepare(df, cap_Ah, soc0=100.0, time_unit="auto",
            dis_positive="auto"):
    """Return dict with t [s], i_dis [A, discharge +], v [V],
    soc [%], plus bookkeeping."""
    cm = map_columns(df)
    t = np.asarray(df[cm["t"]], float)
    i = np.asarray(df[cm["i"]], float)
    v = np.asarray(df[cm["v"]], float)
    ok = np.isfinite(t) & np.isfinite(i) & np.isfinite(v)
    t, i, v = t[ok], i[ok], v[ok]
    order = np.argsort(t)
    t, i, v = t[order], i[order], v[order]
    if time_unit == "auto":
        dt_med = np.median(np.diff(t)) if len(t) > 3 else 1.0
        time_unit = "ms" if dt_med > 20 else "s"
    if time_unit == "ms":
        t = t / 1000.0
    t = t - t[0]
    if dis_positive == "auto":
        sgn = detect_sign(t, i, v)
    else:
        sgn = 1 if dis_positive else -1
    i_dis = sgn * i
    if cm["soc"] is not None:
        soc = np.asarray(df.loc[ok, cm["soc"]], float)[order]
        if np.nanmax(soc) <= 1.5:
            soc = soc * 100.0
    else:
        ah = np.concatenate(
            [[0.0], np.cumsum(0.5 * (i_dis[1:] + i_dis[:-1])
                              * np.diff(t))]) / 3600.0
        soc = soc0 - 100.0 * ah / cap_Ah
    temp = (np.asarray(df.loc[ok, cm["temp"]], float)[order]
            if cm["temp"] is not None else None)
    return dict(t=t, i=i_dis, v=v, soc=soc, temp=temp,
                sign_flipped=(sgn == -1), time_unit=time_unit,
                cap_Ah=cap_Ah)


# ------------------------------------------------------------------ #
#  Pulse segmentation                                                #
# ------------------------------------------------------------------ #
def find_pulses(d, i_thresh=None, min_dur=1.0, min_pre_rest=20.0):
    """Constant-current active segments bounded by rests. Returns a
    list of dicts with sample indices and medians."""
    t, i = d["t"], d["i"]
    if i_thresh is None:
        i_thresh = 0.05 * np.nanmax(np.abs(i))
    act = np.abs(i) > i_thresh
    edges = np.flatnonzero(np.diff(act.astype(int)))
    starts = edges[act[edges + 1]] + 1
    ends = edges[~act[edges + 1]] + 1
    if act[0]:
        starts = np.r_[0, starts]
    if act[-1]:
        ends = np.r_[ends, len(i)]
    pulses = []
    for s, e in zip(starts, ends):
        if t[e - 1] - t[s] < min_dur:
            continue
        pre = t[s] - (t[ends[ends <= s][-1] - 1]
                      if (ends <= s).any() else 0.0)
        i_med = float(np.median(i[s:e]))
        # constant-current check: reject ramps/drive segments
        if np.std(i[s:e]) > 0.1 * abs(i_med) + 1e-9:
            continue
        pulses.append(dict(s=int(s), e=int(e), i=i_med,
                           t0=float(t[s]), t1=float(t[e - 1]),
                           dur=float(t[e - 1] - t[s]),
                           pre_rest=float(pre),
                           soc=float(d["soc"][s])))
    return [p for p in pulses if p["pre_rest"] >= min_pre_rest or
            p["dur"] >= min_dur], float(i_thresh)


# ------------------------------------------------------------------ #
#  Parameter extraction                                              #
# ------------------------------------------------------------------ #
def _fit_relax(tr, vr, order=1):
    """Fit V = Vinf - sum_i A_i exp(-tr/tau_i) with a deterministic
    log-spaced tau grid and linear LSQ amplitudes."""
    taus = np.geomspace(0.5, max(200.0, tr[-1]), 60)
    best = None
    if order == 1:
        for tau in taus:
            X = np.column_stack([np.ones_like(tr),
                                 -np.exp(-tr / tau)])
            coef, res, _, _ = np.linalg.lstsq(X, vr, rcond=None)
            sse = float(res[0]) if len(res) else float(
                np.sum((vr - X @ coef) ** 2))
            if best is None or sse < best[0]:
                best = (sse, coef[0], [(coef[1], tau)])
    else:
        tg = np.geomspace(0.5, max(200.0, tr[-1]), 22)
        for a_ in range(len(tg)):
            for b_ in range(a_ + 1, len(tg)):
                t1, t2 = tg[a_], tg[b_]
                if t2 < 2.5 * t1:
                    continue
                X = np.column_stack([np.ones_like(tr),
                                     -np.exp(-tr / t1),
                                     -np.exp(-tr / t2)])
                coef, res, _, _ = np.linalg.lstsq(X, vr,
                                                  rcond=None)
                sse = float(res[0]) if len(res) else float(
                    np.sum((vr - X @ coef) ** 2))
                if best is None or sse < best[0]:
                    best = (sse, coef[0],
                            [(coef[1], t1), (coef[2], t2)])
    sse, vinf, branches = best
    return vinf, branches, sse


def extract(d, pulses, order=1, relax_window=600.0):
    """Per suitable pulse: R0 from the step edges, RC branches from
    the post-pulse relaxation, OCV points from the rests."""
    t, i, v, soc = d["t"], d["i"], d["v"], d["soc"]
    rows, ocv_pts = [], []
    for k, p in enumerate(pulses):
        s, e = p["s"], p["e"]
        r0s = []
        if s > 0 and abs(i[s] - i[s - 1]) > 1e-6:
            r0s.append((v[s - 1] - v[s]) / (i[s] - i[s - 1]))
        nxt = e if e < len(i) else e - 1
        if nxt < len(i) and abs(i[nxt] - i[e - 1]) > 1e-6:
            r0s.append((v[nxt] - v[e - 1]) /
                       (i[e - 1] - i[nxt]))
        r0 = float(np.mean(r0s)) if r0s else np.nan
        rest_end = pulses[k + 1]["s"] if k + 1 < len(pulses) \
            else len(t)
        tr = t[e:rest_end] - t[e - 1]
        vr = v[e:rest_end]
        keep = tr <= relax_window
        row = dict(soc=p["soc"], i=p["i"], dur=p["dur"],
                   r0=abs(r0), t0=p["t0"])
        if keep.sum() >= 8 and abs(p["i"]) > 1e-6:
            vinf, br, sse = _fit_relax(tr[keep], vr[keep], order)
            fac_ok = True
            for j, (A, tau) in enumerate(br, start=1):
                fac = 1.0 - np.exp(-p["dur"] / tau)
                if fac < 1e-3:
                    fac_ok = False
                    continue
                R = A / (p["i"] * fac)
                row[f"r{j}"] = float(R)
                row[f"tau{j}"] = float(tau)
            if fac_ok or order == 1:
                row["v_inf"] = float(vinf)
                ocv_pts.append((float(soc[min(rest_end - 1,
                                              len(soc) - 1)]),
                                float(vinf)))
        if s > 3 and p["pre_rest"] >= 60.0:
            ocv_pts.append((float(soc[s - 1]), float(v[s - 1])))
        rows.append(row)
    tab = pd.DataFrame(rows)
    # discard unphysical branches (negative R) rather than hide them
    for c in [c for c in tab.columns if c.startswith("r")]:
        tab.loc[tab[c] <= 0, c] = np.nan
    ocv = pd.DataFrame(ocv_pts, columns=["soc", "ocv"])
    ocv = (ocv.groupby(ocv["soc"].round(1)).mean()
              .reset_index(drop=True).sort_values("soc"))
    return tab, ocv


# ------------------------------------------------------------------ #
#  Model tables, simulation, heat                                    #
# ------------------------------------------------------------------ #
def build_model(tab, ocv, order=1):
    """SOC-interpolated parameter functions (edge-clamped)."""
    def mk(col, default):
        sub = tab.dropna(subset=[col]) if col in tab else \
            pd.DataFrame()
        if len(sub) == 0:
            return lambda s: np.full_like(np.asarray(s, float),
                                          default)
        xs = np.asarray(sub["soc"], float)
        ys = np.asarray(sub[col], float)
        o = np.argsort(xs)
        xs, ys = xs[o], ys[o]
        return lambda s: np.interp(np.asarray(s, float), xs, ys)

    ox = np.asarray(ocv["soc"], float)
    oy = np.asarray(ocv["ocv"], float)
    model = dict(order=order,
                 ocv=lambda s: np.interp(np.asarray(s, float),
                                         ox, oy),
                 r0=mk("r0", 0.02), r1=mk("r1", 0.01),
                 tau1=mk("tau1", 30.0))
    if order == 2:
        model["r2"] = mk("r2", 0.005)
        model["tau2"] = mk("tau2", 120.0)
    return model


def simulate(d, model):
    """Exact piecewise-constant-current update of the RC states.
    Returns v_model, per-branch overpotentials, heat series and the
    identity residual."""
    t, i, soc = d["t"], d["i"], d["soc"]
    n = len(t)
    v1 = np.zeros(n)
    v2 = np.zeros(n) if model["order"] == 2 else None
    r0 = model["r0"](soc)
    r1 = model["r1"](soc)
    tau1 = np.maximum(model["tau1"](soc), 1e-3)
    if v2 is not None:
        r2 = model["r2"](soc)
        tau2 = np.maximum(model["tau2"](soc), 1e-3)
    for k in range(1, n):
        dt = t[k] - t[k - 1]
        a = np.exp(-dt / tau1[k])
        v1[k] = v1[k - 1] * a + r1[k] * i[k] * (1 - a)
        if v2 is not None:
            b = np.exp(-dt / tau2[k])
            v2[k] = v2[k - 1] * b + r2[k] * i[k] * (1 - b)
    v_model = model["ocv"](soc) - i * r0 - v1 - \
        (v2 if v2 is not None else 0.0)
    q = i ** 2 * r0 + np.where(r1 > 0, v1 ** 2 / r1, 0.0)
    if v2 is not None:
        q = q + np.where(r2 > 0, v2 ** 2 / r2, 0.0)
    q_alt = i * (model["ocv"](soc) - v_model)
    # The correct bookkeeping: q_alt - q = d/dt(0.5*C*V_rc^2), the
    # capacitor storage rate - so the two INTEGRALS must agree over
    # a profile that starts and ends at rest, while pointwise they
    # differ during transients. q is the instantaneous dissipation
    # (the thermally correct series).
    E_q = float(np.trapezoid(q, t))
    E_alt = float(np.trapezoid(q_alt, t))
    return dict(v_model=v_model, v1=v1, v2=v2, q=q, q_alt=q_alt,
                E_q=E_q, E_alt=E_alt,
                ident_rel=float(abs(E_alt - E_q) /
                                max(abs(E_q), 1e-12)))


def heat_summary(d, sim, i_thresh):
    t, i, q = d["t"], d["i"], sim["q"]
    T = t[-1] - t[0]
    q_mean = float(np.trapezoid(q, t) / max(T, 1e-9))
    act = np.abs(i) > i_thresh
    q_act = float(np.trapezoid(q[act], t[act]) /
                  max(t[act][-1] - t[act][0], 1e-9)) \
        if act.sum() > 3 else q_mean
    return dict(q_mean_profile=q_mean, q_mean_active=q_act,
                q_peak=float(np.nanmax(q)))


def rmse_mv(d, sim, i_thresh=None):
    err = d["v"] - sim["v_model"]
    return float(np.sqrt(np.nanmean(err ** 2)) * 1000.0)


# ------------------------------------------------------------------ #
#  Synthetic HPPC (sample data + ground truth for the tests)         #
# ------------------------------------------------------------------ #
def synth_truth():
    return dict(
        cap_Ah=4.5,
        ocv=lambda s: 2.90 + 1.30 * np.asarray(s / 100.0) ** 0.55
        + 0.05 * (s / 100.0),
        r0=lambda s: 0.018 + 0.006 * (1 - s / 100.0) ** 2,
        r1=lambda s: 0.012 + 0.004 * (1 - s / 100.0),
        tau1=lambda s: 25.0 + 0.0 * np.asarray(s))


def make_synth_hppc(noise_mv=0.0, seed=0):
    """Realistic 21700-style HPPC: at each 10% SOC step a 2C
    discharge pulse and a 1.5C charge pulse (10 s each) with rests,
    then a C/2 leg to the next SOC. Voltage generated by the exact
    same ECM update the extractor must recover."""
    tr = synth_truth()
    cap = tr["cap_Ah"]
    seq = []
    for _ in range(9):
        seq += [(600.0, 0.0), (10.0, 2.0 * cap), (60.0, 0.0),
                (10.0, -1.5 * cap), (60.0, 0.0),
                (720.0, 0.5 * cap)]
    seq += [(600.0, 0.0)]
    t_list, i_list = [0.0], [0.0]
    for dur, cur in seq:
        n = int(dur)
        t0 = t_list[-1]
        t_list += [t0 + k for k in range(1, n + 1)]
        i_list += [cur] * n
    t = np.asarray(t_list)
    i = np.asarray(i_list)
    ah = np.concatenate([[0.0],
                         np.cumsum(0.5 * (i[1:] + i[:-1])
                                   * np.diff(t))]) / 3600.0
    soc = 100.0 - 100.0 * ah / cap
    v1 = np.zeros_like(t)
    for k in range(1, len(t)):
        dt = t[k] - t[k - 1]
        tau = float(tr["tau1"](soc[k]))
        a = np.exp(-dt / tau)
        v1[k] = v1[k - 1] * a + tr["r1"](soc[k]) * i[k] * (1 - a)
    v = tr["ocv"](soc) - i * tr["r0"](soc) - v1
    if noise_mv > 0:
        rng = np.random.default_rng(seed)
        v = v + rng.normal(0, noise_mv / 1000.0, size=len(v))
    return pd.DataFrame({"time_s": t, "current_A": i,
                         "voltage_V": np.round(v, 6)})


def run_pipeline(df, cap_Ah, soc0=100.0, order=1,
                 time_unit="auto", dis_positive="auto",
                 min_pre_rest=20.0):
    """The whole chain; returns everything the tab needs."""
    d = prepare(df, cap_Ah, soc0, time_unit, dis_positive)
    pulses, ithr = find_pulses(d, min_pre_rest=min_pre_rest)
    tab, ocv = extract(d, pulses, order=order)
    if len(ocv) < 2:
        raise ValueError("fewer than two OCV rest points found - "
                         "check rests, threshold and capacity")
    model = build_model(tab, ocv, order=order)
    sim = simulate(d, model)
    heat = heat_summary(d, sim, ithr)
    return dict(d=d, pulses=pulses, tab=tab, ocv=ocv, model=model,
                sim=sim, heat=heat, rmse_mv=rmse_mv(d, sim),
                i_thresh=ithr)
