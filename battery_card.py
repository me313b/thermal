"""Portable battery model card for the Immersion Pack Lab.

Turns a measurement campaign (e.g. the JP50 BioLogic HPPC set) into
ONE tidy CSV - the "card" - that travels by itself:

  * the app re-ingests it without the raw data (upload the card and
    every battery-aware feature uses it; upload nothing and a
    clearly-labelled simplified default is used instead);
  * COMSOL reads it directly: lines starting with '%' are
    COMSOL-native comments, so the file drops straight into an
    Interpolation function (Data source: File; pick the argument
    and value columns);
  * the analytical layer consumes the same numbers (q, R0, OCV).

Card schema v1 (columns): temp_C, rate_C, rep, soc_pct, dir,
i_pulse_A, dur_s, ocv_V, r0_ohm, r1_ohm, tau1_s, r2_ohm, tau2_s,
q_pulse_W. One row per analysed pulse; blanks where a quantity was
not extractable. Metadata and per-condition validation RMSE live in
the '%' header.
"""

import io
import os
import re
import zipfile
import numpy as np
import pandas as pd

CARD_MAGIC = "ipl-battery-card"
CARD_VERSION = 1
CARD_COLS = ["temp_C", "rate_C", "rep", "soc_pct", "dir",
             "i_pulse_A", "dur_s", "ocv_V", "r0_ohm", "r1_ohm",
             "tau1_s", "r2_ohm", "tau2_s", "q_pulse_W"]


# ------------------------------------------------------------------ #
#  Condition names                                                   #
# ------------------------------------------------------------------ #
def parse_condition(name):
    """'JP50_25deg_3C_2nd' -> (25.0, 3.0, 2). Returns None where a
    part is missing."""
    t = re.search(r"(-?\d+(?:\.\d+)?)\s*deg", name, re.I)
    r = re.search(r"_(\d+(?:\.\d+)?)\s*C(?:_|$)", name)
    rep = 2 if re.search(r"2nd|_2($|_)", name) else 1
    return (float(t.group(1)) if t else None,
            float(r.group(1)) if r else None, rep)


# ------------------------------------------------------------------ #
#  Loaders                                                           #
# ------------------------------------------------------------------ #
def load_mpr(path):
    """BioLogic .mpr -> DataFrame(time_s, current_A, voltage_V).
    EC-Lab logs discharge as negative; the HPPC pipeline
    auto-detects and flips."""
    from galvani import BioLogic
    a = BioLogic.MPRfile(path).data
    return pd.DataFrame({
        "time_s": np.asarray(a["time/s"], float),
        "current_A": np.asarray(a["I/mA"], float) / 1000.0,
        "voltage_V": np.asarray(a["Ewe/V"], float)})


def load_split_xlsx(paths):
    """Concatenate EC-Lab split exports (absolute time preserved)
    into one stream; keeps the measured temperature."""
    frames = []
    for p in paths:
        df = pd.read_excel(p)
        cols = {c.lower(): c for c in df.columns}
        v = cols.get("ecell/v") or cols.get("ewe/v")
        f = pd.DataFrame({
            "time_s": pd.to_numeric(df[cols["time/s"]],
                                    errors="coerce"),
            "current_A": pd.to_numeric(df[cols["i/ma"]],
                                       errors="coerce") / 1000.0,
            "voltage_V": pd.to_numeric(df[v], errors="coerce")})
        tc = cols.get("temperature/°c") or cols.get("temperature/c")
        if tc:
            f["temp_C"] = pd.to_numeric(df[tc], errors="coerce")
        frames.append(f)
    out = (pd.concat(frames).dropna(subset=["time_s"])
           .sort_values("time_s"))
    out = out[~out["time_s"].duplicated()].reset_index(drop=True)
    return out


def load_condition_dir(d):
    """Prefer the longest .mpr (the continuous stream); fall back to
    the split Excel files. Returns (df, source_str, temp_series)."""
    mprs = sorted((os.path.join(d, f) for f in os.listdir(d)
                   if f.lower().endswith(".mpr")),
                  key=os.path.getsize)
    temp = None
    xl = []
    for root, _, files in os.walk(d):
        xl += [os.path.join(root, f) for f in files
               if f.lower().endswith(".xlsx")]
    if xl:
        try:
            xdf = load_split_xlsx(sorted(xl))
            if "temp_C" in xdf:
                temp = float(np.nanmedian(xdf["temp_C"]))
        except Exception:
            xdf = None
    else:
        xdf = None
    if mprs:
        try:
            return load_mpr(mprs[-1]), \
                f"mpr:{os.path.basename(mprs[-1])}", temp
        except Exception:
            pass
    if xdf is not None and len(xdf) > 100:
        return xdf[["time_s", "current_A", "voltage_V"]], \
            f"xlsx-splits:{len(xl)}", temp
    raise ValueError(f"no readable BioLogic data in {d}")


def find_condition_dirs(root):
    """Directories that contain .mpr files or Split xlsx exports."""
    out = []
    for r, dirs, files in os.walk(root):
        if any(f.lower().endswith(".mpr") for f in files):
            out.append(r)
    return sorted(set(out))


def extract_zip(uploaded_bytes, workdir):
    os.makedirs(workdir, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(uploaded_bytes)) as z:
        z.extractall(workdir)
    return workdir


# ------------------------------------------------------------------ #
#  Campaign -> card                                                  #
# ------------------------------------------------------------------ #
def campaign_to_card(cond_results, cap_Ah, cell_name="cell"):
    """cond_results: list of dicts with keys name, temp_C, rate_C,
    rep, out (battery_hppc.run_pipeline result), source. Returns
    (card_df, meta)."""
    rows = []
    val = []
    for cr in cond_results:
        out = cr["out"]
        tab = out["tab"]
        d = out["d"]
        model = out["model"]
        for _, r in tab.iterrows():
            i_p = float(r["i"])
            soc = float(r["soc"])
            q_row = np.nan
            if np.isfinite(r.get("r0", np.nan)):
                ocv_here = float(model["ocv"](soc))
                q_row = abs(i_p) * abs(
                    ocv_here - (ocv_here - i_p * r["r0"]
                                - i_p * r.get("r1", 0.0)))
            rows.append(dict(
                temp_C=cr["temp_C"], rate_C=cr["rate_C"],
                rep=cr["rep"], soc_pct=round(soc, 2),
                dir="dis" if i_p > 0 else "chg",
                i_pulse_A=round(i_p, 4),
                dur_s=round(float(r["dur"]), 2),
                ocv_V=round(float(model["ocv"](soc)), 5),
                r0_ohm=r.get("r0", np.nan),
                r1_ohm=r.get("r1", np.nan),
                tau1_s=r.get("tau1", np.nan),
                r2_ohm=r.get("r2", np.nan),
                tau2_s=r.get("tau2", np.nan),
                q_pulse_W=q_row))
        val.append(f"{cr['name']}: rmse_mv={out['rmse_mv']:.2f} "
                   f"pulses={len(out['pulses'])} "
                   f"source={cr['source']}")
    card = pd.DataFrame(rows)
    for c in CARD_COLS:
        if c not in card:
            card[c] = np.nan
    card = card[CARD_COLS].sort_values(
        ["temp_C", "rate_C", "rep", "soc_pct"]).reset_index(
        drop=True)
    meta = dict(cell=cell_name, capacity_Ah=cap_Ah,
                validation=val)
    return card, meta


def card_write(card, meta):
    """Card CSV text with COMSOL-native '%' comment header."""
    hdr = [f"% {CARD_MAGIC} v{CARD_VERSION}",
           f"% cell: {meta.get('cell', '?')}",
           f"% capacity_Ah: {meta.get('capacity_Ah', '?')}",
           "% Columns: " + ", ".join(CARD_COLS),
           "% COMSOL: Global Definitions > Functions > "
           "Interpolation; Data source: File; this file loads "
           "directly ('%' lines are comments). Set the argument "
           "column to soc_pct and the value column to the "
           "quantity you need (filter rows to one temp/rate "
           "first, or use the app's 2-column export).",
           "% Analytical / app: heat per cell at sustained "
           "current I is q = I^2*(r0_ohm + r1_ohm); OCV and "
           "R interpolate in soc_pct."]
    for v in meta.get("validation", []):
        hdr.append(f"% validated: {v}")
    return "\n".join(hdr) + "\n" + card.to_csv(index=False)


def card_read(text_or_buf):
    """Round-trip reader; returns (card_df, meta)."""
    if hasattr(text_or_buf, "read"):
        text = text_or_buf.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8", errors="replace")
    else:
        text = text_or_buf
    lines = text.splitlines()
    meta = dict(validation=[])
    body = []
    magic_ok = False
    for ln in lines:
        if ln.startswith("%"):
            c = ln[1:].strip()
            if c.startswith(CARD_MAGIC):
                magic_ok = True
            elif c.startswith("cell:"):
                meta["cell"] = c.split(":", 1)[1].strip()
            elif c.startswith("capacity_Ah:"):
                try:
                    meta["capacity_Ah"] = float(
                        c.split(":", 1)[1])
                except ValueError:
                    pass
            elif c.startswith("validated:"):
                meta["validation"].append(
                    c.split(":", 1)[1].strip())
        else:
            body.append(ln)
    if not magic_ok:
        raise ValueError("not an ipl-battery-card CSV")
    card = pd.read_csv(io.StringIO("\n".join(body)))
    missing = [c for c in ("soc_pct", "r0_ohm", "ocv_V")
               if c not in card.columns]
    if missing:
        raise ValueError(f"card missing columns {missing}")
    return card, meta


# ------------------------------------------------------------------ #
#  Card -> usable model                                              #
# ------------------------------------------------------------------ #
def conditions_in(card):
    key = card[["temp_C", "rate_C", "rep"]].drop_duplicates()
    return [tuple(x) for x in key.itertuples(index=False)]


def model_from_card(card, temp_C=None, rate_C=None, rep=None,
                    direction="dis"):
    """Interpolants for one condition (nearest available if the
    requested one is absent)."""
    sub = card.copy()
    if direction in ("dis", "chg") and "dir" in sub:
        s2 = sub[sub["dir"] == direction]
        if len(s2) >= 3:
            sub = s2
    conds = sub[["temp_C", "rate_C"]].drop_duplicates()
    if temp_C is not None and len(conds) > 1:
        dt = (conds["temp_C"] - temp_C).abs() + \
             0.01 * (conds["rate_C"] - (rate_C or 0)).abs()
        pick = conds.iloc[int(np.argmin(dt.values))]
        sub = sub[(sub["temp_C"] == pick["temp_C"])
                  & (sub["rate_C"] == pick["rate_C"])]
    if rep is not None and "rep" in sub and \
            (sub["rep"] == rep).sum() >= 3:
        sub = sub[sub["rep"] == rep]
    sub = sub.sort_values("soc_pct")

    def mk(col, default):
        s3 = sub.dropna(subset=[col]) if col in sub else \
            pd.DataFrame()
        if len(s3) == 0:
            return lambda s: np.full_like(
                np.asarray(s, float), default)
        xs = np.asarray(s3["soc_pct"], float)
        ys = np.asarray(s3[col], float)
        return lambda s: np.interp(np.asarray(s, float), xs, ys)

    used = (float(sub["temp_C"].iloc[0]),
            float(sub["rate_C"].iloc[0])) if len(sub) else \
        (np.nan, np.nan)
    return dict(order=1, ocv=mk("ocv_V", 3.7),
                r0=mk("r0_ohm", 0.020), r1=mk("r1_ohm", 0.010),
                tau1=mk("tau1_s", 30.0),
                cond_used=used, n_rows=len(sub))


def default_model():
    """The clearly-labelled fallback when no card or data is
    loaded: constant 20 mOhm ohmic + 10 mOhm polarisation, linear
    3.0-4.2 V OCV."""
    return dict(order=1,
                ocv=lambda s: 3.0 + 1.2 * np.asarray(s, float)
                / 100.0,
                r0=lambda s: np.full_like(np.asarray(s, float),
                                          0.020),
                r1=lambda s: np.full_like(np.asarray(s, float),
                                          0.010),
                tau1=lambda s: np.full_like(np.asarray(s, float),
                                            30.0),
                cond_used=(np.nan, np.nan), n_rows=0)


def q_steady(model, i_A, soc_pct=50.0):
    """Sustained-current heat per cell: I^2 (R0 + R1) at SOC."""
    i = float(i_A)
    return i * i * (float(model["r0"](soc_pct))
                    + float(model["r1"](soc_pct)))


def comsol_two_col(card, quantity, temp_C, rate_C,
                   direction="dis"):
    """A minimal 2-column (soc_pct, value) CSV for a COMSOL
    Interpolation function, '%' header included."""
    sub = card[(card["temp_C"] == temp_C)
               & (card["rate_C"] == rate_C)]
    if "dir" in sub and direction in ("dis", "chg"):
        s2 = sub[sub["dir"] == direction]
        if len(s2) >= 2:
            sub = s2
    sub = sub.dropna(subset=[quantity]).sort_values("soc_pct")
    hdr = (f"% {CARD_MAGIC} 2-col export: {quantity} at "
           f"{temp_C:g} degC, {rate_C:g}C, {direction}\n"
           "% COMSOL: Interpolation function, Data source File; "
           "argument = column 1 (soc_pct), value = column 2.\n")
    return hdr + sub[["soc_pct", quantity]].to_csv(index=False)
