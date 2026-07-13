# Immersion Pack Lab — Red-Team Handover

**Version:** v9.4  **Date:** 13 Jul 2026
*(v9.4 incorporates the first red-team round in full; see REDTEAM_RESPONSE.md for the point-by-point closure. Numbers below are the corrected v9.4 operating point.)*
**Bundle:** `immersion_pack_lab_v9_3_full.zip`
**Primary targets for review:** the **Zones** tab (zonal plate-channel
thermal-hydraulic network + Monte Carlo) and its engine `fea4_channel.py`
/ `zonal.py`. These are the newest and least-exercised subsystems and the
reason for this handover.

---

## 1. What I am asking you to do

Go through the model in detail and try to break it. Specifically:

1. **Find where the physics is wrong or oversold.** Every equation and
   constant is documented in Section 5; check the derivations, the
   correlations, and whether the assumptions in Section 7 actually hold
   for a static/stirred immersion pack.
2. **Attack the numbers.** Section 6 lists every validation gate and
   cross-check. Tell me which ones are real evidence and which are
   circular (a model agreeing with a close relative of itself).
3. **Probe the soft spots I already know about** (Section 8) and, more
   usefully, the ones I have not seen.
4. **Judge the explanations.** The app now shows the governing equation
   beside each result (Zones "Methods" panel; Learn panels; Validate
   FEA4). Tell me where the story is misleading, hand-waved, or wrong.

Good feedback names a specific equation, line, or number, states the
discrepancy, and ideally gives a reference or a reproducing input
(Section 10).

---

## 2. The request that produced this build (context for reviewers)

The supervisor's brief was: *"Make a comprehensive and highly technically
detailed model that addresses the zonal and Monte Carlo case. Make it
highly accurate, state of the art, and a strong comprehensive design
without any simplification. For the explanations, discuss how the results
are achieved and show the equations and the physics beside them."*

So the intent is deliberately **maximal rigour, minimal shortcut**. If a
step looks like a shortcut, that is exactly the thing to challenge. The
one honesty point I want stress-tested: an early conduction-only study of
the slot returned an arc effectiveness of about 0.13 (cell-to-plate
conduction across the gap is *weaker* than the thermosiphon it was meant
to replace), which redirected the whole architecture to **advective
primary cooling** — oil carries the heat, plates are the secondary sink
and structure. Everything downstream depends on that being the right call.

---

## 3. How to run it (about 5 minutes)

```bash
pip install -r requirements.txt
streamlit run app.py            # opens the 12-tab app
SMOKE=1 python app.py           # fast self-test of solver internals
```

Headless render test (exercises every tab, including the zonal solve and
all equation panels):

```python
from streamlit.testing.v1 import AppTest
at = AppTest.from_file("app.py", default_timeout=800)
at.run()
assert len(at.exception) == 0
```

Where the new work lives:
- **Zones** tab: interactive model, the "Methods, physics and how every
  number is reached" panel (equations + schematics + live validation),
  and the tolerance Monte Carlo.
- **Validate** tab: **FEA4** (duct-kernel gates) now sits with FEA1-3.
- `fea4_channel.py` and `zonal.py` both run standalone (`python
  fea4_channel.py`, and the shakedowns in `zonal.py`) and print their own
  gates.

---

## 4. Architecture and file map

| File | Role |
|---|---|
| `app.py` | 12-tab Streamlit app; lumped two-node solver, transient, benchmark, cockpit wiring, all UI and reports |
| `zonal.py` | Plate-channel zonal network solver + Monte Carlo (Section 5 D-H) |
| `fea4_channel.py` | Duct-coefficient engine: velocity Poisson + two-case Graetz kernels on the lens cross-section (Section 5 C) |
| `schematics.py` | Interactive labelled diagrams (unit geometry, lens, solve topology) |
| `cockpit.py` | Browser port of the lumped solver (whole-design-space cockpit) |
| `livepack.py` | Supporting live-pack helpers for the cockpit |
| `coolants.csv` | Fluid property table (ester, HFE, mineral, PAO, water for reference) |
| `pack_benchmark.xlsx` | 58-pack BEV benchmark database |
| `README.md` | Changelog through v9.3 |
| `.streamlit/config.toml` | Brand theming |

The lumped solver (`app.py`) and the zonal solver (`zonal.py`) are
**independent machinery** solving the same pack. Their agreement
(Section 6) is one of the cross-checks — judge how independent it really
is.

---

## 5. The physics stack, with equations and how each result is reached

Notation: `T_c` cell wall, `T_p`/`T_f` plate face, `T_b` oil bulk, `T_w`
water, `q'` line flux [W/m], `a` duct kernel [W/(m·K)], `D_h` hydraulic
diameter, `fRe` friction-factor–Reynolds product, `H` cell height, `s`
oil slot width, `g_x` in-row crevice, `p` pitch.

### A. Layout from the fin rule

The plate is a fin rooted along the tube lines and loaded over its face.
Root pitch and bay efficiency:

$$m=\sqrt{\frac{2\,h_\text{face}}{k_p\,t_p}},\qquad P=\frac{2}{m},\qquad
\eta_\text{bay}=\frac{\tanh(mP/2)}{mP/2}\cdot\frac{\tanh(mH/2)}{mH/2}$$

`P` is snapped to the cell lattice so tubes land on real gaps, never
arbitrary positions. *Representative:* `m ≈ 31 /m`, rule pitch 64 mm →
snapped to 43 mm (every 2 cells) → 17 tubes per plate at the full-pack
default, `η_bay ≈ 0.64`.

**How reached:** closed-form fin theory; the only input is the
face-averaged `h`, taken from the kernel (below). **Challenge:**
`η_bay` is a product-rule lumping of a genuinely 2D plate temperature
field; the plate is not isothermal along its span.

### B. Hydraulics: laminar slot against buoyancy

Each slot column balances available head (pump + thermal buoyancy)
against fully developed laminar friction; viscosity at the local bulk
temperature (Andrade law):

$$\Delta p=\Delta p_\text{pump}+\rho\beta g H(\bar T-T_\text{in}),\qquad
\bar u=\Delta p\,\frac{2D_h^{2}}{fRe\,\mu H},\qquad
\nu(T)=\nu_{25}\,e^{B(1/T-1/298)}$$

**How reached:** `fRe` and `D_h` come from the velocity Poisson solve on
the true lens (Section C), so the friction is geometry-exact, not a
round-duct stand-in. **Result that surprised me:** flow sensitivity to
slot width is nearly linear here (`h ∝ Q''^{~0.7}`-like), **not** the
flat-slot cubic law — see the choke test in Section 6. **Challenge:** the
buoyancy term is a single-`ΔT`, single-`H` lumping; there is no
distributed buoyancy along `z`.

### C. Duct kernels (`fea4_channel.py`): exact developing-flow coupling

On the true lens cross-section (cell crown vs flat plate over one pitch),
the velocity is a Poisson solve and the developing energy problem is
solved twice (one wall heated at a time) and superposed into an exact
2×2 transport matrix:

$$\nabla^2 w=-1\ (\text{no-slip}),\qquad fRe=\frac{2D_h^{2}}{\bar w}$$

$$\begin{bmatrix}q'_c\\ q'_p\end{bmatrix}=\mathbf{a}(z^{*},s)
\begin{bmatrix}T_c-T_b\\ T_p-T_b\end{bmatrix},\qquad
z^{*}=\frac{z}{D_h\,Re\,Pr}$$

This is the influence-coefficient method for doubly connected ducts
(Shah & London). The march is an implicit, splu-prefactorised,
log-spaced `z*` sweep from `3e-4`.

**Validation gates (parallel-plate limits with known answers):**
`fRe = 95.7` vs 96.0; one-wall H2 `Nu = 5.387` vs 5.385; both-walls-T
`Nu = 7.541` vs 7.541 (exact). Two operator bugs were caught *because*
these gates failed first: a masked-cell Dirichlet term baked into the
energy operator (blew a coefficient up by ~1e10), and a wall-flux
extrapolation error. **Representative kernel** at the operating
the operating `z* ≈ 4×10⁻⁴` (mid-height, at the solved velocity): the
stagnant-thermosiphon value, because the whole 70 mm slot sits in the
cell-side film is roughly 3× the stagnant-thermosiphon value, since the whole 70 mm slot sits in the thermal entrance (`Pr ≈ 127`, `Re ≈ 120-160`, ū ≈ 122-162 mm/s after the buoyancy correction).

**Challenge:** the staircase mask of the circular cell overestimates the
cell perimeter (`P_c ≈ 0.042` vs true `0.033`); I use a consistent-flux
formulation so integrated quantities are right, but check the local wall
temperatures. Also check the kernel-bank interpolation (bilinear in `s`,
log-linear in `z*`) between the five tabulated slot widths.

### D. Channel energy march (exact exponential, energy-consistent)

Marching bulk enthalpy up each segment, the walls act through the
conductance-weighted mix `T_eff`. The update is the **exact** solution of
the linear segment (unconditionally stable), and the segment-mean driver
uses the exact factor `f̄` so flux and enthalpy rise agree to machine
precision at any flow, including starved channels:

$$\dot m c_p\frac{dT_b}{dz}=G(T_\text{eff}-T_b),\quad G=\textstyle\sum a,
\quad T_\text{eff}=\frac{G_cT_c+G_pT_p}{G}$$

$$T_b(z{+}dz)=T_\text{eff}+(T_b-T_\text{eff})e^{-x},\quad
\bar f=\frac{1-e^{-x}}{x},\quad x=\frac{G\,dz}{\dot m c_p}$$

**How reached:** this replaced a naive forward step that diverged to
80,000 °C when a channel starved (`x` large). The exponential form is the
integrating factor of the segment ODE. **Challenge:** `T_c` and `T_p` are
frozen within a segment (piecewise-constant wall temperature); with
`nz = 10-16` the z-independence check is 0.02-0.03 °C, but confirm that is
tight enough for your purposes.

### E. Cells: heat that fights back (DCIR feedback)

Each cell is an isothermal balance; generation falls as it warms. Solved
by Newton using the march's own Jacobian, so cell and channels are
consistent rather than lagged:

$$q_\text{gen}(T)=(C\,\text{Ah})^{2}R_{dc}\,e^{-k_{dc}(T-25)}
=\hat q+J(T-T_\text{old})$$

**Challenge (important):** the zonal cell node is the **can**, treated as
isothermal. The jellyroll core-to-can rise (`R_core = 1/(4\pi k_r H)`,
about 2-3 °C at 2C on a 4680, less on a 21700) is **not** added inside the
zonal model — it lives only in the lumped solver. If you care about core
temperature, that term must be superposed.

### F. Plate face and root network (closed form)

The plate face is not iterated. Its heat is linear in its own
temperature, so it is solved directly against the series root resistance,
and the water is marched tube by tube afterwards:

$$R_\text{root}=\frac{1}{h_c A_\text{ct}}+\frac{\ln(d_o/d_i)}{2\pi k_t L}
+\frac{1}{h_w\pi d_i L}$$

$$T_f=\frac{T_w/R_\text{root}-S_0}{W+1/R_\text{root}},\qquad
q_\text{root}=\frac{T_f-T_w}{R_\text{root}}$$

`W` and `S_0` are accumulated from the same march (per bay), so the plate
closure is exact given the bulk field. Water film `h_w` from Gnielinski
(turbulent) or Hausen (laminar developing).

**Representative resistances:** water film `R_film ≈ 6 K/W` per 21.5 mm
crossing dominates; contact (`h_c = 8000 W/m²·K`) is `≈ 0.55 K/W`, a
minor share; tube wall is negligible. **This is the governing resistance
of the whole design** — the design lever is water-side flow and tube
count, not tolerance. **Challenge:** `h_c = 8000 W/m²·K` is an assumption
pending braze data (Section 7); and axial conduction *along the tube
between plates* is neglected (each crossing is an independent segment).

### G. Recirculation plenum (closed-form fixed point)

The mixed exit re-enters as the common inlet. Each column exit is affine
in the inlet, so the recirculation fixed point is solved in one line
rather than ratcheted — this removed the slow mode that stalled
convergence:

$$T_\text{exit}=A+P\,T_\text{plen},\quad P=\prod_z e^{-x}\ \Rightarrow\
T_\text{plen}^{*}=\frac{\langle A\rangle}{1-\langle P\rangle}$$

**Challenge:** this assumes a single well-mixed plenum feeding all
channels equally (a manifold that truly equalises inlet temperature).
For a real serpentine/parallel manifold, check whether inlet
non-uniformity matters.

### H. Tolerance Monte Carlo

$$s_i\sim\mathcal N(s_\text{nom},\sigma_s)\,[\text{trunc.}],\quad
c_k\sim\mathcal N(c,\sigma_c),\quad
P_\text{exceed}=\Pr(T_\text{max}>T_\text{lim})$$

Every sample is a full network solve, warm-started from the converged
nominal case (about 30 s for 60-90 samples). **Result:** the plates and
the mixed plenum homogenise the pack, so manufacturing scatter barely
moves the hottest cell (scatter < 0.1 °C, indistinguishable from solver residue; 0/90 exceedances bound the probability to ≤3.3% at 95%, at
`σ_s = 0.25 mm`, `σ_c = 0.10`). **Challenge:** only slot width and contact
are sampled; plate thickness, tube position, pitch, and fluid property
scatter are not. Add the ones you think matter.

---

## 6. Validation evidence (please classify each as real or circular)

| # | Test | Result | What it proves |
|---|---|---|---|
| V1 | `fRe` parallel-plate limit | 95.7 vs 96.0 | velocity Poisson + `D_h` correct |
| V2 | One-wall H2 Nusselt | 5.387 vs 5.385 | single-wall energy kernel correct |
| V3 | Both-walls-T Nusselt | 7.541 vs 7.541 | two-wall superposition exact |
| V4 | Kernel reciprocity | `a_cp ≈ a_pc` | operator symmetry |
| V5 | Energy closure, full pack | 0.19-0.20% | in = out across the network |
| V6 | z-resolution independence | 0.02-0.03 °C (nz 10 vs 16-18) | march is grid-converged |
| V7 | Uniform-input collapse | spread → 0.00 | no spurious per-channel asymmetry |
| V8 | Lumped cross-check (**matched**) | zonal ~41.9 vs matched-serpentine lumped ~35.8 °C | two different-machinery solvers; the few-°C gap brackets modelling uncertainty (not "agreement") |
| V9 | Choke test | 0.5 mm narrowing → 81% flow share (cubic law predicts 42%) | flow sensitivity is near-linear, decomposed to the lens geometry metric `D_h²A/fRe` ratio 0.81 |
| V10 | Monte Carlo robustness | scatter < 0.1 °C (at/below solver residue); 0/90 exceedances → ≤3.3% (95% UB, rule of three) | tolerance-robust to *scatter*; NOT to the h_c epistemic band |
| V11 | FEA1 (lid, Roark) | +0.0% vs Roark | pressure-vessel mass |
| V12 | FEA2 (core-to-can, exact) | -0.1%; app 20-34% conservative | cell conduction |
| V13 | FEA3 (water-rise spread) | energy closure 0.00% | manifold plumbing |

Numbers are representative at the default operating point (2C, 10 L/min,
21700, ester); they move with inputs. The single result I most want
challenged is **V8**: the two solvers share some correlations (film
laws, fluid properties), so agreement is reassuring but not a fully
independent check.

---

## 7. Assumptions register

| Assumption | Value / form | Basis | Risk if wrong | Enters at |
|---|---|---|---|---|
| Cell can isothermal | `T_c` single node | high-k steel/aluminium can | underestimates core temp | E |
| Core-to-can rise excluded from zonal | — | kept in lumped solver | zonal `T_max` is a can temp, not core | E |
| Root contact conductance | `h_c = 8000 W/m²·K` | braze-grade estimate, **no data yet** | shifts root resistance, but it is the minor share | F |
| Laminar oil in slots | `Re ≈ 120-160`, `fRe` from Poisson | far below transition | entrance/buoyant mixing could raise `h` | B,C,D |
| Water film | Gnielinski/Hausen | standard correlations | modest `h_w` error → level shift | F |
| Plate fin knockdown | product-rule `η_bay` | 1D fin theory | 2D plate field not captured | A,F |
| Single mixed plenum | equal inlet to all channels | idealised manifold | inlet non-uniformity ignored | G |
| Tube axial conduction | neglected between plates | short crossings | small unless tubes are long | F |
| Radiation | excluded | ~1 W at these `ΔT` | negligible | all |
| Staircase cell perimeter | `P_c` overestimated ~27% | FD mask on circle | local wall temps, not integrals | C |
| MC scatter sources | slot width + contact only | dominant tolerances | other tolerances unmodelled | H |

---

## 8. Known soft spots (where I would attack it first)

1. **`h_c = 8000 W/m²·K` is unmeasured.** It is the minor resistance
   share, so the answer is not very sensitive to it, but that claim
   itself should be checked with a sweep.
2. **`η_bay` product-rule fin knockdown** is the crudest step in an
   otherwise exact chain. A proper 2D plate solve (a fea5) would replace
   it. Estimate the error.
3. **V8 is not fully independent.** Both solvers use the same film
   correlations and fluid table. A truly independent check would be a
   CFD run or a physical rig point.
4. **Kernel interpolation** between five slot widths and along `z*`. Check
   the interpolation error against a directly computed kernel at an
   off-grid `s`.
5. **The lens vs the real gap.** The cross-section assumes cell crown vs
   flat plate with a clean lens of oil. Real packs have the in-row
   crevice `g_x`, holder blockage, and end effects the 2D slice omits.
6. **Piecewise-constant wall temperature** within a z-segment. Fine at the
   tested `nz`, but confirm for tall cells or low flow.
7. **The cockpit is a separate JS codebase.** It is audited to +0.00 °C
   against Python on load, but it can drift; do not treat it as a second
   independent model.
8. **Buoyancy model** is single-`ΔT`/single-`H`. In a genuinely
   thermosiphon-dominated case (`Δp_pump → 0`) the distributed buoyancy
   might matter; the choke test held at `dp = 3 Pa` but push it.

---

## 9. Red-team checklist (please work through these)

**Physics**
- [ ] Is the advective-primary reframing correct, or is there a regime
      where conduction-through-slot competes? (Section 2, C)
- [ ] Are the film correlations (Churchill-Chu, Churchill-Bernstein,
      Hausen, Gnielinski) applied in-range? Check `Ra`, `Re`, `Pr` bounds.
- [ ] Does the `η_bay` lumping over- or under-predict plate performance
      vs a 2D plate solve?
- [ ] Is the core-to-can omission in the zonal model acceptable for your
      pass/fail criterion, or must it be superposed everywhere?

**Numerics**
- [ ] Reproduce the three fea4 gates (`python fea4_channel.py`). Do they
      hold on your grid refinement?
- [ ] Push `nz` and `iters` in `zonal.py`; does closure keep falling and
      `T_max` stay put?
- [ ] Test the exponential march at extreme flow ratios (one channel at
      10% of mean); does energy still close?
- [ ] Check the plenum fixed point against a brute-force iterated plenum.

**Code**
- [ ] Read `zonal.solve_zonal` and confirm the plate-face closed form
      (`W_bay`, `S_bay`, `Tf`) matches the derivation in F.
- [ ] Confirm the kernel bank interpolation is not extrapolating outside
      the `s` grid at your inputs.
- [ ] Run the AppTest snippet (Section 3); zero exceptions expected.

**Validation**
- [ ] Classify each V1-V13 as real evidence or circular.
- [ ] Propose the one physical rig point or CFD case that would most
      cheaply break V8.

**Explanations / UX**
- [ ] In the Zones "Methods" panel, does any equation misrepresent what
      the code does?
- [ ] In the Learn tab, is any plain-language claim wrong or oversold?
- [ ] Is anything in the app presented as validated that is actually an
      assumption?

---

## 10. What good feedback looks like

For each issue: **(a)** the exact equation / file+line / number, **(b)**
the discrepancy or objection, **(c)** a reference or a reproducing input
if you have one. "The `η_bay` product rule overpredicts by ~15% vs a 2D
plate solve at `t_p = 1.5 mm`, see [ref]" is worth ten "seems
optimistic"s.

Thank you. The two subsystems that most need your eyes are the fea4
kernel derivation (Section 5C) and the claim that the design is
tolerance-robust because the water film governs (Sections 5F, 6 V9-V10).
If either is wrong, the headline conclusions change.
