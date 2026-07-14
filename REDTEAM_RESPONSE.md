# Response to the red-team review — Immersion Pack Lab v9.4

**From:** the development side. **Date:** 13 Jul 2026.
**Re:** *Red-team review of v9.3, zonal solver and Monte Carlo.*

Thank you. This was exactly the review the tool needed: three real bugs,
two oversold claims, a documentation defect, and a packaging break, every
one reproduced with inputs. All ten suggested fixes are implemented, plus
the core-temperature superposition (F10) and the handover correction
(F8). This note maps each finding to what changed and gives the verifying
number so you can re-run it. Every fix is in the v9.4 bundle; `SMOKE=1
python app.py` and a full headless AppTest both pass.

The headline honesty outcome up front: **the corrected model is less
comfortable than v9.3, and that is the point.** At the default the can
now sits at ~42.4 °C and the **core at ~44.9 °C** (F10; +1.15 °C over v9.3), the buoyancy fix
removes ~9 Pa of spurious head and drops the solved velocity from 166 to
~122-162 mm/s depending on fluid (F1), and the contact conductance is now
a slider with a sweep showing the pack crosses 45 °C node-dependently -
the can below h_c ≈ 2900 W/m²·K, but the core only above ≈ 7500 W/m²·K,
within ~7% of the nominal 8000 (F2). We did **not** engineer the number back to 41.2. The tool now
tells the truth it was hiding: this design is marginal, and its margin
depends on the unmeasured collar contact and on whether the limit is a
can or a core limit.

---

## Bugs

### F1 — Buoyancy referenced to the water inlet [FIXED]

`zonal.py` hydraulics now references the buoyant head to the **loop
return (plenum)** temperature, not the water inlet:
`dp = dp_extra + rho.beta.g.H.clip(T_bar - T_plen, 0, 45)`. In the closed
recirculating loop the net head is `rho.beta.g.H.(T_riser - T_downcomer)`
and the downcomer sits at `T_plen`, so at the default the per-pass rise is
tiny and the net head collapses to ≲0.1 Pa. **Verify:** the shakedown
prints `net head 0.00 Pa of 25 Pa pump (was ~9 Pa)`; ū falls from 166 to
124 mm/s and T_max rises +1.15 °C (F1 and F4 combined). Your dp = 3 Pa case now
correctly exceeds the limit rather than hiding behind phantom head.

### F5 — No laminar-turbulent water bridge [FIXED]

Both solvers now call a single shared function,
`correlations.water_nu(Re, Pr, d_i, L)`, with the continuous 2300-3000
blend. The zonal no longer steps at Re = 2300, and V8's two comparators
can never diverge on this correlation again — which was your explicit
ask. **Verify:** sweep `flow_lpm` through 13.7-13.9 L/min; the −3.3 °C
step is gone.

### F4 — Fin rule on a stale h_face [FIXED]

The layout is now derived from the **kernel-implied** plate film,
`h_face = a11(z*_mid)/pitch` at the design point (`h_face_design()` in
`zonal.py`), so the fin rule is self-consistent with the duct it feeds. At
the default this gives h_face ≈ 265-280 W/m²·K rather than the hardwired
150. The Zones caption now shows the derived h_face and labels the product
rule eta as conservative (your 2D solve put it 3-16% low, i.e. T is an
upper bound in that respect). We took your steer and did **not** build a
live fea5 — the product rule stays, clearly flagged conservative.

---

## Oversold claims

### F2 — h_c = 8000 "minor share" is conditional [FIXED]

`h_contact` and the previously-undocumented `collar_factor` are now named
inputs (`default_d`, and sliders in the Zones tab). A one-click
**sensitivity sweep** plots T_can and T_core against h_c with the 45 °C
line, and reports the h_c at which the can crosses the limit. **Verify:**
the sweep reproduces (can) 20000→41.5, 8000→42.4, 4000→43.9, 2000→46.8 °C;
the caption states both the can and the binding core crossing explicitly. The headline is no longer
presented as unconditional.

### F3 / V8 — Mode-mismatched comparison; omitted oil→tube path [FIXED, two parts]

*The comparison* is now **matched**: the Zones tab solves the lumped
serpentine model at the zonal's own velocity and reports that number
(~35.8 °C against the zonal mean ~42.4 °C), with a caption spelling out
that the two use different machinery (lumped crossflow cell films + area
credit; zonal per-crossing collars) so a few °C gap is expected and
brackets the modelling uncertainty. The misleading "agree within a
fraction of a degree" line is gone, in both the metric and the Methods
panel.

*The omitted path* is now available as an explicit, bounded option:
`wetted_tube_frac > 0` switches on the direct oil→tube bypass. Off by
default (conservative, and its magnitude is header-geometry dependent, so
we expose rather than assume it). **Verify:** at frac 0.5 the shakedown
gives UA ≈ 21 W/K carrying ~350 W (15% of duty), dropping T_max ~2.3 °C —
squarely inside your 12-29% / 3-4 °C estimate. Switching it on *improves*
the margin honestly, exactly as you noted.

### F6 — Monte Carlo statistics [FIXED]

Per-sample tolerance tightened (`iters=60, tol=0.0015`), and the caption
now states plainly that the sample spread is at or below solver residue,
so it reads as *scatter indistinguishable from convergence noise,
< 0.1 °C*. Zero-exceedance is reported as a **rule-of-three** 95% upper
bound (3/M) rather than a bare "P = 0", the core temperature is carried
through the samples, and the caption points at the h_c epistemic band as
the real risk the Monte Carlo does not sample. **Verify:** the exceedance
metric reads "0 / 90, ≤3.3% (95% UB)".

---

## Minor, documentation, packaging

### F7 — Kernel warts and one presentation error [FIXED where it mattered]

The presentation error is corrected: `fea4_channel.py __main__` no longer
prints the 0.40-vs-0.70 cross-coupling as a "reproduction"; it is labelled
a different regime shown for scale. The entrance drift (~2.5%), noisy FD
tail (±7%, starved channels only), and z*-floor (+0.02 °C) are documented
as known first-order/limit behaviour; per your assessment they do not move
the answer, so we did not chase them. V4 is described as a sanity check on
a non-self-adjoint operator, not an exactness gate.

### F8 — Handover described an operating point the app never shipped [FIXED]

`REDTEAM_HANDOVER.md` §5 is corrected to the shipped/solved numbers
(ū ≈ 122-162 mm/s, Re ≈ 120-160, mid-height z* ≈ 4×10⁻⁴, and the
buoyancy-corrected regime), so a reviewer validating against the narrative
no longer chases the slow-point ghosts.

### F9 — Packaging and the Methods preamble [FIXED]

`scipy>=1.10` is now in `requirements.txt` (a clean install no longer
crashes the Zones tab). `zonal.py` has the `__main__` shakedown the
handover promised — it prints the mini/full solves, closure, z-check,
layout, the F1 head, the F2 sweep, and the F3 bypass. The Methods preamble
no longer claims "nothing is asserted"; it names its three engineering
assumptions (h_c, collar factor, conservative fin eta), which was your
answer to our own checklist question turned back on us.

### F10 — Core vs can [FIXED]

The core-to-can rise `R_core = 1/(4 pi k_r H)` is superposed and the peak
**core** temperature reported next to the can (metric, Monte Carlo, and
shakedown). At the default the core is ~44.9 °C against 45 — so if the
limit is a plating/core limit, the margin is ~0.5 °C, and combined with a
sub-nominal h_c it fails. The tool now makes that visible instead of
implying 3.8 °C of headroom.

---

## Your rig point

Adopted as the recommendation: **one tube, two brazed plate crossings, in
a stirred 40 °C oil bath**, run at the pack-representative per-tube flow.
State (b) bare tube isolates the oil→tube path the model omits (F3); state
(a)−(b) isolates `h_c·A_ct` including the collar factor (F2). A plate-face
thermocouple 10 mm from a collar checks the η/spreading picture. This is
in the updated handover as the single highest-value next experiment,
ahead of any CFD, because CFD cannot deliver h_c.

---

## What we deliberately did not change

- **No fea5.** Your 2D solve showed the product-rule eta is conservative
  by 3-16%; we flag it as conservative and leave it. The h_face
  inconsistency (F4), which you rated the bigger issue, is fixed.
- **Kernel entrance drift and FD tail (F7 i, ii).** Documented, not
  chased; they move T_max by ≤0.02 °C.
- **The bypass stays off by default.** It is conservative when off and its
  magnitude is genuinely header-dependent, so we expose it rather than
  bake in a guess. The rig point will fix its magnitude.

Scripts and fixes are all in `immersion_pack_lab_v9_4_full.zip`. If you
want a second pass, the two things most worth re-attacking are the matched
V8 gap (is the truth nearer 36 or 42?) and the core-limit criterion (is
45 a can or a core number?) — both now surfaced, neither yet closed by
measurement.
