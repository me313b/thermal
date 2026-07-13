# Response to the v9.4 verification — final version v9.5

**From:** the development side. **Date:** 13 Jul 2026.
**Re:** *Verification of v9.4 against the red-team findings* (sign-off).

Thank you for the sign-off. You confirmed all ten v9.3 fixes are
implemented as described, every verifying number reproduces, there are no
regressions, `requirements.txt` installs one, and the two claims most in
need of independent confirmation (the buoyancy reference and the shared
water correlation) are correct in code, not just in the note. The v9.3
findings are closed.

You left four minor, explicitly non-blocking observations and answered
the two open questions. v9.5 folds all four in and surfaces both answers
in the tool. Nothing here changes a result; it removes the last stale
numbers and makes the core-node story impossible to miss. `SMOKE=1` and a
full headless AppTest both pass.

---

## The four observations

### Obs 1 — the response memo understated the default by ~0.4 °C [FIXED, doc]

You are right: the code gives **can 42.4 / core 44.9 °C** on the app bank,
a **+1.15 °C** shift over v9.3, not the "+0.4-1.1" the memo quoted.
`REDTEAM_RESPONSE.md` is corrected throughout (default numbers, the F2
verify line, the matched-gap number, the core-margin sentence). The
memo, not the model, was stale.

### Obs 2 — the h_c crossing is node-dependent; "4000-6000" straddled it [FIXED, code + doc]

Corrected everywhere. The **can** crosses 45 at h_c ≈ 2900 W/m²·K; the
**core** crosses at ≈ 7500 - within ~7% of the assumed 8000. The Zones h_c
sweep now interpolates and reports **both** crossings, and its caption
states the consequence plainly: on a can criterion the braze can be quite
imperfect, but on a core criterion the design passes only if the joint is
essentially as good as assumed. The memo's band is replaced by the two
node-specific numbers.

### Obs 3 — the Monte Carlo already logs a core exceedance [FIXED, code]

`monte_carlo` now returns the core-node exceedance count, and the Zones
tab shows it beside the can-based one. When the core clears (0/M) it reads
as a green confirmation with the rule-of-three upper bound; when the core
exceeds - which it does on the default_d design point (worst core sample
45.03 °C, and a meaningful fraction of samples land at 45.0-45.03 because
the hottest-core cell preferentially sits where the contact ran low - it
reads as a warning that states the rule of three no longer applies on that
node and the design is at the limit as drawn. The can-based metric is
unchanged (0/M). So the tab now makes Obs 2 visible in the Monte Carlo,
exactly as you suggested.

### Obs 4 — the layout is formally bank-dependent [FIXED, code]

`h_face_design()` reads a11 from the caller's bank, so the value drifts a
few percent with grid resolution (255 at n=100, 265 at n=70, 278 at the
SMOKE mini-bank). All snap to 17 tubes here, but a future bank change
could move the boundary silently. The `zonal.py` shakedown now **asserts
the tube count** and prints the snap margin: the 43 mm / every-2-cells
snap holds until h_face exceeds ≈ 333 W/m²·K, so the design point has
~26-30% headroom before it would flip to 33 tubes. A regression that
moves the snap now fails the shakedown rather than passing quietly. A
`snap_boundary_h()` helper documents the boundaries. The docstring spells
out the bank-dependence.

---

## Your two open questions, now surfaced in the tool

### (a) Is the truth nearer 36 or 42? → Neither; ~38-40, and the default is ~2-4 °C conservative

Adopted. Your bracket - zonal with the bypass fully on (38.3, across
h_oil_tube 80-250) and the matched lumped with the collar chain added
(38.5 at h_c 8000) landing within 0.2 °C of each other - is now stated in
the matched cross-check caption in the Zones tab: completing each model
with the other's missing physics collapses the gap to **38-40 °C**, so the
shipped default with the bypass off is most likely **~2-4 °C
conservative**, and the rig point buys back exactly that margin. The
residual one-sided bias (the lumped crossflow cell film, optimistic for
axial slot flow) is why we quote 38-40 rather than a point.

### (b) Is 45 a can or a core number? → A specification decision the tool now forces into the open

Agreed it is not decidable in code, and agreed on the framing: cell
datasheet operating limits are conventionally surface/can temperatures,
but plating risk - the phenomenon the CC-CV derate already models - is
core-referenced. v9.5 displays both nodes (metric, sweep, and Monte
Carlo), and the h_c sweep caption ties the two criteria to their h_c
requirements (can → margin to ~2900; core → margin to ~7500). If 45 was
chosen as a longevity/plating proxy, the core is the honest node and this
design is at the limit as drawn. That question goes to the cell datasheet
and the rig, not to more modelling - and the tool no longer lets it hide.

---

## Where this leaves it

The corrected tool is less comfortable than its predecessor and more
truthful for it: the numbers moved in the direction the physics demanded,
not the direction the headline wanted. The remaining uncertainty is no
longer modelling; it is measurement. The single experiment that closes
both open questions is unchanged and still recommended: **one tube, two
brazed plate crossings, in a stirred 40 °C oil bath**, run at the
pack-representative per-tube flow. State (b) bare-tube isolates the
oil-to-tube path; state (a)−(b) isolates h_c·A_ct including the collar
factor; a plate-face thermocouple 10 mm from a collar checks the
spreading picture. It pins the two parameters the 38-40 bracket is
conditional on, and it settles whether the core clears 45.

v9.5 is the final build. Everything you raised across three rounds is
closed in code or surfaced in the UI; what is left is a bath, a tube, and
two thermocouples.
