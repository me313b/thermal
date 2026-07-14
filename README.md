# Immersion Pack Lab v4

Parametric design and teaching app for a **static or stirred immersion-cooled
21700 battery pack** with an **internal water-cooled tube heat exchanger**
(the ICDC architecture of Wang et al. 2023).

## Run
    pip install -r requirements.txt
    streamlit run app.py

Keep `coolants.csv` next to `app.py`. `SMOKE=1 python app.py` runs the physics
self-test (Wang benchmark, thermosiphon sanity, calibration round-trip).

## Model
Two-node thermal network (cells, bulk oil) with water and ambient boundaries;
film coefficients recomputed each step from local dT.

* Cell side: Churchill-Chu vertical plate x confinement penalty below 6 mm gap
  (calibrated to Wang Fig. 9). Tube side: Churchill-Chu horizontal cylinder +
  Schmidt annular fins. Forced/mixed flow: Churchill-Bernstein, Nu^3 blend.
* Water side: Hausen / Gnielinski with laminar-turbulent bridge and axial rise.
* **Thermosiphon (v2)**: buoyant head rho*beta*g*H_loop*dT against laminar bank
  friction + minor losses (K_loop, default 5), solved by bisection; H_loop set
  by the tube-plane position. Feeds the films as the floor on oil velocity.
  Order-of-magnitude validation: Wang measured 0.5-1.8 mm/s; app predicts
  ~4 mm/s for the default 20 kWh pack (bigger head, bigger Q).
* **DCIR(T) (v2)**: R = R25 exp(-k(T-25)), default 1.2 %/K, coupled into heat
  generation each iteration/step. Set k = 0 to disable (recovers v1 numbers).
* **Core temperature (v2)**: quasi-steady radial, R_core = 1/(4 pi k_r H),
  k_r default 0.9 W/mK; peak (centreline) value, mean is half.
* **Worst cell (v2)**: +/- (water rise + stratification)/2 about the mean,
  checked against the 5 K uniformity criterion. An estimate, not a CFD result.

## v2 features
Predicted thermosiphon and tube-plane placement; DCIR(T); core temperature
and core-limit option; worst-cell spread; plan-view Layout tab; Architecture
comparator (static / stirred / bottom cold plate / pumped external HX, with
parasitic power and thermal-system mass); water-pump and stirrer power;
thermal-runaway screening (zone heating, bulk rise, headspace pressure);
design save/load as JSON; A/B pin-and-compare; duty-cycle CSV import (t_s +
C or P_kW); measured-data overlay with one-click calibration-factor fit
(multiplies both oil films, the spray-app pattern); 7x7 two-lever sweep
heatmap with limit contour; student mode (hides Decide).

## Fixed engineering constants (edit in code if you have better numbers)
Cold plate: axial cell path k_ax = 28 W/mK, TIM 0.8 K/W, channel film
3000 W/m2K. Pumped case: 20 cm/s past cells, oil sized for 5 K rise, 30 kPa
at 40% pump efficiency, 25 mL fluid/cell (HPB80 ratio). Runaway: vent gas
heated to 380 K, default 55 kJ/cell with 60% into the local zone. Pump/stirrer
efficiencies 35%/30% (fluid power only; controller overheads excluded).

## Validation
Benchmark tab rebuilds Wang's rig with fixed heater power, cal = 1, no DCIR:
cell 33.0 degC at 1800 s vs 32.3 measured, oil 28.8 vs ~28.5, resistances
within ~30-50% of their Table 6. Treat oil-side h as +/-30% until calibrated.

## Sources
Wang et al. 2023, J. Energy Storage 62, 106821. Zou et al. 2024, J. Energy
Storage 83, 110634. Roe et al. 2022, J. Power Sources 525, 231094.
batterydesign.net (AMG HPB80; Dielectric Immersion Cooling).
coolant_comparison_reviewed.xlsx (reviewed fluid table).

## v3 additions
Electrical: SoC tracking, generic NMC OCV(SoC) and entropic heat, CC-CV
charge with plating-derated current (generic derate map, editable arrays at
the top of app.py), charge/discharge DCIR split, repeated cycling duty.
Drive cycles: WLTP 3b, NEDC, UDDS, HWFET, US06, Artemis Motorway 130 as
coarse breakpoint profiles uniformly speed-scaled to the official distance
(exact by construction; thermally adequate since the pack filters >~0.01 Hz;
use CSV upload for official 1 Hz traces). Vehicle model: mass, CdA, Crr,
drivetrain efficiency, capped regen, accessories. Steady state uses the
duty's RMS C-rate. Structural: cell format presets (18650/21700/4680 from
teardown data; 4680 DCIR is an estimate), busbar heat and mass (auto-sized
at a set current density), cell holders (mass + thermosiphon blockage),
calculated enclosure (stiffened-plate sizing at the burst-disc pressure)
giving Wh/L and an honest total weight, interstitial tube routing option.
Workflow: S x P suggester, one-click HTML design report, +/-30% uncertainty
band, goal-seek on a single lever, production-pack benchmark table
(Tesla 2170 / 4680 / Plaid / LFP, AMG HPB80) with teardown sources.

## v3 constants to challenge
Enclosure: sigma_allow 80 MPa, stiffening knock-down 0.45, design pressure =
burst-disc 0.5 bar g. Busbars: copper at 5 A/mm2, length 1.15 x Ns x pitch.
Holders: 8 g/cell, 20% flow blockage. Interstitial mode: thermosiphon head
8 mm, stratification x0.35 (distributed sinks) - an engineering judgement,
not CFD. OCV/entropic/plating arrays are generic NMC shapes. Vehicle
defaults: 1900 kg, CdA 0.62, Crr 0.009, eta 0.92, regen 0.65 capped 60 kW.

## v4: workbench restructure and 3D view
Layout reorganised around the workflow: 1 Design (all pack/cooling/structure
inputs in three domain columns, with a live 3D or plan view), 2 Duty (load
definition beside the transient response), 3 Results, 4 Improve (diagnosis,
sensitivity, goal-seek, sweep, A/B pin), 5 Safety (runaway screening,
expansion, leak checklist), 6 Compare (architectures, coolant shoot-out,
production packs), 7 Learn, 8 Validate and tune (Wang benchmark, rig
calibration, model-tuning knobs, fluid table). The sidebar is now a status
card plus save/load and report export only. Student version removed.
The 3D view renders every cell (coloured by centre-vs-edge tendency), tube
runs with translucent fin envelopes, the oil fill level, and the enclosure;
it follows format, pitch, tube plane (including interstitial) live.
Old saved designs load unchanged (same widget keys).

## v4.1 fixes
Apply-suggestion / apply-preset / calibration-fit buttons now write widget
values through a pending queue drained at the top of the script, removing the
StreamlitAPIException (state written after widget instantiation). Fluid
scatter charts (Learn panel 3, Coolant shoot-out) auto-thin their labels to
avoid pile-ups; every point keeps its name on hover.

## v5: heat-journey Results, scale-to-C, full report tab
Results is now the heat journey: a Sankey of where the watts go
(cells/busbars -> oil -> water/casing -> chiller -> ambient), the
temperature ladder and resistance ranking, then seven station cards
(source, first film, bulk oil, finned film, wall, water, chiller) each with
its numbers, its share of the ladder, and a quantified "how to improve this
stage" note, closing with the water set-point trade curve (COP vs DCIR).
A chiller model (COP = 0.45 x Carnot on the actual lift) sizes electrical
power. Improve gains "Scale to a target C": solves the design at the target,
tables every subsystem's required change (chiller, single-lever fixes via
goal-seek, stirred combination, busbar resize, water regime, core, plating/
charging), plus a heat/chiller/temperature vs C sweep chart. Tab 9 renders a
full narrative design report (executive summary, design, heat journey, duty,
scaling, uniformity and safety, validation - each with improvement actions
and live figures) and exports it as standalone interactive HTML.
Mass audit moved to Design. Sidebar simplified.

## v6: visual redesign
One design system, applied everywhere: gradient hero header with a live
status chip, KPI cards with large numerals (replacing default metrics),
pill-style tab navigation, soft-shadow white cards for every section,
Inter typography, and a single global plotly template (clean light grid,
modern colorway) so all charts share one look. Chart semantics kept but
refreshed: cell red #EF4444, oil amber #F59E0B, water sky #0EA5E9, brand
indigo/cyan for structure. Sidebar status is a card; station progress bars
use the brand gradient; the HTML report export matches. No physics or
widget-key changes: saved designs load unchanged. All styling lives in the
single CSS block and the "packlab" plotly template at the top of app.py,
so re-theming (including dark mode) is a one-block edit.

## v7: clarity, interactivity, and the serpentine-plate idea
Clearer physics graphics: the heat-flow Sankey is rebuilt with fixed columns,
kW-labelled nodes and stream-coloured ribbons; a new interactive thermal
circuit shows every interface in the cell-to-water chain sized by its share
of resistance, coloured by transport mode (conduction / convection /
advection), with the governing law and live numbers on hover and the weakest
link flagged. Units are typeset properly throughout (m², W/m²·K, °C, mΩ, ν).
Design now has a sticky live panel: the 3D/plan view plus a change summary
(temperatures, margin, spread, mass, deltas vs the pinned design) stays on
screen while you edit; circulation is a first-class design choice
(thermosiphon / open stirring / guided serpentine plates) with a ν(T) curve
for the chosen fluid. Improve opens with a fully interactive Predictor
(what-if sliders, instant re-solve, delta cards). A new Ideas tab compares
concepts against the live baseline, seeded with the serpentine-plate idea:
plates between rows act as conduction fins (fin-efficiency model, contact
factor) and guide parallel manifolded channels driven by a ~0.5 W pump
(slot-flow pressure drop); forced flow reaches every tight gap and the
distributed sink cuts stratification. FEA findings now render as three
charts in Validate. Chart toolbars gain drawing tools; time charts use
unified hover. Learn gains a circulation-options panel.

## v7.1: BEV pack benchmark database
Compare now loads pack_benchmark.xlsx (user-supplied batterydesign.net
export, 58 BEV road-car packs): filterable energy-density map (Wh/kg vs
Wh/L, bubble = kWh, colour = 10 s C-rate), power-vs-energy chart with the
honest 10 s vs continuous caveat, ranked field vs this design, cell-to-pack
ratio comparison, full table with CSV export, percentile cards and a
nearest-neighbour summary. The FEA-honest enclosure variant plots as a
second marker. Keep pack_benchmark.xlsx next to app.py.

## v7.2: clearer heat-flow map and engineering triptych
Sankey rebuilt: single-line kW-in-label nodes with clean typography, short
casing branch (no cross-canvas ribbon), watt-labelled hovers. The live
design panel defaults to a three-view engineering drawing - plan plus two
true-scale sections - showing tubes as circles with their fin annuli,
serpentine plates, oil level and headspace, bottom gap, enclosure wall
thickness and outer dimensions; interstitial routing draws correctly. 3D
stays as an option, resized to fit the panel with a livelier camera and
lighting; plan view height also fits. All three views follow every design
change live.

## v8: the Live Pack and direct-manipulation UX
Results now opens on the Live Pack: a client-side animated cross-section
(pure canvas, no reruns) where oil particles circulate at the solved
velocity - thermosiphon loop, stirred, or guided serpentine lanes with
plates drawn - water beads travel and warm along the tubes, cells and cores
are coloured by their real solved temperatures with heat glow scaled by Q,
and the weakest link is tagged on the picture. Play/speed/layer/exaggerate
controls run in the browser instantly. The station wall is now a stepper:
click a box in the thermal circuit (the chart is selection-enabled) or pick
a chip, and that one station opens with its numbers and improvement moves;
the Sankey, ladder and set-point trade tuck into an expander. The transient
gains a time scrubber with a live readout card (can/core, C, SoC, oil,
heat at the picked instant). Design's three input groups are collapsed
expanders whose titles show the live key values. No API or internet is
required for any of this; livepack.py ships alongside app.py.

## v8.1: units in °C, richer Compare, plain-words Learn
Temperature differences now read in °C everywhere (spread, margins,
core-to-can, water rise, station tolls); compound units like W/m²·K and
kJ/K stay, and the Learn narrative explains that a 1 K difference equals a
1 °C difference. Compare's density map is grouped by energy class (City /
Mid / Large / Flagship) with legend on/off and double-click isolation, plus
toggleable median-crosshair and interquartile zones; a new explorer plots
any two database quantities with per-class colours, optional log axis and a
trend line, with this design overlaid where comparable; a head-to-head
radar duels up to four chosen packs against this design across five
normalised metrics (with the continuous-vs-10 s power caveat stated). Learn
opens with "Where the heat goes - the whole story in plain words": a
seven-step narrative using the live numbers, ending with what matters,
what doesn't, and the current biggest toll.

## v8.2: Live Pack hardened; native theme
The Live Pack component is rewritten to be self-diagnosing: any browser
error paints its own message onto the canvas and status line instead of a
blank panel, and the corner badges are pre-filled server-side so the solved
numbers show even if a browser blocks the script. The canvas is retina-safe
(CSS height + devicePixelRatio buffer scaling, ResizeObserver refit) and
the script passes an automated Node syntax + stub-DOM runtime harness in
the test suite. The toolbar is a floating translucent glass bar on a dark
canvas. A .streamlit/config.toml theme now ships (keep the folder next to
app.py) so native widgets adopt the brand palette on any Streamlit build,
and the tab-pill CSS carries selectors for multiple Streamlit versions.

## v8.3: benchmark comparison always visible; version badge
The header now shows the app version so it is always clear which bundle is
running. The BEV benchmark section no longer hides anything in collapsed
expanders: a sub-navigation (Density map / Power / Rankings / Explorer /
Head-to-head / Full table) shows every comparison full-size, with a green
DATA LOADED chip confirming the source (58 packs) or, if the file is
missing, a clear error plus an in-app upload fallback that works for the
session. Power and Explorer views are grouped by energy class like the map.

## v8.5: readable charts, compact header, gentle white theme
Chart titles and legends no longer collide: the global template puts every
title top-left and every horizontal legend top-right on its own band with
extra headroom, applied to all charts at once. The Ideas comparison puts
pump power on its own right-hand axis so half-watt bars are actually
visible. The oversized banner and the KPI card row are merged into one
compact header band (title, version, status chip, and seven small stats in
two lines). Background is now pure white with a gentler palette: softer
indigo/sky brand, lighter borders and shadows, calmer hero gradient, and a
matching .streamlit theme. The sticky live design panel gains stronger CSS
(column-level :has() selectors, overflow fixes, self-scrolling) so it stays
on screen while scrolling across more Streamlit builds.

## v8.6: thermal X-ray flagship figure; full cooling-architecture trade study
The waterfall is replaced by the Thermal X-ray: a schematic slice from cell
core to water drawn as physical strata (cell, both oil films, bulk oil,
fins and tube, wall, water film, water stream), with the live temperature
profile crossing them, every toll labelled in °C, mechanism and law on
hover, the weakest film flagged, the water inlet-outlet band shown, and the
limit line overlaid. It leads the Results detail and the report. Compare
gains a six-architecture trade study - full immersion static / serpentine /
stirred, partial immersion (60% fill, dry-top axial model), indirect bottom
cold plate (pad + plate + axial-conduction chain), forced air - each with
hotspot, spread vs the 5 °C criterion, max continuous C (bisection), pack
mass, parasitics, chiller electrical and thermal-system cost from editable
assumptions (cells anchored to BNEF Dec 2025: $79/kWh cell, $99/kWh BEV
pack, Europe +56%). Immersion rows use the full solver; dry rows are
closed-form (stated +/-20-30%). The report adds a comprehensive section:
the cooling requirement (duty and 4C-class), weight impact (thermal share
of pack), cost impact, the architecture table and honest verdict.

## v8.7: goal-seek and two-lever sweep actually plot
The two-lever heatmap now runs itself: pick any two different levers and
the 7 x 7 map renders immediately (cached, so revisits are instant), with
the limit contour and a star marking the current design; picking the same
lever twice now says so instead of silently showing nothing. Goal-seek
results persist on screen across reruns instead of vanishing at the next
interaction, and the insufficient-alone case points at the scale-to-C
table. Both expanders open by default.

## v8.8: the heat-flow map explains itself
A new deep-dive generator writes an extensive, structured commentary on
the Sankey from the live solved numbers - encoding and first-law balance
(with the actual ledger residual), a node-by-node analysis (cells with the
C² law and DCIR feedback quantified via real extra solves at 1C and 4C;
busbars; bulk oil as junction, mixer and flywheel with its computed time
constant; the water loop with UA, NTU and effectiveness; casing loss with
its effective UA and the do-not-insulate argument; the chiller as a
vapour-compression cycle with COP at three set points and the dry-cooler
threshold), an efficiency ledger (system cooling COP, overhead fraction,
fan-heater perspective), the key misreading warning (better cooling does
not widen ribbons - it lowers temperatures), and what the map deliberately
cannot show. It renders beneath the map in Results and as its own report
section; every number adapts to the design on screen.

## v8.9: the report rebuilt, in-app and exported
Version is now visible in the app header, the sidebar status card, the
report header and the export footer. The in-app Report opens with its own
header band (spec, version, date, status chip, reading time, key stats)
and a navigable contents bar: read the full report or one numbered section
at a time with previous/next controls. The HTML export is redesigned from
scratch: a proper markdown engine (tables now render as tables - the old
exporter printed raw pipes), a cover with status chips and a stats strip,
a sticky contents sidebar with scrollspy highlighting and smooth scroll,
styled callouts for every Improve action, hover-styled tables, a
back-to-top button, a print stylesheet (contents hidden, sections kept
whole), a responsive single-column layout on narrow screens, and a
provenance footer (Wang 2023, pack_fea_v1, the 58-pack database, BNEF
December 2025).

## v9.0: Live Pack becomes an instrument
The animated cross-section is rebuilt around interaction. Hover any part -
cells, bulk oil, the tube bank, the headspace, or the new chiller glyph
with its spinning fan and live duty/COP/electrical readout - and a tooltip
shows that station's live numbers; click to pin a full inspector card with
the physics, the current tolls in °C, and the levers that move them. The
weakest link is a red pulsing ring drawn on the actual hardware. The oil
flow field is real circulation now: counter-rotating buoyant rolls for
thermosiphon, one large roll when stirred, guided lanes with plates for
serpentine, with an Arrows toggle revealing the field. Water inlet and
outlet temperatures label the tube run; a live temperature legend sits
under the pack and stretches with the Exaggerate toggle; the layout
reserves space so nothing hides under the toolbar. Same self-diagnosing
error surface, retina safety and Node syntax + runtime harness in the
test gate.

## v9.1: the Cockpit - a pilot's panel for the whole system
A dedicated tab where every lever lives on rails around the live pack and
every readout responds instantly. The physics runs in the browser as a
line-faithful port of the app's solver: Andrade viscosity at film
temperature, Churchill-Chu films (vertical cells, horizontal tubes) with
the gap-confinement penalty and the Nu-cubed mixed-convection blend,
Schmidt fin efficiency and the serpentine plate-fin model, the
thermosiphon bisection, Hausen/Gnielinski water side with the 2300-3000
bridge, the DCIR-coupled 60-pass fixed point (busbars resized with
current), pump and stirrer and serpentine-channel laws, and the
45%-Carnot chiller. On load the surrogate is audited against the injected
full-Python solve and the deviation is shown on the panel chip; if it
ever exceeds 0.4 degrees C the chip flags for review.
Left rail (POWER): C-rate 0.2-6, ambient, dielectric fluid selector, cell
pitch with the live gap readout. Right rail (COOLING): water flow, inlet
temperature, tube count, circulation mode (thermosiphon / stirred /
serpentine) with conditional velocity and plate thickness/contact
controls, tube fins switch. Primary flight display: can/core with a
colour margin bar, margin to limit, spread vs the 5 degree criterion,
heat with the casing share, chiller duty/COP/electrical, parasitics, and
added plate mass. Annunciators: flashing LIMIT, LAMINAR, GAP<6, SPREAD>5,
and DRY-COOLER-OK. The canvas keeps the full hover/click station
inspectors, now regenerated from the live in-browser solve, with the
weakest-link ring, flow field, water in/out labels and temperature
legend all live. A copy button exports the settings as JSON; an apply
box under the panel lands them on the real Design widgets (the component
sandbox has no return channel, so the clipboard bridge is the honest
route). Disclosed simplifications: tube length and box size are held
fixed under the pitch lever, and exactness at the current design point
is guaranteed by the audit chip rather than assumed.

## v9.2: Cockpit takes the whole design space
build_geometry, enclosure_calc, build_masses and busbar_props are ported
to the browser alongside the physics, so pack architecture now honestly
reshapes everything live: cell format (18650/21700/4680 presets), series
and parallel counts, pitch and hex/square arrangement, tube plane, tube
OD/wall/material, and the loop fluid all rebuild the box, tube length,
areas, oil volume and pack mass in the panel - v9.1's fixed-geometry
simplification is gone. Rails are grouped and scrollable: POWER, PACK &
CELLS, GEOMETRY, FLUID and LIMITS on the left; WATER LOOP, TUBES & FINS
and CIRCULATION on the right. New instruments: a live temperature ladder
(the X-ray's essence as a stacked strip), a flight-data-recorder trace of
the governing temperature against the limit line, SNAP reference deltas
on the key readouts, a MAX-C autothrottle (in-browser bisection,
debounced), TRIM TURB (closed-form minimum flow for Re >= 3100), DRY SET,
and four scenario presets. The PFD adds max continuous C, energy, pack
mass, and oil tiles. Offline fidelity measurement with the real payload:
the audit chip reads +0.00 degrees C against the full Python solve at the
design point; mass and oil tiles agree with build_masses. The apply
bridge now lands ~24 keys on the real Design widgets including Ns, Np,
capacity, DCIR and its temperature fall, tube spec, tube plane,
arrangement, loop fluid, limit and h_ext (format preset stays
cockpit-only since it drives several linked widgets at once).

## v9.3: Zonal plate-channel model with Monte Carlo (new Zones tab)
The plate-channel architecture solved as a full thermal-hydraulic
network with no hand-waved constants. fea4_channel.py computes exact
developing-flow kernels on the true lens-shaped slot cross-section
(velocity Poisson for fRe; two-case Graetz march with isothermal walls
giving the 2x2 transport matrix a(z*, s)), validated against
parallel-plate limits: fRe 95.7 vs 96.0, one-wall Nu 5.387 vs 5.385,
two-wall Nu 7.541 vs 7.541. An earlier stagnant-slot study returned an
arc effectiveness of ~0.13, which redirected the concept honestly: the
oil is the primary carrier (the whole 70 mm channel operates in the
thermal entrance at Pr ~ 127, cell film ~3x the thermosiphon), and the
plates are the secondary sink plus structure. zonal.py assembles every
slot, cell, plate bay, root and water crossing: exact exponential
z-march per channel, per-cell Newton with the march's own Jacobian and
DCIR feedback, plate faces solved in closed form against the
contact + wall + water-film root chain, water marched tube by tube,
and the recirculation plenum solved from the linear exit map (removing
a slow mode that stalled convergence). Tube count and positions are
derived from the fin rule P = 2/m snapped to the cell lattice (17
tubes at the defaults, one every two cells, bay efficiency 0.64).
Validation: energy closure 0.19% full-size, z-resolution independence
0.02 degrees C, and an independent cross-check against the lumped
serpentine solver (zonal 41.2 vs lumped 40.9 degrees C at 2C).
Findings: the flat-slot cubic law does not apply to the lens duct
(flow sensitivity to slot width is nearly linear, exponent ~0.7), so
with plate homogenisation the design is remarkably tolerance-robust;
the governing resistance is the water film at the roots. The Zones tab
exposes slot width, crevice, pump head, C-rate and three heat maps,
renders the per-cell temperature map with derived tube lines, channel
flow shares and per-tube water outlets, and runs a warm-started
tolerance Monte Carlo with exceedance probability. Smoke gains zonal
gates; all legacy tests green.

## v9.3 explanation layer: equations beside every result
The physics is no longer only in the code. The app went from 2 rendered
equations to 21. The Zones tab gained a "Methods, physics and how every
number is reached" panel that walks the full derivation with LaTeX
(layout fin rule, laminar-vs-buoyancy hydraulics, the fea4 2x2 duct
kernel, the exact exponential energy march, DCIR-coupled cell Newton,
the closed-form plate-face and root network, the closed-form plenum
fixed point, and the Monte Carlo), each equation sitting beside its live
value and a one-line account of how it is obtained, plus a live
validation battery. Three interactive labelled schematics
(schematics.py) give the geometry symbols and the solve topology a
picture: the plan-view plate-channel unit, the lens subchannel fea4
meshes, and the node-and-edge solve graph. The Learn tab panels 4-9 now
each show their governing equation (gap-velocity scaling, fin/Schmidt,
thermosiphon head balance, laminar/turbulent tube Nu, thermal-flywheel
time constant, DCIR and core spreading). FEA4 (the duct-kernel gates:
fRe 95.7 vs 96, one-wall Nu 5.387 vs 5.385, two-wall Nu 7.541 vs 7.541)
now sits with FEA1-3 in the Validate tab instead of being buried in a
caption. The Cockpit gained a panel showing the ported resistance-chain
equations it evaluates. A standalone REDTEAM_HANDOVER.md documents the
whole physics stack with equations, the validation evidence classified
for scrutiny, an assumptions register, the known soft spots, and a
structured red-team checklist.

## v9.4 red-team closure: the corrected, honest model
v9.4 folds in a full external red-team pass on the zonal solver and Monte
Carlo (three bugs, two oversold claims, a doc defect, a packaging break).
The point-by-point response is in REDTEAM_RESPONSE.md. In short:

- **Buoyancy bug (F1):** the buoyant head is now referenced to the loop
  return (plenum) temperature, not the water inlet. At the default this
  removes ~9 Pa of spurious head; the solved slot velocity drops from
  166 to ~122-162 mm/s and the pack runs hotter and more conservatively.
- **Water-film cliff (F5):** the zonal now shares one water-side Nusselt
  function with the lumped solver (correlations.water_nu), with the
  continuous 2300-3000 transition bridge; the -3.3 C step is gone.
- **Contact conductance (F2):** h_contact and the collar engagement
  factor are named inputs with a one-click sensitivity sweep showing
  where the design crosses 45 C (below h_c ~ 4000-6000 W/m2K).
- **Self-consistent layout (F4):** the fin-rule tube pitch is now fed the
  kernel-implied plate film (a11/pitch at the design point), not a
  hardwired 150; the derived h_face and the conservative product-rule eta
  are shown.
- **Matched cross-check (F3/V8):** the lumped comparison is solved in the
  same serpentine architecture at the zonal's own velocity (~35.8 vs
  ~41.9 C), and the misleading "agree within a fraction of a degree" line
  is gone. The direct oil-to-tube bypass the model omitted is now an
  explicit, bounded toggle (off by default, conservative; ~21 W/K and
  ~15% of duty at wetted fraction 0.5).
- **Honest Monte Carlo (F6):** tighter per-sample tolerance, the sample
  spread labelled as convergence noise (< 0.1 C), and a rule-of-three 95%
  upper bound quoted for zero-exceedance runs.
- **Core temperature (F10):** the jellyroll core-to-can rise is superposed
  and the peak CORE temperature reported next to the can (~44.5 vs 45 C at
  the default) - so a core/plating limit is shown to be marginal.
- **Packaging + honesty (F7-F9):** scipy added to requirements; zonal.py
  ships a __main__ shakedown; the fea4 __main__ regime-difference line is
  relabelled; the Methods preamble names its three assumptions instead of
  claiming nothing is asserted; the handover's operating-point narrative
  is corrected.

The net effect is deliberately less flattering than v9.3 and more true:
the design is marginal, and its margin depends on the unmeasured collar
contact and on whether 45 C is a can or a core limit. The recommended
next step is the one-tube two-crossing oil-bath rig point, which measures
both the omitted oil-to-tube path and the contact conductance at once.

## v9.5 final: verification round closed
An independent verification pass signed off every v9.4 fix (all ten
implemented, numbers reproduce, no regressions, requirements installs
one) and raised four minor, non-blocking observations plus answers to the
two open questions. v9.5 folds all four in (see VERIFICATION_RESPONSE.md):

- **Core-node exceedance in the Monte Carlo (obs 3):** monte_carlo now
  returns the core-based exceedance count, and the Zones tab shows it
  beside the can-based one - a green confirmation when the core clears,
  a warning when it does not (the core is the razor-thin node).
- **Node-dependent h_c crossing (obs 2):** the h_c sweep now interpolates
  and reports both crossings - the can clears down to ~2900 W/m2K but the
  core only above ~7500 (within ~7% of the nominal 8000), so on a
  core/plating criterion the braze must be essentially as good as assumed.
- **Layout snap asserted (obs 4):** h_face_design is bank-dependent
  (a11 drifts a few percent with kernel resolution); the shakedown now
  asserts the tube count and prints the snap margin (~26-30% headroom to
  the 33-tube boundary at h_face ~ 333 W/m2K), so a bank change that moves
  the snap fails loudly. A snap_boundary_h() helper documents the edges.
- **Stale memo numbers corrected (obs 1):** the shipped default is
  can ~42.4 / core ~44.9 C (a +1.15 C shift over v9.3, from F1 and F4).

The two open questions are now surfaced in the UI: completing both models
with each other's missing physics converges near 38-40 C (so the
bypass-off default is ~2-4 C conservative), and both temperature nodes are
displayed so the can-vs-core specification decision cannot hide. The
remaining uncertainty is measurement, not modelling: the one-tube,
two-crossing oil-bath rig point pins h_c and the wetted-tube film and
settles whether the core clears 45 C.

## v9.6 Learn tab rebuilt: equations woven into the physics, and it moves
The Learn tab was rebuilt from a stack of collapsed panels (equation in
one place, prose in another) into a single flowing narrative where every
equation sits inside the sentence that explains it, using inline LaTeX.
The through-line is that the heat-transfer coefficient h is a *calculated
result*, not a fixed property: the tab derives h = Nu.k/L from the
boundary-film picture (h = k/delta), shows Ra -> Nu -> h with Churchill-
Chu, and drives the whole chain from live sliders (length L, temperature
difference, and fluid) so the reader watches each number flow through the
equations to the h their surface actually gets - including that a taller
surface has a *lower* h (h ~ L^-1/4).

New interactive pieces:
- An animated in-browser boundary-layer explorer (HTML/JS canvas): warm
  oil rises past a hot cell wall, the no-slip film sets the temperature
  drop, and a flow slider thins the film and updates h = k/delta live.
- A slider-driven film-temperature profile that shows h = k/delta as a
  graph.
- The h-calculator: length/dT/fluid sliders that print Ra, Nu and h as
  fully substituted equations, plus an h-versus-length curve.
- A stirring slider that redraws the whole core-to-water temperature
  stack, showing the two oil films shrink as the flow thins them.
The design levers (viscosity, gap, stirring, fins), the water-side
laminar plateau, the DCIR feedback and core resistance, and the thermal-
flywheel buffering are all retained but reworked to flow with their
equations inline; reference material (moving the oil, safety, the AMG
benchmark, sources) sits in one expander at the end.

Also in v9.6: the three v9.5 sign-off nits are closed - the Monte Carlo
core-exceedance warning no longer hardcodes a false "the can clears"
clause when h_c is low (it now reports the real can count), the
h_face_design docstring says "tube pitch" explicitly, the shakedown's
core-crossing figure is harmonised with the live app, and the response
memos' last stale numbers are corrected. A separate pre-existing render
bug (currency dollar signs in the Report tab parsed as LaTeX math) is
fixed by using USD wording.

## v9.7 System & BOM tab, a 3D serpentine view, cockpit reset, sharper Learn intro
Four things this release.

**New "System & BOM" tab.** Turns the thermal answer into hardware. It
sizes every component from the live design - water pump (flow, head,
hydraulic and electrical power), oil circulation, the dielectric-to-water
heat exchanger (UA = Q/ΔT and area from an editable U), the chiller/dry-
cooler (reject duty, COP, electrical), coolant volume/mass, tubes,
plates, and an expansion vessel sized from the fluid's thermal expansion
over an editable service band. It prints a full bill of materials with
per-item cost, mass and notes, and totals as a percentage of cell cost.
The honest finding it surfaces: for full immersion the flooded coolant
volume dominates the thermal bill (here the dielectric alone is a large
fraction of cell cost), which is single-phase immersion's real economic
penalty versus a cold plate and far worse for two-phase. A dedicated
two-phase section covers the coolant chemistry (fluoroketone FK-5-1-12 as
the low-GWP front-runner, why the PFAS hydrofluoroethers are being retired,
low-GWP refrigerant options), materials compatibility (FKM/FFKM seals,
avoid EPDM/polycarbonate), the extra items a boiling system needs
(condenser, vapour-tight enclosure, pressure management, often no pump),
and when two-phase wins or loses - plus a two-phase BOM delta and an
architecture-comparison table (air / cold plate / single-phase immersion /
two-phase) on surface h, pumping, parts, cost and maturity, and an
applications map. Cost and sizing assumptions are all editable inline. The
Report tab's cost-scale note is corrected to match (it previously
understated the coolant-dominated cost).

**3D serpentine view in the cockpit.** A drag-to-rotate, self-contained 3D
visualization of the architecture the flat instruments cannot show: cells
in a grid, a metal cooling plate in each gap, water tubes bonded along the
plates running through-plane (into the page), and coolant flowing front to
back and warming as it collects heat - the real cell → oil film → plate →
tube → water path. It is built from the live geometry and temperatures,
with an auto-spin toggle and flow/explode sliders, and hand-rolled
orthographic 3D so it needs no external libraries.

**Cockpit reset button.** A Reset control restores the cockpit to its
default state.

**Learn intro.** The opening now frames this specific problem - a high
energy-density 21700 pack rejecting kilowatts into a poor-convecting
dielectric, where two stagnant oil films gate everything and the whole
architecture is an assault on those films - rather than a generic preamble.

Note on approach-specific cockpit parameters: the cockpit already hides
and shows the circulation controls by mode (velocity is hidden for a pure
thermosiphon; plate thickness and plate contact appear only for the
serpentine mode). That behaviour is retained and can be extended to more
parameters on request.
