# Thermal management — interactive tutorials

A small static site of interactive lessons on heat transfer and cooling. No build
step, no frameworks, no dependencies — every page is plain HTML/CSS/JS and works by
opening it in a browser or serving it from GitHub Pages.

## Contents

| File | What it is |
|------|------------|
| `index.html` | Landing page / curriculum map |
| `style.css` | Shared design system (light + dark). Every page links it |
| `coolants.js` | Shared coolant dataset + `prandtl()` and `mouromtseff()` helpers |
| `01-liquid-cooling.html` | Tool: full cold-plate temperature-budget simulator |
| `02-materials.html` | Conduction, material properties, layered-stack calculator |
| `03-interfaces.html` | Boundary conditions, boundary layers, contact resistance |
| `04-coolants.html` | Sortable/filterable coolant catalogue with figure of merit |

## Run it locally

Just open `index.html`. For correct loading of `coolants.js` on tutorial 04, serve
the folder instead of opening from `file://`:

```bash
python3 -m http.server 8000   # then visit http://localhost:8000
```

## Publish on GitHub Pages

1. Create a repo and add every file at its root.
2. `git add . && git commit -m "thermal tutorials" && git push`
3. Repo → Settings → Pages → Source: Deploy from a branch → `main` / root.
4. Live at `https://<you>.github.io/<repo>/`.

## Add a new tutorial

1. Copy `02-materials.html` to e.g. `05-convection.html` — it has the full page
   skeleton (spectrum rule, breadcrumb, header, sections, interactive `.panel`, footer).
2. Keep `<link rel="stylesheet" href="style.css">`. Reuse existing classes:
   `.panel .row .chip .track .seg .legend .grid .card .callout .eq .dtable .pill`.
   Don't invent new visual styles unless a component is genuinely missing.
3. Add a card for it in `index.html` under the right section, changing the `soon`
   tag to `live`.

## Add a coolant

Append one object to the array in `coolants.js`:

```js
{name:"New fluid", cat:"oil", rho:900, cp:2000, k:0.14, mu:0.01,
 bp:250, fp:-40, diel:true, use:"One-line description and tradeoff."}
```

`cat` is one of `water | glycol | oil | dielectric | metal | gas`. `diel` is
`true`, `false`, or `"part"`. Prandtl and the Mouromtseff figure of merit are
computed automatically; tutorial 04 picks it up with no other change.

## Design tokens (in `style.css`)

- Colours are CSS variables that auto-switch for dark mode: `--ink`, `--ink-2`,
  `--bg`, `--surface`, `--surface-2`, `--line`, plus status `--ok/-warn/-bad-t` and
  `-bg`, and the heat ramp `--blue --teal --amber --coral --red`.
- Type: `--sans` for prose, `--mono` for data, equations, and labels.
- The 3 px spectrum rule at the top of each page is the visual signature — keep it.

## Accuracy note

All property values are representative and for teaching. Real coolant and material
properties vary strongly with temperature, grade, and additives; the simulators use
simplified 1-D models. Use them to build intuition, not to size hardware.
