# core/ and ui/ — shared modules

*Status: step 1 of the build order (DESIGN.md §6) landed 2026-09-25.*

## Unit convention

Every tier computes in **metres and years**. Rates are m/yr, so 1 mm/yr is 1e-3 m/yr;
lengths are m; times are yr. Only the UI converts, with `GS.units` (`fromMmyr`, `toKm`,
`fmtTime`, …), and only at the boundary between controls/panels and the model. The two
legacy projects disagreed (uplift_sandbox: km and Myr; terrain_sandbox: m and yr); the
d1 port is verified against the legacy km/Myr series to 1e-12 m in `py/test_core.py`.

Compute stays float32 on grids (d2, d3). `py/test_core.py` reruns the d1 RK4 loop in
float32 and finds a 1 mm difference over 24 Myr against float64.

## Module style

No build step. Each file is an IIFE that attaches to the `window.GS` namespace (and
`module.exports` for node/tests). Load order matters only where noted:

| File | Namespace | Depends on |
|---|---|---|
| `core/units.js` | `GS.units` | — |
| `core/forcing.js` | `GS.forcing` | — |
| `core/noise.js` | `GS.noise` | — |
| `core/lithology.js` | `GS.lithology` | `GS.noise` |
| `core/laws.js` | `GS.laws` | — |
| `ui/urlstate.js` | `GS.ui.urlState` | — |
| `ui/controls.js` | `GS.ui.buildControls`, `GS.ui.toUrlSpecs` | d3 |
| `ui/panels.js` | `GS.ui.panels` | d3 |
| `ui/tooltip.js` | `GS.ui.tip` | d3 (load before controls/panels are built) |
| `ui/styles.css` | shared look, `.badge.asserted` | — |

Python mirrors: `py/gs_core.py` (all of core, function for function) and
`py/d1_reservoir.py`. Validation: `/opt/anaconda3/bin/python3 py/test_core.py`.

## core/forcing.js

- **Time shapes** `shapes[name](t, p)` in [0, 1]: event shapes `bell`, `arc`, `constant`,
  `triangle`, `plateau`, `pulse` on `[0, p.duration]` (ported from uplift_sandbox with the
  peak factored out); cyclic `sine`, `sawtooth` (90 % build, 10 % collapse) on `p.period`
  for ELA glacial cycles.
- `makeSeries({shape, peak, base, duration, period})` → `f(t) = base + peak·shape(t)`.
  The same call builds U(t) (base 0) and ELA(t) (base = interglacial ELA, peak negative).
- `integrate(f, t0, t1, steps)` trapezoid, for cumulative uplift.
- **Spatial patterns** `spatialFactor(u, v, pattern, ramp)` in normalised coordinates:
  `uniform`, `ramp`, `gaussian` (σ = 0.18 of the domain), `tilt`. Ported from
  terrain_sandbox `upliftFactor` (which was in cells). A profile tier passes `v = 0.5`.

Checked: closed-form integrals of every event shape (bell = √(2π)/6·erf(3/√2) ≈ 0.4166 of
peak × duration), boundedness, zero outside the window (bell tail 0.011 at the edge by
construction), periodicity.

## core/noise.js

`mulberry32`, `hash2`/`hash3`, `valueNoise2`/`valueNoise3`, `fbm2`/`fbm3`; `valueNoise` and
`fbm` are aliases of the 2-D forms so terrain_sandbox code ports unchanged. The 3-D forms
exist for pluton blobs. The Python mirror reproduces the 32-bit `Math.imul` arithmetic.

## core/lithology.js — strength in a material frame

Rock moves up through the surface, so bodies are defined on material height
`zm = z − U_cum(x, y)` (`materialZ`). Each tier stores U_cum (scalar / profile / field) and
evaluates the field at the surface every step. Under constant uplift U and erosion E the
surface's material height falls at rate E regardless of U, so a layer buried d below the
initial surface is exposed at t = d/E and stripped at (d + T)/E (tested).

Generators return `r(x, y, zm)` in [0, 1] (0 background, 1 resistant, optional `soft`
transition): `layer {top, thick}`, `slab {top, thick, dipX, dipY}`, `dike {x0, width}`,
`blob {seed, scale, threshold}` (thresholded 3-D fBm), `patch {seed, scale, threshold}`
(2-D, independent of zm: the infinitely thick end-member). `make([...])` unions bodies.

Laws consume r through two multipliers with a per-process contrast c ≥ 1:
`erodibilityFactor(r, c) = 1/(1 + r(c − 1))` for K, K_g, D and
`strengthFactor(r, c) = 1 + r(c − 1)` for cohesion and tan φ.

## core/laws.js — the process laws, once

| Law | Signature | Source |
|---|---|---|
| reservoir | `reservoirK(Upeak, Rlim, n)`, `reservoirErosion(R, k, n)`, `reservoirEquilibrium(U, k, n)` | uplift_sandbox |
| stream power | `streamPower(K, A, S, m, n)` | Whipple & Tucker 1999 |
| erosion–deposition node | `yuanNode(h0, hcur, hr, Kp, dt, dist, dep, nexp, submerged)` | Yuan et al. 2019; Braun & Willett 2013 for G = 0. Extracted verbatim from terrain_sandbox `fluvialStep`; d2_along and d3 both call it |
| creep flux | `creepFluxLinear(D, S)`, `creepFluxRoering(D, S, Sc, fmin)` | Roering et al. 2001 |
| mass balance | `massBalance(zs, ela, bg, accMax, abMax, ablRatio)` | terrain_sandbox (ablation gradient 2.5× accumulation) |
| SIA flux | `siaFlux(H, S, Γ, n)` = Γ H^(n+2) Sⁿ | Glen n = 3, A = 7.57e-17 Pa⁻³ yr⁻¹, Γ = 2A/(n+2)·(ρg)ⁿ |
| SIA thickness inversion | `siaThickness(q, S, Γ, n, Smin)` = (q/(Γ Sⁿ))^(1/(n+2)), slope floored at Smin | MacGregor et al. 2000 / Anderson et al. 2006 flowline form |
| sliding | `slidingSpeed(q, H, fs)` = fs·q/H | |
| abrasion | `abrasion(Kg, Us, l)` = Kg·Us^l | Hallet 1979; l = 1 after Humphrey & Raymond 1994 |
| quarrying | `quarrying(Kq, Us, convexity)` = Kq·Us·max(convexity, 0) | the knob for hypothesis A |
| buzzsaw (asserted ice) | `buzzsawRate(z, ela, width, peak)` Gaussian window round the ELA | d1 only |
| integrator | `rk4(f, t, y, dt)` | |

Landslide laws (Bishop's 2-D method of slices, 3-D method of columns) are added with
d2_across; the d3 scoop code stays in terrain_sandbox until then.

Checked: `yuanNode` satisfies its implicit equation to 1e-12 for n = 1, 1.5, 2, matches
explicit Euler at small dt, approaches the receiver at huge dt, deposits in basins;
the SIA inversion round-trips `siaFlux` to 1e-9 and reproduces the textbook
depth-averaged deformation speed (≈5 m/yr for 300 m ice on 3 %, ≈70 m/yr for 400 m on 5 %).

### A note on the legacy "don't retry" list

`legacy/terrain_sandbox/CLAUDE.md` lists erosion from |q|/H as a failure. It failed there
because the explicit SIA's diffusivity cap flattened q, so q/H carried no channel
information. The steady-discharge rebuild (DESIGN.md §4) routes q directly and has no cap,
so `slidingSpeed` = fs·q/H is the intended sliding law and that entry does not apply to it.
The rest of that list (surface-slope stress, slope cap, rate cap, dt-independent cap for
the reference SIA) still applies.

### Python reference for the glacier bake-off

DESIGN.md §4 says the legacy SIA "stays as the Python reference". No Python mirror of
`glacierStep` was ever checked in, so that reference has to be written from the JS at the
start of step 3 (d2_along), as a flowline SIA first.

## ui/

- `urlState(specs).read()/write(state)`: every control round-trips through the query
  string (clamped to its range / option list). `toUrlSpecs` derives the specs from the
  control specs so a tier declares each control once.
- `buildControls(root, specs, state, onChange)`: sliders, selects, checkboxes grouped in
  fieldsets. `asserted: true` adds the badge (DESIGN.md §2). `log: true` marks a slider
  holding log10 of the quantity (ids end in `Exp`, the terrain_sandbox convention).
- `panels({width, height})`: stacked time-series panels (`make`, `line`, `area`, `clear`,
  `endLabel`, `vline`, `zero`, `axes`), the uplift_sandbox charts generalised with an x
  accessor. `make` takes `method` for the title tooltip.
- `tooltip.js`: one floating methodology box. Every control spec's `method` HTML (equation,
  asserted/computed tag, source) shows on its ⓘ glyph; every panel's on its title. This is
  a standing practice: no control or panel without a `method`.
