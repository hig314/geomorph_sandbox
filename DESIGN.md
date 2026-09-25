# geomorph_sandbox — design and plan

*Status: proposal for discussion (2026-09-25). Nothing below is built yet except the two
legacy projects it grows from.*

## 1. Purpose: illustrate hypotheses, not test them

A browser sandbox whose job is to make geomorphic hypotheses **visible and comparable**.
Output need not be robust enough to *test* a hypothesis, but it must illustrate it clearly:
the same forcing, the same seed, two mechanisms, side by side.

The first target is the origin of **low-angle uplands in glaciated mountain ranges**:

- **Hypothesis A — planation under deep ice.** Uplands are erosional. Deep burial in an ice
  cap erodes high convexities and leaves planed-off surfaces lower down.
- **Hypothesis B — uplifted lowland convexities.** Uplands are inherited. Uplift promotes
  low-angle lowland convexities into low-angle highlands. This requires erosion resistance
  that *persists across uplift-relevant timescales* (a resistant body thick enough not to be
  stripped during uplift). The uplands are ultimately destroyed by **marginal contraction**,
  e.g. landslides eating in from the plateau edge.

Illustrating both needs three ingredients that the legacy projects lack or only partly have:
**strength heterogeneity** (in a material frame that moves with uplift), **approximately
realistic glacial erosion** (valleys carved into uniformly uplifting bedrock as the
foundation), and **landsliding** at plateau margins.

## 2. The dimensional ladder

The project deliberately implements the same physics in **1, 2 and 3 spatial dimensions**
(counting the vertical), as parallel tiers:

| Tier | Coordinates | What it is | Legacy source |
|---|---|---|---|
| **d1** | z(t) | a column / reservoir | `uplift_sandbox` (relief-reservoir model) |
| **d2_across** | z(x,t) | a cross-valley profile | new |
| **d2_along** | z(s,t) | an along-valley (long) profile | new |
| **d3** | z(x,y,t) | a DEM | `terrain_sandbox` |

The organising device is **asserted vs computed**. Every process law is written *once* as a
pure function of local quantities (in `core/`). Each tier supplies those quantities either by
computing them or by **asserting** them as a parameter for the dimension it does not resolve.
The UI labels every asserted quantity as such. The ladder for ice, for example:

| Quantity | d1 | d2_along | d2_across | d3 |
|---|---|---|---|---|
| ice thickness H | asserted | computed (flow-law inversion) | computed | computed |
| ice flux q | asserted | computed (∫ balance × asserted width) | asserted | computed (routed) |
| sliding speed U_s | asserted | computed | computed (transverse distribution) | computed |
| valley width | — | asserted W(s) | resolved | resolved |
| surface slope | asserted | computed | asserted (along-valley) | computed |

Each tier corresponds to a published model class, which keeps the illustrations grounded:
d2_along = MacGregor et al. (2000) / Anderson et al. (2006) flowline glaciers and stream-power
long profiles with Hack's-law area; d2_across = Harbor (1992) U-shape shaping and Bishop's
2D method of slices for landslides; d3 = the existing coupled DEM model (see
`legacy/terrain_sandbox/DESIGN.md`).

**Why the ladder helps the hard problem.** Steady ice flux in a flowline is the integral of
mass balance times width; thickness comes from inverting the flow law locally. The d3 glacier
rebuild (§4) is exactly that, with the flux obtained by *routing* mass balance over the ice
surface instead of integrating along a line. So d2_along is where the thickness inversion
gets validated against a reference SIA, and d3 only adds routing.

## 3. How the two hypotheses appear at each tier

- **d1 column.** Two columns side by side, resistant and weak, under the same uplift.
  Hypothesis B reduces to a timescale statement: (resistant-body thickness ÷ erosion-rate
  contrast) vs uplift duration. Hypothesis A appears as an elevation-dependent glacial
  erosion rate peaking near the ELA (buzzsaw). The relief-reservoir model (dR/dt = U − kRⁿ)
  stays as one d1 scenario.
- **d2_across.** Plateau-margin retreat by landsliding (destruction step of B). Ice-cap burial
  of a ridge crest with a quarrying term on convexity (planation step of A). U-shape carving.
- **d2_along.** Buzzsaw planation at the ELA, overdeepenings, hanging tributaries, and a
  resistant band in the long profile.
- **d3.** The full coupling.

Diagnostic that discriminates A from B in every tier: **slope vs elevation** (where the
low-angle area sits in the hypsometry, and how it moves through time).

## 4. Process laws (shared, `core/`)

- **Forcing.** U(t) and ELA(t) from the d1 forcing library (bell / arc / constant / triangle /
  plateau / pulse), plus spatial uplift patterns (uniform / ramp / gaussian / tilt).
- **Lithology in a material frame.** Rock moves *up through* the surface, so strength is
  defined on material coordinates: `hard = f(x, y, z − U_cum(x, y))`, evaluated per cell each
  step from stored cumulative uplift. `hard` multiplies fluvial K, glacial K_g, creep D,
  cohesion c and tan φ (each with its own sensitivity). Generators: horizontal resistant
  layer, tilted slab, vertical dike, pluton blob (thresholded 3D noise), persistent 2D patch
  map. Rendered as an overlay. This is what makes "persistence across uplift timescales"
  visible: a thin cap is stripped, a thick body survives.
- **Fluvial.** Stream power with erosion–deposition (Yuan et al. 2019, implicit); d2_along
  asserts drainage area by Hack's law.
- **Hillslopes.** Linear + nonlinear (Roering) creep; threshold-slope rockfall as a cheap
  margin-retreat process.
- **Landslides.** Limit-equilibrium failure: 2D method of slices (d2_across), 3D
  method-of-columns scoops (d3, exists). Strength sampled from the lithology field along the
  slip surface. Catastrophic vs slow-slump modes (exist in d3).
- **Glaciers (rebuilt).** Steady ice discharge, not time-stepped SIA:
  1. Route surface mass balance (positive above the ELA, negative below) over the ice surface
     with the existing priority-flood + MFD machinery (d3) or integrate along the flowline
     (d2_along). Flux clamps at zero where ablation exhausts it → terminus.
  2. Invert the SIA flux law locally for thickness from flux and surface slope; iterate a few
     times per frame (surface depends on thickness). Divides need a slope floor and
     under-relaxation — the one numerical risk.
  3. Erosion ∝ sliding speed (E = K_g·U_sˡ, U_s = f_s·q/H), the standard published form.
     Overdeepenings self-limit (filling flattens the surface, thickens the ice, slows it).
     Keep a rate cap as a safety.
  4. Optional **quarrying ∝ U_s × bed convexity** — the knob hypothesis A needs.

  *Why rebuild:* the legacy d3 glacier is an explicit SIA with a fixed diffusivity cap
  (8000 m²/yr on the coarse grid vs 10⁵–10⁷ m²/yr for real valley glaciers). Velocities are
  throttled by 1–3 orders of magnitude, so velocity-based erosion was impossible and the
  fallback (E ∝ basal shear stress, nearly uniform under ice) cannot concentrate erosion in
  trunks or leave tributaries hanging. Erosion was also upsampled from 200 m cells. The
  legacy solver stays as the **Python reference** for the bake-off. The hard-won numerics in
  `legacy/terrain_sandbox/CLAUDE.md` remain the list of things not to retry.

  *Decision gate (Python):* on a cone and a synthetic valley network, thickness within ~30 %
  of the reference SIA; trunk erosion > tributary erosion; U-shaped sections; no runaway at
  K_g × 100. Fallback if it fails: keep the SIA and move it to a worker / WebGPU to drop the
  cap (much bigger job).

## 5. Repository layout

```
geomorph_sandbox/
  core/        forcing functions, process laws, lithology, units   (pure JS, no DOM)
  d1/          column model
  d2_across/   cross-valley profile
  d2_along/    along-valley profile
  d3/          DEM (terrain_sandbox, migrated)
  ui/          shared D3 panels, URL state, presets, A/B lockstep split view
  py/          Python validation mirrors (every scheme, before and after porting)
  docs/        methods (cited references), scenario notes
  legacy/      uplift_sandbox/, terrain_sandbox/ — frozen historic record, recycled from
```

Static bundle, no build step, so any tier drops into a landslidescience.org Django page.

## 6. Build order

1. **Scaffold.** Repo, move legacy in, factor `core/` (forcing, law signatures, lithology) and
   shared `ui/` from the legacy code. Domain scale becomes a parameter (cell size), so d3 can
   run 40 km at 100 m as well as 20 km at 50 m.
2. **d1 upgrades.** Lithology persistence, elevation-dependent glacial rate, two-column A/B.
   Cheap; forces the law signatures to be right.
3. **d2_along.** Stream power with Hack's law; flowline glacier with the thickness inversion
   validated in Python against a reference SIA; sliding-based erosion; quarrying on convexity.
   *De-risks the d3 glacier.*
4. **d2_across.** Harbor-style shaping; creep; threshold slopes; 2D method of slices;
   plateau-edge retreat.
5. **d3.** Swap the glacier to the routed steady-discharge model; per-cell strength fields;
   scoops sampling strength along the slip surface.
6. **Illustration layer (all tiers).** A/B presets, time-forced ELA, lockstep A/B split view
   sharing one seed, slope-vs-elevation diagnostic, "asserted" badges in the UI.

## 7. Open questions / tweaks to discuss

- d2_along before d2_across (as above) or the reverse?
- Domain size and resolution defaults for d3 under an ice cap.
- Whether d1 should be two columns (A/B) or one column with a scenario switch.
- Quarrying-on-convexity: include from the start, or only after the abrasion-only model
  carves convincing valleys?
- Which tiers get the lockstep A/B view first (cheapest in d1 and d2; d3 doubles the cost).
