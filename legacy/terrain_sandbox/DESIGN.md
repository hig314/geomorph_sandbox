# 2D Terrain Evolution Sandbox — Design

A browser-based, interactive landscape-evolution sandbox: a parallel to the 1D
uplift/erosion sandbox, but on a 2D DEM. Same "sandbox feel" — sliders, live
plots, shareable URL state — scaled up to a gridded surface and several coupled
geomorphic processes.

## Goals & scope

- Grid: ~**400×400** cells (160k), interactive framerates in-browser.
- **Compute in float32**, not 1 byte. 8-bit elevation (256 levels) is far too
  coarse for *incremental* erosion — a sub-metre-per-step change would round to
  zero and nothing would evolve. The float field is the source of truth; "1 byte"
  is only a **storage/display/exchange** format (quantize on export, or 16-bit
  PNG packing for better round-trips).
- Not a physically exact model — a *plausible, explorable* one. Correct
  qualitative behavior (relief limits, U-shapes, hanging valleys, episodic
  failures) over perfect physics.
- Reuse from the 1D sandbox: time-forcing functions for `U(t)` and `ELA(t)`,
  D3 plot panels, URL query-string state encoding, the model-core/render/UI
  separation. Static bundle so it drops into a landslidescience.org Django
  template like the 1D one.

## Architecture

```
terrain-model.js   pure core: grid, topo generators, process kernels, stepping.
                   NO DOM, NO D3. Runnable in a worker or node test.
render.js          DEM -> canvas (hillshade + hypsometric colormap).
app.js             controls, animation loop, cross-section (D3), URL state.
index.html / styles.css
```

The process kernels are written as swappable units. The first implementation is
a **CPU reference** (160k cells is tractable in JS at interactive rates and is
far easier to debug); the **hybrid GPU plan** (below) ports the parallel-friendly
kernels to WebGPU while keeping flow routing on the CPU.

## Compute model: hybrid (chosen)

- **WebGPU compute** (ping-pong float32 storage buffers) for the embarrassingly
  parallel / stencil work: uplift, hillslope creep, glacial swath erosion, and
  especially the **SCOOPS3D-style scoop search** (the big GPU win — millions of
  trial surfaces).
- **CPU** for the one inherently sequential step: **flow routing** (steepest
  descent receivers → topological stack → drainage-area accumulation → implicit
  stream-power incision). ~1–2 ms for 160k nodes; run every N steps and upload.
- CPU fallback for browsers without WebGPU.

WebGPU is mild overkill for stencils at 160k cells but pays off for the scoop
search and lets the grid scale later. Start CPU-only for correctness, port hot
kernels to WebGPU once behavior is validated.

## Processes

### 1. Uplift
Uniform `z += U(t)·dt`. `U(t)` reuses the 1D sandbox forcing functions
(bell / plateau / pulse / arc / constant). One compute pass.

### 2. Fluvial erosion — stream power (backbone)
`dz/dt = U − K·Aᵐ·Sⁿ` (m≈0.5, n≈1); A = upslope drainage area, S = local slope.
- Flow routing on CPU: D8 steepest-descent receivers → topological stack →
  drainage-area accumulation → **implicit FastScape incision**
  (Braun & Willett 2013): O(n), unconditionally stable.
- **Priority-flood depression filling** so pits/lakes don't break routing.
- GPU upgrade path: **FastFlow** (Jain et al. 2024) — O(log n) flow accumulation,
  O(log² n) depression routing via parallel pointer-jumping.

### 3. Glacial erosion — approximate U-shaped valleys
Literature: erosion ∝ basal **sliding velocity** (E ∝ U_slide¹⁻²); U-shape arises
because sliding is *suppressed at the narrow valley center*, so peak erosion
migrates toward the walls → widening/flattening of a V into a parabola
(Harbor 1992; MacGregor et al. 2000). Cheap approximation, built as a hybrid:
- **Ice mask via ELA**: erode glacially only above an equilibrium-line altitude
  (time-forced → glacial cycles, glacial-buzzsaw elevation capping).
- **Discharge-weighted, width-scaled lowering**: flow accumulation as ice-flux
  proxy; erode a lateral swath whose **width scales with discharge** (trunk
  glaciers wide, tributaries narrow & shallow → **hanging valleys** emerge).
- **Cross-valley profile shaping**: lowering suppressed at centerline, emphasized
  on lower walls (flat floor + steep walls). Alt: relax flow-perpendicular
  cross-sections toward a discharge-scaled parabolic target.
- **Cirque/headwall erosion**: extra erosion at upstream accumulation tips.

### 4. Mass wasting
**(a) Creep / soil diffusion** — *(first process implemented)*
- Linear: `dz/dt = D·∇²z` (Laplacian stencil → convex hilltops).
- Nonlinear **Roering (2001)**: `q = −D·∇z / (1 − (|∇z|/S_c)²)`. Flux diverges as
  slope → critical `S_c`, producing planar threshold hillslopes and handing off
  to landsliding. Slope-limited dt; auto-substepping for stability.

**(b) Episodic deep failure — SCOOPS3D-like scoop search**
SCOOPS3D (Reid et al. 2015, USGS TM 14-A1): millions of trial **spherical** slip
surfaces, 3D **method of columns** limit-equilibrium (3D Bishop simplified /
Fellenius), min factor-of-safety (FoS) per cell. Browser approximation:
- Candidate scoops `(x, y, r)`, seeded only near over-steepened cells.
- FoS per scoop = one GPU workgroup; columns = cells where the sphere sits below
  ground; cooperatively sum 3D-Bishop numerator/denominator
  (`c·A + (W − u·A)·tanφ` over `W·sinα`); params c, φ, γ, r_u.
- Reduce to per-cell min FoS.
- **Episodic trigger**: when incision/creep/uplift push min FoS < 1 (optional
  stochastic threshold for discreteness), fail the worst scoop: lower z to the
  slip surface, **deposit removed volume at the toe** as runout. Deep failures
  stay discrete, not continuous.

## Initial topography generators (seeded, reproducible)
1. **Tilted plane** — `z = gradient · x`. Control case (diffusion barely acts).
2. **Fractal hills** — summed-octave value-noise FBM (seeded hash lattice,
   smoothstep interpolation, persistence). The interesting case for creep.
3. **Scattered cones** — N random cones, `z = max_k h_k·max(0, 1 − dist/r_k)`.
   Sharp features that visibly round under diffusion.

All keyed by an integer **seed** (in URL state) for reproducibility.

## Coupling & the payoff
Per step: **uplift → (fluvial | glacial) → creep → periodic scoop check**. Uplift
+ incision oversteepen slopes → creep saturates → FoS < 1 → episodic deep failures
feed sediment to valleys → rivers re-excavate. Toggle glacial (drop ELA) → V→U
conversion with hanging tributaries. All driven by the same time-forcing as the
1D model.

## Roadmap
- **P0** — Scaffold: float32 height field, hillshade render, byte-DEM I/O,
  D3 cross-section + time-series, URL state, uplift.
- **P1** — Creep (linear → nonlinear) + the three topo generators.  ✓ *done*
- **P2** — Fluvial: CPU flow routing (priority-flood + MFD area) + implicit
  erosion-deposition (Yuan et al. 2019, Gauss-Seidel; unconditionally stable, no
  sub-stepping) + perimeter base level + spatial uplift.  ✓ *done*
- **P3 / P5** — Glaciers: full ice layer.
  - **P5a** ELA mass balance + SnowSlide-style avalanching (bare peaks).  ✓ *done*
  - **P5b** SIA ice flow (Glen flow law, conservative flux-form, CFL-subcycled
    with a fixed diffusivity cap) + ice thickness/speed panels + white-ice
    rendering.  ✓ *done*
  - **P5b′** Equilibrium tracking: each frame iterates ice to steady state with the
    current bed (windowed volume-convergence; spin-up ~10²–10³ substeps, then ~10²
    to re-converge as the bed barely moves; per-frame budget + safety thickness
    cap). Decouples ice from geomorphic time so glaciers equilibrate "instantly"
    on Myr timescales; erosion accumulates over the full dt. Bedrock + ice-surface
    drawn separately on the cross-sections.  ✓ *done*
  - **P5b″** Coarse ice grid (4×: 100×100): ice runs coarse (smooth field; D cap ∝
    dx² → real flow speeds + ~16× cheaper), bed block-averaged down, ice/speed/
    erosion bilinearly upsampled to carve the fine bed. Big perf win. Avalanche set
    aside for now.  ✓ *done*
  - **P5c** Glacial erosion: abrasion `Ė = K_g·|U|` (∝ basal sliding speed, driven
    by the smooth flux velocity to avoid a slope→erosion feedback); till exported
    as outwash. Validated: bounded, glacier persists, erosion peaks where ice flows
    fastest.  ✓ *done*
  - **P5b‴** Edge outflow boundary (bed drops ~2.5·dxC outside → ice spills off the
    rim): gives a bounded equilibrium instead of domain-filling runaway when
    accumulation ≫ ablation. Rain = 0 above the firn line (no fluvial in the
    accumulation zone); incision skipped under ice. Erosion uses the capped flux
    velocity (uncapped Γ·H⁴·|∇s|³ runs away via |∇s|³); strength set by K_g.  ✓ *done*
  - **P5c′** Meltwater + outwash: ablation releases meltwater (→ fluvial discharge)
    and the glacial sediment load (→ fluvial Qs), both melt-weighted; reuses the
    fluvial solver to deposit tapering outwash downstream (mass-conserving, ~free).
    Tilt-plane ice outflow drains only the down-tilt (low) edge, mirroring the water
    perimeter. Ice surcharge buttresses scoops (debuttressing on deglaciation).  ✓ *done*
  - **P5c″** Lakes: fluvial routing on the surface (bed+ice) so glaciers dam/divert
    flow instead of water vanishing under them; priority-flood fill depth = ponded
    water → ice-dammed + bedrock lakes (render + cross-sections), which trap sediment
    and silt in. Phantom-retreat-patch fix (upwind flux guard); asymmetric ablation
    gradient (2.5×) so tongues don't over-extend below the ELA.  ✓ *done*
  - **Later** cirque/headwall retreat; englacial sediment residence/lag; lake outlet
    incision / level dynamics; subglacial drainage.
  - *Numerics note:* explicit sub-cycled flux-form chosen over semi-implicit ADI —
    the bed-driven SIA flux is advective and ADI was only stable to dt≈2 yr.
    The diffusivity cap must be **dt-independent**; a dt-scaled cap throttles flow
    to zero at large dt and ice runs away.
- **P4** — Deep-seated landslides: stochastic spherical-scoop testing + 3D
  method-of-columns FoS (vector driving) + unified block-rotation + growing
  angle-of-repose runout. Catastrophic, persistent slow rotational slumps (ride
  uplift, creep, reactivate, retire), and vanish modes.  ✓ *done*
  Remaining: WebGPU port of the scoop search.
- **P5** — Coupling, glacial cycles (time-forced ELA), presets, WebGPU port of hot kernels.

## References
- Braun & Willett (2013), *A very efficient O(n), implicit and parallel method to
  solve the stream power equation* — FastScape incision scheme.
- Jain et al. (2024), *FastFlow: GPU Acceleration of Flow and Depression Routing*,
  Computer Graphics Forum — GPU flow routing, O(log n).
- Harbor (1992), *Numerical modeling of the development of U-shaped valleys by
  glacial erosion*, GSA Bulletin; MacGregor et al. (2000) — glacial valley form.
- Roering et al. (2001), *Hillslope evolution by nonlinear creep and landsliding*,
  Geology — nonlinear transport law.
- Reid, Christian, Brien & Henderson (2015), *Scoops3D*, USGS Techniques and
  Methods 14–A1 — 3D method-of-columns slope stability.
- Reference GPU harnesses (hydraulic, not LEM — borrow plumbing only):
  joshbrew/webgpu_hydraulic_thermal_erosion_Jako2011 (WebGPU device + DEM I/O),
  GPU-Gang/WebGPU-Erosion-Simulation, LanLou123/Webgl-Erosion.
- Process-law reference / offline validation: Landlab, fastscapelib (Python).
