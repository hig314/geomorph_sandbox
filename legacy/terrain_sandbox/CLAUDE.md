# Terrain Sandbox — project guide for Claude

A browser-based, interactive **2D landscape-evolution sandbox**: a ~400×400 float32 DEM
evolved by coupled geomorphic processes (uplift, fluvial incision, hillslope creep,
deep-seated landslides, glaciers + lakes). Sibling of the 1D `uplift_sandbox`. Pure
static bundle (no build step) intended to drop into a landslidescience.org Django page.

## How to run / work on it
- **Open it:** `open "file:///Users/Hig/Claude_projects/terrain_sandbox/index.html"` (no server needed).
- **There is NO JS runtime here** (no node). You cannot run or syntax-check the JS. Review
  ported code by hand and open in the browser to test.
- **Only `/opt/anaconda3/bin/python3` is available.** Per the user's hard rule, **validate
  every numerical scheme in Python before and after porting to JS** — mirror the JS logic,
  check mass conservation, stability across `dt`, compare against a reference. This has caught
  almost every bug; do not skip it.
- Keep `methods.html` (cited references, the public "how it works" page) and `DESIGN.md`
  (roadmap P0–P5) current as features land.

## Architecture (keep this separation)
```
terrain-model.js  pure model core: grid, topo generators, all process kernels, step().
                  NO DOM, NO D3. This is the source of truth.
render.js         DEM/fields → canvas (hillshade + hypsometric ramp, ice, lakes, field panels).
app.js            controls, animation loop, D3 cross-sections + history timeline, URL state.
index.html / styles.css / methods.html
```
- Compute is **float32**; 1-byte elevation is export-only (sub-metre per-step changes would
  round to zero). Grid `N=400`, `DX=50` m → 20 km domain.
- UI state round-trips through the **URL query string** (`PARAMS` array in app.js). Adding a
  control = add its id to `PARAMS` + a `syncLabels` line + a `stepParams`/`topoOpts` field.
- **Slider convention:** ids ending `Exp` are log sliders (value = log10 of the quantity).

## Processes (detail in methods.html / DESIGN.md)
- **Uplift:** spatial patterns ramp / **gaussian (default)** / tilt; perimeter base level.
- **Fluvial:** priority-flood depression fill → D8 receivers → topological stack → MFD
  drainage area (Holmgren) → implicit erosion-deposition (Yuan et al. 2019, Gauss-Seidel,
  unconditionally stable). Routing is on the **surface (bed+ice)** so glaciers dam/divert flow.
- **Creep:** linear + nonlinear Roering, auto-substepped.
- **Landslides:** stochastic spherical SCOOPS3D-style scoops, 3D method-of-columns FoS with
  **vector-resolved driving** (symmetric scoops cancel → no flat-ground craters). Catastrophic
  (low AOR runout) vs persistent slow slumps (separate AOR). Ice surcharge buttresses scoops.
- **Glaciers:** SIA flow on a **coarse grid**, equilibrium-tracked; ELA mass balance; abrasion;
  meltwater + outwash; lakes. See the big section below — this is where the subtlety lives.

## ⚠️ GLACIER NUMERICS — hard-won, read before touching glaciers
Glaciers run in `glacierStep()` (terrain-model.js). Many "obvious" approaches blow up; these
are the validated decisions. **Do not re-litigate without a Python test.**

- **Coarse grid:** ice physics runs `ICE_COARSEN=4`× coarser (100×100, dxC=200) than the bed.
  Ice is smooth so sub-grid detail is irrelevant; and because the diffusivity cap scales as
  dx², the coarse grid gives *real* flow speeds (not throttled) and is ~16× cheaper. Bed is
  block-averaged down; ice/speed/erosion bilinearly upsampled back to carve the fine bed.
- **Equilibrium tracking:** ice responds in 10²–10³ yr while terrain evolves over Myr, so each
  frame iterates the SIA flux+balance to **steady state** with the current bed (stop on a
  *windowed total-volume* convergence — a pointwise |∂H/∂t| test never converges due to a
  terminus limit-cycle). Spin-up ~10²–10³ sub-steps; then it re-converges in ~10². Erosion
  uses the equilibrium velocity × full geomorphic dt.
- **Solver = explicit conservative flux-form, CFL-subcycled, fixed D cap.** A semi-implicit ADI
  was prototyped and abandoned (the bed-driven flux is advective → only stable to dt≈2 yr).
  The D cap MUST be **dt-independent** (`0.2·dxC²/ICE_DTSUB_MIN`); a dt-scaled cap throttles
  flow to zero at large dt → runaway. `ICE_DTSUBMAX` caps the sub-step or thin ice (Dmax≈0) →
  huge dt → explosion.
- **Mass balance is referenced to the ice SURFACE (bed+H), and STAYS that way (user's call).**
  This carries a real surface-elevation feedback → small-ice-cap bistability (big long-term
  volume swings). That is an intended feature of the modelled system. **Do NOT switch to a
  bedrock-referenced balance** (less physical; the user rejected it).
- **Ablation gradient is asymmetric:** below the ELA the gradient is `ICE_ABL_RATIO=2.5`× the
  accumulation gradient (standard glaciology). Symmetric melt makes tongues over-extend far
  below the ELA.
- **Edge outflow:** the bed is treated as dropping `ICE_EDGE_DROP=2.5`×dxC just outside the
  domain so ice spills off the rim (bounded equilibrium even when accumulation ≫ ablation). For
  the **tilt** pattern it drains ONLY the down-tilt (low) edge, mirroring the water perimeter.
- **Upwind flux guard (`ICE_FLOWMIN`):** zero the face flux when the upwind (higher-surface)
  cell is ice-free — otherwise an ice-free cell with higher *bedrock* than a neighbour's ice
  *surface* spuriously sources ice (the H≥0 clamp manufactures it), leaving phantom patches far
  below the ELA during retreat.
- **Glacial EROSION (rebuilt to published models — Egholm 2009, MacGregor 2000, Hergarten 2021,
  Schoof 2005):** `E = K_g·τ_b`, `τ_b = ρg·H·min(|∇s_surface|, ICE_SLOPE_CAP=0.25)`, plus a
  per-step ceiling `ICE_ERO_RATECAP=0.02 m/yr·dt`. The key stabiliser: τ_b uses the SURFACE
  slope, so over-deepening flattens the surface → driving stress falls → erosion self-limits at
  ~ice thickness (no runaway). The slope cap = Coulomb regularization (drag can't explode at a
  steep wall); the rate cap stops the bed out-running the ice's re-flattening. **Validated
  bounded and robust to K_g ×100.** Things that all blew up (don't retry): |q|/H (cap-flattened,
  no channelization), Γ·H⁴·|∇s|³ (|∇s|³ explodes), raw τ_b without the slope cap (explodes at
  hole edges), Budd U_s=Cs·τ^m/N with m=2,3 (|∇s|^m explodes), erosion sub-stepping (worse).
  Note: surface-flattening stability fundamentally tensions against sharp narrow incision —
  glaciers make broad U-troughs bounded by ice thickness, not knife-cut slots.
- **Avalanche** (`iceAvalanche`) is defined but **currently unused** (SIA sheds steep ice).

## Meltwater / outwash / lakes coupling
- Above the firn line **rain = 0** (precip is snow → ice); no drainage area, no incision there.
- **Meltwater + glacial sediment** are released melt-weighted in the ablation zone and injected
  into the fluvial system next frame: meltwater → drainage discharge; sediment → fluvial `Qs`,
  which the Yuan scheme spreads downstream as tapering **outwash**. Mass-conserving (1-step lag,
  since fluvial runs before glaciers in `step()`).
- **Lakes:** because routing is on the surface, the priority-flood fill depth (`filled − (bed+ice)`)
  IS the ponded-water level → `state.lake`. Ice-dammed + bedrock lakes render blue (map +
  cross-sections) and trap sediment (silt in via the existing submerged-deposition).
- **No moraine deposition** — considered, but we export till to outwash instead (in-place till
  dumping caused spurious mounds). Adding terminal/lateral moraines is an open option.

## Key constants (terrain-model.js, ~line 846) and notable defaults
- Ice: `ICE_COARSEN 4`, `ICE_DTSUB_MIN 1`, `ICE_DTSUBMAX 5`, `ICE_BUDGET 500`, `ICE_VTOL 4e-4`,
  `ICE_CAP 3000`, `ICE_EDGE_DROP 2.5`, `ICE_FLOWMIN 1e-3`, `ICE_ABL_RATIO 2.5`,
  `ICE_SLOPE_CAP 0.25`, `ICE_ERO_RATECAP 0.02`. ρg=8987, n=3, A(Glen)=7.57e-17.
- Defaults (index.html): glaciers **on**, uplift **gaussian** 2 mm/yr, noise 4 m, ELA 1000 m,
  balGrad 0.007, ice flow ×1, glacial erodibility 1e-8 (log slider −10..−6.5), catastrophic
  AOR 14° / slow AOR 35°, max scoop 4 km, tests/step 5, K 1e-5, dt ~1 kyr.

## State fields of note (model.state)
`h` (bed), `ice`/`iceVel` (fine, upsampled from coarse `iceC`/`nC`/`dxC`), `area`, `sedFlux`,
`lake`, `dz`, `iceMelt`/`glacSed` (fluvial-coupling sources), `iceEroVol`, `iceSubs`, slide/
slump counts, `slideInits` (cumulative for the history timeline). UI panels: Δelevation, sed
flux, ice thickness, ice speed (fixed colour scales), History timeline (relief / ice vol /
cumulative slide initiations).

## Open items / things the user may revisit
- **Lakes sometimes under-fill** a closed basin (user observed; may have been entangled with the
  old erosion blow-up). Now that erosion is stable, re-check; if it persists, get a seed+steps
  and investigate the priority-flood/lake code directly — don't guess-fix.
- Moraine ridges (terminal/lateral) — not built.
- Subglacial hydrology / effective-pressure & reverse-slope shutdown for overdeepenings.
- Cirque/headwall retreat; englacial sediment lag; lake outlet incision / level dynamics.
- Planned but not done: WebGPU port of hot kernels (esp. scoop search).

## Conventions
- Match surrounding code style (ES5-ish `var`, dense numeric loops, no frameworks beyond D3).
- Mass-conservation and stability are non-negotiable; report Python validation results plainly.
- The user is the domain expert (geomorphology/glaciology) — when they say a behaviour is
  physically wrong, trust it and find the mechanism, don't just tune.
