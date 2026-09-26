# d2_along — the along-valley tier, z(s, t)

*Step 3 of the build order (DESIGN.md §6), landed 2026-09-25. Page: `d2_along/index.html`.
Model: `d2_along/model.js` (no DOM). Mirror: `py/d2_along.py` (also holds the reference
SIA). Tests: `py/test_d2_along.py`.*

Every control and panel carries its equations on mouse-over. This file is the same
material in one place, plus what the validation found.

## Geometry and the asserted quantities

Nodes `i = 0 … N−1` at `s = i·ds`, divide at `s = 0`, outlet at `s = L` held at base level
`z = 0` (no uplift, no erosion). The profile cannot resolve its basin or its cross-section,
so two quantities are asserted (DESIGN.md §2 ladder):

- drainage area by Hack's law, `A(s) = k_a (s + s₀)^h`, with `k_a = 5.7`, `h = 1.667`
  (m units; reproduces `L = 1.4 A^0.6` in km) and `s₀ = 200 m`;
- valley width `W(s) = W₀ + k_w s`.

Everything else is computed. Per step: lithology at the surface → uplift → steady glacier and
glacial erosion → fluvial incision on ice-free nodes.

## Glacier: steady ice discharge

Ice equilibrates in 10²–10³ yr while the bed evolves over 10⁵–10⁶ yr, so each step solves
the steady glacier on the current bed (the legacy DEM did the same by iterating an explicit
SIA to convergence every frame; this tier replaces that with a direct solve).

1. **Balance on the ice surface** `z_s = z + H`: `b = β(z_s − ELA)` above the ELA,
   `2.5 β (z_s − ELA)` below, clamped to [−8, 2] m/yr. Surface-referenced by the legacy
   decision, so the surface-elevation feedback is kept.
2. **Routed flux** `Q(s) = ∫₀ˢ b W ds`, clamped at zero: where ablation exhausts the flux the
   glacier ends. `q = Q/W`.
3. **Thickness by upstream march.** At each face between nodes `i` and `i+1` the SIA flux law
   must carry the routed flux:

       Γ′ H_f^(n+2) S_f^n = Q_i / W_f,   H_f = ½(H_i + H_{i+1}),   S_f = (z_i + H_i − z_{i+1} − H_{i+1}) / ds

   with `Γ′ = Γ/(1 − f_s)`, `Γ = 2A/(n+2)·(ρg)ⁿ`, `n = 3`. Starting from `H = 0` at the
   outlet, this is solved for `H_i` node by node upstream. The left side is monotone in
   `H_i`, so the root is unique; bisection to 1 mm in at most 60 steps. Nodes with zero
   outgoing flux carry no ice.
4. **Outer iteration**: surface → balance → flux → march, under-relaxed by 0.5, until the
   thickness changes by less than 1 cm (20–35 iterations from bare bed, a handful once ice
   exists from the previous step).
5. **Sliding and erosion**: `U_s = f_s q / H`; `E_g = K_g U_s^l · f(r, c_g)`, capped at
   2 cm/yr; optional quarrying `K_q U_s · max(−∂²z/∂s², 0)`.

**What was tried first and failed.** The pointwise inversion in DESIGN.md §4,
`H = (q / Γ Sⁿ)^(1/(n+2))` iterated on the surface, diverges: wherever the surface slope
reaches the floor the inverted thickness explodes (3 km of ice where the reference has
245 m), and under-relaxation does not save it. The march is the boundary-value fallback the
design anticipated, and it is exact by construction: at convergence every face satisfies the
flux law to 1e-3 relative.

**Decision gate** (DESIGN.md §4: thickness within ~30 % of a reference SIA). The reference is
`reference_sia()` in `py/d2_along.py`: the legacy explicit conservative flux-form flowline
SIA (CFL sub-cycled, upwind guard, surface-referenced balance, no diffusivity cap) run to
windowed volume convergence on a fixed bed, with sliding entered the same way
(`q_total = q_def/(1 − f_s)`). Results, `f_s ∈ {0, 0.5, 0.8}`:

| bed | max H (ref) | H error (max, nodes with H > 20 m) | volume ratio | terminus |
|---|---|---|---|---|
| linear | 245 m | < 1 % | 1.00 | same node |
| concave | 195 m | < 1 % | 1.00 | same node |
| concave + 150 m bump | 196 m | < 1 % | 1.00 | same node |
| concave with a 200 m step down at 9–14 km | 193 m | < 1 % above the step | 0.97 (f_s = 0), 0.69 (f_s = 0.8) | ref 1–2 km further |

The step bed is the reference's problem, not the march's: its face flux at 9 km is up to
eight times the integrated upstream balance, impossible at steady state. The explicit scheme
surges over the 200 m bed step (an icefall limit cycle) and the windowed convergence test
catches it mid-cycle. The routed solution is a true steady state. This is worth remembering
for d3: the legacy SIA is a reference only on smooth beds.

Also checked: finite thickness at the divide (no blow-up); sliding thins the glacier at
fixed flux; the ice surface never slopes uphill under ice; the depth-averaged deformation
speed matches the textbook Glen-law value.

## Fluvial long profile

Detachment-limited stream power with erosion–deposition (Yuan et al. 2019), one node at a
time from the outlet upstream, using `core/laws.js yuanNode` with the downstream node as
receiver; `K` multiplied by the lithology factor. With `G > 0` the sediment flux is
accumulated downstream and the update is iterated Gauss–Seidel with relaxation 0.6.

**A legacy bug, fixed here and flagged for d3.** The legacy DEM code applied the 0.6
under-relaxation on every pass, including the single pass used when `G = 0`. That pass is
already the exact implicit solve (receivers updated first), so relaxing it applies only 60 %
of the incision each step: steady slopes came out about twice the analytic `U/(K A^m)`, and
erosion depended on dt. Here relaxation applies only inside the deposition iteration. The
d3 migration must carry this fix (its effective K was ~0.6 K, dt-dependent).

Checked: steady-state slope equals `(U / K A^m)^(1/n)` to 1e-6 for `n = 1` and to 1 % for
`n = 1.5`; erosion equals uplift at steady state; eroded minus deposited equals exported to
round-off for `G = 0.5` and `1.5`; steady profiles are dt-independent to 0.01 m, a 200 kyr
transient differs by 15 m between 250 and 1000 yr steps.

## Glacial erosion stability

`K_g = 10⁻⁴` and `10⁻²` (×100) over 300 kyr on the default valley: no NaN, thickness
bounded, erosion at the cap at most. Total glacial erosion at ×100 is within a factor 2 of
×1: the glacier cuts its bed below the ELA and shuts itself off (the buzzsaw), which is the
self-limitation the design relies on. Volume closure of the coupled run
(uplift − fluvial − glacial − Δbed) holds to 1e-6.

## Lithology

`r(s) = body(s, 0, z − U_cum(s))` at the surface, multiplying `K` and `K_g` by
`1/(1 + r(c − 1))`. Checked: a dike at steady fluvial state has slope `c` times the analytic
slope inside and analytic outside, with `E = U` everywhere; a layer's surface pattern is
invariant under pure uplift (material frame).

## Sources

- Flowline glaciers: MacGregor et al. (2000) *Geology*; Anderson, Molnar & Kessler (2006) *JGR*.
- SIA and Glen's law: Cuffey & Paterson (2010); `A = 7.57×10⁻¹⁷ Pa⁻³ yr⁻¹`.
- Sliding-based erosion: Hallet (1979); Humphrey & Raymond (1994); Herman et al. (2015).
- Stream power: Whipple & Tucker (1999); implicit solver Braun & Willett (2013); erosion–deposition Yuan et al. (2019).
- Hack's law: Hack (1957).

## Cell-scale regularisations (2026-09-25)

Exploring extreme scenarios exposed three cell-scale defects, each fixed mechanistically
and exposed as a control with its equation on hover:

1. **Outlet boundary.** The march started from zero ice at the outlet, so a glacier that
   reached it thinned to nothing in one cell and `q/H` gave sliding of 10³ m/yr there. The
   outlet is now free outflow: the ice leaves with the uniform-flow thickness for the local
   bed slope, `Γ′ H^(n+2) S_bed^n = q`. (The Python reference SIA keeps `H = 0` at the outlet;
   the gate beds never reach it.)
2. **Terminus wedge.** `U_s = f_s q/H` diverges as `H → 0` at the last ice cell. Sliding now
   uses `max(H, H_min)` with `H_min = 10 m` (asserted): ice thinner than that is too thin to
   slide erosively. Shared law (`core/laws.js slidingSpeed`), so the d1 node column has it too.
3. **Checkerboard.** In the reported scenario the bed near the outlet grew a cell-scale
   sawtooth of a kilometre in 200 kyr: `E ∝ 1/H` with a face law that averages neighbouring
   thicknesses lets odd and even nodes decouple once thin-ice cells erode at the cap. Fixes 1
   and 2 removed the driver (the thin-ice cap-rate cells); with them in place the raw law no
   longer grows the mode in that scenario. Glacial erosion is nevertheless applied through a
   `[¼ ½ ¼]` footprint (an ice-thickness-scale patch rather than one cell), which annihilates
   a two-cell mode exactly, as a safeguard. Off shows the raw law.

Also: the "steady" initial profile caps its slope at 0.6 (a threshold hillslope). Where
`U/(K A^m)` asks for more, the profile is at threshold, not in fluvial equilibrium, and the
readout says over what length. Without the cap an erodibility of 3×10⁻⁷ at 2.8 mm/yr made a
217 km high starting profile.

Checked (`py/test_d2_along.py`): a glacier reaching the outlet has thickness and sliding
varying by < 30 % over the last four cells; in the reported scenario the interior bed shows
no checkerboard signature (sign of Δ²z flips at < 20 % of nodes, |Δ²z| < 100 m) with the
footprint on or off; the footprint annihilates a two-cell mode; volume still closes; the
capped profile never exceeds slope 0.6. What remains in that scenario is a ~75 m step
400 m from the divide, where thin head ice erodes at the cap: a headwall, not a mode.
