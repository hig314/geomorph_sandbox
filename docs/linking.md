# Linking the tiers: isolate a node

*Landed 2026-09-25. Pages: `d2_along/index.html` (click tool + overlay), `d1/index.html`
scenario "node". Model: `d1/model.js createNodeColumn`, `d2_along/model.js nodeSpec`.
Mirror: `py/d1_node.py`. Tests: `py/test_d1_node.py`.*

## The question it answers

The ladder (DESIGN.md §2) implements one set of process laws at several dimensionalities,
each tier computing what it resolves and asserting the rest. The link makes that concrete
for one point: **how much of what a node does comes from its neighbours?**

Click a node on the along-valley profile. Two curves then grow together as the profile
runs: the node's real elevation, and the same node run as an isolated d1 column with the
same laws, where every quantity that came from adjacency is held at its click-time value.
The two start identical; the growing gap is the behaviour that arises from the extra
dimension.

## What is frozen, what runs, what stays local

| | in the profile | in the isolated column |
|---|---|---|
| downstream bed `z_down` | uplifts and erodes | **keeps doing what it was doing**: `z_down(t) = z_down(t₀) + ΔU_cum − e_down (t − t₀)`, `e_down` its isolation-time erosion rate (frozen, asserted) |
| downstream ice `H_down` | evolves | **frozen, asserted** |
| ice flux `Q` leaving the node | routed from upstream balance | **frozen, asserted** (Q = 0 → fluvial node) |
| area `A`, widths `W`, `W_f`, spacing `ds` | asserted (Hack, geometry) | asserted |
| uplift `U(t)·f(s)` | shared forcing | **runs** |
| ELA(t) | shared forcing, acts through Q | runs, but **inert**: Q is frozen |
| lithology `r = body(s, z − U_cum)` | local | **runs** (material frame) |
| slope `S = (z − z_down)/ds` | local against a moving receiver | local against the frozen receiver |
| thickness from the face flux law | local against a moving downstream surface | local against the frozen surface |

So the column still has feedbacks of its own: ice-free, its slope relaxes toward
`(U f / K A^m)^(1/n)` and it then rises with the block at the uplift rate; under ice, lowering
its bed relative to the downstream surface thickens the ice and slows sliding, a
self-limiting overdeepening. What it cannot do is see its neighbours erode, its glacier
advance or retreat, or a glacial cycle at all.

**What the receiver does, and two rejected variants.** The receiver keeps doing what it was
doing at isolation: uplifting with the shared forcing and eroding at its then rate. For a
profile in steady state those cancel, the receiver holds its height and every isolated
column stays identical to its node (tested: zero divergence over 100 kyr at every seventh
node), so what remains in the divergence is adjacency changing: neighbours speeding up or
slowing down, or the terminus sweeping over a node. Freezing `z_down` itself is the same in
steady state but misses a receiver that changes pace; letting `z_down` uplift without
eroding turned the receiver into an infinitely hard dam, so the node above eroded its slope
away and then rode uplift with no erosion at all. (The observation that prompted this, an
isolated column that seemed to have its uplift switched off, was a steady node correctly
holding its height while the profile's node rose under an advancing terminus.)

Per step the column mirrors `d2_along.step` for one node: lithology → uplift → if `Q > 0`,
`H` from `faceThickness` (the same monotone root the profile's march uses), `U_s = f_s (Q/W)/max(H, H_min)`,
`E_g = K_g f_Kg U_s^l` capped; then the implicit Yuan node against `z_down` at the same under-ice
fluvial efficiency the profile uses (1 on bare ground, ramping to `φ_sub` over `H_f`).

## Validation

- `faceThickness` (now in `core/laws.js` and `gs_core.py`, shared by the march and the column)
  reproduces the profile's marched thickness node by node.
- **Consistency**: with its adjacency refreshed from the true profile every step (the
  downstream bed the march saw, the receiver the fluvial sweep saw, the routed flux), the
  column reproduces the profile node to 0.000 m under ice and at the terminus, and to the
  Gauss–Seidel ordering elsewhere. Same laws, same order of operations.
- Isolated fluvial column with a fixed receiver equilibrates at `z_down + ds (U f / K A^m)^(1/n)`
  with `E = U`; with the default receiver a steady fluvial profile gives zero divergence at
  every node.
- Isolated glacial column: ice thickens and sliding slows as the bed lowers; the face law
  holds at the solve.
- Isolated vs profile diverge, and by hundreds of metres over 100 kyr near a terminus.

## Using it

On the d2_along page every node is isolated as a column at once, in lockstep with the
profile, from the moment of isolation: Reset, a parameter change, or the **Re-isolate**
button, which re-bases all columns on the current state. Move the mouse across the long
profile and the panel to its right shows, for the node under the pointer, the profile's
elevation (solid) and the isolated column's (dashed) since isolation, with ice surfaces;
click to pin a node so the curves stay. The readout gives the divergence and current
erosion and ice in each. "Open in d1" carries the frozen constants and the process
parameters to the d1 page as sliders, each badged asserted with its equation on hover, so
you can ask what the neighbours would have to do to change the column's fate.

Cost: one column step per node per profile step (a face-law root solve and an implicit
fluvial node), about as much again as the profile itself at 301 nodes; histories are four
arrays of N per record, thinned by two beyond 3000 records.

The same device is planned for d2_across (a column from a cross-profile node, where the
frozen quantities are the neighbours' elevations and the transverse ice) and for d3
(a column from a cell, and a long profile from a flow path).
