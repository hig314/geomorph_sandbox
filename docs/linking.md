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
| downstream bed `z_down`, downstream ice `H_down` | evolve | **frozen, asserted** |
| ice flux `Q` leaving the node | routed from upstream balance | **frozen, asserted** (Q = 0 → fluvial node) |
| area `A`, widths `W`, `W_f`, spacing `ds` | asserted (Hack, geometry) | asserted |
| uplift `U(t)·f(s)` | shared forcing | **runs** |
| ELA(t) | shared forcing, acts through Q | runs, but **inert**: Q is frozen |
| lithology `r = body(s, z − U_cum)` | local | **runs** (material frame) |
| slope `S = (z − z_down)/ds` | local against a moving receiver | local against the frozen receiver |
| thickness from the face flux law | local against a moving downstream surface | local against the frozen surface |

So the column still has feedbacks of its own: ice-free, its slope relaxes toward
`(U f / K A^m)^(1/n)`; under ice, lowering its bed against a fixed downstream surface
thickens the ice and slows sliding, a self-limiting overdeepening. What it cannot do is
see its neighbours erode, its glacier advance or retreat, or a glacial cycle at all.

Per step the column mirrors `d2_along.step` for one node: lithology → uplift → if `Q > 0`,
`H` from `faceThickness` (the same monotone root the profile's march uses), `U_s = f_s (Q/W)/max(H, H_min)`,
`E_g = K_g f_Kg U_s^l` capped; else the implicit Yuan node against `z_down`.

## Validation

- `faceThickness` (now in `core/laws.js` and `gs_core.py`, shared by the march and the column)
  reproduces the profile's marched thickness node by node.
- **Consistency**: with its adjacency refreshed from the true profile every step (the
  downstream bed the march saw, the receiver the fluvial sweep saw, the routed flux), the
  column reproduces the profile node to 0.000 m under ice and at the terminus, and to the
  Gauss–Seidel ordering elsewhere. Same laws, same order of operations.
- Isolated fluvial column equilibrates at `z_down + ds (U f / K A^m)^(1/n)` with `E = U`.
- Isolated glacial column: ice thickens and sliding slows as the bed lowers; the face law
  holds at the solve.
- Isolated vs profile diverge, and by hundreds of metres over 100 kyr near a terminus.

## Using it

On the d2_along page, click the long profile; a vertical marker shows the node and the
bottom panel starts the overlay (solid: profile; dashed: isolated). Step or play: both
advance in lockstep. The readout gives the divergence and current erosion and ice in each.
"Open this column in d1" carries the frozen constants and the process parameters to the d1
page as sliders, each badged asserted with its equation on hover, so you can ask what the
neighbours would have to do to change the column's fate. Changing a profile parameter
clears the link (the frozen column would no longer share the profile's laws); click again.

The same device is planned for d2_across (a column from a cross-profile node, where the
frozen quantities are the neighbours' elevations and the transverse ice) and for d3
(a column from a cell, and a long profile from a flow path).
