# geomorph_sandbox — project guide for Claude

Browser sandbox that **illustrates geomorphic hypotheses** by implementing the same process
physics in **1, 2 and 3 spatial dimensions** (column / profile / DEM) as parallel tiers.
Read `DESIGN.md` first: it holds the purpose, the dimensional ladder, the "asserted vs
computed" device, and the build order.

## Layout
- `core/` pure process laws, forcing, lithology (no DOM). `d1/`, `d2_across/`, `d2_along/`,
  `d3/` tiers. `ui/` shared panels + URL state. `py/` Python validation. `docs/` methods.
- `legacy/uplift_sandbox/` (old 1D reservoir model) and `legacy/terrain_sandbox/` (old 2D-DEM
  model) are a **frozen historic record**. Recycle code from them; do not develop them.
  `legacy/terrain_sandbox/CLAUDE.md` holds the hard-won glacier/SIA numerics — read it before
  touching glaciers; its "don't retry" list still applies even though the glacier is being
  rebuilt as a steady-discharge model.

## How to work here
- **No JS runtime** (no node): you cannot run or syntax-check JS. Review ported code by hand,
  open in the browser to test: `open "file:///Users/Hig/Claude_projects/geomorph_sandbox/<tier>/index.html"`.
- **Only `/opt/anaconda3/bin/python3`.** Hard rule: **validate every numerical scheme in
  Python before and after porting to JS** (mirror the JS logic, check mass conservation,
  stability across dt, compare to a reference). Put the mirrors in `py/`.
- Keep the model core of each tier DOM-free; render and UI separate. Compute in float32.
- Keep `DESIGN.md` and `docs/` current as features land. Cite the published model each tier
  is built on.
- Mass conservation and stability are non-negotiable; report validation results plainly.
- The user is the domain expert (geomorphology / glaciology). When they say a behaviour is
  physically wrong, trust it and find the mechanism; don't just tune.
- Match surrounding code style (ES5-ish `var`, dense numeric loops, D3 only).
