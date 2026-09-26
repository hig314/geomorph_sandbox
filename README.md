# geomorph_sandbox

**Live:** <https://hig314.github.io/geomorph_sandbox/> (GitHub Pages, built from `main`).

A browser sandbox for illustrating geomorphic hypotheses, built as parallel models in one,
two and three spatial dimensions. See `DESIGN.md` for the purpose and plan, `docs/core.md`
for the shared modules and unit convention.

## Run

Static bundle, no build step. Open a tier directly:

    open d1/index.html

or serve the repo root (`python3 -m http.server`) and browse to `/d1/`.

## Validate

Every numerical scheme is mirrored and checked in Python before and after porting:

    /opt/anaconda3/bin/python3 py/test_core.py
    /opt/anaconda3/bin/python3 py/test_d1_columns.py
    /opt/anaconda3/bin/python3 py/test_d2_along.py
    /opt/anaconda3/bin/python3 py/test_d1_node.py

## Layout

- `core/` process laws, forcing, lithology, units (pure JS, no DOM)
- `d1/` column · `d2_along/` long profile · `d2_across/` cross profile · `d3/` DEM
- `ui/` shared controls, URL state, D3 panels, stylesheet
- `py/` Python mirrors and tests · `docs/` methods
- `legacy/` the two projects this grew from (`uplift_sandbox`, `terrain_sandbox`), frozen
