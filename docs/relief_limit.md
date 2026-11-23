# Limits to relief — R(t)

*Landed 2026-11-21. Page: `relief_limit/index.html`. Model: `relief_limit/model.js` (no DOM).
Mirror: `py/relief_limit.py`. Tests: `py/test_relief_limit.py`.*

A presentation figure more than a model: relief against time, in seven steps that each add
one idea. It sits below d1 on the ladder. d1's relief reservoir (`dR/dt = U − kRⁿ`) with
large n already behaves as a hard ceiling; this page replaces that smooth ceiling with
discrete failures and asks what sets the ceiling.

## The rule

Dimensionless: relief in units of the functional limit L, time in units of L/g₀, where g is
the constant growth rate of relief (uplift and valley incision) and g₀ its reference value.
The "Uplift rate" control sets g/g₀ from 0 to 1 (default 1, the maximum); the earthquake
sequence does not depend on it. At g = 0 relief never grows: a failure needs an earthquake
that reaches the standing relief or, in step 7, enough damage to bring the limit down to it. With ceiling C (L in steps
1–3, the static limit Ls from step 4):

    between events   dR/dt = g; arriving at R = C is a failure
    earthquake       D = f·C, a fraction of the limit as it stands; L_eff = C − D;
                     failure if R ≥ L_eff
    failure          R ← max(L_eff − d, R_floor)      (L_eff = C for a failure at the ceiling)
    floor            R_floor = 10 % of the graph height (1.1 Ls) by default: a landslide never
                     takes relief to zero; a slope at the floor does not fail, and an
                     earthquake that would have failed it leaves no damage
    step 1           no landslide, relief holds at the ceiling
    landslide size   d constant in steps 2–5; in steps 6–7 each failure draws its own d
    damage (step 7)  an earthquake that causes no failure: C ← C − k·D
                     any failure (at an earthquake or at the damaged ceiling): C ← Ls

The solution is piecewise linear, so the model steps from event to event and is exact.
There is no time step and nothing to tune for stability.

## The seven steps

| Step | Adds | Terms on |
|---|---|---|
| 1 | constant limit, relief flatlines | d = 0, no earthquakes |
| 2 | small landslides | d = d_small |
| 3 | large landslides | d = d_large |
| 4 | a static limit above; constant earthquakes explain L | C = Ls, f = (Ls − L)/Ls at each instant R reaches L |
| 5 | the earthquakes drawn: dashed drop-line from Ls to a dot at Ls·(1 − f) | as 4 |
| 6 | a spectrum of earthquakes | Poisson arrivals (rate λ), Gutenberg–Richter magnitudes, seeded |
| 7 | damage, removed by the next landslide | k > 0 |

(The request numbered two steps "4"; they are 4 and 5 here.)

Steps 3, 4 and 5 are one relief history, so the page keeps the time cursor and fades the
new elements in. In steps 4 and 5 the coincidence is stipulated (user direction,
2026-11-21): an earthquake strikes at each instant relief reaches L, so the recurrence is
the regrowth time of one landslide, d/g, and follows the uplift rate and landslide size
instead of being a control. At g = 0 relief never reaches L and no earthquake is drawn.

## Earthquake sizes (steps 6–7)

Magnitudes are drawn from the doubly truncated Gutenberg–Richter law on
[Mmax − ΔM, Mmax], `N(>M) ∝ 10^(−bM)`, b = 1 by default (Gutenberg & Richter 1944;
truncated form as in Cornell & Vanmarcke 1969). Shaking amplitude at a fixed distance grows
as `10^(cM)`; the fraction of the limit an earthquake removes is taken proportional to it:

    f = f_max · 10^(−c (Mmax − M)),      D = f · Ls(t)

so f follows a truncated power law, `N(>f) ∝ f^(−b/c)` (exponent 2 at the defaults), between
`f_max·10^(−cΔM)` and `f_max`. Only magnitude differences enter, so no absolute magnitude is
needed. **c = 0.5 is asserted** (ground-motion models give roughly 0.3–0.5 for peak
acceleration) and **distance to the source is not sampled**: every earthquake is equally
far away. A site hazard curve, which folds distance in, is the natural replacement.

The depression is relative to the limit as it stands (user direction, 2026-11-21): a
damaged slope is pulled down from its damaged limit by the same fraction, and damage is
`Ls ← Ls − k·D = Ls(1 − k f)`.

## Landslide sizes (steps 6–7)

Each failure draws its increment from a truncated power law (log-log linear in
size–frequency space), by inversion:

    P(d > x) = (x^−α − d_max^−α) / (d_min^−α − d_max^−α),   d_min ≤ d ≤ d_max

Defaults d_min = 0.05 L, d_max = 0.6 L, α = 1. The sizes use their own random stream
(seed + 1000003), consumed in order of occurrence, so they never alter the earthquake
sequence. Power-law size statistics are established for landslide areas and volumes (e.g.
Hovius et al. 1997; Stark & Hovius 2001; Malamud et al. 2004); applying one to the relief
increment is an assertion of this figure. "Constant" restores the step 3–5 size.

## Choices made, open to change

- **A failure leaves relief d below the limit in force**, not d below where it stood. In
  steps 2–5 the two are the same. From step 6 it means a large earthquake removes all the
  relief above its dot and d more, so failure size varies with earthquake size as well as
  timing. The alternative (a fixed drop d from the current relief) is a one-line change.
- **A landslide removes the damaged rock and restores the static limit to its original
  value** (added 2026-11-21 at the user's direction). Damage accumulates over the
  earthquakes a slope survives and lasts until its next failure, of either kind. The
  failure itself is measured from the damaged limit (R ← C_damaged − D − d); the limit is
  restored after it.
- **k ≤ 1**, so damage alone never brings the static limit down to the relief (the
  surviving slope has R < Ls − D < Ls − kD).
- Defaults (rate 30, ΔM 2, f_max 0.6) draw 113 earthquakes of which 9 cause landslides in
  step 6, with mean relief near L.

## The figure

One SVG, viewBox 1920 × 1080, attributes only (the saved SVG is self-contained). Relief is
the only heavy line; limits are grey dashes, labelled directly at the right; earthquakes
are the only colour, filled dots where they fail the slope and open where they do not. No
ticks or grid: the axes are dimensionless. Two line weights a factor of two apart: full
weight for axes, relief, limits and circle outlines, half weight for the earthquake
drop-lines; a double weight is reserved for faded reference lines behind the others. The vertical frame is the same in every step so
the static limit arrives in empty space. In step 7 the static limit is drawn stepped up to
the time cursor and level beyond it.

**In frame and out.** On the last step two sliders sit inside the 16:9 frame, under the
time axis: weakening (the damage fraction k) and uplift rate (g). They are HTML laid over the SVG inside
`#slide`, sized in container units so they scale with the figure and go full screen with
it, and they mirror the controls of the same name below the frame. Every other control is
outside the frame as presenter's notes. "Sliders in frame" switches them to every step or
off. The saved SVG does not include them.

**Pace.** Two drawing times, in seconds for the line to cross the slide: `drawEarly` for
steps 1–5 and `drawTime` for steps 6–7 (both 12 s by default, 1–40 s). Steps 4 and 5 keep
the line of step 3, so `drawEarly` matters for them only when one is opened directly. Both
are URL parameters, so a slide can pin its own pace (`…?embed=1&drawEarly=5`).

**Embedding.** `relief_limit/index.html?embed=1` shows only the 16:9 frame, for use as a
slide inside another deck (the `Presentation` project). The host sets the step with
`postMessage({type: "relief-step", step: n, replay: bool})`; step keys pressed inside the
figure are sent back as `{type: "deck-key", key}`. Other URL parameters still apply.

Keys: ← → (or a clicker's PageUp / PageDown) step, space play / pause, R replay, F full
screen, 1–7 jump. Every parameter is in the URL.

## Validation (`py/test_relief_limit.py`, 80 checks)

- every step: time monotone; each segment is growth at exactly g, a hold, or a drop;
  relief budget `R0 + g·t − Σdrops = R_end` to round-off; relief never above the ceiling;
  after each earthquake `R ≤ Ls − D`
- step 1 reaches L at (L − R0)/g and holds; steps 2–3 have period d/g, drop d, the expected
  count; long-run mean relief L − d/2
- steps 3 = 4 = 5 vertex for vertex; in steps 4–5 every earthquake finds relief exactly at
  L, for a range of uplift rates, landslide sizes, starting reliefs and floors
- negligible earthquakes: failures at the static limit
- step 7: each earthquake pulls down from Ls − k·ΣD over the earthquakes survived since the
  last landslide; every landslide (including one at the damaged ceiling) restores Ls exactly;
  k = 0 is step 6; steps 6 and 7 share the earthquake sequence
- sample rate of a long sequence; magnitudes against the exact truncated Gutenberg–Richter
  exceedance at three levels; f against its definition and bounds; the b/c power-law ratio;
  D = f × the limit in force for every earthquake, including damaged limits
- an independent time-stepped integration of the same rules agrees at dt 1e-3 and 2.5e-4
  (at 4e-3 a step can hold a ceiling failure and an earthquake and misorder them)
- g = 0 in every step: relief non-increasing, no failure without an earthquake, no division
  by zero; in step 7 damage alone brings the limit down to the relief
- landslide sizes: constant in steps 3–5; in steps 6–7 each failure leaves relief its own d
  below the limit in force, within the cut-offs; sampled exceedance against the exact
  truncated power law at α = 0.5, 1, 2; independent of the earthquake stream
- the floor: no landslide leaves relief below it in any step under four stress settings;
  an earthquake is marked as failing the slope only if a landslide happened; floor = 0
  reproduces the unfloored rule
- g = 0.2, 0.6, 1.7: period d/g, and the step-7 budget, growth rate, earthquake bound and
  reference; the earthquake sequence is independent of g
- fingerprints per step; `GS.relief.fingerprint(n)` in the browser console prints the same
  line (all seven equal to nine decimals, 2026-11-21)
