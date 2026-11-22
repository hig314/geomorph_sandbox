/*
 * relief_limit/model.js — "Limits to relief": relief against time, R(t). Pure model, no DOM.
 *
 * Relief grows at a constant rate g until a limit fails it. Seven stages build the idea:
 *   1  constant limit L, relief flatlines at it
 *   2  each arrival at L is a landslide: relief drops by a small increment d and regrows
 *   3  the same with a large d
 *   4  a static limit Ls stands above; L is where constant earthquakes (depression
 *      D = Ls − L) put the limit at the moment they strike. Their timing is stipulated to
 *      coincide exactly with relief reaching L, so the recurrence is d/g, not a parameter
 *   5  the same model, with each earthquake drawn (a drop-line from Ls to Ls − D)
 *   6  earthquakes are a Poisson process with Gutenberg–Richter magnitudes; each depresses
 *      the limit by a fraction f(M) of the limit as it stands; one fails the slope only if
 *      Ls·(1 − f) reaches the relief
 *   7  an earthquake that does not fail the slope damages it: Ls drops by k·D; the next
 *      landslide removes the damaged rock and restores Ls
 *
 * One rule covers every stage. With ceiling C (L in stages 1–3, Ls after):
 *   between events   dR/dt = g, and R = C is a failure
 *   earthquake       D = f·C (a fraction of the current limit); Leff = C − D; failure if R ≥ Leff
 *   failure          R ← max(Leff − d, 0)   (Leff = C for a failure at the ceiling)
 *   landslide size   d constant in stages 2–5; in stages 6–7 each failure draws its own d from
 *                    a truncated power law, N(>d) ∝ d^(−α) on [dMin, dMax] (log-log linear
 *                    in size–frequency space), from a stream separate from the earthquakes
 *   floor            no landslide takes relief below the floor: R ← max(Leff − d, floor); a slope
 *                    already at or below it does not fail (and an earthquake that would have
 *                    failed it leaves no damage either)
 *   g = 0            relief does not grow; only earthquakes (and damage) can fail the slope
 *   d = 0            no landslide: relief holds at the ceiling (stage 1)
 *   damage (7)       earthquake without failure: C ← C − k·D;  any failure: C ← Ls
 *
 * Earthquake sizes (stages 6–7): magnitudes from the doubly truncated Gutenberg–Richter law
 * on [Mmax − dM, Mmax], N(>M) ∝ 10^(−b·M); shaking amplitude at a fixed distance scales as
 * 10^(c·M), and the depression fraction is taken proportional to it:
 *   f = fmax · 10^(−c·(Mmax − M))
 * so f is a truncated power law, N(>f) ∝ f^(−b/c). Distance to the source is not sampled.
 * Stages 4–5 use the constant fraction f = (Ls − L) / Ls.
 *
 * The solution is piecewise linear, so it is advanced event to event and is exact: there
 * is no time step. Dimensionless: relief in units of L, time in units of L/g.
 *
 * Mirrored in py/relief_limit.py; validated in py/test_relief_limit.py (budget, ceiling,
 * saw-tooth period and count, stages 3 = 4 = 5, damage and its removal by landslides, a time-stepped
 * reference, and fingerprints compared with GS.relief.fingerprint in the browser).
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  var EPS = 1e-9;

  var DEFAULTS = {
    stage: 3,
    g: 1,          // growth rate of relief (uplift + incision), L per time unit
    L: 1,          // functional limit
    R0: 0.3,       // initial relief
    tmax: 4,       // length of the record
    dSmall: 0.05,  // landslide increment, stage 2
    dLarge: 0.25,  // landslide increment, stages 3–5 (and 6–7 when dDist = "constant")
    dDist: "power", // landslide sizes in stages 6–7: "power" (truncated power law) or "constant"
    dMin: 0.05,    // smallest landslide increment, stages 6–7
    dMax: 0.6,     // largest landslide increment, stages 6–7
    dAlpha: 1,     // exponent of the cumulative size distribution, N(>d) ∝ d^(−α)
    floor: 0.176,  // relief no landslide can cut below (the page sets it to 10 % of the graph height)
    Ls: 1.6,       // static limit
    rate: 30,      // earthquakes per time unit, stages 6–7
    b: 1,          // Gutenberg–Richter b-value, stages 6–7
    dM: 2,         // magnitude range drawn, Mmax − Mmin, stages 6–7
    fmax: 0.6,     // fraction of the current limit removed by the largest earthquake (Mmax)
    c: 0.5,        // shaking amplitude ∝ 10^(c·M) at a fixed distance
    damage: 0.15,  // k: drop of Ls per unit depression, until the next landslide, stage 7
    seed: 82
  };

  // [[t, f, m]] in time order, t ≤ tmax: f = fraction of the current limit removed,
  // m = Mmax − M (magnitude units below the largest earthquake; 0 in stages 4–5).
  function makeQuakes(p) {
    var st = p.stage, out = [], t, f, m, k;
    if (st <= 3) return out;
    if (st <= 5) {
      // stipulated coincidence: an earthquake strikes at each instant relief reaches L, so
      // the recurrence is the regrowth time of one landslide, (L − max(L − d, floor)) / g
      f = (p.Ls - p.L) / p.Ls;
      if (p.g <= 0 || p.R0 > p.L) return out; // relief never arrives at L: no coincidence to draw
      var T = (p.L - Math.max(p.L - p.dLarge, Math.max(p.floor, 0))) / p.g;
      t = (p.L - p.R0) / p.g;
      for (k = 0; t + k * T <= p.tmax; k++) { out.push([t + k * T, f, 0]); if (T <= 0) break; }
      return out;
    }
    var rnd = GS.noise.mulberry32(p.seed | 0);
    var beta = p.b * Math.log(10), span = 1 - Math.exp(-beta * p.dM);
    t = 0;
    for (;;) {
      t += -Math.log(1 - rnd()) / p.rate;
      m = p.dM + Math.log(1 - rnd() * span) / beta; // Mmax − M, in [0, dM]
      if (t > p.tmax) return out;
      out.push([t, p.fmax * Math.pow(10, -p.c * m), m]);
    }
  }

  // Returns nextD(): the increment of the next landslide. Constant in stages 1–5 (0 in
  // stage 1) or when dDist = "constant"; otherwise a truncated power law on [dMin, dMax],
  // P(d > x) = (x^−a − hi^−a) / (lo^−a − hi^−a), drawn by inversion from its own stream.
  function makeSizes(p) {
    var st = p.stage, cons = st === 1 ? 0 : (st === 2 ? p.dSmall : p.dLarge);
    var lo = Math.min(p.dMin, p.dMax), hi = Math.max(p.dMin, p.dMax), a = p.dAlpha;
    if (st <= 5 || p.dDist !== "power") return function () { return cons; };
    if (hi <= lo * (1 + 1e-9) || a <= 0) return function () { return lo; };
    var rnd = GS.noise.mulberry32((p.seed | 0) + 1000003), span = 1 - Math.pow(lo / hi, a);
    return function () { return lo * Math.pow(1 - rnd() * span, -1 / a); };
  }

  function simulate(params) {
    var p = Object.assign({}, DEFAULTS, params || {});
    p.stage = Math.round(+p.stage);
    var st = p.stage, g = p.g, tmax = p.tmax;
    var nextD = makeSizes(p), hold = st === 1, floor = Math.max(p.floor, 0), d;
    var C = st <= 3 ? p.L : p.Ls;
    var k = st === 7 ? p.damage : 0, C0 = C;
    var qs = makeQuakes(p);
    var t = 0, R = Math.min(p.R0, C);
    var path = [[t, R]], lim = [[0, C]], fails = [], quakes = [];
    var qi = 0, tq, tn, th, D, Leff, trig, Rn, top;
    for (;;) {
      tq = qi < qs.length ? qs[qi][0] : Infinity;
      tn = Math.min(tq, tmax);
      // grow to tn; every arrival at the ceiling on the way is a failure
      for (;;) {
        th = g > 0 ? t + (C - R) / g : Infinity;
        if (th >= tn - EPS) break;
        path.push([th, C]); t = th;
        if (hold) { R = C; t = tn; break; } // no landslide: hold at the ceiling
        d = nextD();
        R = Math.min(Math.max(C - d, floor), C);
        fails.push({ t: t, from: C, to: R, cause: "limit", d: d, lim: C });
        path.push([t, R]);
        if (C < C0) { lim.push([t, C]); lim.push([t, C0]); C = C0; } // the landslide removes the damaged rock
      }
      R = Math.min(R + g * (tn - t), C); t = tn;
      path.push([t, R]);
      if (tq > tmax) break;
      top = C;
      D = qs[qi][1] * C; // relative to the limit as it stands
      Leff = Math.max(C - D, 0);
      trig = false;
      if (R >= Leff - EPS) { // the depressed limit reaches the relief
        d = nextD();
        Rn = Math.max(Leff - d, floor);
        if (Rn < R) {
          trig = true;
          fails.push({ t: t, from: R, to: Rn, cause: "quake", d: d, lim: Leff });
          R = Rn; path.push([t, R]);
          if (C < C0) { lim.push([t, C]); lim.push([t, C0]); C = C0; } // the landslide removes the damaged rock
        }
      } else if (k > 0) {
        lim.push([t, C]); C = C - k * D; lim.push([t, C]); // still > R, since k ≤ 1 and R < C − D
      }
      quakes.push({ t: t, D: D, f: qs[qi][1], m: qs[qi][2], top: top, y: Leff, trig: trig });
      qi++;
    }
    lim.push([tmax, C]);
    return { path: path, lim: lim, quakes: quakes, fails: fails, p: p };
  }

  // Ceiling in force at time t (after any step at exactly t) on the stepped lim polyline.
  function ceilingAt(lim, t) {
    var c = lim[0][1];
    for (var i = 0; i < lim.length; i++) { if (lim[i][0] <= t) c = lim[i][1]; else break; }
    return c;
  }

  // Time-averaged relief over the record (trapezoids; drops have zero width).
  function meanRelief(path) {
    var a = 0, i, T = path[path.length - 1][0] - path[0][0];
    for (i = 0; i + 1 < path.length; i++) a += 0.5 * (path[i][1] + path[i + 1][1]) * (path[i + 1][0] - path[i][0]);
    return T > 0 ? a / T : path[0][1];
  }

  // Same line as py/test_relief_limit.py prints, for checking the port in the console.
  function fingerprint(stage, extra) {
    var s = simulate(Object.assign({ stage: stage }, extra || {})), sum = 0, tr = 0, i;
    for (i = 0; i < s.path.length; i++) sum += s.path[i][1];
    for (i = 0; i < s.quakes.length; i++) if (s.quakes[i].trig) tr++;
    return stage + ": " + s.path.length + " " + s.fails.length + " " + s.quakes.length + " " + tr + " " +
      sum.toFixed(9) + " " + s.lim[s.lim.length - 1][1].toFixed(9) + " " + meanRelief(s.path).toFixed(9);
  }

  GS.relief = { DEFAULTS: DEFAULTS, simulate: simulate, makeQuakes: makeQuakes, ceilingAt: ceilingAt, meanRelief: meanRelief, fingerprint: fingerprint };
})(typeof window !== "undefined" ? window : this);
