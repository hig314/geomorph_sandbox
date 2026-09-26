/*
 * core/forcing.js — time forcing shapes and spatial uplift patterns.
 *
 * Pure functions, no DOM. Units everywhere in geomorph_sandbox: metres, years
 * (1 mm/yr = 1e-3 m/yr; the UI converts for display, see core/units.js).
 *
 * Ported from legacy/uplift_sandbox/model.js (time shapes, with the peak factored
 * out so the same shapes drive U(t) and ELA(t)) and legacy/terrain_sandbox/
 * terrain-model.js upliftFactor (spatial patterns, in normalised coordinates).
 *
 * Exposed as GS.forcing (browser) or module.exports (node). Mirrored in py/gs_core.py.
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});

  // ---- Time shapes s(t, p) in [0, 1] --------------------------------------
  // Event shapes use p.duration (the event window [0, duration]); cyclic shapes
  // use p.period. A series is base + peak * s(t).
  var shapes = {
    // Gaussian bell peaking at T/2; sigma = T/6 so ±3σ fits the window.
    bell: function (t, p) {
      var T = p.duration, mu = T / 2, sigma = T / 6, z = (t - mu) / sigma;
      return Math.exp(-0.5 * z * z);
    },
    // Upper half of a circle mapped onto the window.
    arc: function (t, p) {
      var x = (t - p.duration / 2) / (p.duration / 2);
      if (x <= -1 || x >= 1) return 0;
      return Math.sqrt(1 - x * x);
    },
    // Square pulse: on during the event, then off.
    constant: function (t, p) {
      return t > p.duration ? 0 : 1;
    },
    // Linear ramp up then down.
    triangle: function (t, p) {
      var x = Math.abs((t - p.duration / 2) / (p.duration / 2));
      return Math.max(0, 1 - x);
    },
    // Trapezoid: 25 % ramp up, 50 % hold, 25 % ramp down.
    plateau: function (t, p) {
      var T = p.duration, r = T * 0.25;
      if (t <= 0 || t >= T) return 0;
      if (t < r) return t / r;
      if (t > T - r) return (T - t) / r;
      return 1;
    },
    // Tectonic pulse: sharp rise over the first 10 %, exponential decay (τ = 0.3 T).
    pulse: function (t, p) {
      var T = p.duration;
      if (t <= 0) return 0;
      var rise = T * 0.1;
      if (t < rise) return t / rise;
      return Math.exp(-(t - rise) / (T * 0.3));
    },
    // Cyclic shapes (for ELA glacial cycles). Phase 0 at t = 0.
    sine: function (t, p) {
      return 0.5 * (1 - Math.cos((2 * Math.PI * t) / p.period));
    },
    // Asymmetric sawtooth: slow build over 90 % of the period, fast collapse.
    sawtooth: function (t, p) {
      var ph = t / p.period; ph -= Math.floor(ph);
      return ph < 0.9 ? ph / 0.9 : (1 - ph) / 0.1;
    }
  };
  var SHAPE_NAMES = ["bell", "plateau", "pulse", "arc", "triangle", "constant", "sine", "sawtooth"];

  // Build f(t) = base + peak * shape(t). spec: {shape, peak, base, duration, period}.
  function makeSeries(spec) {
    var p = {
      shape: spec.shape || "bell",
      peak: spec.peak != null ? spec.peak : 1,
      base: spec.base != null ? spec.base : 0,
      duration: spec.duration != null ? spec.duration : 1e6,
      period: spec.period != null ? spec.period : 1e5
    };
    var s = shapes[p.shape] || shapes.bell;
    var f = function (t) { return p.base + p.peak * s(t, p); };
    f.spec = p;
    return f;
  }

  // Trapezoid integral of f over [t0, t1] with `steps` panels (cumulative forcing).
  function integrate(f, t0, t1, steps) {
    steps = steps || 1000;
    var dt = (t1 - t0) / steps, sum = 0.5 * (f(t0) + f(t1));
    for (var i = 1; i < steps; i++) sum += f(t0 + i * dt);
    return sum * dt;
  }

  // ---- Spatial patterns: factor in [0, 1] at normalised (u, v) in [0, 1]² ----
  // "uniform" 1 everywhere; "ramp" 0 at the rim rising to 1 over `ramp` (fraction
  // of the domain); "gaussian" peak at the centre, σ = 0.18 of the domain;
  // "tilt" linear in u. For a 1-D profile pass v = 0.5.
  function spatialFactor(u, v, pattern, ramp) {
    if (pattern === "gaussian") {
      var ru = u - 0.5, rv = v - 0.5, sig = 0.18;
      return Math.exp(-(ru * ru + rv * rv) / (2 * sig * sig));
    }
    if (pattern === "tilt") return u;
    if (pattern === "ramp") {
      var dE = Math.min(u, v, 1 - u, 1 - v);
      return ramp > 0 ? (dE < ramp ? dE / ramp : 1) : 1;
    }
    return 1; // uniform
  }
  var PATTERN_NAMES = ["uniform", "ramp", "gaussian", "tilt"];

  var api = {
    shapes: shapes,
    SHAPE_NAMES: SHAPE_NAMES,
    makeSeries: makeSeries,
    integrate: integrate,
    spatialFactor: spatialFactor,
    PATTERN_NAMES: PATTERN_NAMES
  };
  GS.forcing = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : this);
