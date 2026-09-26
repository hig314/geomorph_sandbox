/*
 * core/lithology.js — erosion-resistance field in a MATERIAL frame.
 *
 * Rock moves up through the surface, so strength is defined on material
 * coordinates: the material height zm of a surface point is its current elevation
 * minus the cumulative uplift it has experienced,
 *
 *     zm = z − U_cum(x, y).
 *
 * A body defined in (x, y, zm) therefore rides up with uplift and is exposed when
 * erosion has removed everything above it. Each tier stores U_cum (a scalar in d1,
 * a profile in d2, a field in d3) and evaluates the field at the surface each step.
 *
 * A generator returns r(x, y, zm) in [0, 1]: 0 = background rock, 1 = fully
 * resistant body, with an optional soft transition. The laws consume r through
 *     erodibilityFactor(r, c) = 1 / (1 + r (c − 1))   → multiplies K, K_g, D
 *     strengthFactor(r, c)    = 1 + r (c − 1)         → multiplies cohesion, tan φ
 * so a contrast c means the resistant body erodes c× slower / is c× stronger,
 * and every process has its own contrast.
 *
 * Coordinates x, y, zm are in metres. Pure functions, no DOM.
 * Exposed as GS.lithology (browser) or module.exports (node). Mirrored in py/gs_core.py.
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  var noise = GS.noise || (typeof require !== "undefined" ? require("./noise.js") : null);

  function materialZ(z, ucum) { return z - ucum; }

  // Smooth window: 1 inside [lo, hi], 0 outside, linear over `soft` m at each edge.
  function window1(v, lo, hi, soft) {
    if (soft > 0) {
      var a = (v - lo) / soft, b = (hi - v) / soft;
      if (a < 0) a = 0; else if (a > 1) a = 1;
      if (b < 0) b = 0; else if (b > 1) b = 1;
      return a < b ? a : b;
    }
    return v >= lo && v <= hi ? 1 : 0;
  }
  function ramp01(v) { return v < 0 ? 0 : (v > 1 ? 1 : v); }

  // Generators. Each takes a spec and returns r(x, y, zm).
  var generators = {
    none: function () { return function () { return 0; }; },

    // Horizontal resistant layer: top at zm = top, thickness `thick`.
    layer: function (s) {
      var top = s.top || 0, thick = s.thick || 200, soft = s.soft || 0;
      return function (x, y, zm) { return window1(zm, top - thick, top, soft); };
    },

    // Tilted slab: layer whose top rises with x and y at gradients dipX, dipY (m/m).
    slab: function (s) {
      var top = s.top || 0, thick = s.thick || 200, soft = s.soft || 0;
      var dipX = s.dipX || 0, dipY = s.dipY || 0;
      return function (x, y, zm) {
        var t = top + dipX * x + dipY * y;
        return window1(zm, t - thick, t, soft);
      };
    },

    // Vertical dike: band of width `width` centred on x = x0 (independent of zm, y).
    dike: function (s) {
      var x0 = s.x0 || 0, w = s.width || 100, soft = s.soft || 0;
      return function (x) { return window1(x, x0 - w / 2, x0 + w / 2, soft); };
    },

    // Pluton blob: thresholded 3-D fBm on (x, y, zm) / scale. `threshold` in (0,1)
    // sets the volume fraction; `soft` (in noise units) widens the boundary.
    blob: function (s) {
      var seed = s.seed | 0, scale = s.scale || 2000, thr = s.threshold != null ? s.threshold : 0.6;
      var oct = s.octaves || 3, pers = s.persistence != null ? s.persistence : 0.5, soft = s.soft != null ? s.soft : 0.05;
      return function (x, y, zm) {
        var v = noise.fbm3(x / scale, y / scale, zm / scale, seed, oct, pers, 2);
        return soft > 0 ? ramp01((v - thr) / soft) : (v >= thr ? 1 : 0);
      };
    },

    // Persistent 2-D patch map: thresholded 2-D fBm on (x, y), independent of zm —
    // a resistant body of unlimited depth (the "infinitely thick" end-member).
    patch: function (s) {
      var seed = s.seed | 0, scale = s.scale || 2000, thr = s.threshold != null ? s.threshold : 0.6;
      var oct = s.octaves || 3, pers = s.persistence != null ? s.persistence : 0.5, soft = s.soft != null ? s.soft : 0.05;
      return function (x, y) {
        var v = noise.fbm2(x / scale, y / scale, seed, oct, pers, 2);
        return soft > 0 ? ramp01((v - thr) / soft) : (v >= thr ? 1 : 0);
      };
    }
  };
  var GENERATOR_NAMES = ["none", "layer", "slab", "dike", "blob", "patch"];

  // make(spec) → r(x, y, zm). spec.type selects the generator. An array of specs
  // combines by max (union of bodies).
  function make(spec) {
    if (Array.isArray(spec)) {
      var fns = spec.map(make);
      return function (x, y, zm) {
        var r = 0;
        for (var i = 0; i < fns.length; i++) { var v = fns[i](x, y, zm); if (v > r) r = v; }
        return r;
      };
    }
    var g = generators[(spec && spec.type) || "none"] || generators.none;
    return g(spec || {});
  }

  // Multipliers applied by the process laws. contrast c ≥ 1.
  function erodibilityFactor(r, c) { return 1 / (1 + r * (c - 1)); }
  function strengthFactor(r, c) { return 1 + r * (c - 1); }

  var api = {
    materialZ: materialZ,
    generators: generators,
    GENERATOR_NAMES: GENERATOR_NAMES,
    make: make,
    erodibilityFactor: erodibilityFactor,
    strengthFactor: strengthFactor
  };
  GS.lithology = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : this);
