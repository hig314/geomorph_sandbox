/*
 * render.js — DEM -> canvas. Hillshade (directional lighting) blended with a
 * hypsometric (elevation) colormap. Pure rendering; depends only on a model
 * with {state:{h,n,dx}} and a min/max range. window.TerrainRender.
 */
(function (root) {
  "use strict";

  // Simple hypsometric ramp: low = green, mid = brown, high = white.
  var STOPS = [
    { t: 0.0, c: [60, 110, 70] },
    { t: 0.35, c: [120, 150, 80] },
    { t: 0.6, c: [150, 120, 80] },
    { t: 0.8, c: [120, 90, 70] },
    { t: 1.0, c: [245, 245, 245] }
  ];

  function rampLinear(t) {
    if (t <= 0) return STOPS[0].c;
    if (t >= 1) return STOPS[STOPS.length - 1].c;
    for (var i = 1; i < STOPS.length; i++) {
      if (t <= STOPS[i].t) {
        var a = STOPS[i - 1], b = STOPS[i];
        var f = (t - a.t) / (b.t - a.t);
        return [
          a.c[0] + (b.c[0] - a.c[0]) * f,
          a.c[1] + (b.c[1] - a.c[1]) * f,
          a.c[2] + (b.c[2] - a.c[2]) * f
        ];
      }
    }
    return STOPS[STOPS.length - 1].c;
  }

  // Smoothed 256-entry lookup table. The piecewise-linear ramp has slope
  // discontinuities at the stops, which show as Mach-band "terraces" at fixed
  // elevations on smooth slopes (false horizontal beds). Box-blurring the LUT
  // removes those discontinuities while keeping the green→brown→white character.
  var LUT = (function () {
    var N = 256, cur = [], i, k;
    for (i = 0; i < N; i++) cur.push(rampLinear(i / (N - 1)));
    for (var pass = 0; pass < 220; pass++) {
      var prev = cur.map(function (c) { return c.slice(); });
      for (i = 1; i < N - 1; i++) for (k = 0; k < 3; k++) cur[i][k] = (prev[i - 1][k] + prev[i][k] + prev[i + 1][k]) / 3;
    }
    return cur;
  })();
  function ramp(t) {
    if (t < 0) t = 0; else if (t > 1) t = 1;
    return LUT[(t * 255) | 0];
  }

  // Light from the upper-left (image coords: -x = left, -y = up), 45° above horizon.
  var alt = (45 * Math.PI) / 180;
  var Lx = -Math.cos(alt) * Math.SQRT1_2;
  var Ly = -Math.cos(alt) * Math.SQRT1_2;
  var Lz = Math.sin(alt);

  // Render model height field into a canvas.
  // opts: {min, max, vertExag, mode}. mode ∈ "shaded" (color×hillshade, default),
  // "hillshade" (grayscale relief shading only), "color" (hypsometric, no
  // shading), "flow" (drainage-area / channel network over a hillshade base).
  function render(canvas, model, opts) {
    opts = opts || {};
    var mode = opts.mode || "shaded";
    var area = model.state.area; // drainage area (m^2), may be all zero
    var logMax = Math.log10(model.state.n * model.state.n); // for flow normalization
    var cellArea = model.state.dx * model.state.dx;
    var h = model.state.h, n = model.state.n, dx = model.state.dx;
    var showIce = opts.showIce !== false; // false → render bare substrate
    var iw = showIce ? 1 : 0;             // ice weight in the surface used for shading
    var iceArr = model.state.ice; // glacier ice thickness (may be all zero)
    var ICE_MIN = 2; // m of ice to render as glacier
    var lakeArr = model.state.lake; // ponded-water depth (may be undefined/zero)
    var LAKE_MIN = 1.5; // m of water to render as a lake (filters epsilon-fill)
    var min = opts.min, max = opts.max;
    if (min == null || max == null) {
      min = Infinity; max = -Infinity;
      for (var q = 0; q < h.length; q++) {
        if (h[q] < min) min = h[q];
        if (h[q] > max) max = h[q];
      }
    }
    var range = max - min || 1;
    var ve = opts.vertExag != null ? opts.vertExag : 2;

    if (canvas.width !== n || canvas.height !== n) {
      canvas.width = n;
      canvas.height = n;
    }
    var ctx = canvas.getContext("2d");
    var img = ctx.createImageData(n, n);
    var data = img.data;

    for (var j = 0; j < n; j++) {
      var jm = j > 0 ? j - 1 : 0;
      var jp = j < n - 1 ? j + 1 : n - 1;
      for (var i = 0; i < n; i++) {
        var im = i > 0 ? i - 1 : 0;
        var ip = i < n - 1 ? i + 1 : n - 1;
        var idx = j * n + i;

        // Surface normal from central differences of the TOP surface (bed + ice),
        // so glacier topography is shaded too. (vertical exaggeration ve)
        var dzdx = ((h[j * n + ip] + iw * iceArr[j * n + ip]) - (h[j * n + im] + iw * iceArr[j * n + im])) / (2 * dx) * ve;
        var dzdy = ((h[jp * n + i] + iw * iceArr[jp * n + i]) - (h[jm * n + i] + iw * iceArr[jm * n + i])) / (2 * dx) * ve;
        var nx = -dzdx, ny = -dzdy, nz = 1;
        var len = Math.sqrt(nx * nx + ny * ny + nz * nz);
        var shade = (nx * Lx + ny * Ly + nz * Lz) / len; // [-1,1]
        shade = 0.25 + 0.75 * Math.max(0, shade); // ambient + diffuse

        var o = idx * 4, r0, g0, b0;
        if (mode === "flow") {
          // Grayscale hillshade base; channels painted blue, intensity by log area.
          var base = shade * 230;
          var cells = (area ? area[idx] : 0) / cellArea;
          if (cells > 4) {
            var tt = Math.min(1, Math.log10(cells) / logMax);
            r0 = base * (1 - tt);
            g0 = base * (1 - 0.55 * tt);
            b0 = base * (1 - tt) + 255 * tt;
          } else {
            r0 = g0 = b0 = base;
          }
        } else if (mode === "hillshade") {
          var g = shade * 255;
          r0 = g0 = b0 = g;
        } else if (mode === "color") {
          var c1 = ramp((h[idx] - min) / range);
          r0 = c1[0]; g0 = c1[1]; b0 = c1[2];
        } else {
          var c2 = ramp((h[idx] - min) / range);
          r0 = c2[0] * shade; g0 = c2[1] * shade; b0 = c2[2] * shade;
        }
        // Lakes: ponded water, blue and darkening with depth (drawn over terrain,
        // under ice). Skipped in flow mode, where blue means channels.
        if (lakeArr && mode !== "flow" && lakeArr[idx] > LAKE_MIN) {
          var ld = lakeArr[idx], dm = ld > 40 ? 1 : ld / 40; // depth mix
          r0 = (70 - 45 * dm) * shade;
          g0 = (120 - 50 * dm) * shade;
          b0 = (190 - 30 * dm) * shade;
        }
        // Glacier ice: realistic white with a faint hillshade (subtle relief).
        if (showIce && iceArr[idx] > ICE_MIN) {
          var w = 255 * (0.82 + 0.18 * Math.max(0, (nx * Lx + ny * Ly + nz * Lz) / len));
          r0 = g0 = b0 = w;
        }
        data[o] = r0;
        data[o + 1] = g0;
        data[o + 2] = b0;
        data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }

  // Render a scalar field into a small canvas.
  //   type "diverging": red (negative) – white (0) – blue (positive). For dz
  //     (erosion vs deposition). Symmetric scale from max |value|.
  //   type "log": white -> purple by log magnitude. For sediment flux.
  function renderField(canvas, field, n, opts) {
    opts = opts || {};
    var type = opts.type || "diverging";
    var maxAbs = 0, maxv = 0;
    for (var q = 0; q < field.length; q++) {
      var v = field[q];
      if (v > maxv) maxv = v;
      var a = v < 0 ? -v : v;
      if (a > maxAbs) maxAbs = a;
    }
    if (maxAbs === 0) maxAbs = 1;
    if (opts.scale) maxAbs = opts.scale; // fixed diverging scale (frame-to-frame stable)
    // Log normalization reference. opts.scale fixes it for the "ice"/"log" ramps too,
    // so the colour mapping doesn't flicker as the field's own max wobbles each frame.
    var logRef = Math.log10((opts.scale != null ? opts.scale : maxv) + 1) || 1;

    if (canvas.width !== n) { canvas.width = n; canvas.height = n; }
    var ctx = canvas.getContext("2d");
    var img = ctx.createImageData(n, n);
    var d = img.data;
    for (var i = 0; i < n * n; i++) {
      var o = i * 4, val = field[i], r, g, b;
      if (type === "ice") {
        var ti = val > 0 ? Math.log10(val + 1) / logRef : 0;
        r = 255 - ti * 210; g = 255 - ti * 120; b = 255; // white -> blue
      } else if (type === "log") {
        var t = val > 0 ? Math.log10(val + 1) / logRef : 0;
        r = 255 - t * 175; g = 255 - t * 235; b = 255 - t * 60; // white -> purple
      } else {
        var s = val / maxAbs; // [-1,1]
        if (s >= 0) { r = 255 - s * 215; g = 255 - s * 140; b = 255; } // white->blue (deposition)
        else { r = 255; g = 255 + s * 175; b = 255 + s * 215; }        // white->red (erosion)
      }
      d[o] = r; d[o + 1] = g; d[o + 2] = b; d[o + 3] = 255;
    }
    ctx.putImageData(img, 0, 0);
  }

  root.TerrainRender = { render: render, renderField: renderField, ramp: ramp };
})(typeof window !== "undefined" ? window : this);
