/*
 * terrain-model.js — Pure 2D terrain-evolution core (no DOM, no D3).
 *
 * Phase 1: initial-topography generators + hillslope creep (linear and the
 * nonlinear Roering 2001 law) + uniform uplift. Float32 height field is the
 * source of truth; quantize to bytes only for I/O. Later phases add fluvial,
 * glacial, and SCOOPS3D-style deep failure as additional kernels.
 *
 * Units: elevation in metres, distance in metres, time in years.
 *
 * Exposed as window.TerrainModel (browser) or module.exports (node).
 */
(function (root) {
  "use strict";

  // ---- Seeded PRNG + integer hash (reproducible topography) --------------
  function mulberry32(a) {
    return function () {
      a |= 0;
      a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Deterministic [0,1) value at integer lattice point (ix, iy) for a seed.
  function hash2(ix, iy, seed) {
    var h = Math.imul(ix, 374761393) + Math.imul(iy, 668265263) + Math.imul(seed, 982451653);
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    h = h ^ (h >>> 16);
    return (h >>> 0) / 4294967296;
  }

  function smoothstep(t) {
    return t * t * (3 - 2 * t);
  }

  // Single-octave value noise sampled at (x, y) in lattice units.
  function valueNoise(x, y, seed) {
    var x0 = Math.floor(x), y0 = Math.floor(y);
    var fx = smoothstep(x - x0), fy = smoothstep(y - y0);
    var v00 = hash2(x0, y0, seed), v10 = hash2(x0 + 1, y0, seed);
    var v01 = hash2(x0, y0 + 1, seed), v11 = hash2(x0 + 1, y0 + 1, seed);
    var a = v00 + (v10 - v00) * fx;
    var b = v01 + (v11 - v01) * fx;
    return a + (b - a) * fy; // [0,1]
  }

  // Fractal Brownian motion (summed octaves of value noise) -> [0,1].
  function fbm(x, y, seed, octaves, persistence, lacunarity) {
    var amp = 1, freq = 1, sum = 0, norm = 0;
    for (var o = 0; o < octaves; o++) {
      sum += amp * valueNoise(x * freq, y * freq, seed + o * 1013);
      norm += amp;
      amp *= persistence;
      freq *= lacunarity;
    }
    return sum / norm;
  }

  // ---- Topography generators --------------------------------------------
  // Each fills the height array h (length n*n) in place. opts.seed makes them
  // reproducible. Returned heights are in metres.

  var generators = {
    // Uniform tilted plane: rises across x. Control case — diffusion barely acts.
    plane: function (h, n, dx, opts) {
      var grad = opts.gradient != null ? opts.gradient : 0.05; // rise/run
      for (var j = 0; j < n; j++) {
        for (var i = 0; i < n; i++) {
          h[j * n + i] = grad * (i * dx);
        }
      }
    },

    // Fractal hills: multi-octave value noise. The interesting case for creep.
    hills: function (h, n, dx, opts) {
      var seed = opts.seed | 0;
      var relief = opts.relief != null ? opts.relief : 800; // peak-to-trough (m)
      var baseFreq = opts.baseFreq != null ? opts.baseFreq : 4; // cycles across domain
      var octaves = opts.octaves || 6;
      var persistence = opts.persistence != null ? opts.persistence : 0.5;
      for (var j = 0; j < n; j++) {
        for (var i = 0; i < n; i++) {
          var x = (i / n) * baseFreq;
          var y = (j / n) * baseFreq;
          h[j * n + i] = fbm(x, y, seed, octaves, persistence, 2) * relief;
        }
      }
    },

    // Scattered cones: max of N random cones. Sharp features that round under creep.
    cones: function (h, n, dx, opts) {
      var seed = opts.seed | 0;
      var count = opts.count || 25;
      var maxHeight = opts.maxHeight != null ? opts.maxHeight : 600; // m
      var rMin = opts.rMin != null ? opts.rMin : 0.05; // fraction of domain
      var rMax = opts.rMax != null ? opts.rMax : 0.18;
      var rand = mulberry32(seed);
      var cones = [];
      for (var c = 0; c < count; c++) {
        cones.push({
          cx: rand() * n,
          cy: rand() * n,
          r: (rMin + rand() * (rMax - rMin)) * n,
          h: maxHeight * (0.4 + 0.6 * rand())
        });
      }
      for (var j = 0; j < n; j++) {
        for (var i = 0; i < n; i++) {
          var z = 0;
          for (var k = 0; k < cones.length; k++) {
            var ck = cones[k];
            var d = Math.sqrt((i - ck.cx) * (i - ck.cx) + (j - ck.cy) * (j - ck.cy));
            var zc = ck.h * Math.max(0, 1 - d / ck.r);
            if (zc > z) z = zc;
          }
          h[j * n + i] = z;
        }
      }
    }
  };

  // ---- Hillslope creep kernels ------------------------------------------
  // Linear diffusion: dz/dt = D * laplacian(z). Explicit 5-point stencil.
  // Writes into `out`. Boundaries: zero-flux (reflect) via clamped neighbours.
  function creepLinearStep(h, out, n, dx, D, dt) {
    var inv = 1 / (dx * dx);
    for (var j = 0; j < n; j++) {
      var jm = j > 0 ? j - 1 : 0;
      var jp = j < n - 1 ? j + 1 : n - 1;
      for (var i = 0; i < n; i++) {
        var im = i > 0 ? i - 1 : 0;
        var ip = i < n - 1 ? i + 1 : n - 1;
        var c = h[j * n + i];
        var lap = (h[j * n + im] + h[j * n + ip] + h[jm * n + i] + h[jp * n + i] - 4 * c) * inv;
        out[j * n + i] = c + D * dt * lap;
      }
    }
  }

  // Nonlinear Roering (2001): flux q = -D * grad / (1 - (|grad|/Sc)^2), diverging
  // as slope -> Sc. Finite-volume: compute fluxes on the +x and +y faces, then
  // dz/dt = -div(q). |grad|/Sc is clamped just below 1 for numerical safety.
  function creepNonlinearStep(h, out, n, dx, D, Sc, dt) {
    var qx = creepNonlinearStep._qx;
    if (!qx || qx.length !== n * n) qx = creepNonlinearStep._qx = new Float32Array(n * n);
    var qy = creepNonlinearStep._qy;
    if (!qy || qy.length !== n * n) qy = creepNonlinearStep._qy = new Float32Array(n * n);
    var Sc2 = Sc * Sc;

    // Flux across the face between cell (i,j) and (i+1,j), stored at (i,j).
    for (var j = 0; j < n; j++) {
      for (var i = 0; i < n; i++) {
        var idx = j * n + i;
        // +x face
        if (i < n - 1) {
          var sx = (h[idx + 1] - h[idx]) / dx; // slope component
          var fx = 1 - (sx * sx) / Sc2;
          if (fx < 0.02) fx = 0.02; // clamp near-critical to keep finite
          qx[idx] = -D * sx / fx;
        } else qx[idx] = 0;
        // +y face
        if (j < n - 1) {
          var sy = (h[idx + n] - h[idx]) / dx;
          var fy = 1 - (sy * sy) / Sc2;
          if (fy < 0.02) fy = 0.02;
          qy[idx] = -D * sy / fy;
        } else qy[idx] = 0;
      }
    }

    // dz/dt = -(dqx/dx + dqy/dy), using fluxes on the cell's two -faces too.
    var invdx = 1 / dx;
    for (var j2 = 0; j2 < n; j2++) {
      for (var i2 = 0; i2 < n; i2++) {
        var id = j2 * n + i2;
        var qxIn = i2 > 0 ? qx[id - 1] : 0;
        var qyIn = j2 > 0 ? qy[id - n] : 0;
        var div = (qx[id] - qxIn) * invdx + (qy[id] - qyIn) * invdx;
        out[id] = h[id] - dt * div;
      }
    }
  }

  // ---- Fluvial incision (stream power) ----------------------------------
  // CPU flow routing + implicit FastScape incision. Pipeline per call:
  //   1. priority-flood depression fill (so every cell drains to a boundary)
  //   2. D8 steepest-descent receivers on the filled surface
  //   3. topological "stack" ordering (receivers before donors)
  //   4. drainage-area accumulation up the stack
  //   5. implicit stream-power incision  dz/dt = -K A^m S^n  on the real DEM
  // Boundary cells are fixed base level (outlets). Buffers are cached on the
  // function object and reused between calls.
  //
  // Refs: O'Callaghan & Mark (1984) D8; Barnes et al. (2014) priority-flood;
  //       Braun & Willett (2013) O(n) implicit stream-power solver.
  function FluvialState(n) {
    var size = n * n;
    this.filled = new Float32Array(size);
    this.rs = new Float32Array(size);   // routing surface (bed + ice): water flows on this
    this.lake = new Float32Array(size); // ponded-water depth (filled − routing surface)
    this.receiver = new Int32Array(size);
    this.dist = new Float32Array(size);
    this.area = new Float32Array(size);
    this.stack = new Int32Array(size);
    this.ndon = new Int32Array(size);
    this.delta = new Int32Array(size + 1);
    this.donors = new Int32Array(size);
    this.closed = new Uint8Array(size);
    this.floodOrder = new Int32Array(size); // cells in increasing filled order
    this.floodCount = 0;
    this.h0 = new Float32Array(size); // elevation snapshot at start of fluvial step
    this.Qs = new Float64Array(size); // sediment flux (m^3/yr), Gauss-Seidel scratch
    this.sedFlux = new Float32Array(size); // sediment flux per cell (for display)
    this.lastIters = 1; // Gauss-Seidel iterations used last call
    // binary min-heap (parallel arrays of key/value), capacity = size
    this.heapKey = new Float64Array(size);
    this.heapVal = new Int32Array(size);
    this.heapLen = 0;
  }

  // 8-neighbour offsets and their planimetric distances (in cells).
  var NB_DX = [-1, 1, 0, 0, -1, -1, 1, 1];
  var NB_DY = [0, 0, -1, 1, -1, 1, -1, 1];
  var NB_DIAG = [1, 1, 1, 1, Math.SQRT2, Math.SQRT2, Math.SQRT2, Math.SQRT2];
  var MFD_W = new Float64Array(8); // scratch: per-neighbour flow weights

  function heapPush(fs, key, val) {
    var i = fs.heapLen++;
    fs.heapKey[i] = key;
    fs.heapVal[i] = val;
    while (i > 0) {
      var par = (i - 1) >> 1;
      if (fs.heapKey[par] <= fs.heapKey[i]) break;
      var tk = fs.heapKey[par]; fs.heapKey[par] = fs.heapKey[i]; fs.heapKey[i] = tk;
      var tv = fs.heapVal[par]; fs.heapVal[par] = fs.heapVal[i]; fs.heapVal[i] = tv;
      i = par;
    }
  }
  function heapPop(fs) {
    var top = fs.heapVal[0];
    var last = --fs.heapLen;
    fs.heapKey[0] = fs.heapKey[last];
    fs.heapVal[0] = fs.heapVal[last];
    var i = 0;
    for (;;) {
      var l = 2 * i + 1, r = l + 1, sm = i;
      if (l < fs.heapLen && fs.heapKey[l] < fs.heapKey[sm]) sm = l;
      if (r < fs.heapLen && fs.heapKey[r] < fs.heapKey[sm]) sm = r;
      if (sm === i) break;
      var tk = fs.heapKey[sm]; fs.heapKey[sm] = fs.heapKey[i]; fs.heapKey[i] = tk;
      var tv = fs.heapVal[sm]; fs.heapVal[sm] = fs.heapVal[i]; fs.heapVal[i] = tv;
      i = sm;
    }
    return top;
  }

  // Domain-edge boundary, set relative to the inward neighbour (stable; no
  // feedback). Default: a low rim 1 m below the interior all around — every edge
  // is a base-level outlet, sediment exits wherever the interior slopes outward.
  // "tilt": extrapolate the regional tilt one cell out — the down-tilt edge (−x)
  // sits 1 m lower (outflow), the up-tilt edge (+x) sits 1 m higher (no outflow,
  // flow turns inward), and the cross edges (±y) have no off-grid slope.
  function setPerimeter(h, n, pattern) {
    if (pattern === "tilt") {
      for (var j = 0; j < n; j++) {
        h[j * n] = h[j * n + 1] - 1;                 // left  (down-tilt): outflow
        h[j * n + (n - 1)] = h[j * n + (n - 2)] + 1;  // right (up-tilt):  no outflow
      }
      for (var i = 1; i < n - 1; i++) {
        h[i] = h[n + i];                              // top    (cross): no off-grid slope
        h[(n - 1) * n + i] = h[(n - 2) * n + i];      // bottom (cross)
      }
      return;
    }
    var drop = 1;
    for (var i2 = 1; i2 < n - 1; i2++) {
      h[i2] = h[n + i2] - drop;
      h[(n - 1) * n + i2] = h[(n - 2) * n + i2] - drop;
    }
    for (var j2 = 1; j2 < n - 1; j2++) {
      h[j2 * n] = h[j2 * n + 1] - drop;
      h[j2 * n + (n - 1)] = h[j2 * n + (n - 2)] - drop;
    }
    h[0] = h[n + 1] - drop;
    h[n - 1] = h[n + (n - 2)] - drop;
    h[(n - 1) * n] = h[(n - 2) * n + 1] - drop;
    h[(n - 1) * n + (n - 1)] = h[(n - 2) * n + (n - 2)] - drop;
  }

  // Spatial uplift factor in [0,1] for cell (ix,iy). "gaussian": peak at centre,
  // smoothly → 0 at the edges. "tilt": linear plane increasing left → right.
  // "ramp" (default): 0 at the rim, ramping to 1 over `ramp` cells.
  function upliftFactor(ix, iy, n, pattern, ramp) {
    if (pattern === "gaussian") {
      var c = (n - 1) / 2, rx = ix - c, ry = iy - c, sig = n * 0.18;
      return Math.exp(-(rx * rx + ry * ry) / (2 * sig * sig));
    }
    if (pattern === "tilt") return ix / (n - 1);
    var dE = ix;
    if (iy < dE) dE = iy;
    if (n - 1 - ix < dE) dE = n - 1 - ix;
    if (n - 1 - iy < dE) dE = n - 1 - iy;
    return ramp > 0 ? (dE < ramp ? dE / ramp : 1) : 1;
  }

  function fluvialStep(h, n, dx, K, m, nexp, dt, G, mfdExp, perimPattern, ice, ela, meltw, glacSed) {
    G = G != null ? G : 1; // deposition coefficient (0 = pure detachment-limited)
    var P = mfdExp != null ? mfdExp : 4; // MFD flow-concentration exponent (1 = max spread)
    var fs = fluvialStep._fs;
    if (!fs || fs.receiver.length !== n * n) fs = fluvialStep._fs = new FluvialState(n);
    var size = n * n;
    var filled = fs.filled, closed = fs.closed;
    var EPS = 1e-3; // m, enforces a tiny gradient across filled flats

    setPerimeter(h, n, perimPattern); // low-rim outlet, or tilt-extrapolated edges

    // Water flows on the SURFACE (bed + ice), not the bare bedrock — so a glacier is
    // a real barrier: rivers route around it or pond against it instead of vanishing
    // under it. The priority-flood fill of this surface IS the ponded-water level, so
    // lakes (ice-dammed and bedrock basins alike) fall out for free as filled − rs.
    var rs = fs.rs;
    for (var z0 = 0; z0 < size; z0++) rs[z0] = ice ? h[z0] + ice[z0] : h[z0];

    // --- 1. priority-flood fill (+epsilon) ---
    filled.set(rs);
    closed.fill(0);
    fs.heapLen = 0;
    for (var i = 0; i < n; i++) {
      pushBoundary(fs, filled, closed, i, n);                 // top row
      pushBoundary(fs, filled, closed, (n - 1) * n + i, n);    // bottom row
    }
    for (var j = 1; j < n - 1; j++) {
      pushBoundary(fs, filled, closed, j * n, n);              // left col
      pushBoundary(fs, filled, closed, j * n + (n - 1), n);    // right col
    }
    fs.floodCount = 0;
    while (fs.heapLen > 0) {
      var c = heapPop(fs);
      fs.floodOrder[fs.floodCount++] = c; // increasing filled order
      var cx = c % n, cy = (c / n) | 0;
      for (var d = 0; d < 8; d++) {
        var nx = cx + NB_DX[d], ny = cy + NB_DY[d];
        if (nx < 0 || ny < 0 || nx >= n || ny >= n) continue;
        var nb = ny * n + nx;
        if (closed[nb]) continue;
        closed[nb] = 1;
        var lim = filled[c] + EPS;
        if (filled[nb] < lim) filled[nb] = lim; // raise pits to keep descent
        heapPush(fs, filled[nb], nb);
      }
    }

    // Ponded water = how much the flood had to raise each cell above the surface.
    var lake = fs.lake;
    for (var z1 = 0; z1 < size; z1++) { var dl = filled[z1] - rs[z1]; lake[z1] = dl > 0 ? dl : 0; }

    // --- 2. D8 receivers on the filled surface (boundary drains to self) ---
    var receiver = fs.receiver, dist = fs.dist;
    for (var idx = 0; idx < size; idx++) {
      var ix = idx % n, iy = (idx / n) | 0;
      if (ix === 0 || iy === 0 || ix === n - 1 || iy === n - 1) {
        receiver[idx] = idx; dist[idx] = 0; continue; // base level
      }
      var best = idx, bestSlope = 0, bestDist = dx;
      for (var dd = 0; dd < 8; dd++) {
        var nb2 = (iy + NB_DY[dd]) * n + (ix + NB_DX[dd]);
        var L = NB_DIAG[dd] * dx;
        var slope = (filled[idx] - filled[nb2]) / L;
        if (slope > bestSlope) { bestSlope = slope; best = nb2; bestDist = L; }
      }
      receiver[idx] = best; dist[idx] = bestDist;
    }

    // --- 3. topological stack (receivers before donors) ---
    var ndon = fs.ndon, delta = fs.delta, donors = fs.donors, stack = fs.stack;
    ndon.fill(0);
    for (var a = 0; a < size; a++) if (receiver[a] !== a) ndon[receiver[a]]++;
    delta[size] = size;
    for (var b = size - 1; b >= 0; b--) delta[b] = delta[b + 1] - ndon[b];
    // fill donor CSR using a moving cursor (reuse area as scratch cursor)
    var cursor = fs.area; // temporarily
    for (var q = 0; q < size; q++) cursor[q] = delta[q];
    for (var e = 0; e < size; e++) {
      var rr = receiver[e];
      if (rr !== e) { donors[cursor[rr]] = e; cursor[rr]++; }
    }
    var nstack = 0;
    // iterative DFS from each base-level node
    var work = fs.heapVal; // reuse as work stack (heap is idle now)
    var wp = 0;
    for (var s0 = 0; s0 < size; s0++) {
      if (receiver[s0] !== s0) continue;
      work[wp++] = s0;
      while (wp > 0) {
        var node = work[--wp];
        stack[nstack++] = node;
        for (var dptr = delta[node]; dptr < delta[node + 1]; dptr++) work[wp++] = donors[dptr];
      }
    }

    // --- 4. drainage area, multiple-flow-direction (Holmgren) ---
    // Distribute each cell's accumulated area to ALL lower neighbours, weighted
    // by slope^P. Processed in decreasing filled order (highest first) so each
    // cell's upstream area is complete before it is shared out. P controls
    // concentration: P=1 spreads widely (breaks the parallel-stripe artifact of
    // single-flow routing), large P approaches steepest-descent. This is the
    // drainage area used by stream power; the bed update still follows the single
    // steepest receiver for the stable implicit solver.
    // Each cell seeds its own rainfall (area = cell). Above the firn line precip
    // falls as snow (feeds the ice, no runoff), so rain = 0 there: A→0 upslope →
    // no discharge and no incision above the ELA. (Meltwater below the ELA is not
    // yet routed — a deferred refinement.)
    var area = fs.area, cell = dx * dx;
    var rainMask = (ice && ela != null);
    for (var z = 0; z < size; z++) {
      area[z] = (rainMask && (h[z] + ice[z]) > ela) ? 0 : cell;
      if (meltw) area[z] += meltw[z] * cell; // glacial meltwater adds discharge (rain-equiv)
    }
    var fo = fs.floodOrder, fcount = fs.floodCount;
    for (var fi = fcount - 1; fi >= 0; fi--) {
      var c2 = fo[fi];
      var Acur = area[c2];
      var cx2 = c2 % n, cy2 = (c2 / n) | 0;
      var fc = filled[c2], sumw = 0;
      for (var d2 = 0; d2 < 8; d2++) {
        var nx2 = cx2 + NB_DX[d2], ny2 = cy2 + NB_DY[d2];
        if (nx2 < 0 || ny2 < 0 || nx2 >= n || ny2 >= n) { MFD_W[d2] = 0; continue; }
        var s2 = (fc - filled[ny2 * n + nx2]) / (NB_DIAG[d2] * dx);
        if (s2 > 0) {
          var w = s2;
          for (var kk = 1; kk < P; kk++) w *= s2; // s2^P (P integer)
          MFD_W[d2] = w; sumw += w;
        } else MFD_W[d2] = 0;
      }
      if (sumw > 0) {
        var share = Acur / sumw;
        for (var d3 = 0; d3 < 8; d3++) {
          if (MFD_W[d3] > 0) area[(cy2 + NB_DY[d3]) * n + (cx2 + NB_DX[d3])] += MFD_W[d3] * share;
        }
      }
    }


    // --- 5. implicit erosion + deposition (Yuan et al. 2019), Gauss-Seidel ---
    // Solves dh/dt = -K A^m S^n + G*Qs/A implicitly in one framework. Per
    // iteration: (a) accumulate sediment flux Qs up the stack from the current
    // estimate, (b) update every cell implicitly against its (already-updated)
    // receiver plus the deposition source G*Qs/A. Unconditionally stable, so a
    // large dt needs no sub-stepping; iterations scale with G (G=0 -> 1 pass =
    // plain FastScape). Under-relaxation (OMEGA) tames the deposition feedback.
    // Mass-conserving by construction (eroded = deposited + exported).
    var h0 = fs.h0, Qs = fs.Qs, sedFlux = fs.sedFlux;
    h0.set(h);
    var linear = Math.abs(nexp - 1) < 1e-9;
    var OMEGA = 0.6;
    var maxIter = G <= 0 ? 1 : Math.min(25, Math.ceil(2 + G * 20));
    var iters = 0;
    for (var it = 0; it < maxIter; it++) {
      iters = it + 1;
      // (a) sediment flux: Qs_out = Qs_in + (eroded-deposited) volume rate
      for (var z2 = 0; z2 < size; z2++) Qs[z2] = 0;
      for (var sp = nstack - 1; sp >= 0; sp--) {
        var nd = stack[sp];
        // qout = inflow + glacial-sediment source + (eroded − deposited) this cell
        var qout = Qs[nd] + (glacSed ? glacSed[nd] : 0) + ((h0[nd] - h[nd]) * cell) / dt;
        var rA = receiver[nd];
        if (rA !== nd) Qs[rA] += qout;
      }
      // (b) implicit elevation update down the stack
      var maxd = 0;
      for (var t = 0; t < nstack; t++) {
        var node = stack[t];
        var r = receiver[node];
        if (r === node) continue; // base level
        if (area[node] <= 0) continue; // no rain (above firn line) → no fluvial work
        if (ice && ice[node] > 2) continue; // ice-covered: glacial, not subaerial fluvial
        var hr = h[r];
        var dep = (dt * G * Qs[node]) / area[node];
        if (dep < 0) dep = 0;
        // no bedrock incision into/within a lake; deposition fills those cells
        var submerged = filled[node] - h0[node] > 1 || filled[r] - hr > 1;
        var nh;
        if (h0[node] > hr && !submerged) {
          var Kp = K * Math.pow(area[node], m);
          if (linear) {
            var f = (Kp * dt) / dist[node];
            nh = (h0[node] + f * hr + dep) / (1 + f);
          } else {
            var z = h[node];
            for (var nit = 0; nit < 20; nit++) {
              var s = (z - hr) / dist[node];
              if (s < 0) s = 0;
              var g = z - h0[node] + Kp * dt * Math.pow(s, nexp) - dep;
              var dg = 1 + (Kp * dt * nexp * Math.pow(s, nexp - 1)) / dist[node];
              var dzz = g / dg;
              z -= dzz;
              if (Math.abs(dzz) < 1e-4) break;
            }
            nh = z;
          }
          if (nh < hr) nh = hr;
        } else {
          nh = h0[node] + dep; // basin/flat: deposit only
        }
        nh = h[node] + OMEGA * (nh - h[node]); // under-relaxation
        var dd = nh - h[node];
        if (dd < 0) dd = -dd;
        if (dd > maxd) maxd = dd;
        h[node] = nh;
      }
      if (G <= 0 || maxd < 0.02) break;
    }
    fs.lastIters = iters;
    for (var z3 = 0; z3 < size; z3++) sedFlux[z3] = Qs[z3]; // expose flux for display
  }
  function pushBoundary(fs, filled, closed, idx, n) {
    if (closed[idx]) return;
    closed[idx] = 1;
    heapPush(fs, filled[idx], idx);
  }

  // ---- Deep-seated landslides (SCOOPS3D-style, stochastic) ---------------
  // Scratch buffers, generation-tagged so each failure gets a fresh "floor"
  // (pre-deposit surface) and active set without clearing the whole grid.
  var LS_floor = null, LS_fgen = null, LS_inGen = null, LS_srcIdx = null, LS_srcT = null;
  var LS_active = [], LS_failId = 0;
  var LS_fosOut = { valid: false, fos: 0, sx: 0, sy: 0, blockVol: 0 };
  var SLOW_DROP = 2; // slow-slump creep: max head-drop per active step (m)
  var SLUMP_MAX = 200; // cap on simultaneously-tracked persistent slumps

  // floor of a cell for the current failure = its surface just before debris
  // arrived (captured lazily). Talus may never lower a cell below this, so only
  // the debris flows — the underlying terrain is never carved.
  function lsFloorOf(h, c) {
    if (LS_fgen[c] !== LS_failId) { LS_fgen[c] = LS_failId; LS_floor[c] = h[c]; }
    return LS_floor[c];
  }

  // Growing angle-of-repose relaxation. Debris in `active` cascades to downslope
  // neighbours until no slope exceeds repose; the active region GROWS to wherever
  // material lands (no fixed box → no clipped downslope margin). Full-excess to
  // the steepest neighbour in high-to-low order makes it cascade in few sweeps.
  function talusGrow(h, n, dx, reposeDeg, active, maxsweep) {
    var tanr = Math.tan((reposeDeg * Math.PI) / 180);
    for (var sweep = 0; sweep < maxsweep; sweep++) {
      active.sort(function (a, b) { return h[b] - h[a]; }); // high -> low
      var moved = 0;
      for (var k = 0; k < active.length; k++) {
        var c = active[k];
        var avail = h[c] - lsFloorOf(h, c);
        if (avail <= 0) continue;
        var cx = c % n, cy = (c / n) | 0;
        var best = -1, bestEx = 0;
        for (var d = 0; d < 8; d++) {
          var ni = cx + NB_DX[d], nj = cy + NB_DY[d];
          if (ni < 0 || nj < 0 || ni >= n || nj >= n) continue;
          var nb = nj * n + ni;
          var ex = h[c] - h[nb] - tanr * NB_DIAG[d] * dx;
          if (ex > bestEx) { bestEx = ex; best = nb; }
        }
        if (best >= 0) {
          var amt = bestEx; if (amt > avail) amt = avail;
          lsFloorOf(h, best);
          h[c] -= amt; h[best] += amt; moved += amt;
          if (LS_inGen[best] !== LS_failId) { LS_inGen[best] = LS_failId; active.push(best); }
        }
      }
      if (moved < 1e-2) break;
    }
  }

  // One pass of stochastic failure testing. Each test is a random trial scoop
  // (the "pre-existing weakness"); if its factor of safety is below threshold
  // (the "geometry that supports failure"), it fails. Catastrophic mode removes
  // the mass and spreads it downslope to the angle of repose; vanish mode exports
  // it. Returns {count, vol}. rng is a persistent seeded PRNG.
  // Factor of safety of one trial scoop (3D method of columns, Ordinary).
  // Driving is resolved as a VECTOR (net slide direction): each column's driving
  // W·sinα points down-dip (horizontally toward the sphere axis), so a symmetric
  // scoop on flat ground cancels (FoS huge → stable, no spurious "craters") and
  // only an oversteepened, asymmetric mass drives a net downslope failure.
  // Results are returned in the shared LS_fosOut: {valid, fos, sx, sy, blockVol}.
  function scoopFoS(h, n, dx, Cx, Cy, Cz, R, c, phi, gamma, ru, ice, gammaIce) {
    var cellA = dx * dx, tanphi = Math.tan((phi * Math.PI) / 180), R2 = R * R;
    var ic = (Cx / dx) | 0, jc = (Cy / dx) | 0, rc = ((R / dx) | 0) + 1;
    var i0 = ic - rc < 0 ? 0 : ic - rc, i1 = ic + rc > n - 1 ? n - 1 : ic + rc;
    var j0 = jc - rc < 0 ? 0 : jc - rc, j1 = jc + rc > n - 1 ? n - 1 : jc + rc;
    var Rsum = 0, dvx = 0, dvy = 0, nc = 0, blockVol = 0;
    for (var j = j0; j <= j1; j++) {
      var ddy = j * dx - Cy;
      for (var i = i0; i <= i1; i++) {
        var ddx = i * dx - Cx, d2 = ddx * ddx + ddy * ddy;
        if (d2 >= R2) continue;
        var root = Math.sqrt(R2 - d2), zslip = Cz - root, zg = h[j * n + i];
        if (zslip >= zg) continue;
        var dz = zg - zslip, cosA = root / R;
        if (cosA < 1e-3) continue;
        var sinA = Math.sqrt(1 - cosA * cosA);
        // Column weight = rock above the slip surface + overlying-ice surcharge.
        // Ice (γ_ice ≈ 9 kN/m³) is ~0.45× rock (≈ 20), so it loads slopes — adding
        // both driving (W·sinα) and frictional resistance (W·cosα·tanφ) — but with
        // less than half the punch of rock. On the gentle toe/valley (cosα≈1,
        // sinα≈0) it is almost pure resistance → buttressing; on steep upper slopes
        // (where ice is thin anyway) it would add driving. Net: valley ice stabilises
        // the slope above it, and losing the ice (deglaciation) drops the FoS
        // (paraglacial debuttressing). Pore pressure is from rock depth only.
        var W = gamma * dz * cellA;
        if (ice) W += gammaIce * ice[j * n + i] * cellA;
        var A = cellA / cosA, u = ru * gamma * dz;
        Rsum += c * A + (W * cosA - u * A) * tanphi;
        var d = Math.sqrt(d2);
        if (d > 1e-6) { var ws = W * sinA; dvx += (ws * -ddx) / d; dvy += (ws * -ddy) / d; }
        blockVol += dz * cellA;
        nc++;
      }
    }
    var den = Math.sqrt(dvx * dvx + dvy * dvy);
    LS_fosOut.blockVol = blockVol;
    if (den <= 0 || nc < 5) { LS_fosOut.valid = false; return LS_fosOut; }
    LS_fosOut.valid = true; LS_fosOut.fos = Rsum / den;
    LS_fosOut.sx = dvx / den; LS_fosOut.sy = dvy / den;
    return LS_fosOut;
  }

  // Apply a failure: remove the block above the slip surface, displace it
  // downslope by `disp` cells along the slide direction (the rotation — large
  // for catastrophic, small for a slow creep increment), then resolve with the
  // growing angle-of-repose relaxation. vanish = export instead of depositing.
  // Mass-conserving. Returns the removed volume.
  function applyScoop(h, n, dx, Cx, Cy, Cz, R, sx, sy, disp, repose, vanish) {
    var cellA = dx * dx, R2 = R * R;
    var ic = (Cx / dx) | 0, jc = (Cy / dx) | 0, rc = ((R / dx) | 0) + 1;
    var i0 = ic - rc < 0 ? 0 : ic - rc, i1 = ic + rc > n - 1 ? n - 1 : ic + rc;
    var j0 = jc - rc < 0 ? 0 : jc - rc, j1 = jc + rc > n - 1 ? n - 1 : jc + rc;
    LS_failId++;
    var nsrc = 0, removed = 0;
    for (var j = j0; j <= j1; j++) {
      var ddy = j * dx - Cy;
      for (var i = i0; i <= i1; i++) {
        var ddx = i * dx - Cx, d2 = ddx * ddx + ddy * ddy;
        if (d2 >= R2) continue;
        var zs = Cz - Math.sqrt(R2 - d2), idx = j * n + i, zg = h[idx];
        if (zs >= zg) continue;
        LS_srcIdx[nsrc] = idx; LS_srcT[nsrc] = zg - zs; nsrc++;
        removed += (zg - zs) * cellA;
        LS_fgen[idx] = LS_failId; LS_floor[idx] = zs;
        h[idx] = zs;
      }
    }
    if (removed <= 0) return 0;
    if (vanish) return removed;
    LS_active.length = 0;
    for (var k = 0; k < nsrc; k++) {
      var si = LS_srcIdx[k];
      var sti = Math.round((si % n) + disp * sx), stj = Math.round(((si / n) | 0) + disp * sy);
      if (sti < 0 || stj < 0 || sti >= n || stj >= n) continue; // ran off domain (export)
      var tc = stj * n + sti;
      lsFloorOf(h, tc);
      h[tc] += LS_srcT[k];
      if (LS_inGen[tc] !== LS_failId) { LS_inGen[tc] = LS_failId; LS_active.push(tc); }
    }
    talusGrow(h, n, dx, repose, LS_active, 150);
    return removed;
  }

  // Slow rotational creep: an incremental rotation of the block about a horizontal
  // axis through its centroid (perpendicular to the slide direction) — the head
  // drops and the toe rises by k·(u − ū), k = drop/R, conserving volume (raise is
  // scaled to match the clamped lowering). Only a thin slab moves per step, so a
  // slow slump creeps gradually instead of mobilizing the whole block like a
  // catastrophic failure. The toe bulge is then resolved with the growing AOR.
  function applySlowRotation(h, n, dx, Cx, Cy, Cz, R, sx, sy, drop, repose) {
    var cellA = dx * dx, R2 = R * R;
    var ic = (Cx / dx) | 0, jc = (Cy / dx) | 0, rc = ((R / dx) | 0) + 1;
    var i0 = ic - rc < 0 ? 0 : ic - rc, i1 = ic + rc > n - 1 ? n - 1 : ic + rc;
    var j0 = jc - rc < 0 ? 0 : jc - rc, j1 = jc + rc > n - 1 ? n - 1 : jc + rc;
    var nsrc = 0, su = 0;
    for (var j = j0; j <= j1; j++) {
      var ddy = j * dx - Cy;
      for (var i = i0; i <= i1; i++) {
        var ddx = i * dx - Cx, d2 = ddx * ddx + ddy * ddy;
        if (d2 >= R2) continue;
        var idx = j * n + i;
        if (Cz - Math.sqrt(R2 - d2) >= h[idx]) continue; // not in the block
        var u = ddx * sx + ddy * sy;
        LS_srcIdx[nsrc] = idx; LS_srcT[nsrc] = u; nsrc++; su += u;
      }
    }
    if (nsrc < 5) return 0;
    var ubar = su / nsrc, k = drop / R;
    var loweredTot = 0, raiseTot = 0;
    for (var a = 0; a < nsrc; a++) {
      var ia = LS_srcIdx[a], dz = k * (LS_srcT[a] - ubar);
      if (dz < 0) {
        var ddxa = (ia % n) * dx - Cx, ddya = ((ia / n) | 0) * dx - Cy;
        var zslip = Cz - Math.sqrt(R2 - (ddxa * ddxa + ddya * ddya));
        var room = h[ia] - zslip;
        loweredTot += -dz < room ? -dz : room;
      } else raiseTot += dz;
    }
    if (raiseTot <= 0 || loweredTot <= 0) return 0;
    var scale = loweredTot / raiseTot;
    LS_failId++;
    LS_active.length = 0;
    for (var b = 0; b < nsrc; b++) {
      var ib = LS_srcIdx[b], dzb = k * (LS_srcT[b] - ubar);
      if (dzb < 0) {
        var ddxb = (ib % n) * dx - Cx, ddyb = ((ib / n) | 0) * dx - Cy;
        var zsb = Cz - Math.sqrt(R2 - (ddxb * ddxb + ddyb * ddyb));
        var lo = -dzb < h[ib] - zsb ? -dzb : h[ib] - zsb;
        h[ib] -= lo; // head drops (leaves a headscarp)
      } else {
        lsFloorOf(h, ib);
        h[ib] += dzb * scale; // toe rises
        if (LS_inGen[ib] !== LS_failId) { LS_inGen[ib] = LS_failId; LS_active.push(ib); }
      }
    }
    talusGrow(h, n, dx, repose, LS_active, 80);
    return loweredTot * cellA;
  }

  // Stochastic failure testing + persistent slow-slump processing. `slumps` is a
  // mutable list of persistent rotational slumps that ride uplift, re-check their
  // stability each step, creep when unstable, and are retired once eroded away.
  function landslideStep(h, n, dx, p, rng, slumps) {
    var size = n * n;
    if (!LS_floor || LS_floor.length !== size) {
      LS_floor = new Float32Array(size); LS_fgen = new Int32Array(size);
      LS_inGen = new Int32Array(size); LS_srcIdx = new Int32Array(size); LS_srcT = new Float32Array(size);
    }
    var count = 0, vol = 0, inits = 0; // inits = NEW failures (excludes slump reactivation)
    var du = p.upliftDu || 0, ramp = p.upliftRamp || 0;

    // --- process persistent slow slumps ---
    for (var s = slumps.length - 1; s >= 0; s--) {
      var sl = slumps[s];
      if (du !== 0) { // ride uplift (the slip surface rises with the surface)
        sl.Cz += du * upliftFactor((sl.Cx / dx) | 0, (sl.Cy / dx) | 0, n, p.upliftPattern, ramp);
      }
      var fr = scoopFoS(h, n, dx, sl.Cx, sl.Cy, sl.Cz, sl.R, p.c, p.phi, p.gamma, p.ru, p.ice, p.gammaIce);
      if (!fr.valid || fr.blockVol < sl.vol0 * 0.05) { slumps.splice(s, 1); continue; } // eroded away
      if (fr.fos < p.fos) { // still unstable -> one small rotational creep increment
        vol += applySlowRotation(h, n, dx, sl.Cx, sl.Cy, sl.Cz, sl.R, fr.sx, fr.sy, SLOW_DROP, p.reposeSlow);
        count++;
      }
    }

    // --- stochastic failure tests (the "pre-existing weakness") ---
    for (var t = 0; t < p.tests; t++) {
      var R = p.Rmin + rng() * (p.Rmax - p.Rmin);
      // Keep the whole scoop footprint inside the grid: a clipped footprint
      // (e.g. near a corner) has an asymmetric column set whose vector-driving
      // doesn't cancel, spuriously failing. Restrict centres so the bbox fits.
      var rcT = ((R / dx) | 0) + 1;
      var lo = rcT, hi = n - 1 - rcT;
      if (hi < lo) continue; // scoop too large to fit anywhere
      var ic = lo + ((rng() * (hi - lo + 1)) | 0), jc = lo + ((rng() * (hi - lo + 1)) | 0);
      var D = R * (0.1 + rng() * 0.4);
      var Cx = ic * dx, Cy = jc * dx, Cz = h[jc * n + ic] + R - D;
      var fr2 = scoopFoS(h, n, dx, Cx, Cy, Cz, R, p.c, p.phi, p.gamma, p.ru, p.ice, p.gammaIce);
      if (!fr2.valid || fr2.fos >= p.fos) continue;
      if (p.mode === "slow") {
        if (slumps.length < SLUMP_MAX) {
          slumps.push({ Cx: Cx, Cy: Cy, Cz: Cz, R: R, sx: fr2.sx, sy: fr2.sy, vol0: fr2.blockVol });
          vol += applySlowRotation(h, n, dx, Cx, Cy, Cz, R, fr2.sx, fr2.sy, SLOW_DROP, p.reposeSlow);
          count++; inits++;
        }
      } else {
        var rf = Math.sqrt(Math.max(0, 2 * R * D - D * D));
        vol += applyScoop(h, n, dx, Cx, Cy, Cz, R, fr2.sx, fr2.sy, rf / dx, p.repose, p.mode === "vanish");
        count++; inits++;
      }
    }
    return { count: count, vol: vol, inits: inits };
  }

  // ---- Glaciers (P5a: ice layer + ELA mass balance + avalanche) ----------
  // Snow/ice avalanching = debris-limited talus on the ice surface (floor = bed):
  // ice cascades downslope until no ice-surface slope exceeds the repose angle.
  // Bedrock steeper than repose cannot hold ice (surface can't be relaxed with
  // ice on it) → it goes bare, exposing high peaks; ice pools on gentler ground.
  function iceAvalanche(h, ice, n, dx, reposeDeg, maxsweep) {
    var tanr = Math.tan((reposeDeg * Math.PI) / 180);
    for (var sweep = 0; sweep < maxsweep; sweep++) {
      var active = [];
      for (var k = 0; k < n * n; k++) if (ice[k] > 1e-3) active.push(k);
      active.sort(function (a, b) { return (h[b] + ice[b]) - (h[a] + ice[a]); });
      var moved = 0;
      for (var q = 0; q < active.length; q++) {
        var c = active[q], avail = ice[c];
        if (avail <= 0) continue;
        var cx = c % n, cy = (c / n) | 0, sc = h[c] + ice[c], best = -1, bestEx = 0;
        for (var d = 0; d < 8; d++) {
          var ni = cx + NB_DX[d], nj = cy + NB_DY[d];
          if (ni < 0 || nj < 0 || ni >= n || nj >= n) continue;
          var nb = nj * n + ni;
          var ex = sc - (h[nb] + ice[nb]) - tanr * NB_DIAG[d] * dx;
          if (ex > bestEx) { bestEx = ex; best = nb; }
        }
        if (best >= 0) {
          var amt = bestEx; if (amt > avail) amt = avail;
          ice[c] -= amt; ice[best] += amt; moved += amt;
        }
      }
      if (moved < 1e-2) break;
    }
  }

  // SIA ice flow (P5b/c). Shallow-Ice Approximation, depth-integrated flux form:
  //   q = −D ∇s,  D = Γ·H^(n+2)·|∇s|^(n−1),  Γ = 2A/(n+2)·(ρg)^n,  n = 3
  //   ∂H/∂t = ḃ − ∇·q     (ḃ = ELA mass balance)
  // Solved EXPLICITLY in conservative finite-volume form (exactly mass-conserving;
  // D = 0 on ice-free faces keeps H ≥ 0).
  //
  // EQUILIBRIUM MODE: on a landscape that evolves over Myr, ice (which responds in
  // 10²–10³ yr) is effectively always at steady state with the current bed. So each
  // frame we iterate the flux+balance to EQUILIBRIUM rather than advancing a fixed
  // time: sub-cycle on the CFL limit dt_sub ≤ 0.2·dx²/D_max (D capped at a fixed
  // 0.2·dx²/dt_sub_min, dt_sub itself capped at DTSUBMAX), and stop when the total
  // ice volume stops changing (windowed test — robust to the terminus limit-cycle
  // that a pointwise |∂H/∂t| test trips on). First spin-up costs ~10²–10³ sub-steps;
  // afterwards the bed barely moves frame-to-frame so the ice is already near
  // equilibrium and re-converges in ~10²; a per-frame BUDGET caps the worst case and
  // spreads big changes (e.g. an ELA move) over a few frames. A safety ICE_CAP
  // bounds runaway when accumulation ≫ ablation capacity (no bounded equilibrium).
  //
  // EROSION (over the full geomorphic dt): ΔE = K_g·τ_b·dt, abrasion ∝ basal shear
  // τ_b = ρg·H·min(|∇s_surface|, slope_cap) (Hallet 1979; MacGregor 2000; Egholm 2009).
  // Built from the SURFACE slope → over-deepening flattens the surface → self-limits at
  // ~ice thickness (Hergarten 2021); Coulomb-regularized (Schoof 2005) + a cm/yr per-step
  // ceiling → bounded, no runaway holes. Till is handed to fluvial as outwash. (A
  // semi-implicit ADI solver was prototyped but the bed-driven flux is advective and
  // only stable to dt≈2 yr — no gain; see methods.html.)
  // ρg = ρ_ice·g ≈ 917·9.8; A = temperate-ice Glen coefficient (Pa⁻³ yr⁻¹).
  // COARSE GRID: the ice physics runs on a grid ICE_COARSEN× coarser than the
  // bedrock. Ice is smooth, so sub-grid bed detail is irrelevant; and because the
  // SIA D cap scales as dx², a coarser grid gives far faster flow (real velocities,
  // not throttled) and ~ICE_COARSEN² fewer cells, so equilibration is cheap. The
  // bed is block-averaged down each step; ice thickness, speed and erosion are
  // bilinearly interpolated back up to the fine grid for rendering and to carve the
  // fine bedrock. Avalanching is omitted for now (steep ice mostly flows off via the
  // SIA term anyway). ρg = ρ_ice·g ≈ 917·9.8; A = temperate-ice Glen coeff (Pa⁻³ yr⁻¹).
  var ICE_RHOG = 8987.0, ICE_NEXP = 3, ICE_A = 7.57e-17;
  var ICE_GAMMA = (2 * ICE_A / (ICE_NEXP + 2)) * Math.pow(ICE_RHOG, ICE_NEXP);
  var ICE_COARSEN = 4;     // bed cells per ice cell (400→100 ice grid)
  var ICE_DTSUB_MIN = 1.0; // yr — sets the fixed D cap (D_cap = 0.2·dxC²/dt_sub_min)
  var ICE_DTSUBMAX = 5.0;  // yr — max sub-step (bounds thin-ice/source overshoot)
  var ICE_BUDGET = 500;    // max equilibration sub-steps per frame
  var ICE_CHECK = 20;      // volume-convergence check interval (sub-steps)
  var ICE_VTOL = 4e-4;     // relative volume-change tolerance for equilibrium
  var ICE_CAP = 3000;      // m — safety thickness cap (rarely hit now there's an outlet)
  var ICE_EDGE_DROP = 2.5; // boundary outflow: bed drops this×dxC outside the domain
  var ICE_FLOWMIN = 1e-3;  // m — ice below this can't source flux (upwind guard)
  var ICE_ABL_RATIO = 2.5; // ablation gradient ÷ accumulation gradient (real glaciers ~2–3)
  var ICE_SLOPE_CAP = 0.25; // Coulomb-regularized basal drag: surface slope in τ_b saturates here
  var ICE_ERO_RATECAP = 0.02; // m/yr — max bedrock-lowering rate (safety; real glacial erosion ≲ cm/yr)
  var _gbedC = null, _gqxC = null, _gqyC = null, _gsC = null, _gvelC = null, _geroC = null, _gNC = 0;
  function glacierScratch(sizeC) {
    if (_gNC !== sizeC) {
      _gbedC = new Float32Array(sizeC); _gqxC = new Float32Array(sizeC);
      _gqyC = new Float32Array(sizeC); _gsC = new Float32Array(sizeC);
      _gvelC = new Float32Array(sizeC); _geroC = new Float32Array(sizeC); _gNC = sizeC;
    }
  }
  // Block-average a fine field down to the coarse grid.
  function downsampleBlock(fine, n, coarse, nC) {
    var rx = n / nC;
    for (var J = 0; J < nC; J++) {
      var j0 = (J * rx) | 0, j1 = ((J + 1) * rx) | 0; if (j1 <= j0) j1 = j0 + 1; if (j1 > n) j1 = n;
      for (var I = 0; I < nC; I++) {
        var i0 = (I * rx) | 0, i1 = ((I + 1) * rx) | 0; if (i1 <= i0) i1 = i0 + 1; if (i1 > n) i1 = n;
        var sum = 0, cnt = 0;
        for (var jj = j0; jj < j1; jj++) for (var ii = i0; ii < i1; ii++) { sum += fine[jj * n + ii]; cnt++; }
        coarse[J * nC + I] = sum / cnt;
      }
    }
  }
  // Bilinearly interpolate a coarse field up to the fine grid (into `fine`).
  function upsampleBilinear(coarse, nC, fine, n) {
    var rx = n / nC;
    for (var j = 0; j < n; j++) {
      var fcy = (j + 0.5) / rx - 0.5, J0 = Math.floor(fcy), ty = fcy - J0;
      var J0c = J0 < 0 ? 0 : (J0 > nC - 1 ? nC - 1 : J0);
      var J1c = J0 + 1 < 0 ? 0 : (J0 + 1 > nC - 1 ? nC - 1 : J0 + 1);
      for (var i = 0; i < n; i++) {
        var fcx = (i + 0.5) / rx - 0.5, I0 = Math.floor(fcx), tx = fcx - I0;
        var I0c = I0 < 0 ? 0 : (I0 > nC - 1 ? nC - 1 : I0);
        var I1c = I0 + 1 < 0 ? 0 : (I0 + 1 > nC - 1 ? nC - 1 : I0 + 1);
        var c00 = coarse[J0c * nC + I0c], c10 = coarse[J0c * nC + I1c];
        var c01 = coarse[J1c * nC + I0c], c11 = coarse[J1c * nC + I1c];
        var top = c00 + (c10 - c00) * tx, bot = c01 + (c11 - c01) * tx;
        fine[j * n + i] = top + (bot - top) * ty;
      }
    }
  }
  // Equilibrate the ice (coarse) to the current bed, carve the fine bedrock over the
  // geomorphic dt p.dt, and upsample ice/speed to the fine grid. Returns {sub,eroVol}.
  function glacierStep(h, n, dx, iceC, nC, dxC, fineIce, fineVel, fineMelt, fineSed, p) {
    var sizeC = nC * nC, frame = p.dt;
    glacierScratch(sizeC);
    var bedC = _gbedC, qx = _gqxC, qy = _gqyC, s = _gsC, velC = _gvelC, eroC = _geroC;
    downsampleBlock(h, n, bedC, nC);
    var GAM = ICE_GAMMA * (p.flow != null ? p.flow : 1);
    var dx2 = dxC * dxC;
    var DCAP = (0.2 * dx2) / ICE_DTSUB_MIN; // fixed, dt-independent (no flow throttling)
    var ela = p.ela, bg = p.balGrad, accMax = p.accMax, abMax = p.abMax;
    // Outflow at domain edges: the bed is treated as dropping ICE_EDGE_DROP·dxC just
    // outside, so ice spills over the rim and leaves the domain (an outlet → there is
    // a bounded equilibrium even when accumulation exceeds the domain's ablation).
    // A moderate drop drains without violently vacuuming the edge (the H⁵ flux is
    // self-limiting). edgeOut(H) = flux magnitude leaving across one boundary face.
    var BDROP = ICE_EDGE_DROP * dxC;
    function edgeOut(H) {
      if (H <= 0) return 0;
      var sl = (H + BDROP) / dxC, Hf = 0.5 * H;
      var D = GAM * Hf * Hf * Hf * Hf * Hf * sl * sl;
      if (D > DCAP) D = DCAP;
      return D * sl;
    }
    // Which edges drain. For a tilting plane, mirror the fluvial perimeter logic:
    // the bed tilts up toward +x, so ice flows DOWN-tilt and leaves only the LEFT
    // (low) edge; the right edge is up-tilt (uphill outside → no outflow) and the
    // top/bottom are cross-slope (no off-grid flow). Other patterns: drain all
    // edges (the perimeter sits below the interior, so ice spills off any side).
    var tilt = p.pattern === "tilt";
    var outL = true, outR = !tilt, outT = !tilt, outB = !tilt;
    var i, j, k, sub = 0, volPrev = 0;
    for (k = 0; k < sizeC; k++) volPrev += iceC[k];
    // ---- equilibration loop (coarse grid) ----
    while (sub < ICE_BUDGET) {
      for (k = 0; k < sizeC; k++) s[k] = bedC[k] + iceC[k];
      var Dmax = 0;
      for (j = 0; j < nC; j++) { // x-faces
        var jm = j > 0 ? j - 1 : 0, jp = j < nC - 1 ? j + 1 : nC - 1, jd = jp - jm;
        for (i = 0; i < nC - 1; i++) {
          k = j * nC + i;
          var sx = (s[k + 1] - s[k]) / dxC;
          // Upwind guard: ice flows from the higher-surface (source) cell. If that
          // cell is ice-free, there is no ice to transport — without this, a centered
          // Hf lets the H≥0 clamp manufacture ice from an ice-free, higher-bedrock
          // neighbour, leaving phantom patches that persist far below the ELA during
          // retreat. (No effect on normal flow, where the source is the icy uphill.)
          if (iceC[sx > 0 ? k + 1 : k] <= ICE_FLOWMIN) { qx[k] = 0; continue; }
          var Hf = 0.5 * (iceC[k] + iceC[k + 1]);
          var syx = jd ? ((s[jp * nC + i] + s[jp * nC + i + 1]) - (s[jm * nC + i] + s[jm * nC + i + 1])) / (2 * jd * dxC) : 0;
          var slx = Math.sqrt(sx * sx + syx * syx);
          var Dx = GAM * Hf * Hf * Hf * Hf * Hf * slx * slx;
          if (Dx > DCAP) Dx = DCAP;
          if (Dx > Dmax) Dmax = Dx;
          qx[k] = -Dx * sx;
        }
        qx[j * nC + nC - 1] = 0;
      }
      for (j = 0; j < nC - 1; j++) { // y-faces
        for (i = 0; i < nC; i++) {
          k = j * nC + i;
          var sy = (s[k + nC] - s[k]) / dxC;
          if (iceC[sy > 0 ? k + nC : k] <= ICE_FLOWMIN) { qy[k] = 0; continue; } // upwind guard (see x-faces)
          var Hfy = 0.5 * (iceC[k] + iceC[k + nC]);
          var im = i > 0 ? i - 1 : 0, ip = i < nC - 1 ? i + 1 : nC - 1, id = ip - im;
          var sxy = id ? ((s[j * nC + ip] + s[(j + 1) * nC + ip]) - (s[j * nC + im] + s[(j + 1) * nC + im])) / (2 * id * dxC) : 0;
          var sly = Math.sqrt(sy * sy + sxy * sxy);
          var Dy = GAM * Hfy * Hfy * Hfy * Hfy * Hfy * sly * sly;
          if (Dy > DCAP) Dy = DCAP;
          if (Dy > Dmax) Dmax = Dy;
          qy[k] = -Dy * sy;
        }
      }
      for (i = 0; i < nC; i++) qy[(nC - 1) * nC + i] = 0;
      var dtsub = Dmax > 0 ? 0.2 * dx2 / Dmax : ICE_DTSUBMAX;
      if (dtsub > ICE_DTSUBMAX) dtsub = ICE_DTSUBMAX;
      for (j = 0; j < nC; j++) { // update H: ∂H/∂t = balance − div(q); clamp [0, ICE_CAP]
        for (i = 0; i < nC; i++) {
          k = j * nC + i;
          // Mass balance is referenced to the ice SURFACE (bed + H), the physical
          // firn line. This carries a real surface-elevation feedback (thicker ice →
          // higher, colder surface → more accumulation) and the resulting small-ice-
          // cap bistability — an intended feature of the modelled system, not a bug.
          // The ablation gradient (below the ELA) is steeper than the accumulation
          // gradient — a standard glaciological asymmetry. A symmetric gradient melts
          // too weakly below the ELA and tongues stretch unrealistically far.
          var ds = s[k] - ela;
          var b = ds > 0 ? ds * bg : ds * (bg * ICE_ABL_RATIO);
          if (b > accMax) b = accMax; else if (b < -abMax) b = -abMax;
          // boundary faces drain out of the domain (edgeOut), interior faces use q
          var qE = i < nC - 1 ? qx[k] : (outR ? edgeOut(iceC[k]) : 0);
          var qW = i > 0 ? qx[k - 1] : (outL ? -edgeOut(iceC[k]) : 0);
          var qN = j < nC - 1 ? qy[k] : (outB ? edgeOut(iceC[k]) : 0);
          var qS = j > 0 ? qy[k - nC] : (outT ? -edgeOut(iceC[k]) : 0);
          var v = iceC[k] + dtsub * (b - ((qE - qW) + (qN - qS)) / dxC);
          iceC[k] = v < 0 ? 0 : (v > ICE_CAP ? ICE_CAP : v);
        }
      }
      sub++;
      if (sub % ICE_CHECK === 0) {
        var vol = 0;
        for (k = 0; k < sizeC; k++) vol += iceC[k];
        if (Math.abs(vol - volPrev) <= ICE_VTOL * (vol + 1)) break;
        volPrev = vol;
      }
    }
    // ---- display speed (|q|/H) + glacial abrasion over dt ----
    // Published-model recipe (MacGregor 2000; Egholm 2009; Herman; Hergarten/Deal &
    // Prasicek 2021). Abrasion ∝ basal shear stress (a proxy for basal sliding speed;
    // Hallet 1979; l=1 per Humphrey & Raymond 1994 / Cook 2020):
    //     E = K_g · τ_b,  τ_b = ρg·H·min(|∇s|, slope_cap)
    // τ_b is built from the ICE-SURFACE slope, so the over-deepening it creates flattens
    // the surface → driving stress falls → erosion shuts off; the overdeepening depth
    // self-limits to ~the ice thickness rather than running away (Hergarten 2021). The
    // surface slope is Coulomb-regularized (saturates at slope_cap; basal drag cannot
    // exceed a Coulomb yield — Schoof 2005), so a steep wall can't explode τ_b. A
    // realistic per-frame ceiling (≲ cm/yr) keeps the bed from out-running the ice's
    // ability to re-flatten its surface — the combination is bounded and robust even
    // when K_g is cranked hard (validated: no runaway holes; carves U-troughs).
    var Kg = p.iceErode != null ? p.iceErode : 0;
    var eroPerCap = ICE_ERO_RATECAP * frame; // m of bedrock per frame ceiling
    for (k = 0; k < sizeC; k++) s[k] = bedC[k] + iceC[k]; // fresh equilibrium surface
    for (j = 0; j < nC; j++) {
      for (i = 0; i < nC; i++) {
        k = j * nC + i;
        var H = iceC[k];
        if (H <= 1) { velC[k] = 0; eroC[k] = 0; continue; }
        var ux = 0.5 * ((i < nC - 1 ? qx[k] : 0) + (i > 0 ? qx[k - 1] : 0));
        var uy = 0.5 * ((j < nC - 1 ? qy[k] : 0) + (j > 0 ? qy[k - nC] : 0));
        velC[k] = Math.sqrt(ux * ux + uy * uy) / H; // actual flow speed, for display
        if (Kg > 0 && frame > 0) {
          var ip = i < nC - 1 ? i + 1 : i, im = i > 0 ? i - 1 : i;
          var jp = j < nC - 1 ? j + 1 : j, jm = j > 0 ? j - 1 : j;
          var sx2 = ip !== im ? (s[j * nC + ip] - s[j * nC + im]) / ((ip - im) * dxC) : 0;
          var sy2 = jp !== jm ? (s[jp * nC + i] - s[jm * nC + i]) / ((jp - jm) * dxC) : 0;
          var sl = Math.sqrt(sx2 * sx2 + sy2 * sy2);
          if (sl > ICE_SLOPE_CAP) sl = ICE_SLOPE_CAP; // Coulomb regularization
          var e = Kg * (ICE_RHOG * H * sl) * frame;
          eroC[k] = e > eroPerCap ? eroPerCap : e;
        } else eroC[k] = 0;
      }
    }
    // ---- upsample to fine grid: ice + speed for display, erosion to carve bedrock ----
    upsampleBilinear(iceC, nC, fineIce, n);
    upsampleBilinear(velC, nC, fineVel, n);
    var eroVol = 0;
    if (Kg > 0) {
      var rx = n / nC; // inline bilinear upsample of eroC, subtracting from bedrock h
      for (j = 0; j < n; j++) {
        var fcy = (j + 0.5) / rx - 0.5, J0 = Math.floor(fcy), ty = fcy - J0;
        var J0c = J0 < 0 ? 0 : (J0 > nC - 1 ? nC - 1 : J0), J1c = J0 + 1 < 0 ? 0 : (J0 + 1 > nC - 1 ? nC - 1 : J0 + 1);
        for (i = 0; i < n; i++) {
          var fcx = (i + 0.5) / rx - 0.5, I0 = Math.floor(fcx), tx = fcx - I0;
          var I0c = I0 < 0 ? 0 : (I0 > nC - 1 ? nC - 1 : I0), I1c = I0 + 1 < 0 ? 0 : (I0 + 1 > nC - 1 ? nC - 1 : I0 + 1);
          var c00 = eroC[J0c * nC + I0c], c10 = eroC[J0c * nC + I1c], c01 = eroC[J1c * nC + I0c], c11 = eroC[J1c * nC + I1c];
          var top = c00 + (c10 - c00) * tx, bot = c01 + (c11 - c01) * tx, e = top + (bot - top) * ty;
          if (e > 0) { h[j * n + i] -= e; eroVol += e; }
        }
      }
    }
    // ---- meltwater + glacial-sediment release (for fluvial outwash coupling) ----
    // Where ice melts (ablation), release meltwater (discharge) and the glacially
    // eroded sediment, both proportional to local melt — i.e. the glacier delivers
    // its water and ground-up load to where it wastes. fineMelt is a melt rate
    // (m/yr, added to fluvial discharge); fineSed is a sediment rate (m³/yr,
    // injected into the fluvial Qs). Σ fineSed·dt = eroVol → mass conserved across
    // the glacial-erosion → fluvial-deposition hand-off. Consumed next frame.
    if (fineMelt && fineSed) {
      var ela2 = p.ela, bg2 = p.balGrad, abMax2 = p.abMax, totMelt = 0, fk;
      var ablG = bg2 * ICE_ABL_RATIO; // ablation gradient (matches the balance)
      for (fk = 0; fk < n * n; fk++) {
        var surf = h[fk] + fineIce[fk], mlt = 0; // melt where the ice SURFACE is below the ELA
        if (fineIce[fk] > 1 && surf < ela2) {
          mlt = (ela2 - surf) * ablG;
          if (mlt > abMax2) mlt = abMax2;
        }
        fineMelt[fk] = mlt; totMelt += mlt;
      }
      var sedScale = (totMelt > 0 && frame > 0) ? eroVol / (totMelt * frame) : 0; // m²
      for (fk = 0; fk < n * n; fk++) fineSed[fk] = fineMelt[fk] * sedScale; // m³/yr
    }
    return { sub: sub, eroVol: eroVol };
  }

  // ---- Model object ------------------------------------------------------
  function createModel(opts) {
    opts = opts || {};
    var n = opts.n || 400;
    var dx = opts.dx || 50; // metres per cell (400 * 50 = 20 km domain)
    var size = n * n;
    var h = new Float32Array(size);
    var tmp = new Float32Array(size);
    var initial = new Float32Array(size); // snapshot at last generate()

    var area = new Float32Array(size); // drainage area (display + stream power)
    var dz = new Float32Array(size); // elevation change this step (erosion/deposition/creep)
    var sedFlux = new Float32Array(size); // sediment throughput this step (m^3)
    var lake = new Float32Array(size); // ponded-water depth (lakes), from fluvial routing
    var dzBuf = new Float32Array(size); // scratch: elevation snapshot
    var slideRng = mulberry32(12345); // persistent PRNG for stochastic failures
    var slumps = []; // persistent slow rotational slumps
    // Glaciers run on a COARSER grid (ice is smooth; sub-grid bed detail doesn't
    // matter, and a bigger dx makes the stiff SIA flow far cheaper to equilibrate).
    var nC = Math.max(8, Math.round(n / ICE_COARSEN));
    var dxC = (n * dx) / nC; // keep the same physical domain
    var iceC = new Float32Array(nC * nC); // authoritative ice thickness (m), coarse
    var ice = new Float32Array(size); // ice thickness upsampled to fine grid (display/erosion)
    var iceVel = new Float32Array(size); // ice speed (m/yr) upsampled to fine grid
    var iceMelt = new Float32Array(size); // meltwater release rate (m/yr) → fluvial discharge
    var glacSed = new Float32Array(size); // glacial sediment release (m³/yr) → fluvial Qs
    var state = { n: n, dx: dx, h: h, initial: initial, area: area, dz: dz, sedFlux: sedFlux, lake: lake, ice: ice, iceVel: iceVel, iceC: iceC, nC: nC, dxC: dxC, time: 0, steps: 0, fluvIters: 0, slideCount: 0, slideVol: 0, slideInits: 0, slumpCount: 0, iceEroVol: 0, iceSubs: 0, iceMelt: iceMelt, glacSed: glacSed };

    function generate(g) {
      g = g || {};
      var type = g.type || "hills";
      (generators[type] || generators.hills)(h, n, dx, g);
      area.fill(0);
      dz.fill(0);
      sedFlux.fill(0);
      lake.fill(0);
      ice.fill(0);
      iceVel.fill(0);
      iceC.fill(0);
      iceMelt.fill(0);
      glacSed.fill(0);
      slumps.length = 0;
      slideRng = mulberry32((((g.seed | 0) ^ 0x5f3759df) >>> 0) || 1);
      // Uncorrelated white-noise perturbation. A little roughness breaks the
      // symmetry of smooth surfaces (planes, cones) so dendritic drainage can
      // self-organize instead of forming parallel/degenerate channels.
      var noise = g.noise != null ? g.noise : 8; // metres (peak amplitude)
      if (noise > 0) {
        var rand = mulberry32(((g.seed | 0) ^ 0x9e3779b9) >>> 0);
        for (var k = 0; k < size; k++) h[k] += (rand() * 2 - 1) * noise;
      }
      initial.set(h);
      state.time = 0;
      state.steps = 0;
    }

    // Advance one frame of model time. dt in years. Auto-substeps creep so the
    // explicit diffusion stays stable regardless of slider values.
    function step(p) {
      p = p || {};
      var dt = p.dt != null ? p.dt : 1000; // yr per frame
      var D = p.D != null ? p.D : 1; // m^2/yr
      var uplift = p.uplift != null ? p.uplift : 0; // m/yr (uniform)
      var nonlinear = !!p.nonlinear;
      var Sc = p.Sc != null ? p.Sc : 0.6;
      var fluvial = !!p.fluvial;
      var K = p.K != null ? p.K : 2e-5; // stream-power erodibility
      var m = p.m != null ? p.m : 0.5; // drainage-area exponent
      var nexp = p.nexp != null ? p.nexp : 1; // slope exponent
      var G = p.depo != null ? p.depo : 1; // deposition coefficient
      var mfdExp = p.mfd != null ? p.mfd : 4; // MFD flow-concentration exponent

      // Spatial uplift: zero at the rim, ramping linearly up to full rate over
      // `upliftRamp` cells from the edge, then flat in the interior. This keeps
      // the rim a low, stationary base level while the centre rises, so relief
      // (and the interior-to-edge gradient that drives fluvial incision) grows
      // with time. upliftRamp = 0 recovers spatially uniform uplift.
      var du = uplift * dt; // per-step uplift increment (also used by slow slumps)
      var ramp = p.upliftRamp != null ? p.upliftRamp : 50;
      var pattern = p.upliftPattern || "ramp";
      if (uplift !== 0) {
        for (var jj = 0; jj < n; jj++) {
          for (var ii = 0; ii < n; ii++) {
            h[jj * n + ii] += du * upliftFactor(ii, jj, n, pattern, ramp);
          }
        }
      }

      // Snapshot (post-uplift) so dz captures only the geomorphic change.
      dzBuf.set(h);

      // Fluvial incision + deposition: one implicit (Yuan 2019) solve per frame,
      // unconditionally stable at any dt — no sub-stepping needed.
      if (fluvial) {
        // When glaciers are on, suppress rain above the firn line (precip there is
        // snow → ice, no runoff) and skip incision under ice. Uses last frame's ice.
        var fluvEla = p.glaciers ? (p.ela != null ? p.ela : 1000) : null;
        fluvialStep(h, n, dx, K, m, nexp, dt, G, mfdExp, pattern, p.glaciers ? ice : null, fluvEla,
          p.glaciers ? iceMelt : null, p.glaciers ? glacSed : null);
        if (fluvialStep._fs) {
          area.set(fluvialStep._fs.area);
          sedFlux.set(fluvialStep._fs.sedFlux);
          lake.set(fluvialStep._fs.lake);
          state.fluvIters = fluvialStep._fs.lastIters;
        }
      } else {
        sedFlux.fill(0);
        lake.fill(0);
        state.fluvIters = 0;
      }

      // Hillslope creep. Auto-substep so explicit diffusion stays stable.
      var safety = nonlinear ? 0.1 : 0.2;
      var sub = Math.max(1, Math.ceil((D * dt) / (dx * dx) / safety));
      var sdt = dt / sub;
      for (var s = 0; s < sub; s++) {
        if (nonlinear) creepNonlinearStep(h, tmp, n, dx, D, Sc, sdt);
        else creepLinearStep(h, tmp, n, dx, D, sdt);
        h.set(tmp);
      }

      // Deep-seated landslides (stochastic failure testing).
      if (p.landslides) {
        var lr = landslideStep(h, n, dx, {
          tests: p.lsTests != null ? p.lsTests : 50,
          fos: p.lsFos != null ? p.lsFos : 1,
          c: p.lsC != null ? p.lsC : 10,
          phi: p.lsPhi != null ? p.lsPhi : 32,
          gamma: 20,
          ru: p.lsRu != null ? p.lsRu : 0.2,
          Rmin: 100,
          Rmax: p.lsRmax != null ? p.lsRmax : 1200,
          repose: p.lsRepose != null ? p.lsRepose : 14,
          reposeSlow: p.lsReposeSlow != null ? p.lsReposeSlow : 35,
          mode: p.lsMode || "catastrophic",
          // ice surcharge for buttressing (γ_ice ≈ 9 kN/m³); only when glaciers on
          ice: p.glaciers ? ice : null,
          gammaIce: ICE_RHOG / 1000,
          upliftDu: du,
          upliftRamp: ramp,
          upliftPattern: pattern
        }, slideRng, slumps);
        state.slideCount = lr.count;
        state.slideVol = lr.vol;
        state.slideInits = lr.inits;
        state.slumpCount = slumps.length;
      } else {
        state.slideCount = 0;
        state.slideVol = 0;
        state.slideInits = 0;
      }

      // Glaciers: equilibrate SIA ice (coarse grid) then carve the fine bed (P5b/c).
      if (p.glaciers) {
        var gl = glacierStep(h, n, dx, iceC, nC, dxC, ice, iceVel, iceMelt, glacSed, {
          dt: dt,
          ela: p.ela != null ? p.ela : 1000,
          balGrad: p.balGrad != null ? p.balGrad : 0.007,
          accMax: p.accMax != null ? p.accMax : 1.5,
          abMax: p.abMax != null ? p.abMax : 10,
          flow: p.iceFlow != null ? p.iceFlow : 1,
          iceErode: p.iceErode != null ? p.iceErode : 0,
          pattern: pattern
        });
        state.iceEroVol = gl.eroVol;
        state.iceSubs = gl.sub;
      } else {
        ice.fill(0); iceVel.fill(0); iceC.fill(0); iceMelt.fill(0); glacSed.fill(0); // glaciers off
        state.iceEroVol = 0; state.iceSubs = 0;
      }

      // dz = geomorphic elevation change this step (excludes uplift).
      for (var q2 = 0; q2 < size; q2++) dz[q2] = h[q2] - dzBuf[q2];

      state.time += dt;
      state.steps += 1;
      return state;
    }

    function reset() {
      h.set(initial);
      area.fill(0);
      dz.fill(0);
      sedFlux.fill(0);
      lake.fill(0);
      ice.fill(0);
      iceVel.fill(0);
      iceC.fill(0);
      iceMelt.fill(0);
      glacSed.fill(0);
      slumps.length = 0;
      slideRng = mulberry32(12345);
      state.time = 0;
      state.steps = 0;
    }

    function stats() {
      var min = Infinity, max = -Infinity, sum = 0;
      for (var k = 0; k < size; k++) {
        var v = h[k];
        if (v < min) min = v;
        if (v > max) max = v;
        sum += v;
      }
      return { min: min, max: max, mean: sum / size, relief: max - min };
    }

    return {
      state: state,
      generate: generate,
      step: step,
      reset: reset,
      stats: stats
    };
  }

  var api = {
    createModel: createModel,
    generators: generators,
    creepLinearStep: creepLinearStep,
    creepNonlinearStep: creepNonlinearStep,
    fluvialStep: fluvialStep,
    fbm: fbm
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.TerrainModel = api;
  }
})(typeof window !== "undefined" ? window : this);
