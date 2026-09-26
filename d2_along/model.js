/*
 * d2_along/model.js — the along-valley tier, z(s, t). Pure model, no DOM.
 *
 * s along the valley: node i at s = i·ds, i = 0 (divide) … N−1 (outlet, fixed base
 * level z = 0). Units m, yr. State in Float32Array; arithmetic in doubles.
 *
 * Asserted: drainage area A(s) = ka (s + s0)^h (Hack), valley width W(s) = W0 + kw s.
 * Computed: ice flux Q (∫ b W ds), thickness H (upstream march on the SIA flux law),
 * sliding U_s = f_s q/H, erosion, fluvial incision (Yuan implicit along the line).
 *
 * Mirrors py/d2_along.py loop for loop (validated there: py/test_d2_along.py — the
 * DESIGN.md §4 gate against a reference explicit SIA passes within 2 %, fluvial
 * steady state is the analytic slope–area law, mass closes to round-off, K_g × 100
 * self-limits). The only step that differs is the reference SIA, which is Python-only.
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  var forcing = GS.forcing, laws = GS.laws, litho = GS.lithology, ICE = laws.ICE;

  var DEFAULTS = {
    L: 30000, N: 301,
    initProfile: "steady", zHead: 2500, noise: 0, seed: 1,
    shape: "constant", peakUplift: 1e-3, duration: 1e7, pattern: "uniform", ramp: 0.1,
    K: 3.2e-6, m: 0.5, nexp: 1.0, G: 0, ka: 5.7, hack: 1.667, s0: 200, W0: 300, kw: 0.05,
    glacierOn: true, elaBase: 1500, elaAmp: 0, elaPeriod: 1e5, balGrad: 0.007,
    accMax: 2.0, abMax: 8.0, fs: 0.5, flow: 1.0, Kg: 1e-4, lexp: 1.0, eroCap: 0.02,
    Kq: 0, iceIters: 30, relax: 0.5, smin: 1e-3, iceTol: 0.01, Hmin: 10, eroSmooth: true, Hf: 100, phiSub: 0, Sc: 0.8,
    litho: null, contrastK: 5, contrastKg: 5,
    dt: 500
  };

  function hackArea(s, ka, h, s0) { return ka * Math.pow(s + s0, h); }
  function valleyWidth(s, W0, kw) { return W0 + kw * s; }

  // Q(s) = ∫ b W ds downstream, clamped at zero (terminus).
  function routeFlux(b, W, ds, Q, N) {
    var acc = 0;
    for (var i = 0; i < N; i++) {
      acc += b[i] * W[i] * ds;
      if (acc < 0) acc = 0;
      Q[i] = acc;
    }
  }

  // Given Q_i leaving node i, march upstream from the outlet solving the face flux law for
  // H_i given H_{i+1} (laws.faceThickness: monotone → unique root).
  // Outlet: free outflow — ice leaves with the uniform-flow thickness for the local bed slope
  // (a zero-thickness outlet forces a one-cell wedge and a sliding spike there).
  function marchThickness(z, W, Q, ds, Gam, n, N, H, smin) {
    H[N - 1] = 0;
    if (Q[N - 1] > 0) {
      var Sb = (z[N - 2] - z[N - 1]) / ds; if (Sb < smin) Sb = smin;
      H[N - 1] = Math.pow((Q[N - 1] / W[N - 1]) / (Gam * Math.pow(Sb, n)), 1 / (n + 2));
    }
    for (var i = N - 2; i >= 0; i--) {
      if (Q[i] <= 0) { H[i] = 0; continue; }
      var Wf = 0.5 * (W[i] + W[i + 1]);
      H[i] = laws.faceThickness(z[i], z[i + 1] + H[i + 1], H[i + 1], Q[i] / Wf, ds, Gam, n);
    }
  }

  function createModel(params) {
    var p = Object.assign({}, DEFAULTS, params || {});
    var N = p.N, L = p.L, ds = L / (N - 1);
    var st = {
      p: p, t: 0, N: N, ds: ds,
      s: new Float32Array(N), z: new Float32Array(N), ucum: new Float32Array(N),
      H: new Float32Array(N), Us: new Float32Array(N), q: new Float32Array(N), b: new Float32Array(N),
      Ef: new Float32Array(N), Eg: new Float32Array(N), Er: new Float32Array(N), r: new Float32Array(N),
      A: new Float32Array(N), W: new Float32Array(N),
      exported: 0, eroFluv: 0, eroGlac: 0, eroRock: 0, upliftVol: 0, iceIters: 0,
      history: []
    };
    // scratch (doubles)
    var Hn = new Float64Array(N), Q = new Float64Array(N), bb = new Float64Array(N), zz = new Float64Array(N);
    var h0 = new Float64Array(N), Qs = new Float64Array(N), fK = new Float64Array(N), fKg = new Float64Array(N), fac = new Float64Array(N);
    var iceMask = new Uint8Array(N), wfl = new Float64Array(N);

    var rng = GS.noise.mulberry32(p.seed | 0);
    var i;
    for (i = 0; i < N; i++) {
      var s = i * ds; st.s[i] = s;
      var f = 1 - s / L;
      st.z[i] = p.initProfile === "linear" ? p.zHead * f : p.zHead * Math.pow(f, 1.5);
      if (p.noise > 0) { // Box–Muller
        var u1 = rng() || 1e-12, u2 = rng();
        st.z[i] += p.noise * Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
      }
      st.A[i] = hackArea(s, p.ka, p.hack, p.s0);
      st.W[i] = valleyWidth(s, p.W0, p.kw);
      fac[i] = forcing.spatialFactor(s / L, 0.5, p.pattern, p.ramp);
    }
    st.z[N - 1] = 0;
    st.steadyCapped = 0; // nodes where the steady profile hit the threshold-slope cap
    if (p.initProfile === "steady" && p.K > 0) { // analytic fluvial steady state for U_peak·f(s): S = (U f / K A^m)^(1/n)
      var zz0 = 0;
      for (i = N - 2; i >= 0; i--) {
        var Ss = Math.pow(p.peakUplift * fac[i] / (p.K * Math.pow(st.A[i], p.m)), 1 / p.nexp);
        if (Ss > p.Sc) { Ss = p.Sc; st.steadyCapped++; } // threshold-hillslope cap
        zz0 += ds * Ss;
        st.z[i] = zz0 + (p.noise > 0 ? st.z[i] - (p.zHead * Math.pow(1 - st.s[i] / L, 1.5)) : 0);
      }
    }
    st.z0 = new Float32Array(st.z);

    var U = forcing.makeSeries({ shape: p.shape, peak: p.peakUplift, duration: p.duration });
    var ELA = forcing.makeSeries({ shape: "sine", peak: -p.elaAmp, base: p.elaBase, period: p.elaPeriod });
    var lithoFn = p.litho ? litho.make(p.litho) : null;
    var n = ICE.N, Gam = ICE.GAMMA * p.flow / (1 - p.fs);

    // Steady ice on the current bed: surface → balance → routed flux → march, under-relaxed.
    function steadyIce(ela) {
      var H = st.H, z = st.z, W = st.W, it, k;
      for (k = 0; k < N; k++) zz[k] = z[k];
      for (it = 1; it <= p.iceIters; it++) {
        for (k = 0; k < N; k++) bb[k] = laws.massBalance(zz[k] + H[k], ela, p.balGrad, p.accMax, p.abMax);
        routeFlux(bb, W, ds, Q, N);
        marchThickness(zz, W, Q, ds, Gam, n, N, Hn, p.smin);
        var maxd = 0;
        for (k = 0; k < N; k++) {
          var d = Hn[k] - H[k];
          H[k] = H[k] + p.relax * d;
          if (d < 0) d = -d;
          if (d > maxd) maxd = d;
        }
        if (maxd < p.iceTol) break;
      }
      st.iceIters = it;
      for (k = 0; k < N; k++) bb[k] = laws.massBalance(zz[k] + H[k], ela, p.balGrad, p.accMax, p.abMax);
      routeFlux(bb, W, ds, Q, N);
      for (k = 0; k < N; k++) {
        st.q[k] = Q[k] / W[k]; st.b[k] = bb[k];
        st.Us[k] = laws.slidingSpeed(st.q[k], H[k], p.fs, p.Hmin);
      }
    }

    // Fluvial efficiency under ice: 1 on bare ground, ramping down over the first Hf metres of
    // ice to a floor phiSub (subglacial meltwater). A binary on/off at the margin made each cell
    // the terminus vacated take a burst of incision and left a cell-scale sawtooth in bed and ice.
    function fluvialWeight(H) {
      var r = 1 - H / p.Hf; if (r < 0) r = 0; if (r > 1) r = 1;
      return p.phiSub + (1 - p.phiSub) * r;
    }
    // Implicit Yuan erosion–deposition along the line (receiver = downstream node).
    function fluvialStep() {
      var z = st.z, A = st.A, W = st.W, dt = p.dt, G = p.G, k;
      for (k = 0; k < N; k++) { h0[k] = z[k]; wfl[k] = fluvialWeight(st.H[k]); }
      var OMEGA = G > 0 ? 0.6 : 1.0; // relax only while iterating the deposition feedback
      var maxIter = G <= 0 ? 1 : Math.min(25, Math.ceil(2 + G * 20));
      for (var it = 0; it < maxIter; it++) {
        for (k = 0; k < N; k++) Qs[k] = 0;
        for (k = 0; k < N - 1; k++) Qs[k + 1] += Qs[k] + (h0[k] - z[k]) * ds * W[k] / dt;
        var maxd = 0;
        for (k = N - 2; k >= 0; k--) {
          if (wfl[k] <= 0) continue;
          var hr = z[k + 1];
          var dep = dt * G * Qs[k] / A[k];
          var Kp = p.K * fK[k] * wfl[k] * Math.pow(A[k], p.m);
          var nh = laws.yuanNode(h0[k], z[k], hr, Kp, dt, ds, dep, p.nexp, false);
          nh = z[k] + OMEGA * (nh - z[k]);
          var dd = nh - z[k]; if (dd < 0) dd = -dd; if (dd > maxd) maxd = dd;
          z[k] = nh;
        }
        if (G <= 0 || maxd < 0.02) break;
      }
      for (k = 0; k < N; k++) Qs[k] = 0;
      var vol = 0;
      for (k = 0; k < N - 1; k++) {
        Qs[k + 1] += Qs[k] + (h0[k] - z[k]) * ds * W[k] / dt;
        st.Ef[k] = (h0[k] - z[k]) / dt;
        vol += (h0[k] - z[k]) * ds * W[k];
      }
      st.Ef[N - 1] = 0;
      st.exported += Qs[N - 1] * dt;
      st.eroFluv += vol;
    }

    function step() {
      var z = st.z, dt = p.dt, k;
      // lithology at the surface (material frame)
      for (k = 0; k < N; k++) {
        var r = lithoFn ? lithoFn(st.s[k], 0, litho.materialZ(z[k], st.ucum[k])) : 0;
        st.r[k] = r;
        fK[k] = litho.erodibilityFactor(r, p.contrastK);
        fKg[k] = litho.erodibilityFactor(r, p.contrastKg);
      }
      // 1. uplift (outlet fixed)
      var u = U(st.t), uv = 0;
      for (k = 0; k < N - 1; k++) {
        var du = u * fac[k] * dt;
        z[k] += du; st.ucum[k] += du; uv += du * ds * st.W[k];
      }
      st.upliftVol += uv;
      // 2. glacier (steady) → mask + erosion
      if (p.glacierOn) {
        steadyIce(ELA(st.t));
        var H = st.H, gv = 0, Eraw = h0; // reuse the fluvial scratch as raw-erosion buffer
        for (k = 0; k < N; k++) {
          iceMask[k] = H[k] > 1 ? 1 : 0;
          var E = laws.abrasion(p.Kg * fKg[k], st.Us[k], p.lexp);
          if (p.Kq > 0 && k > 0 && k < N - 1) {
            var conv = -(z[k - 1] - 2 * z[k] + z[k + 1]) / (ds * ds); // convex-up: z'' < 0
            E += laws.quarrying(p.Kq * fKg[k], st.Us[k], conv);
          }
          Eraw[k] = iceMask[k] ? E : 0;
        }
        for (k = 0; k < N; k++) {
          // erosion footprint [¼ ½ ¼]: glacial erosion acts over an ice-thickness-scale patch,
          // not one cell; it also removes the two-cell (checkerboard) mode of the E ∝ 1/H feedback.
          var Es = !p.eroSmooth ? Eraw[k] : (k > 0 && k < N - 1) ? 0.25 * Eraw[k - 1] + 0.5 * Eraw[k] + 0.25 * Eraw[k + 1]
                 : (k === 0 ? 0.5 * (Eraw[0] + Eraw[1]) : 0.5 * (Eraw[N - 2] + Eraw[N - 1]));
          Es = laws.capRate(Es, p.eroCap);
          if (k === N - 1) Es = 0;
          st.Eg[k] = Es;
          z[k] -= Es * dt; gv += Es * dt * ds * st.W[k];
        }
        st.eroGlac += gv;
      } else {
        for (k = 0; k < N; k++) { iceMask[k] = 0; st.H[k] = 0; st.Us[k] = 0; st.Eg[k] = 0; st.q[k] = 0; }
      }
      // 3. fluvial on ice-free nodes
      fluvialStep();
      // 4. threshold-slope failure (headwall retreat / step collapse): no cell may stand steeper than
      //    Sc above its downstream neighbour, ice or no ice (a 100 m cell step steeper than Sc is
      //    treated as unstable: rockfall subaerially, block failure of the step under ice); the
      //    excess is removed and exported (debris evacuated by the glacier / river below —
      //    asserted). Sweep from the outlet up so each receiver is final.
      var lim = p.Sc * ds, rock = 0;
      for (k = 0; k < N; k++) st.Er[k] = 0;
      for (k = N - 2; k >= 0; k--) {
        var ex = z[k] - z[k + 1] - lim;
        if (ex > 0) { z[k] -= ex; rock += ex * ds * st.W[k]; st.Er[k] = ex / dt; }
      }
      st.eroRock += rock;
      st.t += dt;
    }

    // Solve the steady ice on the initial bed without eroding, so the t = 0 state (and a
    // node isolated from it) carries the glacier the first step would see.
    function primeIce() {
      if (!p.glacierOn) return;
      steadyIce(ELA(st.t));
    }

    function diagnostics() {
      var z = st.z, H = st.H, zmax = -1e9, zmin = 1e9, vol = 0, term = -1, hmax = 0, k;
      for (k = 0; k < N; k++) {
        if (z[k] > zmax) zmax = z[k];
        if (z[k] < zmin) zmin = z[k];
        vol += H[k] * st.W[k] * ds;
        if (H[k] > 1) term = st.s[k];
        if (H[k] > hmax) hmax = H[k];
      }
      var dV = 0;
      for (k = 0; k < N; k++) dV += (z[k] - st.z0[k]) * ds * st.W[k];
      return { t: st.t, relief: zmax - zmin, zmax: zmax, iceVol: vol, terminus: term, maxH: hmax,
               closure: st.upliftVol - st.eroFluv - st.eroGlac - st.eroRock - dV, upliftVol: st.upliftVol };
    }

    function record() { st.history.push(diagnostics()); }

    // Everything a d1 node column needs to be born from node i, frozen at this instant.
    function nodeSpec(i) {
      if (i < 0) i = 0; if (i > N - 2) i = N - 2;
      return {
        i: i, s: st.s[i], z0: st.z[i], ucum0: st.ucum[i], zDown: st.z[i + 1], HDown: st.H[i + 1],
        Q: st.q[i] * st.W[i], A: st.A[i], W: st.W[i], Wf: 0.5 * (st.W[i] + st.W[i + 1]), ds: ds, fac: fac[i],
        eDown: st.Ef[i + 1] + st.Eg[i + 1] + st.Er[i + 1],
        t0: st.t, dt: p.dt, params: p, litho: p.litho
      };
    }

    return { state: st, p: p, step: step, primeIce: primeIce, diagnostics: diagnostics, record: record, nodeSpec: nodeSpec, U: U, ELA: ELA };
  }

  // Lithology spec from the shared control values (SI units). Used by the d2_along page
  // and by the d1 node scenario so a column opened from a node carries the same body.
  function lithoFromControls(type, topM, thickM, pos, Lm) {
    if (type === "layer") return { type: "layer", top: topM, thick: thickM, soft: 5 };
    if (type === "slab") return { type: "slab", top: topM, thick: thickM, dipX: -pos * 0.2, soft: 5 };
    if (type === "dike") return { type: "dike", x0: pos * Lm, width: thickM, soft: 5 };
    return null;
  }

  var api = { DEFAULTS: DEFAULTS, createModel: createModel, lithoFromControls: lithoFromControls, hackArea: hackArea, valleyWidth: valleyWidth, marchThickness: marchThickness, routeFlux: routeFlux };
  GS.d2along = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : this);
