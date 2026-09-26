/*
 * d1/model.js — the column tier (z(t)). Pure model, no DOM.
 *
 * Scenario "reservoir": the legacy uplift_sandbox relief-reservoir model,
 *     dR/dt = U(t) − k Rⁿ,   k = U_peak / R_limⁿ,
 * integrated with RK4 (core/laws.js). Units: metres, years (1 mm/yr = 1e-3 m/yr).
 * Everything the column does not resolve (the drainage network that sets k and n)
 * is asserted, and the UI says so.
 *
 * Mirrored in py/d1_reservoir.py; validated in py/test_core.py (steady state,
 * dt-independence, scipy reference, float32, legacy km/Myr equivalence).
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  var forcing = GS.forcing, laws = GS.laws;

  function runReservoir(params) {
    var p = Object.assign({
      duration: 12e6,      // yr, length of the uplift event
      tailFactor: 2,       // simulate to duration × tailFactor (relaxation tail)
      peakUplift: 3e-3,    // m/yr
      exponent: 4,         // erosion nonlinearity n
      reliefLimit: 3000,   // m, equilibrium relief at peak uplift
      initialRelief: 0,    // m
      steps: 2000,
      shape: "bell"
    }, params || {});

    var U = forcing.makeSeries({ shape: p.shape, peak: p.peakUplift, duration: p.duration });
    var k = laws.reservoirK(p.peakUplift, p.reliefLimit, p.exponent);
    var n = p.exponent;
    var simTime = p.duration * p.tailFactor;
    var dt = simTime / p.steps;
    var dRdt = function (t, R) { return U(t) - laws.reservoirErosion(R, k, n); };

    var series = [];
    var R = p.initialRelief;
    for (var i = 0; i <= p.steps; i++) {
      var t = i * dt;
      var u = U(t), E = laws.reservoirErosion(R, k, n);
      series.push({ t: t, uplift: u, relief: R, erosion: E, net: u - E, reliefEq: laws.reservoirEquilibrium(u, k, n) });
      R = laws.rk4(dRdt, t, R, dt);
      if (R < 0) R = 0;
    }
    return { params: p, k: k, series: series };
  }

  var api = { runReservoir: runReservoir };
  GS.d1 = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : this);

/*
 * Scenario "columns": two lockstep columns A (weak) and B (with a resistant body)
 * under one U(t) and one ELA(t). Mirrored in py/d1_columns.py, validated in
 * py/test_d1_columns.py (reservoir regression, bookkeeping, persistence timescale
 * T·c/U, exposure at depth/U, analytic buzzsaw cap, ELA cycle, dt-independence).
 *
 *   dz/dt = U(t) − E_sub − E_gl                     per column, z above base level
 *   E_ref = k z_Aⁿ                                  A is the landscape's relief reservoir
 *   E_sub,B = D_B · f(r_B, c),  f = 1/(1 + r(c−1))   D_B = asserted demand (p.coupling):
 *       "landscape" k z_Aⁿ · "own" k z_Bⁿ · "max" k max(z_A, z_B)ⁿ
 *   E_gl = P g(z) · f(r, c_g),  g = exp(−½u²) below the ELA, φ + (1−φ) exp(−½u²) above
 *   r_B = layer(z_B − U_cum)                        body in the material frame
 */
(function (root) {
  "use strict";
  var GS = root.GS, forcing = GS.forcing, laws = GS.laws, litho = GS.lithology;

  var COLUMN_DEFAULTS = {
    duration: 12e6, tailFactor: 2, peakUplift: 3e-3, shape: "bell",
    exponent: 4, reliefLimit: 3000, initialRelief: 0,
    bodyDepth: 0, bodyThick: 500, bodySoft: 10, contrast: 5, contrastG: 5,
    glacierOn: false, elaBase: 2000, elaAmp: 0, elaPeriod: 1e5,
    buzzPeak: 1e-3, buzzWidth: 300, buzzAbove: 0,
    coupling: "landscape",
    steps: 2000
  };

  function runColumns(params) {
    var p = Object.assign({}, COLUMN_DEFAULTS, params || {});
    var U = forcing.makeSeries({ shape: p.shape, peak: p.peakUplift, duration: p.duration });
    var ELA = forcing.makeSeries({ shape: "sine", peak: -p.elaAmp, base: p.elaBase, period: p.elaPeriod });
    var k = laws.reservoirK(p.peakUplift, p.reliefLimit, p.exponent);
    var n = p.exponent, z0 = p.initialRelief;
    var body = litho.make({ type: "layer", top: z0 - p.bodyDepth, thick: p.bodyThick, soft: p.bodySoft });
    var gl = p.glacierOn, P = p.buzzPeak, w = p.buzzWidth, phi = p.buzzAbove, c = p.contrast, cg = p.contrastG;
    var coupling = p.coupling;
    function buzz(z, ela) {
      var e = laws.buzzsawRate(z, ela, w, P);
      return z <= ela ? e : phi * P + (1 - phi) * e;
    }

    // y = [zA, zB, Ucum, EcumA, EcumB]; returns rates and the diagnostic pieces
    var aux = {};
    function rates(t, y, out) {
      var zA = y[0], zB = y[1], ucum = y[2];
      var u = U(t), ela = ELA(t);
      var Eref = laws.reservoirErosion(zA, k, n);
      var DB = coupling === "own" ? laws.reservoirErosion(zB, k, n)
             : coupling === "max" ? laws.reservoirErosion(zA > zB ? zA : zB, k, n) : Eref;
      var rB = body(0, 0, litho.materialZ(zB, ucum));
      var EgA = gl ? buzz(zA, ela) : 0;
      var EgB = gl ? buzz(zB, ela) * litho.erodibilityFactor(rB, cg) : 0;
      var EsubB = DB * litho.erodibilityFactor(rB, c);
      var EA = Eref + EgA, EB = EsubB + EgB;
      out[0] = u - EA; out[1] = u - EB; out[2] = u; out[3] = EA; out[4] = EB;
      aux.u = u; aux.ela = ela; aux.Eref = Eref; aux.EsubB = EsubB; aux.EgA = EgA; aux.EgB = EgB; aux.rB = rB;
      return out;
    }

    var sim = p.duration * p.tailFactor, dt = sim / p.steps;
    var y = [z0, z0, 0, 0, 0], k1 = [0, 0, 0, 0, 0], k2 = [0, 0, 0, 0, 0], k3 = [0, 0, 0, 0, 0], k4 = [0, 0, 0, 0, 0], tmp = [0, 0, 0, 0, 0];
    var series = [], j;
    for (var i = 0; i <= p.steps; i++) {
      var t = i * dt;
      rates(t, y, k1);
      var top = z0 - p.bodyDepth + y[2];
      series.push({
        t: t, uplift: aux.u, ela: aux.ela, zA: y[0], zB: y[1], ucum: y[2], ecumA: y[3], ecumB: y[4],
        EA: aux.Eref + aux.EgA, EB: aux.EsubB + aux.EgB, EsubA: aux.Eref, EsubB: aux.EsubB, EgA: aux.EgA, EgB: aux.EgB,
        rB: aux.rB, bodyTop: top, bodyBot: top - p.bodyThick, diff: y[1] - y[0]
      });
      for (j = 0; j < 5; j++) tmp[j] = y[j] + dt / 2 * k1[j];
      rates(t + dt / 2, tmp, k2);
      for (j = 0; j < 5; j++) tmp[j] = y[j] + dt / 2 * k2[j];
      rates(t + dt / 2, tmp, k3);
      for (j = 0; j < 5; j++) tmp[j] = y[j] + dt * k3[j];
      rates(t + dt, tmp, k4);
      for (j = 0; j < 5; j++) y[j] = y[j] + dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
      if (y[0] < 0) y[0] = 0;
      if (y[1] < 0) y[1] = 0;
    }
    return { params: p, k: k, series: series };
  }

  GS.d1.runColumns = runColumns;
  GS.d1.COLUMN_DEFAULTS = COLUMN_DEFAULTS;
})(typeof window !== "undefined" ? window : this);
