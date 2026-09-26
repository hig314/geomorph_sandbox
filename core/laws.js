/*
 * core/laws.js — the process laws, each written ONCE as a pure function of local
 * quantities. Tiers supply those quantities by computing them or by asserting them
 * (DESIGN.md §2). Nothing here knows about grids, profiles or the DOM.
 *
 * Units: metres, years (core/units.js). Mirrored in py/gs_core.py; every scheme
 * that uses these laws is validated there before and after porting.
 *
 * Sources:
 *   fluvial   Yuan et al. (2019) implicit erosion–deposition node update, extracted
 *             from legacy/terrain_sandbox fluvialStep; Braun & Willett (2013) for G = 0.
 *   creep     linear diffusion; Roering et al. (2001) nonlinear flux.
 *   glacier   SIA flux (Glen n = 3), MacGregor et al. (2000) / Anderson et al. (2006)
 *             flowline glaciers; abrasion E = K_g U_s^l (Hallet 1979; l = 1 after
 *             Humphrey & Raymond 1994); mass balance with an asymmetric ablation
 *             gradient (legacy terrain_sandbox, ICE_ABL_RATIO).
 *   reservoir the legacy uplift_sandbox relief-reservoir law dR/dt = U − k Rⁿ.
 * Landslide (Bishop slices / method of columns) laws are added with d2_across.
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});

  // ---- Glacier constants (temperate ice, Glen n = 3) ------------------------
  var ICE = {
    RHOG: 8987.0,        // ρ_ice g  (Pa/m)
    N: 3,                // Glen exponent
    A: 7.57e-17,         // Glen rate factor (Pa⁻³ yr⁻¹)
    ABL_RATIO: 2.5,      // ablation gradient ÷ accumulation gradient
    SLOPE_MIN: 1e-3      // slope floor for the thickness inversion (divides)
  };
  ICE.GAMMA = (2 * ICE.A / (ICE.N + 2)) * Math.pow(ICE.RHOG, ICE.N); // SIA flux coefficient

  // ---- Reservoir (d1) --------------------------------------------------------
  // E(R) = k Rⁿ. k from an intuitive "relief limit": equilibrium relief at U_peak.
  function reservoirK(Upeak, Rlim, n) { return Upeak / Math.pow(Rlim, n); }
  function reservoirErosion(R, k, n) { return k * Math.pow(R > 0 ? R : 0, n); }
  function reservoirEquilibrium(U, k, n) { return Math.pow(U / k, 1 / n); }

  // ---- Fluvial ---------------------------------------------------------------
  // Detachment-limited stream power rate (m/yr).
  function streamPower(K, A, S, m, n) { return K * Math.pow(A, m) * Math.pow(S > 0 ? S : 0, n); }

  // One node of the implicit Yuan et al. (2019) erosion–deposition update.
  //   h0    elevation at the start of the step
  //   hcur  current iterate (Newton start for n ≠ 1)
  //   hr    receiver elevation (already updated this sweep)
  //   Kp    K·Aᵐ (× lithology factor) at the node
  //   dt    step (yr);  dist  distance to receiver (m)
  //   dep   deposition this step, dt·G·Qs/A (m, ≥ 0)
  //   nexp  slope exponent
  //   submerged  true → no incision (lake), deposit only
  // Returns the new elevation (caller applies any under-relaxation).
  function yuanNode(h0, hcur, hr, Kp, dt, dist, dep, nexp, submerged) {
    if (dep < 0) dep = 0;
    var nh;
    if (h0 > hr && !submerged) {
      if (Math.abs(nexp - 1) < 1e-9) {
        var f = (Kp * dt) / dist;
        nh = (h0 + f * hr + dep) / (1 + f);
      } else {
        var z = hcur;
        for (var it = 0; it < 20; it++) {
          var s = (z - hr) / dist;
          if (s < 0) s = 0;
          var g = z - h0 + Kp * dt * Math.pow(s, nexp) - dep;
          var dg = 1 + (Kp * dt * nexp * Math.pow(s, nexp - 1)) / dist;
          var dz = g / dg;
          z -= dz;
          if (Math.abs(dz) < 1e-4) break;
        }
        nh = z;
      }
      if (nh < hr) nh = hr;
    } else {
      nh = h0 + dep;
    }
    return nh;
  }

  // ---- Hillslope creep (flux per unit width, m²/yr, positive downslope of S<0) --
  function creepFluxLinear(D, S) { return -D * S; }
  // Roering: q = −D S / (1 − (S/Sc)²), the denominator clamped at fmin near Sc.
  function creepFluxRoering(D, S, Sc, fmin) {
    var f = 1 - (S * S) / (Sc * Sc);
    if (f < (fmin != null ? fmin : 0.02)) f = fmin != null ? fmin : 0.02;
    return -D * S / f;
  }

  // ---- Glacier ---------------------------------------------------------------
  // Surface mass balance (m/yr ice) at surface elevation zs. Positive above the
  // ELA with gradient bg, negative below with gradient bg·ablRatio, clamped.
  function massBalance(zs, ela, bg, accMax, abMax, ablRatio) {
    var ds = zs - ela;
    var b = ds > 0 ? ds * bg : ds * bg * (ablRatio != null ? ablRatio : ICE.ABL_RATIO);
    if (accMax != null && b > accMax) b = accMax;
    if (abMax != null && b < -abMax) b = -abMax;
    return b;
  }
  // SIA deformational flux per unit width: q = Γ H^(n+2) S^n  (m²/yr).
  function siaFlux(H, S, Gamma, n) {
    Gamma = Gamma != null ? Gamma : ICE.GAMMA; n = n != null ? n : ICE.N;
    return Gamma * Math.pow(H, n + 2) * Math.pow(Math.abs(S), n);
  }
  // Inverse: thickness carrying flux q down surface slope S. Slope floored at Smin
  // (divides). H = (q / (Γ S^n))^(1/(n+2)).
  function siaThickness(q, S, Gamma, n, Smin) {
    Gamma = Gamma != null ? Gamma : ICE.GAMMA; n = n != null ? n : ICE.N;
    Smin = Smin != null ? Smin : ICE.SLOPE_MIN;
    var s = Math.abs(S); if (s < Smin) s = Smin;
    if (q <= 0) return 0;
    return Math.pow(q / (Gamma * Math.pow(s, n)), 1 / (n + 2));
  }
  // Thickness H_i at a node whose downstream neighbour has surface zsDown and thickness
  // HDown, such that the SIA face flux carries qf (m²/yr):
  //     Gamma · (½(H_i + HDown))^(n+2) · ((zi + H_i − zsDown)/ds)^n = qf
  // Monotone in H_i → unique root; bisection to 1e-3 m. qf ≤ 0 → 0. Shared by the
  // d2_along upstream march and the d1 node column. Mirrored in py/gs_core.py.
  function faceThickness(zi, zsDown, HDown, qf, ds, Gamma, n) {
    if (qf <= 0) return 0;
    var lo = zsDown - zi; if (lo < 0) lo = 0;
    var hi = lo + 10;
    var f = function (x) {
      var Hf = 0.5 * (x + HDown), Sf = (zi + x - zsDown) / ds;
      return Gamma * Math.pow(Hf, n + 2) * Math.pow(Sf, n) - qf;
    };
    while (f(hi) < 0 && hi < 1e5) hi *= 2;
    var x = 0.5 * (lo + hi);
    for (var it = 0; it < 60; it++) {
      if (f(x) > 0) hi = x; else lo = x;
      x = 0.5 * (lo + hi);
      if (hi - lo < 1e-3) break;
    }
    return x;
  }
  // Sliding speed as a fraction fs of the depth-averaged speed q/H (m/yr).
  function slidingSpeed(q, H, fs) { return H > 0 ? fs * q / H : 0; }
  // Abrasion E = K_g U_s^l (m/yr).
  function abrasion(Kg, Us, l) { return Kg * Math.pow(Us > 0 ? Us : 0, l != null ? l : 1); }
  // Quarrying ∝ sliding × bed convexity (curvature > 0 = convex up), the knob for
  // hypothesis A. Zero on concave bed.
  function quarrying(Kq, Us, convexity) { return convexity > 0 ? Kq * Us * convexity : 0; }
  // Elevation-dependent glacial rate for tiers that assert the ice (d1): a Gaussian
  // "buzzsaw" window of half-width `width` around the ELA.
  function buzzsawRate(z, ela, width, peak) {
    var u = (z - ela) / width;
    return peak * Math.exp(-0.5 * u * u);
  }
  function capRate(E, cap) { return E > cap ? cap : E; }

  // ---- Integrators ------------------------------------------------------------
  // One RK4 step of dy/dt = f(t, y) (scalar).
  function rk4(f, t, y, dt) {
    var k1 = f(t, y);
    var k2 = f(t + dt / 2, y + (dt / 2) * k1);
    var k3 = f(t + dt / 2, y + (dt / 2) * k2);
    var k4 = f(t + dt, y + dt * k3);
    return y + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
  }

  var api = {
    ICE: ICE,
    reservoirK: reservoirK, reservoirErosion: reservoirErosion, reservoirEquilibrium: reservoirEquilibrium,
    streamPower: streamPower, yuanNode: yuanNode,
    creepFluxLinear: creepFluxLinear, creepFluxRoering: creepFluxRoering,
    massBalance: massBalance, siaFlux: siaFlux, siaThickness: siaThickness, faceThickness: faceThickness,
    slidingSpeed: slidingSpeed, abrasion: abrasion, quarrying: quarrying,
    buzzsawRate: buzzsawRate, capRate: capRate,
    rk4: rk4
  };
  GS.laws = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : this);
