/*
 * model.js — Pure 1D uplift / relief / erosion reservoir model.
 *
 * No DOM, no D3. Just the physics, so the same code can run in a Django
 * template, a worker, or a test. Exposed as window.UpliftModel (browser)
 * or module.exports (node).
 *
 * Conceptual model
 * ----------------
 *   Relief R is a reservoir (units: km).
 *   Uplift U(t) adds to the reservoir.
 *   Erosion E(R) removes from it, and is *nonlinear* in relief:
 *       E(R) = k * R^n          (n defaults to 4)
 *   so erosion becomes overwhelmingly strong as relief grows. This creates a
 *   "relief limit": the equilibrium relief barely changes even when uplift
 *   rate changes a lot, because R sits under an n-th root.
 *
 *       dR/dt = U(t) - k * R^n
 *
 *   Equilibrium (steady) relief for a sustained uplift U:
 *       R_eq = (U / k)^(1/n)
 *
 * Units
 * -----
 *   time   : Myr
 *   relief : km
 *   uplift : km/Myr   (note: 1 mm/yr = 1 km/Myr exactly, so these read the
 *            same as the familiar mm/yr rock-uplift rates)
 *   erosion: km/Myr
 *   net    : km/Myr
 *   All rates share km/Myr, so the rate axes are directly comparable and no
 *   unit conversion is needed anywhere.
 */
(function (root) {
  "use strict";

  // ---- Uplift forcing functions U(t) -------------------------------------
  // Each returns uplift rate in km/Myr for a given time t (Myr) over [0, T].

  var upliftFunctions = {
    // Gaussian "bell curve": smooth rise and fall peaking at T/2. sigma is set
    // so ~99.7% (±3 sigma) of the bell fits inside the event window, leaving
    // the edges essentially zero.
    bell: function (t, p) {
      var mu = p.duration / 2;
      var sigma = p.duration / 6;
      var z = (t - mu) / sigma;
      return p.peakUplift * Math.exp(-0.5 * z * z);
    },

    // Circular / semicircular arc: rises and falls smoothly, peaking at T/2.
    // Shape of the upper half of a circle mapped onto the time window.
    arc: function (t, p) {
      var x = (t - p.duration / 2) / (p.duration / 2); // -1 .. 1
      if (x <= -1 || x >= 1) return 0;
      return p.peakUplift * Math.sqrt(1 - x * x);
    },

    // Constant uplift during the event, then off — a square pulse.
    constant: function (t, p) {
      return t > p.duration ? 0 : p.peakUplift;
    },

    // Linear ramp up then down (triangle) — for comparison with the arc.
    triangle: function (t, p) {
      var x = Math.abs((t - p.duration / 2) / (p.duration / 2)); // 0..1
      return p.peakUplift * Math.max(0, 1 - x);
    },

    // Plateau / trapezoid: ramp up, hold steady, then ramp down. Rise and fall
    // each take 25% of the event; the steady hold occupies the middle 50%.
    plateau: function (t, p) {
      var r = p.duration * 0.25; // ramp length
      if (t <= 0 || t >= p.duration) return 0;
      if (t < r) return p.peakUplift * (t / r); // ramp up
      if (t > p.duration - r) return p.peakUplift * (p.duration - t) / r; // ramp down
      return p.peakUplift; // hold
    },

    // Tectonic pulse: quick onset then exponential decay within the event.
    // Sharp rise over the first 10%, then relaxes with a 0.3*duration timescale.
    pulse: function (t, p) {
      if (t <= 0) return 0;
      var rise = p.duration * 0.1;
      if (t < rise) return p.peakUplift * (t / rise);
      var tau = p.duration * 0.3;
      return p.peakUplift * Math.exp(-(t - rise) / tau);
    }
  };

  // ---- Erosion law -------------------------------------------------------
  // E(R) = k * R^n  (km/Myr). k is derived from an intuitive "relief limit":
  // the equilibrium relief that *would* be reached under sustained peak uplift.
  //     R_lim = (U_peak / k)^(1/n)  =>  k = U_peak / R_lim^n
  function erodibility(p) {
    return p.peakUplift / Math.pow(p.reliefLimit, p.exponent); // km/Myr, km
  }

  function erosionRate(R, k, n) {
    return k * Math.pow(Math.max(R, 0), n); // km/Myr
  }

  // Instantaneous equilibrium relief (km) the current uplift is "aiming" at.
  function equilibriumRelief(U, k, n) {
    return Math.pow(U / k, 1 / n);
  }

  // ---- Integration -------------------------------------------------------
  // RK4 on dR/dt = U(t) - k*R^n (all km/Myr). RK4 keeps the stiff R^n stable.
  function run(params) {
    var p = Object.assign(
      {
        duration: 10, // Myr, length of the uplift event
        tailFactor: 2, // simulate out to duration * tailFactor (relaxation tail)
        peakUplift: 2, // km/Myr
        exponent: 4, // erosion nonlinearity
        reliefLimit: 3, // km, equilibrium relief at peak uplift
        initialRelief: 0, // km
        steps: 2000, // integration / output resolution
        upliftType: "bell"
      },
      params || {}
    );

    var Ufn = upliftFunctions[p.upliftType] || upliftFunctions.arc;
    var k = erodibility(p);
    var n = p.exponent;
    var simTime = p.duration * p.tailFactor; // total time including relaxation tail
    var dt = simTime / p.steps;

    var U_at = function (t) {
      return Ufn(t, p); // km/Myr
    };
    var dRdt = function (t, R) {
      return U_at(t) - erosionRate(R, k, n);
    };

    var series = [];
    var R = p.initialRelief;

    for (var i = 0; i <= p.steps; i++) {
      var t = i * dt;
      var U = U_at(t); // km/Myr
      var E = erosionRate(R, k, n); // km/Myr
      series.push({
        t: t,
        uplift: U, // km/Myr
        upliftFlux: U, // km/Myr (same units now; overlaid vs erosion)
        relief: R, // km
        erosion: E, // km/Myr
        net: U - E, // km/Myr, net rate = dR/dt
        reliefEq: equilibriumRelief(U, k, n) // km, instantaneous target
      });

      // advance with RK4
      var k1 = dRdt(t, R);
      var k2 = dRdt(t + dt / 2, R + (dt / 2) * k1);
      var k3 = dRdt(t + dt / 2, R + (dt / 2) * k2);
      var k4 = dRdt(t + dt, R + dt * k3);
      R = R + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
      if (R < 0) R = 0;
    }

    return { params: p, k: k, series: series };
  }

  var api = {
    run: run,
    upliftFunctions: upliftFunctions,
    erodibility: erodibility,
    erosionRate: erosionRate,
    equilibriumRelief: equilibriumRelief
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.UpliftModel = api;
  }
})(typeof window !== "undefined" ? window : this);
