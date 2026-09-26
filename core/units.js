/*
 * core/units.js — the one unit convention, and display formatting.
 *
 * Internal units in every tier: metres and years. Rates are m/yr (so 1 mm/yr is
 * 1e-3 m/yr), lengths m, times yr. The UI converts at the edge with these helpers;
 * the model core never sees mm, km, kyr or Myr.
 *
 * Exposed as GS.units (browser) or module.exports (node).
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});

  var KYR = 1e3, MYR = 1e6, KM = 1e3, MM = 1e-3;

  var api = {
    KYR: KYR, MYR: MYR, KM: KM, MM: MM,
    // to internal
    fromMmyr: function (v) { return v * MM; },
    fromKm: function (v) { return v * KM; },
    fromMyr: function (v) { return v * MYR; },
    fromKyr: function (v) { return v * KYR; },
    // to display
    toMmyr: function (v) { return v / MM; },
    toKm: function (v) { return v / KM; },
    toMyr: function (v) { return v / MYR; },
    toKyr: function (v) { return v / KYR; },
    // formatted strings for readouts
    fmtRate: function (mPerYr, digits) {
      return (mPerYr / MM).toFixed(digits != null ? digits : 2) + " mm/yr";
    },
    fmtLen: function (m, digits) {
      if (Math.abs(m) >= KM) return (m / KM).toFixed(digits != null ? digits : 2) + " km";
      return m.toFixed(0) + " m";
    },
    fmtTime: function (yr, digits) {
      if (Math.abs(yr) >= MYR) return (yr / MYR).toFixed(digits != null ? digits : 2) + " Myr";
      if (Math.abs(yr) >= KYR) return (yr / KYR).toFixed(digits != null ? digits : 0) + " kyr";
      return yr.toFixed(0) + " yr";
    }
  };
  GS.units = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : this);
