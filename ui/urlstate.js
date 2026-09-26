/*
 * ui/urlstate.js — round-trip control state through the URL query string so any
 * configured view can be bookmarked or shared. Shared by every tier.
 *
 * specs: [{id, type: "num" | "select" | "bool", val, min, max, options}]
 *   read()       → state object (URL values clamped to [min, max] / valid options,
 *                  defaults otherwise)
 *   write(state) → history.replaceState (no history entries; silently ignored on
 *                  file:// where some browsers block it)
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  GS.ui = GS.ui || {};

  GS.ui.urlState = function (specs) {
    function read() {
      var q = new URLSearchParams(root.location ? root.location.search : "");
      var state = {};
      specs.forEach(function (c) {
        var raw = q.get(c.id), v = c.val;
        if (raw !== null && raw !== "") {
          if (c.type === "bool") v = raw === "1" || raw === "true";
          else if (c.type === "select") v = c.options.indexOf(raw) !== -1 ? raw : c.val;
          else {
            var num = +raw;
            if (!isNaN(num)) {
              if (c.min != null && num < c.min) num = c.min;
              if (c.max != null && num > c.max) num = c.max;
              v = num;
            }
          }
        }
        state[c.id] = v;
      });
      return state;
    }
    function write(state) {
      var q = new URLSearchParams();
      specs.forEach(function (c) {
        var v = state[c.id];
        q.set(c.id, c.type === "bool" ? (v ? 1 : 0) : v);
      });
      var qs = "?" + q.toString();
      if (root.location && qs !== root.location.search) {
        try { root.history.replaceState(null, "", root.location.pathname + qs); } catch (e) { /* nicety only */ }
      }
    }
    return { read: read, write: write, specs: specs };
  };
})(typeof window !== "undefined" ? window : this);
