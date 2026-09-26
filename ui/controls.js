/*
 * ui/controls.js — build slider / select / checkbox rows from a spec list (D3).
 *
 * spec fields: id, label, type ("range" | "select" | "check"), val, min, max, step,
 *   unit, fmt(v) → string, options [{v, t}], group (fieldset legend), asserted
 *   (true → an "asserted" badge: the tier does not resolve this quantity and takes
 *   it as a parameter, DESIGN.md §2), note (small text under the control),
 *   log (true → the slider holds log10 of the quantity; ids conventionally end in
 *   "Exp"; the displayed value is 10^v), method (HTML: equation + method + source,
 *   shown on hover of the ⓘ glyph via ui/tooltip.js — every control should have one).
 *
 * buildControls(root, specs, state, onChange) mutates state in place, calls
 * onChange(id) after each change, and returns {refresh()} to push state → DOM.
 * The urlState specs are derived with GS.ui.toUrlSpecs(specs).
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  GS.ui = GS.ui || {};

  function display(c, v) {
    if (c.fmt) return c.fmt(v);
    var x = c.log ? Math.pow(10, v) : v;
    var s = c.log ? x.toExponential(1) : String(x);
    return s + (c.unit ? " " + c.unit : "");
  }

  GS.ui.buildControls = function (rootSel, specs, state, onChange) {
    var groups = {}, order = [];
    specs.forEach(function (c) {
      if (c.hidden) return; // URL-only state, no control
      var g = c.group || "";
      if (!groups[g]) { groups[g] = []; order.push(g); }
      groups[g].push(c);
    });
    var inputs = {};
    order.forEach(function (g) {
      var host = rootSel;
      if (g) {
        var fs = rootSel.append("fieldset");
        fs.append("legend").text(g);
        host = fs;
      }
      groups[g].forEach(function (c) {
        var row = host.append("div").attr("class", "control" + (c.type === "check" ? " check" : ""));
        var lab = row.append("label").attr("for", c.id);
        if (c.type === "check") {
          var cb = lab.append("input").attr("type", "checkbox").attr("id", c.id)
            .property("checked", !!state[c.id])
            .on("change", function () { state[c.id] = this.checked; onChange(c.id); });
          lab.append("span").text(" " + c.label);
          inputs[c.id] = cb;
        } else {
          lab.append("span").text(c.label + " ");
          if (c.type === "range") lab.append("span").attr("class", "val").attr("id", c.id + "-val").text(display(c, state[c.id]));
        }
        if (c.asserted) lab.append("span").attr("class", "badge asserted").attr("title", "Asserted: this tier does not resolve it; it is taken as a parameter.").text("asserted");
        if (c.method && GS.ui.tip) GS.ui.tip.glyph(lab, c.method);
        if (c.type === "range") {
          inputs[c.id] = row.append("input").attr("type", "range").attr("id", c.id)
            .attr("min", c.min).attr("max", c.max).attr("step", c.step).property("value", state[c.id])
            .on("input", function () {
              state[c.id] = +this.value;
              d3.select("#" + c.id + "-val").text(display(c, state[c.id]));
              onChange(c.id);
            });
        } else if (c.type === "select") {
          var sel = row.append("select").attr("id", c.id).on("change", function () { state[c.id] = this.value; onChange(c.id); });
          c.options.forEach(function (o) {
            sel.append("option").attr("value", o.v).property("selected", o.v === state[c.id]).text(o.t);
          });
          inputs[c.id] = sel;
        }
        if (c.note) row.append("div").attr("class", "note").text(c.note);
      });
    });
    function refresh() {
      specs.forEach(function (c) {
        var inp = inputs[c.id];
        if (!inp) return;
        if (c.type === "check") inp.property("checked", !!state[c.id]);
        else inp.property("value", state[c.id]);
        if (c.type === "range") d3.select("#" + c.id + "-val").text(display(c, state[c.id]));
      });
    }
    return { refresh: refresh, inputs: inputs };
  };

  // Hold a box's height at the tallest it has been (grow-only), so a readout whose text
  // length changes step to step does not bounce the controls below it. holdHeight(sel, true)
  // releases the hold (e.g. on Reset).
  GS.ui.holdHeight = function (sel, release) {
    sel.each(function () {
      if (release) { this.style.minHeight = ""; this._held = 0; return; }
      var h = this.offsetHeight;
      if (!this._held || h > this._held) { this._held = h; this.style.minHeight = h + "px"; }
    });
  };

  // Control specs → urlState specs.
  GS.ui.toUrlSpecs = function (specs) {
    return specs.map(function (c) {
      return {
        id: c.id,
        type: c.type === "check" ? "bool" : (c.type === "select" ? "select" : "num"),
        val: c.val, min: c.min, max: c.max,
        options: c.options ? c.options.map(function (o) { return o.v; }) : null
      };
    });
  };
})(typeof window !== "undefined" ? window : this);
