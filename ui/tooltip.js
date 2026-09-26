/*
 * ui/tooltip.js — one floating methodology tooltip for the whole page.
 *
 * Every control, panel and readout carries a `method` HTML string: the governing
 * equation(s), what the tier computes vs asserts, and the source. Hovering the
 * ⓘ glyph (controls) or a panel title shows it.
 *
 *   GS.ui.tip.attach(selection, htmlOrFn)   show on mouseenter, follow, hide on leave
 *   GS.ui.tip.glyph(parentSel, html)        append an ⓘ span with the tooltip attached
 * Equations are plain HTML (Unicode, <sup>, <sub>, <i>); no maths library.
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  GS.ui = GS.ui || {};

  var box = null;
  function ensure() {
    if (!box) box = d3.select("body").append("div").attr("class", "gs-tip").style("display", "none");
    return box;
  }
  function place(ev) {
    var b = ensure(), pad = 14, w = 360;
    var x = ev.clientX + pad, y = ev.clientY + pad;
    var vw = root.innerWidth, vh = root.innerHeight, bh = b.node().offsetHeight;
    if (x + w > vw - 8) x = ev.clientX - w - pad;
    if (y + bh > vh - 8) y = Math.max(8, ev.clientY - bh - pad);
    b.style("left", x + "px").style("top", y + "px").style("width", w + "px");
  }
  function attach(sel, html) {
    sel.classed("has-tip", true)
      .on("mouseenter.tip", function (ev) {
        var h = typeof html === "function" ? html() : html;
        if (!h) return;
        ensure().html(h).style("display", "block");
        place(ev);
      })
      .on("mousemove.tip", function (ev) { if (box) place(ev); })
      .on("mouseleave.tip", function () { if (box) box.style("display", "none"); });
    return sel;
  }
  function glyph(parent, html) {
    var g = parent.append("span").attr("class", "tip-glyph").attr("aria-label", "methodology").text("ⓘ");
    return attach(g, html);
  }
  GS.ui.tip = { attach: attach, glyph: glyph };
})(typeof window !== "undefined" ? window : this);
