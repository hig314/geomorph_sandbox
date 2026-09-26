/*
 * ui/panels.js — stacked D3 time-series panels sharing an x axis.
 * Ported from legacy/uplift_sandbox/app.js and parameterised by geometry.
 *
 *   var P = GS.ui.panels({width: 820, height: 150});
 *   var pn = P.make(parentSel, {title, yLabel, color, showX, xLabel, method});
 *   (method: HTML shown on hovering the panel title — the equations behind the lines)
 *   P.line(pn, "name", data, yAcc, color, dashed, xAcc, tip)   // xAcc defaults to d.t
 *   P.area(pn, "name", data, y0Acc, y1Acc, color, xAcc, tip)        // filled band
 *   tip: HTML shown on hovering the curve (a wide invisible hit stroke catches the pointer);
 *   build it with GS.ui.curveTip(label, color, dashed, text). Visible curves themselves
 *   ignore the pointer so the hit strokes and the panel's own handlers see it.
 *   P.endLabel(pn, "name", lastDatum, yAcc, color, text)
 *   P.vline(pn, "name", xValue)   P.zero(pn)   P.axes(pn)
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  GS.ui = GS.ui || {};

  GS.ui.panels = function (opts) {
    opts = opts || {};
    var margin = opts.margin || { top: 18, right: 70, bottom: 8, left: 64 };
    var fullW = opts.width || 820, panelH = opts.height || 150;
    var innerW = fullW - margin.left - margin.right;
    var tAcc = function (d) { return d.t; };

    function make(parent, o) {
      var h = panelH, extra = o.showX ? 30 : 0;
      var svg = parent.append("svg").attr("class", "panel")
        .attr("width", fullW).attr("height", h + margin.top + margin.bottom + extra)
        .attr("viewBox", "0 0 " + fullW + " " + (h + margin.top + margin.bottom + extra));
      var g = svg.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
      var x = d3.scaleLinear().range([0, innerW]);
      var y = d3.scaleLinear().range([h, 0]);
      var xAxisG = g.append("g").attr("class", "axis x-axis").attr("transform", "translate(0," + h + ")");
      var yAxisG = g.append("g").attr("class", "axis y-axis");
      g.append("text").attr("class", "axis-label y-label").attr("transform", "rotate(-90)")
        .attr("x", -h / 2).attr("y", -margin.left + 14).attr("text-anchor", "middle").text(o.yLabel || "");
      var title = g.append("text").attr("class", "panel-title").attr("x", 2).attr("y", -4).text(o.title || "").style("fill", o.color || null);
      if (o.method && GS.ui.tip) { title.text((o.title || "") + " ⓘ"); GS.ui.tip.attach(title, o.method); }
      if (o.showX && o.xLabel) {
        g.append("text").attr("class", "axis-label x-label").attr("x", innerW / 2).attr("y", h + 30).attr("text-anchor", "middle").text(o.xLabel);
      }
      return { svg: svg, g: g, x: x, y: y, h: h, w: innerW, xAxisG: xAxisG, yAxisG: yAxisG, showX: !!o.showX };
    }

    function hitPath(panel, name, d, tip, width, fillHit) {
      if (!tip || !GS.ui.tip) return;
      var hit = panel.g.selectAll("path.hit-" + name).data([0]);
      var ent = hit.enter().append("path").attr("class", "hit hit-" + name)
        .attr("fill", fillHit ? "transparent" : "none").attr("stroke", "transparent").attr("stroke-width", width || 12)
        .style("pointer-events", fillHit ? "all" : "stroke").style("cursor", "help");
      GS.ui.tip.attach(ent, tip);
      ent.merge(hit).attr("d", d).raise();
    }

    function line(panel, name, data, yAcc, color, dashed, xAcc, tip) {
      xAcc = xAcc || tAcc;
      var ln = d3.line().x(function (d) { return panel.x(xAcc(d)); }).y(function (d) { return panel.y(yAcc(d)); }).curve(d3.curveMonotoneX);
      var path = panel.g.selectAll("path.line-" + name).data([data]);
      path.enter().append("path").attr("class", "line line-" + name).attr("fill", "none").style("pointer-events", "none")
        .attr("stroke", color).attr("stroke-width", dashed ? 1.5 : 2).attr("stroke-dasharray", dashed ? "4 3" : null)
        .merge(path).attr("d", ln);
      hitPath(panel, name, ln(data), tip, 12, false);
    }

    function area(panel, name, data, y0Acc, y1Acc, color, xAcc, tip) {
      xAcc = xAcc || tAcc;
      var ar = d3.area().x(function (d) { return panel.x(xAcc(d)); })
        .y0(function (d) { return panel.y(y0Acc(d)); }).y1(function (d) { return panel.y(y1Acc(d)); });
      var path = panel.g.selectAll("path.area-" + name).data([data]);
      path.enter().append("path").attr("class", "area area-" + name).attr("fill", color).attr("stroke", "none").style("pointer-events", "none")
        .merge(path).attr("d", ar);
      hitPath(panel, name, ar(data), tip, 0, true);
    }
    function clear(panel) { panel.g.selectAll("path.line, path.area, path.hit, text.endlabel, line.vline, line.zero-line").remove(); }

    function endLabel(panel, name, lastDatum, yAcc, color, text) {
      var sel = panel.g.selectAll("text.endlabel-" + name).data([lastDatum]);
      sel.enter().append("text").attr("class", "endlabel endlabel-" + name)
        .attr("x", innerW + 4).attr("dy", "0.32em").style("fill", color).text(text)
        .merge(sel).attr("y", function (d) { return panel.y(yAcc(d)); });
    }

    function vline(panel, name, xv) {
      var sel = panel.g.selectAll("line.vline-" + name).data([xv]);
      sel.enter().append("line").attr("class", "vline vline-" + name).attr("y1", 0).attr("y2", panel.h)
        .attr("stroke", "#bbb").attr("stroke-width", 1).attr("stroke-dasharray", "2 3")
        .merge(sel).attr("x1", function (d) { return panel.x(d); }).attr("x2", function (d) { return panel.x(d); });
    }

    function zero(panel) {
      var sel = panel.g.selectAll("line.zero-line").data([0]);
      sel.enter().append("line").attr("class", "zero-line").attr("x1", 0).attr("stroke", "#ccc").attr("stroke-width", 1)
        .merge(sel).attr("x2", innerW).attr("y1", panel.y(0)).attr("y2", panel.y(0));
    }

    function axes(panel, xTicks, yTicks) {
      panel.yAxisG.call(d3.axisLeft(panel.y).ticks(yTicks || 5));
      panel.xAxisG.call(d3.axisBottom(panel.x).ticks(xTicks || 8));
      if (!panel.showX) panel.xAxisG.selectAll("text").style("display", "none");
    }

    return { make: make, line: line, area: area, hitPath: hitPath, clear: clear, endLabel: endLabel, vline: vline, zero: zero, axes: axes, innerW: innerW, panelH: panelH };
  };

  // Tooltip HTML for a curve: swatch (solid/dashed/fill), name, and what it is.
  GS.ui.curveTip = function (label, color, style, text, eq) {
    var sw = style === "fill"
      ? '<span class="sw" style="background:' + color + '"></span>'
      : '<span class="sw line' + (style === "dashed" ? " dashed" : "") + '" style="border-color:' + color + '"></span>';
    return "<h4>" + sw + label + "</h4>" + (eq ? '<span class="eq">' + eq + "</span>" : "") + (text ? "<div>" + text + "</div>" : "");
  };
})(typeof window !== "undefined" ? window : this);
