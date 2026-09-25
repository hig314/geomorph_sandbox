/*
 * app.js — D3 rendering + interactive controls for the uplift/relief/erosion
 * reservoir model. Depends on model.js (window.UpliftModel) and d3 v7.
 *
 * Renders three vertically stacked line charts that share a time axis:
 *   1. Uplift  (mm/yr)   — the forcing
 *   2. Relief  (m)       — the reservoir state, with the chased equilibrium
 *   3. Erosion (m/Myr)   — the nonlinear response, vs uplift flux
 */
(function () {
  "use strict";

  var COLORS = {
    uplift: "#2c7fb8",
    relief: "#238b45",
    erosion: "#d94801",
    flux: "#888",
    eq: "#999",
    net: "#6a51a3",
    zero: "#ccc"
  };

  // Chart geometry
  var margin = { top: 18, right: 70, bottom: 8, left: 64 };
  var fullW = 820;
  var panelH = 150;
  var gap = 10;
  var innerW = fullW - margin.left - margin.right;

  // Build one panel <svg> with axes; returns helpers to update it.
  function makePanel(parent, opts) {
    var h = panelH;
    var svg = parent
      .append("svg")
      .attr("class", "panel")
      .attr("width", fullW)
      .attr("height", h + margin.top + margin.bottom)
      .attr("viewBox", "0 0 " + fullW + " " + (h + margin.top + margin.bottom));

    var g = svg
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var x = d3.scaleLinear().range([0, innerW]);
    var y = d3.scaleLinear().range([h, 0]);

    var xAxisG = g
      .append("g")
      .attr("class", "axis x-axis")
      .attr("transform", "translate(0," + h + ")");
    var yAxisG = g.append("g").attr("class", "axis y-axis");

    // y label
    g.append("text")
      .attr("class", "axis-label y-label")
      .attr("transform", "rotate(-90)")
      .attr("x", -h / 2)
      .attr("y", -margin.left + 14)
      .attr("text-anchor", "middle")
      .text(opts.yLabel);

    // panel title (top-left, inside)
    g.append("text")
      .attr("class", "panel-title")
      .attr("x", 2)
      .attr("y", -4)
      .text(opts.title)
      .style("fill", opts.color);

    return {
      svg: svg,
      g: g,
      x: x,
      y: y,
      h: h,
      xAxisG: xAxisG,
      yAxisG: yAxisG,
      showX: !!opts.showX
    };
  }

  // Append/refresh a named line path in a panel.
  function drawLine(panel, name, data, accessor, color, dashed) {
    var line = d3
      .line()
      .x(function (d) {
        return panel.x(d.t);
      })
      .y(function (d) {
        return panel.y(accessor(d));
      })
      .curve(d3.curveMonotoneX);

    var path = panel.g.selectAll("path.line-" + name).data([data]);
    path
      .enter()
      .append("path")
      .attr("class", "line line-" + name)
      .attr("fill", "none")
      .attr("stroke", color)
      .attr("stroke-width", dashed ? 1.5 : 2)
      .attr("stroke-dasharray", dashed ? "4 3" : null)
      .merge(path)
      .attr("d", line);
  }

  // Vertical guide marking the end of the uplift event.
  function drawEventEnd(panel, tEnd) {
    var sel = panel.g.selectAll("line.event-end").data([tEnd]);
    sel
      .enter()
      .append("line")
      .attr("class", "event-end")
      .attr("y1", 0)
      .attr("y2", panel.h)
      .attr("stroke", "#bbb")
      .attr("stroke-width", 1)
      .attr("stroke-dasharray", "2 3")
      .merge(sel)
      .attr("x1", function (d) { return panel.x(d); })
      .attr("x2", function (d) { return panel.x(d); });
  }

  // Small text label riding at the right end of a line.
  function drawEndLabel(panel, name, lastDatum, accessor, color, text) {
    var sel = panel.g.selectAll("text.endlabel-" + name).data([lastDatum]);
    sel
      .enter()
      .append("text")
      .attr("class", "endlabel endlabel-" + name)
      .attr("x", innerW + 4)
      .attr("dy", "0.32em")
      .style("fill", color)
      .text(text)
      .merge(sel)
      .attr("y", function (d) { return panel.y(accessor(d)); });
  }

  function updateAxes(panel) {
    var xAxis = d3.axisBottom(panel.x).ticks(8);
    var yAxis = d3.axisLeft(panel.y).ticks(5);
    panel.yAxisG.call(yAxis);
    if (panel.showX) {
      panel.xAxisG.call(xAxis).style("display", null);
    } else {
      panel.xAxisG.call(xAxis);
      panel.xAxisG.selectAll("text").style("display", "none");
    }
  }

  // ---- Set up DOM --------------------------------------------------------
  var chartsRoot = d3.select("#charts");
  var upliftPanel = makePanel(chartsRoot, {
    title: "Uplift",
    yLabel: "U (km/Myr)",
    color: COLORS.uplift,
    showX: false
  });
  var reliefPanel = makePanel(chartsRoot, {
    title: "Relief",
    yLabel: "R (km)",
    color: COLORS.relief,
    showX: false
  });
  var erosionPanel = makePanel(chartsRoot, {
    title: "Erosion  vs  uplift flux",
    yLabel: "rate (km/Myr)",
    color: COLORS.erosion,
    showX: false
  });
  var netPanel = makePanel(chartsRoot, {
    title: "Net rate  (uplift − erosion = dR/dt)",
    yLabel: "net (km/Myr)",
    color: COLORS.net,
    showX: true
  });
  // x-axis title under the last panel
  netPanel.g
    .append("text")
    .attr("class", "axis-label x-label")
    .attr("x", innerW / 2)
    .attr("y", panelH + 34)
    .attr("text-anchor", "middle")
    .text("time (Myr)");

  // ---- Controls ----------------------------------------------------------
  var controlSpecs = [
    { id: "peakUplift", label: "Peak uplift", min: 0.1, max: 10, step: 0.1, val: 3, unit: "km/Myr" },
    { id: "duration", label: "Event duration", min: 1, max: 40, step: 0.5, val: 12, unit: "Myr" },
    { id: "reliefLimit", label: "Relief limit @ peak U", min: 0.5, max: 8, step: 0.1, val: 3, unit: "km" },
    { id: "exponent", label: "Erosion exponent n", min: 1, max: 6, step: 0.5, val: 4, unit: "" },
    { id: "initialRelief", label: "Initial relief", min: 0, max: 4, step: 0.05, val: 0, unit: "km" }
  ];

  var UPLIFT_TYPES = ["bell", "plateau", "pulse", "arc", "triangle", "constant"];

  // ---- State, seeded from the URL query string ---------------------------
  // Every setting is encodable as ?peakUplift=..&duration=..&upliftType=..
  // so a configured view can be bookmarked or shared.
  var urlParams = new URLSearchParams(window.location.search);

  function seedNumeric(c) {
    var raw = urlParams.get(c.id);
    if (raw === null || raw === "") return c.val;
    var v = +raw;
    if (isNaN(v)) return c.val;
    return Math.min(c.max, Math.max(c.min, v)); // clamp to slider range
  }

  var state = {};
  controlSpecs.forEach(function (c) {
    c.init = seedNumeric(c); // resolved initial value (URL or default)
    state[c.id] = c.init;
  });
  var urlType = urlParams.get("upliftType");
  state.upliftType = UPLIFT_TYPES.indexOf(urlType) !== -1 ? urlType : "bell";

  // Write current state back into the URL without adding history entries.
  function syncURL() {
    var p = new URLSearchParams();
    controlSpecs.forEach(function (c) {
      p.set(c.id, state[c.id]);
    });
    p.set("upliftType", state.upliftType);
    var qs = "?" + p.toString();
    if (qs !== window.location.search) {
      // Some browsers block replaceState on file:// URLs; don't let that break the app.
      try {
        history.replaceState(null, "", window.location.pathname + qs);
      } catch (e) {
        /* ignore — URL syncing is a nicety, not required for the model */
      }
    }
  }

  var controlsRoot = d3.select("#controls");
  controlSpecs.forEach(function (c) {
    var row = controlsRoot.append("div").attr("class", "control");
    row
      .append("label")
      .attr("for", c.id)
      .html(c.label + ' <span class="val" id="' + c.id + '-val">' + c.init + (c.unit ? " " + c.unit : "") + "</span>");
    row
      .append("input")
      .attr("type", "range")
      .attr("id", c.id)
      .attr("min", c.min)
      .attr("max", c.max)
      .attr("step", c.step)
      .attr("value", c.init)
      .on("input", function () {
        state[c.id] = +this.value;
        d3.select("#" + c.id + "-val").text(this.value + (c.unit ? " " + c.unit : ""));
        render();
      });
  });

  // Uplift shape selector
  var shapeRow = controlsRoot.append("div").attr("class", "control");
  shapeRow.append("label").text("Uplift shape");
  var sel = shapeRow
    .append("select")
    .attr("id", "upliftType")
    .on("change", function () {
      state.upliftType = this.value;
      render();
    });
  [
    { v: "bell", t: "Bell curve" },
    { v: "plateau", t: "Ramp · hold · decay" },
    { v: "pulse", t: "Pulse (sharp rise, decay)" },
    { v: "arc", t: "Circular arc" },
    { v: "triangle", t: "Triangle" },
    { v: "constant", t: "Constant (square pulse)" }
  ].forEach(function (o) {
    sel
      .append("option")
      .attr("value", o.v)
      .attr("selected", o.v === state.upliftType ? "selected" : null)
      .text(o.t);
  });

  // ---- Render ------------------------------------------------------------
  function render() {
    var result = window.UpliftModel.run({
      duration: state.duration,
      peakUplift: state.peakUplift,
      exponent: state.exponent,
      reliefLimit: state.reliefLimit,
      initialRelief: state.initialRelief,
      upliftType: state.upliftType,
      steps: 2000
    });
    var data = result.series;

    var tMax = data[data.length - 1].t; // full simulated span (event + relaxation tail)
    var uMax = d3.max(data, function (d) { return d.uplift; }) || 1;
    var rMax = d3.max(data, function (d) { return Math.max(d.relief, d.reliefEq); }) || 1;
    var eMax = d3.max(data, function (d) { return Math.max(d.erosion, d.upliftFlux); }) || 1;
    var netMin = d3.min(data, function (d) { return d.net; });
    var netMax = d3.max(data, function (d) { return d.net; });
    var netAbs = Math.max(Math.abs(netMin), Math.abs(netMax), 1);

    [upliftPanel, reliefPanel, erosionPanel, netPanel].forEach(function (pn) {
      pn.x.domain([0, tMax]);
    });
    upliftPanel.y.domain([0, uMax * 1.1]);
    reliefPanel.y.domain([0, rMax * 1.1]);
    erosionPanel.y.domain([0, eMax * 1.1]);
    netPanel.y.domain([-netAbs * 1.1, netAbs * 1.1]);

    // Uplift panel
    drawLine(upliftPanel, "uplift", data, function (d) { return d.uplift; }, COLORS.uplift);

    // Relief panel: actual relief + equilibrium it chases
    drawLine(reliefPanel, "reliefEq", data, function (d) { return d.reliefEq; }, COLORS.eq, true);
    drawLine(reliefPanel, "relief", data, function (d) { return d.relief; }, COLORS.relief);
    var last = data[data.length - 1];
    drawEndLabel(reliefPanel, "relief", last, function (d) { return d.relief; }, COLORS.relief, "relief");
    drawEndLabel(reliefPanel, "reliefEq", last, function (d) { return d.reliefEq; }, COLORS.eq, "equilibrium");

    // Erosion panel: erosion + uplift flux (both m/Myr) to show flux balance
    drawLine(erosionPanel, "flux", data, function (d) { return d.upliftFlux; }, COLORS.flux, true);
    drawLine(erosionPanel, "erosion", data, function (d) { return d.erosion; }, COLORS.erosion);

    // Net panel: dR/dt, with a zero baseline (above = building, below = decaying)
    var zeroY = netPanel.y(0);
    var zline = netPanel.g.selectAll("line.zero-line").data([0]);
    zline
      .enter()
      .append("line")
      .attr("class", "zero-line")
      .attr("x1", 0)
      .attr("stroke", COLORS.zero)
      .attr("stroke-width", 1)
      .merge(zline)
      .attr("x2", innerW)
      .attr("y1", zeroY)
      .attr("y2", zeroY);
    drawLine(netPanel, "net", data, function (d) { return d.net; }, COLORS.net);

    [upliftPanel, reliefPanel, erosionPanel, netPanel].forEach(function (pn) {
      drawEventEnd(pn, state.duration);
      updateAxes(pn);
    });

    // Readout: peak relief reached vs the relief limit
    var peakR = d3.max(data, function (d) { return d.relief; });
    d3.select("#readout").html(
      "Peak relief reached: <b>" + peakR.toFixed(2) + " km</b>" +
      "  ·  relief limit @ peak U: <b>" + state.reliefLimit.toFixed(2) + " km</b>" +
      "  ·  erodibility k = " + result.k.toExponential(2)
    );

    syncURL();
  }

  render();
})();
