/*
 * app.js — UI, animation loop, cross-section plot, and URL state for the 2D
 * terrain sandbox. Depends on terrain-model.js, render.js, and d3 v7.
 */
(function () {
  "use strict";

  var N = 400, DX = 50; // 400 cells * 50 m = 20 km domain
  var model = window.TerrainModel.createModel({ n: N, dx: DX });

  var mapCanvas = document.getElementById("map");
  var overlay = document.getElementById("mapOverlay");

  // 10-step running average of dz for a less flickery Δelevation panel.
  var DZ_WINDOW = 10;
  var dzRing = [];
  for (var _i = 0; _i < DZ_WINDOW; _i++) dzRing.push(new Float32Array(N * N));
  var dzRingIdx = 0, dzCount = 0;
  var dzSum = new Float32Array(N * N);
  var dzAvg = new Float32Array(N * N);
  function accumulateDz() {
    var dz = model.state.dz, slot = dzRing[dzRingIdx], size = N * N;
    for (var k = 0; k < size; k++) { dzSum[k] += dz[k] - slot[k]; slot[k] = dz[k]; }
    dzRingIdx = (dzRingIdx + 1) % DZ_WINDOW;
    if (dzCount < DZ_WINDOW) dzCount++;
    var inv = 1 / dzCount;
    for (var q = 0; q < size; q++) dzAvg[q] = dzSum[q] * inv;
  }
  function clearDzAvg() {
    for (var i = 0; i < DZ_WINDOW; i++) dzRing[i].fill(0);
    dzRingIdx = 0; dzCount = 0; dzSum.fill(0); dzAvg.fill(0);
  }

  // Controls that round-trip through the URL query string.
  var PARAMS = ["topo", "seed", "relief", "noise", "Dexp", "nonlinear", "Sc", "fluvial", "Kexp", "m", "nexp", "depo", "mfd",
    "landslides", "lsMode", "lsStrength", "lsTestsExp", "lsC", "lsPhi", "lsRu", "lsRmaxExp", "lsRepose", "lsReposeSlow",
    "glaciers", "ela", "balGrad", "iceFlowExp", "iceErodeExp",
    "upliftExp", "upliftPattern", "upliftRamp", "dtExp", "shading", "vertExag", "showIce"];

  // Cross-sections at the three lines dividing the square into quarters.
  var ROWS = [Math.round(N * 0.25), Math.round(N * 0.5), Math.round(N * 0.75)];
  var ROW_LABELS = ["¼", "½", "¾"];
  var el = {};
  PARAMS.forEach(function (id) { el[id] = document.getElementById(id); });

  // ---- URL state ---------------------------------------------------------
  function readURL() {
    var p = new URLSearchParams(window.location.search);
    PARAMS.forEach(function (id) {
      if (!p.has(id)) return;
      var v = p.get(id);
      if (el[id].type === "checkbox") el[id].checked = v === "1" || v === "true";
      else el[id].value = v;
    });
  }
  function writeURL() {
    var p = new URLSearchParams();
    PARAMS.forEach(function (id) {
      p.set(id, el[id].type === "checkbox" ? (el[id].checked ? 1 : 0) : el[id].value);
    });
    var qs = "?" + p.toString();
    if (qs !== window.location.search) {
      try { history.replaceState(null, "", window.location.pathname + qs); } catch (e) {}
    }
  }

  // Reflect a control's value into its little readout span.
  function syncLabels() {
    setVal("relief", el.relief.value + " m");
    setVal("noise", el.noise.value + " m");
    setVal("Dexp", Math.pow(10, +el.Dexp.value).toFixed(3) + " m²/yr");
    setVal("Sc", (+el.Sc.value).toFixed(2));
    setVal("upliftExp", Math.pow(10, +el.upliftExp.value).toFixed(2) + " mm/yr");
    setVal("upliftRamp", (+el.upliftRamp.value === 0 ? "uniform" : el.upliftRamp.value + " cells"));
    var dtv = Math.pow(10, +el.dtExp.value);
    setVal("dtExp", dtv >= 1000 ? (dtv / 1000).toFixed(dtv >= 10000 ? 0 : 1) + " kyr/step" : Math.round(dtv) + " yr/step");
    setVal("vertExag", (+el.vertExag.value).toFixed(1) + "×");
    setVal("seed", el.seed.value.trim() === "" ? "clock" : el.seed.value);
    setVal("Kexp", Math.pow(10, +el.Kexp.value).toExponential(1));
    setVal("m", (+el.m.value).toFixed(2));
    setVal("nexp", (+el.nexp.value).toFixed(1));
    setVal("depo", (+el.depo.value).toFixed(2));
    setVal("mfd", el.mfd.value + (+el.mfd.value === 1 ? " (most spread)" : +el.mfd.value >= 8 ? " (steepest)" : ""));
    setVal("lsStrength", (+el.lsStrength.value).toFixed(2) + "×");
    setVal("lsTestsExp", Math.round(Math.pow(10, +el.lsTestsExp.value)) + "");
    setVal("lsC", el.lsC.value + " kPa");
    setVal("lsPhi", el.lsPhi.value + "°");
    setVal("lsRu", (+el.lsRu.value).toFixed(2));
    setVal("lsRmaxExp", Math.round(Math.pow(10, +el.lsRmaxExp.value)) + " m");
    setVal("lsRepose", el.lsRepose.value + "°");
    setVal("lsReposeSlow", el.lsReposeSlow.value + "°");
    setVal("ela", el.ela.value + " m");
    setVal("balGrad", (+el.balGrad.value).toFixed(3) + " /yr");
    setVal("iceFlowExp", "×" + Math.pow(10, +el.iceFlowExp.value).toFixed(2));
    setVal("iceErodeExp", Math.pow(10, +el.iceErodeExp.value).toExponential(1));
  }
  function setVal(id, text) {
    var s = document.getElementById(id + "-val");
    if (s) s.textContent = text;
  }

  // ---- Reading control state into model params --------------------------
  function stepParams() {
    return {
      dt: Math.pow(10, +el.dtExp.value),
      D: Math.pow(10, +el.Dexp.value),
      uplift: Math.pow(10, +el.upliftExp.value) * 0.001, // mm/yr -> m/yr (log slider)
      upliftPattern: el.upliftPattern.value,
      upliftRamp: +el.upliftRamp.value, // cells; 0 = uniform
      nonlinear: el.nonlinear.checked,
      Sc: +el.Sc.value,
      fluvial: el.fluvial.checked,
      K: Math.pow(10, +el.Kexp.value),
      m: +el.m.value,
      nexp: +el.nexp.value,
      depo: +el.depo.value,
      mfd: +el.mfd.value,
      landslides: el.landslides.checked,
      lsMode: el.lsMode.value,
      lsFos: 1 / (+el.lsStrength.value), // strength 1.0 = physical FoS<1; higher = fewer slides
      lsTests: Math.round(Math.pow(10, +el.lsTestsExp.value)),
      lsC: +el.lsC.value,
      lsPhi: +el.lsPhi.value,
      lsRu: +el.lsRu.value,
      lsRmax: Math.pow(10, +el.lsRmaxExp.value),
      lsRepose: +el.lsRepose.value,
      lsReposeSlow: +el.lsReposeSlow.value,
      glaciers: el.glaciers.checked,
      ela: +el.ela.value,
      balGrad: +el.balGrad.value,
      iceFlow: Math.pow(10, +el.iceFlowExp.value),
      iceErode: Math.pow(10, +el.iceErodeExp.value)
    };
  }
  function topoOpts() {
    var relief = +el.relief.value;
    // Blank seed -> use the clock so each (re)generation is a fresh terrain.
    var blank = el.seed.value.trim() === "";
    var seed = blank ? (Date.now() >>> 0) % 1000000 : (+el.seed.value | 0);
    return {
      type: el.topo.value,
      seed: seed,
      relief: relief,
      maxHeight: relief,
      noise: +el.noise.value,
      fromClock: blank
    };
  }

  // ---- Fixed color/plot range (set at generate, so flattening is visible) -
  var range = { min: 0, max: 1 };
  function recomputeRange() {
    var s = model.stats();
    var pad = (s.relief || 1) * 0.05;
    range.min = s.min - pad;
    range.max = s.max + pad;
  }

  // ---- Cross-section plot ------------------------------------------------
  function rowProfile(arr, r) {
    var out = [];
    for (var i = 0; i < N; i++) out.push({ x: (i * DX) / 1000, z: arr[r * N + i] });
    return out;
  }
  // Bed + ice-surface + lake profile for a row (z = bedrock, zi = ice surface,
  // zw = water surface, ice = thickness, lk = water depth).
  function rowProfileIce(r) {
    var h = model.state.h, ice = model.state.ice, lk = model.state.lake, out = [];
    for (var i = 0; i < N; i++) {
      var k = r * N + i;
      var lake = lk ? lk[k] : 0;
      out.push({ x: (i * DX) / 1000, z: h[k], zi: h[k] + ice[k], ice: ice[k], lk: lake, zw: h[k] + lake });
    }
    return out;
  }

  // Build one cross-section panel bound to a fixed row; returns {update}.
  function makeXSec(row, label, showX) {
    var wrap = d3.select("#xsections").append("div").attr("class", "xsection-wrap");
    wrap.append("div").attr("class", "xsection-title").text("Cross-section " + label + "  (row " + row + ")");
    var svg = wrap.append("svg");
    var m = { top: 6, right: 10, bottom: showX ? 26 : 8, left: 46 };
    var W = 360, H = showX ? 122 : 104;
    var iw = W - m.left - m.right, ih = H - m.top - m.bottom;
    svg.attr("width", W).attr("height", H).attr("viewBox", "0 0 " + W + " " + H);
    var g = svg.append("g").attr("transform", "translate(" + m.left + "," + m.top + ")");
    var xs = d3.scaleLinear().range([0, iw]);
    var ys = d3.scaleLinear().range([ih, 0]);
    var xAxisG = g.append("g").attr("class", "axis").attr("transform", "translate(0," + ih + ")");
    var yAxisG = g.append("g").attr("class", "axis");
    g.append("text").attr("class", "axis-label").attr("transform", "rotate(-90)")
      .attr("x", -ih / 2).attr("y", -36).attr("text-anchor", "middle").text("elev (m)");
    if (showX) {
      g.append("text").attr("class", "axis-label").attr("x", iw / 2).attr("y", ih + 22)
        .attr("text-anchor", "middle").text("distance (km)");
    }
    function join(cls, data, color, w, line) {
      var p = g.selectAll("path.l-" + cls).data([data]);
      p.enter().append("path").attr("class", "l-" + cls).attr("fill", "none")
        .attr("stroke", color).attr("stroke-width", w)
        .merge(p).attr("d", line);
    }
    function update() {
      xs.domain([0, (N * DX) / 1000]);
      var ip = rowProfileIce(row);
      var glOn = el.glaciers.checked;
      ys.domain([range.min, Math.max(range.max, iceYMax)]);
      var line = d3.line().x(function (d) { return xs(d.x); }).y(function (d) { return ys(d.z); });
      // ice body: filled area between bed and ice surface, plus an ice-surface line
      var iceArea = d3.area().defined(function (d) { return d.ice > 0.5; })
        .x(function (d) { return xs(d.x); })
        .y0(function (d) { return ys(d.z); })
        .y1(function (d) { return ys(d.zi); });
      var iceLine = d3.line().defined(function (d) { return d.ice > 0.5; })
        .x(function (d) { return xs(d.x); }).y(function (d) { return ys(d.zi); });
      var ap = g.selectAll("path.l-icefill").data([glOn ? ip : []]);
      ap.enter().append("path").attr("class", "l-icefill").attr("fill", "rgba(120,170,225,0.45)")
        .attr("stroke", "none").merge(ap).attr("d", iceArea);
      // lake body: filled water between bed and water surface
      var lakeArea = d3.area().defined(function (d) { return d.lk > 0.5; })
        .x(function (d) { return xs(d.x); })
        .y0(function (d) { return ys(d.z); })
        .y1(function (d) { return ys(d.zw); });
      var lp = g.selectAll("path.l-lakefill").data([ip]);
      lp.enter().append("path").attr("class", "l-lakefill").attr("fill", "rgba(40,110,200,0.55)")
        .attr("stroke", "none").merge(lp).attr("d", lakeArea);
      join("init", rowProfile(model.state.initial, row), "#bbb", 1, line);
      join("cur", rowProfile(model.state.h, row), "#8a5a2b", 1.6, line); // bedrock (brown)
      join("ice", glOn ? ip : [], "#3a78c2", 1.4, iceLine);              // ice surface (blue)
      xAxisG.call(d3.axisBottom(xs).ticks(6));
      if (!showX) xAxisG.selectAll("text").style("display", "none");
      yAxisG.call(d3.axisLeft(ys).ticks(4));
    }
    return { update: update };
  }

  var xSecs = ROWS.map(function (r, i) {
    return makeXSec(r, ROW_LABELS[i], i === ROWS.length - 1); // x-axis on bottom panel only
  });
  var iceYMax = 0; // shared ice-surface max across the cross-section rows (consistent y-scale)
  function updatePlots() {
    iceYMax = 0;
    if (el.glaciers.checked) {
      var h = model.state.h, ice = model.state.ice;
      ROWS.forEach(function (r) {
        for (var i = 0; i < N; i++) { var k = r * N + i; if (ice[k] > 0.5) { var z = h[k] + ice[k]; if (z > iceYMax) iceYMax = z; } }
      });
    }
    xSecs.forEach(function (s) { s.update(); });
  }

  // ---- Overlay: the three cross-section lines ---------------------------
  function drawOverlay() {
    if (overlay.width !== N) { overlay.width = N; overlay.height = N; }
    var ctx = overlay.getContext("2d");
    ctx.clearRect(0, 0, N, N);
    ctx.strokeStyle = "rgba(44,127,184,0.9)";
    ctx.lineWidth = 1.5;
    ROWS.forEach(function (r) {
      ctx.beginPath();
      ctx.moveTo(0, r + 0.5);
      ctx.lineTo(N, r + 0.5);
      ctx.stroke();
    });
  }

  // ---- Long-term history timeline ---------------------------------------
  // Records a few scalars each step and plots their trends over model time. Each
  // series is normalized to its own range (absolute current values in the legend).
  var HIST_MAX = 1200; // decimate beyond this to bound memory
  var hist = { t: [], relief: [], ice: [], inits: [] };
  var cumInits = 0;
  function clearHistory() { hist = { t: [], relief: [], ice: [], inits: [] }; cumInits = 0; }
  function recordHistory() {
    cumInits += model.state.slideInits || 0;
    var s = model.stats();
    var iceC = model.state.iceC, dxC = model.state.dxC, vol = 0;
    for (var k = 0; k < iceC.length; k++) vol += iceC[k];
    vol *= dxC * dxC; // m^3
    hist.t.push(model.state.time);
    hist.relief.push(s.relief);
    hist.ice.push(vol);
    hist.inits.push(cumInits);
    if (hist.t.length > HIST_MAX) {
      ["t", "relief", "ice", "inits"].forEach(function (key) {
        var a = hist[key], b = []; for (var i = 0; i < a.length; i += 2) b.push(a[i]); hist[key] = b;
      });
    }
  }
  var TL_SERIES = [
    { key: "relief", color: "#8a5a2b", label: "relief", unit: "m", fmt: function (v) { return Math.round(v) + " m"; } },
    { key: "ice", color: "#3a78c2", label: "ice vol", unit: "km³", fmt: function (v) { return (v / 1e9).toFixed(1) + " km³"; } },
    { key: "inits", color: "#c0392b", label: "slides", unit: "", fmt: function (v) { return Math.round(v) + ""; } }
  ];
  var tlInit = false, tlScales = null;
  function setupTimeline() {
    var svg = d3.select("#timeline");
    var W = 560, H = 132, m = { top: 8, right: 10, bottom: 22, left: 10 };
    svg.attr("viewBox", "0 0 " + W + " " + H).attr("width", W).attr("height", H);
    var g = svg.append("g").attr("transform", "translate(" + m.left + "," + m.top + ")");
    var iw = W - m.left - m.right, ih = H - m.top - m.bottom;
    var xs = d3.scaleLinear().range([0, iw]);
    var ys = d3.scaleLinear().domain([0, 1]).range([ih, 0]);
    var xAxisG = g.append("g").attr("class", "axis").attr("transform", "translate(0," + ih + ")");
    g.append("text").attr("class", "axis-label").attr("x", iw / 2).attr("y", ih + 20).attr("text-anchor", "middle").text("time (kyr)");
    tlScales = { g: g, xs: xs, ys: ys, iw: iw, ih: ih, xAxisG: xAxisG };
    tlInit = true;
  }
  function drawTimeline() {
    if (!tlInit) setupTimeline();
    var ts = tlScales, n = hist.t.length;
    if (n < 2) { ts.g.selectAll("path.tl").remove(); d3.select("#timeline-legend").text(""); return; }
    ts.xs.domain([hist.t[0] / 1000, hist.t[n - 1] / 1000]);
    var line = d3.line().x(function (d) { return ts.xs(d.x); }).y(function (d) { return ts.ys(d.y); });
    var legend = [];
    TL_SERIES.forEach(function (ser) {
      var arr = hist[ser.key], mx = 0;
      for (var i = 0; i < n; i++) if (arr[i] > mx) mx = arr[i];
      var norm = mx > 0 ? mx : 1, data = [];
      for (var q = 0; q < n; q++) data.push({ x: hist.t[q] / 1000, y: arr[q] / norm });
      var p = ts.g.selectAll("path.tl-" + ser.key).data([data]);
      p.enter().append("path").attr("class", "tl tl-" + ser.key).attr("fill", "none")
        .attr("stroke", ser.color).attr("stroke-width", 1.6).merge(p).attr("d", line);
      legend.push('<b style="color:' + ser.color + '">' + ser.label + "</b> " + ser.fmt(arr[n - 1]));
    });
    ts.xAxisG.call(d3.axisBottom(ts.xs).ticks(6));
    document.getElementById("timeline-legend").innerHTML = "&nbsp;&nbsp;" + legend.join(" &nbsp; ");
  }

  // ---- Render + readout --------------------------------------------------
  function draw() {
    recomputeRange(); // auto-rescale to current elevations (uplift raises the dome)
    window.TerrainRender.render(mapCanvas, model, {
      min: range.min, max: range.max, vertExag: +el.vertExag.value, mode: el.shading.value,
      showIce: el.showIce.checked
    });
    drawOverlay();
    window.TerrainRender.renderField(document.getElementById("dzMap"), dzAvg, N, { type: "diverging", scale: 2 });
    window.TerrainRender.renderField(document.getElementById("sedMap"), model.state.sedFlux, N, { type: "log" });
    // Fixed scales so the ice panels don't flicker as the field max wobbles frame-to-frame.
    window.TerrainRender.renderField(document.getElementById("iceMap"), model.state.ice, N, { type: "ice", scale: 800 });
    window.TerrainRender.renderField(document.getElementById("iceVelMap"), model.state.iceVel, N, { type: "ice", scale: 150 });
    updatePlots();
    drawTimeline();
    updateReadout();
  }
  var _ice = { pct: 0, maxH: 0, maxU: "0" };
  function updateReadout() {
    if (el.glaciers.checked) _ice = iceStats();
    var s = model.stats();
    document.getElementById("readout").innerHTML =
      "time <b>" + (model.state.time / 1000).toFixed(1) + " kyr</b>" +
      " · steps <b>" + model.state.steps + "</b>" +
      " · relief <b>" + Math.round(s.relief) + " m</b>" +
      " · mean <b>" + Math.round(s.mean) + " m</b>" +
      " · range <b>" + Math.round(s.min) + "–" + Math.round(s.max) + " m</b>" +
      (el.fluvial.checked ? " · fluvial iters <b>" + model.state.fluvIters + "</b>" : "") +
      (el.landslides.checked ? " · slides <b>" + model.state.slideCount + "</b> · slumps <b>" + model.state.slumpCount + "</b>" : "") +
      (el.glaciers.checked ? " · glaciated <b>" + _ice.pct + "%</b>" +
        " · max ice <b>" + _ice.maxH + " m</b>" +
        " · max flow <b>" + _ice.maxU + " m/yr</b>" +
        " · ice iters <b>" + model.state.iceSubs + "</b>" : "");
  }
  function iceStats() {
    var ice = model.state.ice, vel = model.state.iceVel, c = 0, mh = 0, mu = 0;
    for (var k = 0; k < ice.length; k++) {
      if (ice[k] > 2) c++;
      if (ice[k] > mh) mh = ice[k];
      if (vel[k] > mu) mu = vel[k];
    }
    return { pct: Math.round((100 * c) / ice.length), maxH: Math.round(mh), maxU: mu.toFixed(1) };
  }

  // ---- Generate / reset --------------------------------------------------
  function regenerate() {
    var opts = topoOpts();
    model.generate(opts);
    clearDzAvg();
    clearHistory();
    recomputeRange();
    // Show which seed was actually used (esp. the clock-derived one).
    setVal("seed", opts.fromClock ? "clock → " + opts.seed : String(opts.seed));
    draw();
    writeURL();
  }

  // ---- Animation loop ----------------------------------------------------
  var playing = false;
  function loop() {
    if (playing) {
      model.step(stepParams());
      accumulateDz();
      recordHistory();
      draw();
    }
    requestAnimationFrame(loop);
  }

  // ---- Wire up controls --------------------------------------------------
  function onControlChange() { syncLabels(); writeURL(); }
  PARAMS.forEach(function (id) {
    el[id].addEventListener("input", function () {
      onControlChange();
      // Live view-only changes redraw without stepping.
      if (id === "vertExag" || id === "shading" || id === "showIce") draw();
    });
  });

  document.getElementById("regen").addEventListener("click", regenerate);
  // Changing topo type / seed / feature-height / noise implies a regenerate.
  ["topo", "seed", "relief", "noise"].forEach(function (id) {
    el[id].addEventListener("change", regenerate);
  });

  var ppBtn = document.getElementById("playPause");
  ppBtn.addEventListener("click", function () {
    playing = !playing;
    ppBtn.textContent = playing ? "❚❚ Pause" : "▶ Play";
  });
  document.getElementById("stepBtn").addEventListener("click", function () {
    model.step(stepParams());
    accumulateDz();
    recordHistory();
    draw();
  });
  document.getElementById("resetBtn").addEventListener("click", function () {
    model.reset();
    clearDzAvg();
    clearHistory();
    draw();
  });

  // ---- Init --------------------------------------------------------------
  readURL();
  syncLabels();
  regenerate();
  requestAnimationFrame(loop);
})();
