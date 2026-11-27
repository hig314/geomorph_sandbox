/*
 * relief_limit/glacier_app.js — the stand-alone 16:9 slide "Glacier buttressing": relief and
 * threshold over the last 200 kyr, with ice thickness in a thin panel below. Depends on
 * relief_limit/glacier_model.js, core/noise.js, ui/*.js (window.GS) and d3 v7.
 *
 * Same conventions as app.js: one SVG with a fixed 1920 × 1080 viewBox styled by attributes
 * only (Save SVG is self-contained); the model is solved once per parameter change
 * (GS.glacier.simulate, exact) and the animation moves a time cursor tc across it, revealing
 * relief, threshold and ice together. Two sliders (weakening, erosion rate) sit in the frame.
 */
(function () {
  "use strict";
  var GS = window.GS, GL = GS.glacier;

  // The frame of app.js (X0, X1, Y1, labels, line weights, slider row) with the relief axis
  // shortened to make room for the ice panel between it and the slider row.
  var W = 1920, H = 1080, X0 = 130, X1 = 1590, Y0 = 700, Y1 = 90, IY0 = 900, IY1 = 760, HMAX = 1.15;
  var THEMES = {
    light: { bg: "#ffffff", ink: "#1a1a1a", axis: "#9a9a9a", text: "#555555", limit: "#7a7a7a", ice: "#3f7cb0", iceFill: "#d5e5f2" },
    dark:  { bg: "#111111", ink: "#f2f2f2", axis: "#777777", text: "#b8b8b8", limit: "#a0a0a0", ice: "#86b8e0", iceFill: "#1f3547" }
  };
  var WT = 4, FONT = "Helvetica Neue, Helvetica, Arial, sans-serif";
  var W_REF = GL.DEFAULTS.W; // weakening at "× 1"
  // glacial stages, labelled under the ice panel's time axis: [age at the label, text, anchor]
  var PHASES = [[155, "glacial", "middle"], [123, "interglacial", "middle"], [62, "glacial", "middle"], [0, "interglacial", "end"]];

  // ---- Methodology text (same helper as app.js) -----------------------------
  function M(title, eqs, text, src, tag) {
    var h = "<h4>" + (tag ? '<span class="tag ' + tag + '">' + tag + "</span>" : "") + title + "</h4>";
    (eqs || []).forEach(function (e) { h += '<span class="eq">' + e + "</span>"; });
    if (text) h += "<div>" + text + "</div>";
    if (src) h += '<div class="src">' + src + "</div>";
    return h;
  }
  var SRC = "relief_limit/glacier_model.js · py/test_glacier_buttress.py";
  var m = {
    weak: M("Weakening", ["stress change since the last fracture = change in buttressing:  S = ∫ |d(B·b)/dt| dt", "fracture when S reaches q,  q = 0.025 × (1 ± 0.25, uniform):  D ← D + W·q", "W = (slider) × " + W_REF + ",   T = T<sub>0</sub> − D + B·b(H)"],
      "How much the rock beside the glacier is weakened by changes in the support the ice gives it. Every change in buttressing, up or down, changes the stress in the rock; the change accumulates until the rock fractures, after a nearly fixed amount, and the threshold then steps down by W times that amount. So the steps are small and regular, and they follow the buttressing: frequent where the ice is thin and its support changes fast (the terminations, the interglacials, the final plunge), rare under thick ice, whose fluctuations hardly change the support. The damage stays until a landslide removes the damaged rock, after which the threshold is back to its buttressed value. At × 1 the threshold loses as much to damage as the buttressing has changed since the last landslide. At zero the threshold is set by buttressing alone and failures come only where erosion reaches it or the ice thins under the relief; raise it and failures come earlier, at lower relief, and more of them are set off by a fracture.", SRC, "asserted"),
    erosion: M("Erosion rate", ["dR/dt = E = (slider) × " + GL.DEFAULTS.Eref + " T<sub>0</sub> per kyr"],
      "Glacial erosion deepening the valley and so raising the relief of its walls, at a constant rate wherever there is ice: the whole record for the thick glacier, and not while the retreating glacier is gone (relief then stands still). At the default, relief would grow by 2.4 times the ice-free threshold over the 200 kyr if nothing failed. Faster erosion brings relief up to the threshold sooner after each landslide, so failures are more frequent and relief stays close to the threshold; at zero relief only falls, and only the thinning ice or fracturing can fail the slope.", SRC, "asserted"),
    drawTime: M("Drawing time", ["t<sub>cursor</sub> advances 200 kyr / (drawing time) per second"],
      "How many seconds the record takes to draw across the slide. It changes nothing in the model, only the pace of the animation.", null),
    theme: M("Slide background", [], "Light or dark, to match the deck. Colours only; the saved SVG uses the theme shown.", null),
    frame: M("Sliders in the frame", [], "The two sliders, weakening and erosion rate, sit inside the 16:9 frame under the ice panel, where they stay visible in full screen. They are the same quantities as the controls of the same name below and move together with them. The saved SVG does not include them.", null),
    ref: M("Threshold without damage", ["T<sub>intact</sub> = T<sub>0</sub> + B·b(H),   b = (1 − e<sup>−H/H<sub>s</sub></sup>) / (1 − e<sup>−1/H<sub>s</sub></sup>),   H<sub>s</sub> = 0.3"],
      "A faded, wide line behind the threshold: what the threshold would be with buttressing alone and no fracturing. The gap between it and the dashed threshold is the damage standing in the rock at that moment; it closes at every landslide.", SRC, null),
    R0: M("Initial relief R<sub>0</sub>", ["R(200 ka) = R<sub>0</sub>"],
      "Relief at the start of the record, 200 kyr ago, as a multiple of the ice-free threshold T<sub>0</sub>. The default, 1.3, is above the mean relief of the record and above what the rock could hold without ice: the slope starts out held up by the glacier. With weakening at zero and a slow erosion rate (about 0.25) relief then builds slowly under the buttressed threshold through the long glacials and fails when the ice thins. A low value gives a long ramp of growing relief before the first failure.", SRC, "asserted"),
    ice: M("Glacier history", ["thick:  H ≥ 0.05 throughout;  interglacials keep ~40 % of the glacial maximum; plunge in the last 600 years", "retreating:  H = 0 from 128 to 118 ka and from 10 to 4.5 ka;  re-advance to a maximum at 0.3 ka;  H = 0 today"],
      "Which ice history the slope lives through. Both share the two glacials and the glacier's own few-thousand-year fluctuations. <b>Thick</b>: the glacier never leaves; it thins in the interglacials and plunges with climate change at the end. <b>Retreating</b>: the glacier retreats completely early in the last interglacial and early in the Holocene, re-advances in the late Holocene (a Little Ice Age maximum about 300 years ago) and then retreats completely. With no ice there is no buttressing (the threshold is the ice-free one, less any damage), no glacial erosion (relief stands still) and no change in buttressing to fracture the rock.", SRC, "asserted"),
    seed: M("Random sequence", ["mulberry32(seed): ice fluctuations; seed + 2000003: fracture quotas (± 25 %); seed + 1000003: landslide sizes"],
      "Which realisation is drawn: the glacier's internal fluctuations, the small spread in the buttressing change each fracture takes, and the landslide sizes, each from its own stream. The same seed always gives the same figure. The glacial–interglacial envelope is fixed.", "core/noise.js", null)
  };

  var specs = [
    { id: "weak", label: "Weakening", type: "range", val: 1, min: 0, max: 3, step: 0.05, fmt: function (v) { return "× " + v.toFixed(2); }, group: "Glacier", method: m.weak },
    { id: "erosion", label: "Erosion rate", type: "range", val: 1, min: 0, max: 3, step: 0.05, fmt: function (v) { return "× " + v.toFixed(2); }, group: "Glacier", method: m.erosion },
    { id: "R0", label: "Initial relief", type: "range", val: 1.3, min: 0.2, max: 1.5, step: 0.05, fmt: function (v) { return v.toFixed(2) + " T₀"; }, group: "Glacier", method: m.R0 },
    { id: "ice", label: "Glacier history", type: "select", val: "thick", group: "Glacier", method: m.ice,
      options: [{ v: "thick", t: "thick" }, { v: "retreat", t: "retreating" }] },
    { id: "seed", label: "Random sequence", type: "range", val: 7, min: 1, max: 99, step: 1, group: "Glacier", method: m.seed },
    { id: "drawTime", label: "Drawing time", type: "range", val: 15, min: 1, max: 60, step: 1, unit: "s", group: "Figure", method: m.drawTime },
    { id: "theme", label: "Background", type: "select", val: "light", group: "Figure", method: m.theme,
      options: [{ v: "light", t: "light" }, { v: "dark", t: "dark" }] },
    { id: "frame", label: "Sliders in frame", type: "select", val: "on", group: "Figure", method: m.frame,
      options: [{ v: "on", t: "shown" }, { v: "off", t: "hidden" }] },
    { id: "ref", label: "Threshold without damage", type: "select", val: "on", group: "Figure", method: m.ref,
      options: [{ v: "on", t: "shown" }, { v: "off", t: "hidden" }] }
  ];
  var url = GS.ui.urlState(GS.ui.toUrlSpecs(specs));
  var state = url.read();
  // Keep only what differs from the defaults in the address bar, so a reload picks up a changed
  // default instead of the value the page wrote there itself.
  var urlWriteAll = url.write;
  url.write = function (st) {
    urlWriteAll(st);
    var q = new URLSearchParams(window.location.search);
    specs.forEach(function (c) { if (st[c.id] === c.val) q["delete"](c.id); });
    var qs = q.toString();
    try { window.history.replaceState(null, "", window.location.pathname + (qs ? "?" + qs : "")); } catch (e) { /* nicety only */ }
  };

  // Embed mode (?embed=1), as app.js: only the 16:9 frame; the host sends
  // {type: "relief-step" | "deck:step", step, replay, dur} and gets {type: "deck-key", key} back.
  var EMBED = /[?&]embed=1\b/.test(window.location.search);
  if (EMBED) { document.body.classList.add("embed"); url.write = function () {}; }

  // ---- Figure ---------------------------------------------------------------
  var slide = d3.select("#stage").append("div").attr("id", "slide");
  var svg = slide.append("svg").attr("xmlns", "http://www.w3.org/2000/svg")
    .attr("viewBox", "0 0 " + W + " " + H).attr("font-family", FONT);
  var bg = svg.append("rect").attr("width", W).attr("height", H);
  var clip = svg.append("defs").append("clipPath").attr("id", "gb-clip").append("rect").attr("x", 0).attr("y", 0).attr("height", H);
  var iceArea = svg.append("path").attr("clip-path", "url(#gb-clip)");
  var iceLine = svg.append("path").attr("fill", "none").attr("stroke-width", WT).attr("stroke-linejoin", "round").attr("clip-path", "url(#gb-clip)");
  var axis = svg.append("path").attr("fill", "none").attr("stroke-width", WT)
    .attr("d", "M" + X0 + "," + Y1 + "V" + Y0 + "H" + X1 + "M" + X0 + "," + IY1 + "V" + IY0 + "H" + X1);
  var yLab = svg.append("text").attr("x", X0).attr("y", Y1 - 28).attr("font-size", 38).text("Relief");
  var iceLab = svg.append("text").attr("x", X1 + 18).attr("font-size", 34);
  iceLab.append("tspan").attr("x", X1 + 18).attr("y", (IY0 + IY1) / 2 - 8).text("Ice");
  iceLab.append("tspan").attr("x", X1 + 18).attr("dy", 40).text("thickness");
  var tLeft = svg.append("text").attr("x", X0).attr("y", IY0 + 46).attr("font-size", 34).text("200 kyr ago");
  var xLab = svg.append("text").attr("x", X1).attr("y", 1018).attr("text-anchor", "end").attr("font-size", 38).text("Time"); // where app.js has it
  var phaseG = svg.append("g").attr("font-size", 34);
  var refLine = svg.append("path").attr("fill", "none").attr("stroke-width", 2 * WT).attr("stroke-opacity", 0.22).attr("stroke-linejoin", "round").attr("clip-path", "url(#gb-clip)");
  var refLab = svg.append("text").attr("x", X1 + 18).attr("font-size", 34).attr("dy", "0.34em").attr("fill-opacity", 0.6).text("Undamaged");
  // the threshold, drawn as app.js draws the step-7 limit: up to the cursor and level beyond it
  // (as it stands now), dashed; its steps down at fractures drawn solid on top
  var thrLine = svg.append("path").attr("fill", "none").attr("stroke-width", WT).attr("stroke-dasharray", "16 12");
  var fracG = svg.append("g");
  var thrLab = svg.append("text").attr("x", X1 + 18).attr("font-size", 34).attr("dy", "0.34em").text("Threshold");
  var relief = svg.append("path").attr("fill", "none").attr("stroke-width", WT).attr("stroke-linejoin", "round").attr("clip-path", "url(#gb-clip)");
  function hit(node) {
    return svg.append(node).attr("fill", "none").attr("stroke", "transparent").attr("stroke-width", 28)
      .attr("clip-path", "url(#gb-clip)").style("pointer-events", "stroke").style("cursor", "help");
  }
  var hitIce = hit("path"), hitRef = hit("path"), hitThr = hit("path"), hitRelief = hit("path");

  var sim = null, pal = THEMES.light, tc = 0, playing = false, lastNow = 0, yMax = 1.65, pace = 0;
  function x(t) { return X0 + (X1 - X0) * t / sim.p.tmax; }
  function y(r) { return Y0 - (Y0 - Y1) * r / yMax; }
  function yi(h) { return IY0 - (IY0 - IY1) * h / HMAX; }
  function poly(pts, fy) {
    var s = "", i;
    for (i = 0; i < pts.length; i++) s += (i ? "L" : "M") + x(pts[i][0]).toFixed(2) + "," + fy(pts[i][1]).toFixed(2);
    return s;
  }
  function bAt(t) { // the buttressing at t: node values, linear between (as glacier_model.js)
    var i = Math.max(0, Math.min(Math.floor(t / sim.p.dtH), sim.ice.length - 2)), a = sim.ice[i], b = sim.ice[i + 1];
    var ba = sim.p.B * GL.butt(a[1], sim.p), bb = sim.p.B * GL.butt(b[1], sim.p);
    return ba + (bb - ba) * (t - a[0]) / (b[0] - a[0]);
  }
  function modelParams() {
    return { W: state.weak * W_REF, E: state.erosion, R0: state.R0, seed: state.seed, ice: state.ice };
  }

  function rebuild() {
    sim = GL.simulate(modelParams());
    var p = sim.p, ice = sim.ice;
    pal = THEMES[state.theme] || THEMES.light;
    yMax = 1.1 * (p.T0 + p.B);
    d3.select("#stage").style("background", pal.bg);
    bg.attr("fill", pal.bg);
    axis.attr("stroke", pal.axis);
    [yLab, iceLab, tLeft, xLab, thrLab].forEach(function (s) { s.attr("fill", pal.text); });
    refLab.attr("fill", pal.text);

    var dIce = poly(ice, yi);
    iceLine.attr("stroke", pal.ice).attr("d", dIce); hitIce.attr("d", dIce);
    iceArea.attr("fill", pal.iceFill).attr("d", dIce + "L" + x(p.tmax) + "," + IY0 + "L" + x(0) + "," + IY0 + "Z");
    phaseG.selectAll("text").remove();
    PHASES.forEach(function (ph) {
      phaseG.append("text").attr("x", ph[2] === "end" ? X1 : x(p.tmax - ph[0])).attr("y", IY0 + 46)
        .attr("text-anchor", ph[2]).attr("fill", pal.text).text(ph[1]);
    });

    var dRef = poly(ice.map(function (v) { return [v[0], p.T0 + p.B * GL.butt(v[1], p)]; }), y); // node values, linear between, as the model
    refLine.attr("stroke", pal.limit).attr("d", dRef).style("display", state.ref === "off" ? "none" : null);
    hitRef.attr("d", dRef).style("display", state.ref === "off" ? "none" : null);
    thrLine.attr("stroke", pal.limit);
    hitThr.attr("d", poly(sim.thr, y));
    fracG.selectAll("line").remove();
    fracG.selectAll("line").data(sim.fracs).enter().append("line")
      .attr("x1", function (d) { return x(d.t); }).attr("x2", function (d) { return x(d.t); })
      .attr("y1", function (d) { return y(d.from); }).attr("y2", function (d) { return y(d.to); })
      .attr("stroke", pal.ink).attr("stroke-width", 2.5 * WT).attr("stroke-linecap", "butt");
    var dRel = poly(sim.path, y);
    relief.attr("stroke", pal.ink).attr("d", dRel); hitRelief.attr("d", dRel);
  }

  // Everything that follows the time cursor.
  function draw() {
    var p = sim.p, T = GL.valueAt(sim.thr, tc), Ti = p.T0 + bAt(tc);
    clip.attr("width", tc >= p.tmax ? W : x(tc));
    var pts = [], i;
    for (i = 0; i < sim.thr.length && sim.thr[i][0] <= tc; i++) pts.push(sim.thr[i]);
    pts.push([tc, T]); pts.push([p.tmax, T]);
    thrLine.attr("d", poly(pts, y));
    fracG.selectAll("line").style("display", function (d) { return d.t <= tc ? null : "none"; });
    thrLab.attr("y", y(T));
    refLab.attr("y", y(Ti)).style("display", state.ref === "off" || Math.abs(y(Ti) - y(T)) < 40 ? "none" : null);
  }

  // ---- In-frame sliders (as app.js) ------------------------------------------
  var frame = slide.append("div").attr("id", "frame-controls");
  var frameCtl = [
    { id: "weak", label: "Weakening", fmt: function (v) { return "× " + v.toFixed(2); } },
    { id: "erosion", label: "Erosion rate", fmt: function (v) { return "× " + v.toFixed(2); } }
  ];
  frameCtl.forEach(function (fc) {
    var c = specs.filter(function (s) { return s.id === fc.id; })[0];
    var box = frame.append("label").attr("class", "fc");
    var head = box.append("span").attr("class", "fc-head");
    head.append("span").text(fc.label);
    fc.val = head.append("span").attr("class", "fc-val");
    fc.inp = box.append("input").attr("type", "range").attr("min", c.min).attr("max", c.max).attr("step", c.step)
      .on("input", function () { state[fc.id] = +this.value; controls.refresh(); onChange(fc.id); })
      .on("change", function () { this.blur(); }); // hand the keys back to the slide
  });
  function syncFrame() {
    frame.style("display", state.frame === "off" ? "none" : null)
      .style("--fc-text", pal.text).style("--fc-line", pal.axis).style("--fc-ink", pal.ink);
    frameCtl.forEach(function (fc) { fc.inp.property("value", state[fc.id]); fc.val.text(fc.fmt(state[fc.id])); });
  }

  // ---- Curve tooltips, readout ---------------------------------------------
  function kyr(t) { var a = sim.p.tmax - t; return a < 0.05 ? "today" : a.toFixed(1) + " ka"; }
  GS.ui.tip.attach(hitRelief, function () {
    return GS.ui.curveTip("Relief", pal.ink, "solid",
      "Grows at a constant rate by glacial erosion wherever there is ice (it stands still while the glacier is gone) and drops at each landslide, to d below the threshold in force (d from a truncated power law, as in the relief figure). " + sim.fails.length + " landslides in this record.",
      "dR/dt = E with ice, 0 without;  failure where R = T:  R ← max(T − d, floor)");
  });
  GS.ui.tip.attach(hitThr, function () {
    return GS.ui.curveTip("Threshold", pal.limit, "dashed",
      "The relief the valley wall can hold. The ice buttresses it, on a saturating curve: the response is stretched where the ice is thin and compressed where it is thick, so the threshold hardly moves through the glacials and drops fast in the terminations and when the ice is nearly gone. Each change in buttressing stresses the rock, which fractures after a nearly fixed amount of change and steps the threshold down (the black marks). A landslide removes the damaged rock and the threshold returns to its buttressed value. Right of the cursor it is drawn level, as it stands now.",
      "T = T<sub>0</sub> + B·b(H) − D,   B = " + sim.p.B + ",   H<sub>s</sub> = " + sim.p.Hs + ";   fracture: D ← D + W·q;   landslide: D ← 0");
  });
  GS.ui.tip.attach(hitRef, function () {
    return GS.ui.curveTip("Threshold without damage", pal.limit, "solid", "Buttressing alone: what the threshold would be without fracturing. The gap below it is the damage standing in the rock.", "T<sub>0</sub> + B·b(H)");
  });
  GS.ui.tip.attach(hitIce, function () {
    return GS.ui.curveTip("Ice thickness", pal.ice, "solid",
      "Two glacial cycles: the penultimate glacial, the last interglacial, the last glacial to its maximum about 24 ka, and the Holocene, on a fixed schematic envelope; on top of it, the glacier thickens and thins over a few thousand years by its own dynamics (random, seeded). " +
      (sim.p.ice === "retreat" ? "A retreating glacier: it is gone early in the last interglacial (128–118 ka) and early in the Holocene (10–4.5 ka), re-advances to a Little Ice Age maximum about 300 years ago, and then retreats completely." :
        "A thick glacier: it keeps about 40 % of its glacial maximum through the interglacials, until climate change takes it, in the last 600 years, lower than in any interglacial.") + " In units of the largest glacial thickness.",
      sim.p.ice === "retreat" ? "H = max(0, envelope + σ·(0.4 + 0.6·envelope)·min(1, envelope/0.15)·ξ(t))" : "H = max(H<sub>min</sub>, envelope + σ·(0.4 + 0.6·envelope)·ξ(t))");
  });

  function readout() {
    var n = sim.fails.length, nf = 0, nd = 0;
    sim.fails.forEach(function (f) {
      if (f.cause === "fracture") nf++;
      else if (GL.valueAt(sim.ice, Math.min(f.t + 0.25, sim.p.tmax)) < GL.valueAt(sim.ice, Math.max(f.t - 0.25, 0))) nd++;
    });
    var h = "<b>Glacier buttressing" + (sim.p.ice === "retreat" ? ", a retreating glacier" : "") + ".</b> The ice raises the threshold, most sensitively where it is thin (saturating); each change in buttressing fractures the rock and steps the threshold down (black marks); erosion raises relief steadily.<br>";
    h += "Landslides: <b>" + n + "</b> (" + nf + " set off by a fracture, " + nd + " as the ice thinned, " + (n - nf - nd) + " reached by erosion)  ·  fractures: <b>" + sim.fracs.length + "</b>";
    if (n) h += "  ·  first landslide " + kyr(sim.fails[0].t);
    d3.select("#readout").html(h);
  }

  // ---- Animation and transport ----------------------------------------------
  function tick(now) {
    if (!playing) return;
    // the first frame's timestamp can precede the performance.now() taken in play(): never step back
    var dt = Math.max(0, Math.min((now - lastNow) / 1000, 0.1)); lastNow = Math.max(lastNow, now);
    tc = Math.min(sim.p.tmax, tc + dt * sim.p.tmax / (pace || state.drawTime));
    draw();
    if (tc >= sim.p.tmax) { playing = false; buttons(); } else requestAnimationFrame(tick);
  }
  function play() {
    if (tc >= sim.p.tmax) tc = 0;
    if (!playing) { playing = true; lastNow = performance.now(); requestAnimationFrame(tick); }
    buttons();
  }
  function pause() { playing = false; buttons(); }
  function replay() { tc = 0; draw(); playing = false; play(); }
  function fullscreen() {
    var el = document.getElementById("stage");
    if (document.fullscreenElement) document.exitFullscreen(); else if (el.requestFullscreen) el.requestFullscreen();
  }
  function saveSvg() {
    var src = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg.node());
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([src], { type: "image/svg+xml" }));
    a.download = "glacier_buttressing.svg";
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 1000);
  }

  var tp = d3.select("#transport"), playBtn;
  playBtn = tp.append("button").on("click", function () { if (playing) pause(); else play(); this.blur(); });
  tp.append("button").text("Replay").on("click", function () { replay(); this.blur(); });
  tp.append("button").text("Full screen").on("click", function () { fullscreen(); this.blur(); });
  tp.append("button").text("Save SVG").attr("title", "The frame as shown, 1920 × 1080, self-contained").on("click", function () { saveSvg(); this.blur(); });
  tp.append("span").attr("class", "keys").html("<kbd>space</kbd> play / pause  ·  <kbd>R</kbd> replay  ·  <kbd>F</kbd> full screen");
  function buttons() { playBtn.text(playing ? "Pause" : "Play"); }

  d3.select(window).on("keydown", function (ev) {
    var tag = ev.target && ev.target.tagName;
    if (tag === "INPUT" && ev.target.type !== "range" || tag === "SELECT" || tag === "TEXTAREA" || ev.metaKey || ev.ctrlKey || ev.altKey) return;
    var k = ev.key;
    if (EMBED && (k === "ArrowRight" || k === "ArrowLeft" || k === "PageDown" || k === "PageUp" || k === " " || k === "e" || k === "E" || k === "m" || k === "M")) {
      window.parent.postMessage({ type: "deck-key", key: k }, "*");
      ev.preventDefault();
      return;
    }
    if (tag === "INPUT") return;
    if (k === " ") { if (playing) pause(); else play(); }
    else if (k === "r" || k === "R") replay();
    else if (k === "f" || k === "F") fullscreen();
    else return;
    ev.preventDefault();
  });

  // ---- Controls -------------------------------------------------------------
  function onChange(id) {
    rebuild();
    if (id !== "theme" && id !== "drawTime" && id !== "frame" && id !== "ref") { playing = false; tc = sim.p.tmax; } // tuning: show the whole record
    draw(); syncFrame(); readout(); buttons(); url.write(state);
  }
  var controls = GS.ui.buildControls(d3.select("#controls"), specs, state, onChange);

  rebuild(); draw(); syncFrame(); readout(); buttons(); url.write(state);
  play();

  // A deck drives the slide: arriving going forwards replays it (in dur seconds if given),
  // arriving from a later slide shows the whole record.
  window.addEventListener("message", function (ev) {
    var d = ev.data;
    if (!EMBED || !d || (d.type !== "relief-step" && d.type !== "deck:step")) return;
    pace = +d.dur > 0 ? +d.dur : 0;
    if (d.replay) replay(); else { playing = false; tc = sim.p.tmax; draw(); buttons(); }
  });
})();
