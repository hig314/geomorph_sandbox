/*
 * relief_limit/app.js — the animated 16:9 figure, transport and controls for "Limits to
 * relief". Depends on relief_limit/model.js, core/noise.js, ui/*.js (window.GS) and d3 v7.
 *
 * The figure is one SVG with a fixed 1920 × 1080 viewBox, styled by attributes only so the
 * "Save SVG" frame is self-contained. The model is solved once per parameter change
 * (GS.relief.simulate, exact); the animation only moves a time cursor tc across it: a clip
 * rectangle reveals the relief line, earthquakes appear as the cursor passes them, and the
 * static limit is drawn stepped up to tc and level beyond it (the limit as it stands now).
 * Moving between stages that share a relief history (3 → 4 → 5) keeps the cursor and fades
 * the new elements in; any other change of stage redraws from t = 0.
 */
(function () {
  "use strict";
  var GS = window.GS, RL = GS.relief;

  var W = 1920, H = 1080, X0 = 130, X1 = 1590, Y0 = 960, Y1 = 90;
  var THEMES = {
    light: { bg: "#ffffff", ink: "#1a1a1a", axis: "#9a9a9a", text: "#555555", limit: "#7a7a7a", quake: "#c0392b" },
    dark:  { bg: "#111111", ink: "#f2f2f2", axis: "#777777", text: "#b8b8b8", limit: "#a0a0a0", quake: "#ff7b66" }
  };
  // Line weights. Two weights a factor of two apart carry the figure: WT for the axes, relief,
  // limits and circle outlines; WT / 2 for the momentary dips of the limit at earthquakes.
  // A third, 2 × WT, is reserved for faded reference lines drawn behind the others (not used yet).
  var WT = 4;
  var FONT = "Helvetica Neue, Helvetica, Arial, sans-serif";
  var STAGES = [
    ["A limit to relief", "Relief grows until it meets a limit, and stays there."],
    ["Small landslides", "At the limit a landslide removes a small increment of relief, which regrows and fails again."],
    ["Large landslides", "The same, with larger failures: bigger saw-teeth, lower mean relief."],
    ["A static limit above", "Slopes fail below their static limit because something perturbs them: earthquakes of constant size, each arriving just as relief reaches the limit."],
    ["Earthquakes lower the limit", "Each earthquake pulls the limit down from the static limit for a moment; the dot is the limit at that moment."],
    ["A spectrum of earthquakes", "Earthquakes of many sizes and intervals. One fails the slope only if it lowers the limit to the relief."],
    ["Earthquakes damage the slope", "An earthquake that does not fail the slope weakens it and the static limit steps down; the next landslide removes the damaged rock and restores it."]
  ];

  // ---- Methodology text (same helper as d1/app.js) --------------------------
  function M(title, eqs, text, src, tag) {
    var h = "<h4>" + (tag ? '<span class="tag ' + tag + '">' + tag + "</span>" : "") + title + "</h4>";
    (eqs || []).forEach(function (e) { h += '<span class="eq">' + e + "</span>"; });
    if (text) h += "<div>" + text + "</div>";
    if (src) h += '<div class="src">' + src + "</div>";
    return h;
  }
  var SRC = "relief_limit/model.js · py/test_relief_limit.py";
  var m = {
    stage: M("Step of the argument", ["dR/dt = g", "failure if R ≥ L<sub>s</sub>(1 − f),  then R ← L<sub>s</sub>(1 − f) − d"],
      "Seven figures that add one idea at a time: a flat limit (1), landslides small (2) and large (3), a static limit above with constant earthquakes (4), the earthquakes drawn (5), a spectrum of earthquakes (6), and damage (7). One rule runs all of them; each step only switches on another term. Steps 3, 4 and 5 are the same relief history seen three ways, so moving between them does not redraw the line.", SRC),
    drawEarly: M("Drawing time, steps 1–5", ["t<sub>cursor</sub> advances t<sub>max</sub> / (drawing time) per second"],
      "How many seconds the line takes to draw across the slide in the early steps. It changes nothing in the model, only the pace of the animation. Steps 1, 2 and 3 each draw a new line; steps 4 and 5 keep the line of step 3, so this only applies to them when one is opened directly.", null),
    drawTime: M("Drawing time, steps 6–7", ["t<sub>cursor</sub> advances t<sub>max</sub> / (drawing time) per second"],
      "How many seconds the line takes to draw across the slide in the last two steps. It changes nothing in the model, only the pace of the animation: slow it down to talk over it, speed it up to get to the point.", null),
    frame: M("Sliders in the frame", [],
      "Two sliders, weakening and uplift rate, can sit inside the 16:9 frame under the time axis, where they stay visible in full screen. By default they appear on the last step only. They are the same quantities as the controls of the same name below and move together with them; every other control stays out here as presenter's notes. The saved SVG does not include them.", null),
    g: M("Uplift rate g", ["dR/dt = g", "saw-tooth period = d / g"],
      "How fast relief grows between failures, as a fraction of the reference rate g<sub>0</sub> that sets the time unit (L/g<sub>0</sub>); the default is the maximum. At zero relief never grows: nothing regenerates what a landslide removes, and the only way to a failure is an earthquake large enough to reach the standing relief, or in the last step enough accumulated damage to bring the limit down to it. It stands for rock uplift and the valley incision that keeps pace with it. Faster uplift rebuilds relief sooner after each landslide, so the slope spends more of its time high: landslides are more frequent, mean relief is higher, and in the last step less damage accumulates before the next failure. Slower uplift leaves the slope low for long stretches, where only the rare large earthquakes reach it. The earthquake sequence does not change with g.", SRC, "asserted"),
    theme: M("Slide background", [], "Light or dark, to match the deck. Colours only; the saved SVG uses the theme shown.", null),
    R0: M("Initial relief R<sub>0</sub>", ["R(0) = R<sub>0</sub>"],
      "Where the record starts, as a fraction of the functional limit L. A low value shows a long ramp of growing relief before the first failure; 1 starts at the limit.", SRC, "asserted"),
    tmax: M("Length of the record t<sub>max</sub>", ["time unit = L / g"],
      "How much time the slide spans, in units of the time relief takes to grow from nothing to the limit (L/g). With L = 3 km and g = 1 mm/yr one unit is 3 Myr. Longer records show more saw-teeth and more earthquakes, each drawn smaller.", null, "asserted"),
    dSmall: M("Landslide increment d, step 2", ["R ← L − d", "period = d / g,   mean relief = L − d/2"],
      "How much relief one landslide removes in step 2, as a fraction of L. Small failures keep relief pinned just under the limit: the saw-tooth is fine and the mean relief is almost L.", SRC, "asserted"),
    dLarge: M("Landslide increment d, steps 3–5", ["R ← L<sub>eff</sub> − d,   L<sub>eff</sub> = L (step 3),  L<sub>s</sub> − D (steps 4–5)"],
      "How far below the limit in force a failure leaves the slope, the same for every landslide in steps 3 to 5. Larger failures are rarer (period d/g) and hold the mean relief further below the limit (L − d/2). Steps 6 and 7 use it only if their landslide sizes are set to constant.", SRC, "asserted"),
    floorPct: M("Relief floor", ["R ← max(L<sub>eff</sub> − d,  R<sub>floor</sub>)", "R<sub>floor</sub> = (floor %) × graph height,   graph height = 1.1 L<sub>s</sub>"],
      "The lowest relief a landslide can leave, as a share of the height of the graph above its base. A landslide removes part of a slope, never the whole of it, so however large the failure and however slow the uplift, relief does not drop to zero. A slope already at the floor does not fail. Zero removes the floor.", SRC, "asserted"),
    dDist: M("Landslide sizes, steps 6–7", ["power law:  N(&gt;d) ∝ d<sup>−α</sup>,   d<sub>min</sub> ≤ d ≤ d<sub>max</sub>", "constant:  d = the step 3–5 size"],
      "Whether every landslide removes the same increment of relief or each draws its own. With the power law, landslide sizes plot as a straight line in log size against log frequency between the two cut-offs: many small failures that barely notch the slope and a few that take a large part of it. The sizes come from their own random stream, in order of occurrence, so changing them never changes the earthquakes.", "Truncated Pareto, drawn by inversion. Power-law landslide size statistics: e.g. Hovius et al. 1997; Stark & Hovius 2001; Malamud et al. 2004 (for areas and volumes; applying one to the relief increment is an assertion of this figure).", "asserted"),
    dMin: M("Smallest landslide d<sub>min</sub>, steps 6–7", ["d ≥ d<sub>min</sub>"],
      "The lower cut-off of the landslide size distribution, as a fraction of L. Because small landslides are the most common, this sets the typical failure: lower it and most failures become notches that relief regrows through almost at once.", SRC, "asserted"),
    dMax: M("Largest landslide d<sub>max</sub>, steps 6–7", ["d ≤ d<sub>max</sub>"],
      "The upper cut-off, as a fraction of L: the largest increment of relief a single failure can remove below the limit in force at that moment. Raise it and the rare large failures reset the slope further, after which it spends a long time out of reach of all but the largest earthquakes.", SRC, "asserted"),
    dAlpha: M("Size exponent α, steps 6–7", ["P(d &gt; x) = (x<sup>−α</sup> − d<sub>max</sub><sup>−α</sup>) / (d<sub>min</sub><sup>−α</sup> − d<sub>max</sub><sup>−α</sup>)"],
      "The slope of the size–frequency line in log-log space. A large α concentrates landslides near the smallest size, so relief rides high with a fine irregular saw-tooth; a small α spreads them toward the largest size, so large failures are common and mean relief is lower.", SRC, "asserted"),
    Ls: M("Static limit L<sub>s</sub>", ["L<sub>eff</sub>(t) = L<sub>s</sub> · (1 − f(t))", "steps 4–5:  f = (L<sub>s</sub> − L) / L<sub>s</sub>"],
      "The relief the slope could hold if nothing ever shook it, as a multiple of L. Relief that reaches it fails with no trigger at all. In steps 4 and 5 the constant earthquake depresses the limit by exactly L<sub>s</sub> − L, which is what makes L the functional limit; raising L<sub>s</sub> there lengthens the drop-lines and changes nothing else. Those earthquakes are stipulated to strike at the exact instants relief reaches L, whatever the uplift rate and landslide size, so their recurrence is d/g and is not a control.", SRC, "asserted"),
    rate: M("Earthquake rate λ, steps 6–7", ["Δt ~ Exponential(mean 1/λ)  (Poisson arrivals)"],
      "How many earthquakes of at least the smallest magnitude drawn strike per unit time, at random intervals. Most are small. More earthquakes sample the rare large ones more often, so relief is held lower; in step 7 they also wear the static limit down faster between landslides.", SRC, "asserted"),
    b: M("Gutenberg–Richter b-value, steps 6–7", ["N(&gt;M) ∝ 10<sup>−bM</sup>,   M<sub>max</sub> − ΔM ≤ M ≤ M<sub>max</sub>", "N(&gt;f) ∝ f<sup>−b/c</sup>,   c = 0.5"],
      "The size distribution of earthquakes. With b = 1, the value observed almost everywhere, each unit of magnitude down is ten times as frequent. A higher b gives a population dominated by small earthquakes that rarely reach the relief; a lower b gives more large ones, more landslides and lower relief. Magnitudes are drawn from the law truncated at both ends.", "Gutenberg & Richter 1944; doubly truncated form as in Cornell & Vanmarcke 1969", "asserted"),
    dM: M("Magnitude range ΔM, steps 6–7", ["M<sub>min</sub> = M<sub>max</sub> − ΔM", "f<sub>min</sub> = f<sub>max</sub> · 10<sup>−cΔM</sup>"],
      "How far below the largest earthquake the population extends. Each extra unit adds ten times as many earthquakes (at b = 1), each depressing the limit about a third as far, so a wide range fills the figure with short lines that never fail the slope. The rate counts all earthquakes in the range, so widening the range at a fixed rate makes large ones rarer.", SRC, "asserted"),
    fmax: M("Depression by the largest earthquake f<sub>max</sub>, steps 6–7", ["D = f · L<sub>s</sub>(t),   L<sub>eff</sub> = L<sub>s</sub>(t) · (1 − f)", "f = f<sub>max</sub> · 10<sup>−c (M<sub>max</sub> − M)</sup>,   c = 0.5"],
      "How far the largest earthquake pulls the limit down, as a fraction of the static limit as it stands at that moment (so a damaged slope is pulled down from its damaged limit, by the same fraction). Smaller earthquakes depress it in proportion to their shaking amplitude, which at a fixed distance grows about 10<sup>0.5</sup> ≈ 3.2 times per unit of magnitude. Distance to the source is not sampled: every earthquake is treated as equally far away, which is the main thing a real site would add.", "Amplitude scaling c ≈ 0.5 is an assertion of this figure (ground-motion models give roughly 0.3–0.5 for peak acceleration).", "asserted"),
    damage: M("Damage per earthquake k, step 7", ["no failure:  L<sub>s</sub> ← L<sub>s</sub> − k·D = L<sub>s</sub> (1 − k f)", "landslide:  L<sub>s</sub> ← L<sub>s,0</sub>"],
      "The degree of weakening by an earthquake that does not fail the slope: how much of its momentary depression of the limit is left behind as damage. In the frame it is the slider called Weakening. Zero is step 6. The static limit ratchets down with every earthquake that does not cause a landslide, large ones most, so each later earthquake reaches lower and the slope fails sooner than its intact strength would allow. The landslide, when it comes, removes the damaged rock and the static limit returns to its original value. Damage therefore lives for one saw-tooth: the longer a slope goes without failing, the weaker it has become.", SRC, "asserted"),
    seed: M("Random sequence, steps 6–7", ["mulberry32(seed)"],
      "Which sequence of earthquakes is drawn. The same seed always gives the same figure, so a slide can be reproduced; steps 6 and 7 share the sequence, so the only difference between them is the damage.", "core/noise.js", null)
  };

  function fL(v) { return v.toFixed(2) + " L"; }
  function fT(v) { return v.toFixed(2) + " L/g₀"; }
  var specs = [
    { id: "stage", label: "Step", type: "select", val: "1", group: "Figure", method: m.stage,
      options: STAGES.map(function (s, i) { return { v: String(i + 1), t: (i + 1) + " · " + s[0] }; }) },
    { id: "drawEarly", label: "Drawing time, steps 1–5", type: "range", val: 12, min: 1, max: 40, step: 1, unit: "s", group: "Figure", method: m.drawEarly },
    { id: "drawTime", label: "Drawing time, steps 6–7", type: "range", val: 12, min: 1, max: 40, step: 1, unit: "s", group: "Figure", method: m.drawTime },
    { id: "theme", label: "Background", type: "select", val: "light", group: "Figure", method: m.theme,
      options: [{ v: "light", t: "light" }, { v: "dark", t: "dark" }] },
    { id: "frame", label: "Sliders in frame", type: "select", val: "last", group: "Figure", method: m.frame,
      options: [{ v: "last", t: "last step only" }, { v: "all", t: "every step" }, { v: "off", t: "hidden" }] },
    { id: "g", label: "Uplift rate", type: "range", val: 1, min: 0, max: 1, step: 0.02, fmt: function (v) { return v.toFixed(2) + " × g₀"; }, group: "Relief and landslides", method: m.g },
    { id: "R0", label: "Initial relief", type: "range", val: 0.3, min: 0, max: 1, step: 0.05, fmt: fL, group: "Relief and landslides", method: m.R0 },
    { id: "tmax", label: "Length of record", type: "range", val: 4, min: 2, max: 12, step: 0.5, fmt: fT, group: "Relief and landslides", method: m.tmax },
    { id: "dSmall", label: "Landslide size, step 2", type: "range", val: 0.05, min: 0.01, max: 0.2, step: 0.01, fmt: fL, group: "Relief and landslides", method: m.dSmall },
    { id: "dLarge", label: "Landslide size, steps 3–5", type: "range", val: 0.25, min: 0.05, max: 0.6, step: 0.01, fmt: fL, group: "Relief and landslides", method: m.dLarge },
    { id: "floorPct", label: "Relief floor", type: "range", val: 10, min: 0, max: 30, step: 1, fmt: function (v) { return v + " % of graph height"; }, group: "Relief and landslides", method: m.floorPct },
    { id: "dDist", label: "Landslide sizes, steps 6–7", type: "select", val: "power", group: "Relief and landslides", method: m.dDist,
      options: [{ v: "power", t: "random, power law" }, { v: "constant", t: "constant (step 3–5 size)" }] },
    { id: "dMin", label: "Smallest landslide, steps 6–7", type: "range", val: 0.05, min: 0.01, max: 0.3, step: 0.01, fmt: fL, group: "Relief and landslides", method: m.dMin },
    { id: "dMax", label: "Largest landslide, steps 6–7", type: "range", val: 0.6, min: 0.1, max: 1.2, step: 0.01, fmt: fL, group: "Relief and landslides", method: m.dMax },
    { id: "dAlpha", label: "Size exponent α, steps 6–7", type: "range", val: 1, min: 0.2, max: 3, step: 0.1, fmt: function (v) { return v.toFixed(1); }, group: "Relief and landslides", method: m.dAlpha },
    { id: "Ls", label: "Static limit", type: "range", val: 1.6, min: 1.1, max: 2.5, step: 0.05, fmt: fL, group: "Earthquakes", method: m.Ls },
    { id: "rate", label: "Rate, steps 6–7", type: "range", val: 30, min: 2, max: 80, step: 1, fmt: function (v) { return v + " per L/g₀"; }, group: "Earthquakes", method: m.rate },
    { id: "b", label: "b-value, steps 6–7", type: "range", val: 1, min: 0.5, max: 1.5, step: 0.05, fmt: function (v) { return v.toFixed(2); }, group: "Earthquakes", method: m.b },
    { id: "dM", label: "Magnitude range, steps 6–7", type: "range", val: 2, min: 0.5, max: 3, step: 0.1, fmt: function (v) { return v.toFixed(1) + " units"; }, group: "Earthquakes", method: m.dM },
    { id: "fmax", label: "Largest depression, steps 6–7", type: "range", val: 0.6, min: 0.1, max: 0.95, step: 0.01, fmt: function (v) { return Math.round(v * 100) + " % of limit"; }, group: "Earthquakes", method: m.fmax },
    { id: "damage", label: "Damage, step 7", type: "range", val: 0.15, min: 0, max: 1, step: 0.01, fmt: function (v) { return v.toFixed(2) + " of D"; }, group: "Earthquakes", method: m.damage },
    { id: "seed", label: "Random sequence", type: "range", val: 82, min: 1, max: 99, step: 1, group: "Earthquakes", method: m.seed }
  ];
  var url = GS.ui.urlState(GS.ui.toUrlSpecs(specs));
  var state = url.read();

  // Embed mode (?embed=1): the page is a slide inside another deck (an iframe). Only the
  // 16:9 frame is shown; the host drives the steps with postMessage({type: "relief-step",
  // step: n, replay: bool}) and the step keys are handed back to it as {type: "deck-key",
  // key}, so a clicker keeps working when the figure has the focus. The URL is left alone.
  var EMBED = /[?&]embed=1\b/.test(window.location.search);
  if (EMBED) { document.body.classList.add("embed"); url.write = function () {}; }

  // ---- Figure ---------------------------------------------------------------
  var slide = d3.select("#stage").append("div").attr("id", "slide");
  var svg = slide.append("svg").attr("xmlns", "http://www.w3.org/2000/svg")
    .attr("viewBox", "0 0 " + W + " " + H).attr("font-family", FONT);
  var bg = svg.append("rect").attr("width", W).attr("height", H);
  var clip = svg.append("defs").append("clipPath").attr("id", "rl-clip").append("rect").attr("x", 0).attr("y", 0).attr("height", H);
  var axis = svg.append("path").attr("fill", "none").attr("stroke-width", WT).attr("d", "M" + X0 + "," + Y1 + "V" + Y0 + "H" + X1);
  var yLab = svg.append("text").attr("x", X0).attr("y", Y1 - 28).attr("font-size", 38).text("Relief");
  var xLab = svg.append("text").attr("x", X1).attr("y", Y0 + 58).attr("text-anchor", "end").attr("font-size", 38).text("Time");

  function fader(sel, sec) { return sel.style("transition", "opacity " + sec + "s").style("opacity", 0); }
  var funcG = fader(svg.append("g"), 0.6);
  var funcLine = funcG.append("line").attr("x1", X0).attr("x2", X1).attr("stroke-width", WT).attr("stroke-dasharray", "16 12");
  var funcLab = funcG.append("text").attr("x", X1 + 18).attr("font-size", 34);
  var statG = fader(svg.append("g"), 0.6);
  var statLine = statG.append("path").attr("fill", "none").attr("stroke-width", WT).attr("stroke-dasharray", "16 12");
  var statLab = statG.append("text").attr("x", X1 + 18).attr("font-size", 34).attr("dy", "0.34em").text("Static limit");
  var relief = svg.append("path").attr("fill", "none").attr("stroke-width", WT).attr("stroke-linejoin", "round")
    .attr("clip-path", "url(#rl-clip)").style("pointer-events", "none");
  var quakeG = svg.append("g");
  // wide invisible strokes that catch the pointer for the curve tooltips (as ui/panels.js)
  function hit(node) {
    return svg.append(node).attr("fill", "none").attr("stroke", "transparent").attr("stroke-width", 28)
      .style("pointer-events", "stroke").style("cursor", "help");
  }
  var hitFunc = hit("line").attr("x1", X0).attr("x2", X1), hitStat = hit("path"), hitRelief = hit("path").attr("clip-path", "url(#rl-clip)");

  var sim = null, pal = THEMES.light, tc = 0, playing = false, lastNow = 0, yMax = 1;
  function x(t) { return X0 + (X1 - X0) * t / sim.p.tmax; }
  function y(r) { return Y0 - (Y0 - Y1) * r / yMax; }
  function poly(pts) {
    var s = "", i;
    for (i = 0; i < pts.length; i++) s += (i ? "L" : "M") + x(pts[i][0]).toFixed(2) + "," + y(pts[i][1]).toFixed(2);
    return s;
  }
  function modelParams() {
    return { stage: +state.stage, g: state.g, floor: state.floorPct / 100 * 1.1 * Math.max(state.Ls, 1), R0: state.R0, tmax: state.tmax, dSmall: state.dSmall, dLarge: state.dLarge,
             dDist: state.dDist, dMin: state.dMin, dMax: state.dMax, dAlpha: state.dAlpha, Ls: state.Ls,
             rate: state.rate, b: state.b, dM: state.dM, fmax: state.fmax, damage: state.damage, seed: state.seed };
  }
  function samePath(a, b) {
    if (!a || !b || a.length !== b.length) return false;
    for (var i = 0; i < a.length; i++) if (Math.abs(a[i][0] - b[i][0]) > 1e-9 || Math.abs(a[i][1] - b[i][1]) > 1e-9) return false;
    return true;
  }

  // Solve the model and lay out everything that does not depend on the cursor.
  function rebuild() {
    var st = +state.stage, p;
    sim = RL.simulate(modelParams()); p = sim.p;
    pal = THEMES[state.theme] || THEMES.light;
    yMax = 1.1 * Math.max(p.Ls, p.L); // the same frame in every step, so the static limit arrives in empty space
    d3.select("#stage").style("background", pal.bg);
    bg.attr("fill", pal.bg);
    axis.attr("stroke", pal.axis);
    yLab.attr("fill", pal.text); xLab.attr("fill", pal.text);

    funcLine.attr("y1", y(p.L)).attr("y2", y(p.L)).attr("stroke", pal.limit);
    hitFunc.attr("y1", y(p.L)).attr("y2", y(p.L)).style("display", st <= 5 ? null : "none");
    var lines = st <= 3 ? ["Limit"] : ["Limit during", "earthquakes"];
    funcLab.attr("fill", pal.text).attr("y", y(p.L) + (lines.length > 1 ? -8 : 12)).selectAll("tspan").remove();
    lines.forEach(function (s, i) { funcLab.append("tspan").attr("x", X1 + 18).attr("dy", i ? 40 : 0).text(s); });
    funcG.style("opacity", st <= 5 ? 1 : 0);

    statLine.attr("stroke", pal.limit); statLab.attr("fill", pal.text);
    statG.style("opacity", st >= 4 ? 1 : 0);
    hitStat.style("display", st >= 4 ? null : "none");

    relief.attr("stroke", pal.ink).attr("d", poly(sim.path));
    hitRelief.attr("d", poly(sim.path));

    quakeG.selectAll("g").remove();
    if (st >= 5) {
      var q = quakeG.selectAll("g").data(sim.quakes).enter().append("g").style("transition", "opacity 0.3s").style("opacity", 0).style("cursor", "help");
      q.append("line").attr("x1", function (d) { return x(d.t); }).attr("x2", function (d) { return x(d.t); })
        .attr("y1", function (d) { return y(d.top); }).attr("y2", function (d) { return y(d.y); })
        .attr("stroke", pal.quake).attr("stroke-width", WT / 2).attr("stroke-dasharray", "6 6")
        .attr("stroke-opacity", function (d) { return d.trig ? 1 : 0.6; });
      q.append("circle").attr("cx", function (d) { return x(d.t); }).attr("cy", function (d) { return y(d.y); }).attr("r", 9)
        .attr("stroke", pal.quake).attr("stroke-width", WT).attr("fill", function (d) { return d.trig ? pal.quake : pal.bg; });
      q.each(function (d) {
        GS.ui.tip.attach(d3.select(this), GS.ui.curveTip(d.trig ? "Earthquake that fails the slope" : "Earthquake the slope survives", pal.quake, "dashed",
          "At t = " + fT(d.t) + " the limit is pulled from the static limit " + fL(d.top) + " down to the dot at " + fL(d.y) + ": " + (100 * d.f).toFixed(0) + " % of the limit as it stands" + (st >= 6 ? ", an earthquake " + d.m.toFixed(1) + " magnitude units below the largest" : "") + ". " +
          (d.trig ? "The dot is at or below the relief, so the slope fails and relief drops to d below the dot."
                  : "The dot stays above the relief: no landslide." + (st === 7 ? " The slope is damaged instead; the static limit drops by k·D = " + fL(sim.p.damage * d.D) + " until the next landslide." : "")),
          "L<sub>eff</sub> = L<sub>s</sub> · (1 − f)"));
      });
      window.getComputedStyle(svg.node()).opacity; // flush, so the quakes fade in from 0
    }
  }

  // Everything that follows the time cursor.
  function draw() {
    var p = sim.p, lim = sim.lim, pts = [], i, c = lim[0][1];
    clip.attr("width", tc >= p.tmax ? W : x(tc));
    for (i = 0; i < lim.length && lim[i][0] <= tc; i++) { pts.push(lim[i]); c = lim[i][1]; }
    pts.push([tc, c]); pts.push([p.tmax, c]);
    var dStat = poly(pts);
    statLine.attr("d", dStat); hitStat.attr("d", dStat);
    statLab.attr("y", y(c));
    quakeG.selectAll("g").style("opacity", function (d) { return d.t <= tc ? 1 : 0; })
      .style("pointer-events", function (d) { return d.t <= tc ? null : "none"; });
  }

  // ---- In-frame sliders ------------------------------------------------------
  // HTML laid over the SVG inside #slide, so they go full screen with it. Each mirrors a
  // control of the presenter panel (same state id) and the two stay in step.
  var frame = slide.append("div").attr("id", "frame-controls");
  var frameCtl = [
    { id: "damage", label: "Weakening", fmt: function (v) { return Math.round(v * 100) + " % of each dip"; } },
    { id: "g", label: "Uplift rate", fmt: function (v) { return "× " + v.toFixed(2); } }
  ];
  frameCtl.forEach(function (fc) {
    var c = specs.filter(function (s) { return s.id === fc.id; })[0];
    var box = frame.append("label").attr("class", "fc");
    var head = box.append("span").attr("class", "fc-head");
    head.append("span").text(fc.label);
    fc.val = head.append("span").attr("class", "fc-val");
    fc.inp = box.append("input").attr("type", "range").attr("min", c.min).attr("max", c.max).attr("step", c.step)
      .on("input", function () { state[fc.id] = +this.value; controls.refresh(); onChange(fc.id); })
      .on("change", function () { this.blur(); }); // hand the arrow keys back to the step control
  });
  function syncFrame() {
    var on = state.frame === "all" || (state.frame === "last" && +state.stage === STAGES.length);
    frame.style("display", on ? null : "none")
      .style("--fc-text", pal.text).style("--fc-line", pal.axis).style("--fc-ink", pal.ink);
    frameCtl.forEach(function (fc) { fc.inp.property("value", state[fc.id]); fc.val.text(fc.fmt(state[fc.id])); });
  }

  // ---- Curve tooltips, readout ---------------------------------------------
  GS.ui.tip.attach(hitRelief, function () {
    return GS.ui.curveTip("Relief", pal.ink, "solid",
      "Grows at the constant rate g (uplift and valley incision) and drops instantly at each landslide" + (+state.stage === 1 ? "; in this step it simply holds at the limit." : sim.fails.length ? ", to d below the limit in force (d from " + fL(d3.min(sim.fails, function (f) { return f.d; })) + " to " + fL(d3.max(sim.fails, function (f) { return f.d; })) + " in this record)." : "."),
      "dR/dt = g;  failure: R ← L<sub>eff</sub> − d");
  });
  GS.ui.tip.attach(hitFunc, function () {
    return +state.stage <= 3
      ? GS.ui.curveTip("Limit to relief", pal.limit, "dashed", "The relief at which the slope fails. Constant, and for now unexplained.", "R ≤ L")
      : GS.ui.curveTip("Limit during earthquakes", pal.limit, "dashed", "The static limit less the depression of the constant earthquake. It is in force only at the instants an earthquake strikes; between them the slope could stand higher.", "L = L<sub>s</sub> · (1 − f)");
  });
  GS.ui.tip.attach(hitStat, function () {
    return GS.ui.curveTip("Static limit", pal.limit, "dashed",
      "The relief the unshaken slope can hold." + (+state.stage === 7 ? " Each earthquake that does not cause a landslide lowers it by k·D; the next landslide removes the damaged rock and restores it. Right of the cursor it is drawn level, as it stands now." : " Failures happen below it because earthquakes lower the limit for a moment."),
      +state.stage === 7 ? "L<sub>s</sub> ← L<sub>s</sub> − k·D;  landslide: L<sub>s</sub> ← L<sub>s,0</sub>" : "L<sub>eff</sub> = L<sub>s</sub> · (1 − f)");
  });

  function readout() {
    var st = +state.stage, nq = sim.quakes.length, tr = sim.quakes.filter(function (q) { return q.trig; }).length;
    var h = "<b>" + st + " · " + STAGES[st - 1][0] + ".</b> " + STAGES[st - 1][1] + "<br>";
    h += "Landslides: <b>" + sim.fails.length + "</b>  ·  mean relief <b>" + fL(RL.meanRelief(sim.path)) + "</b>";
    if (nq) h += "  ·  earthquakes: <b>" + nq + "</b>, of which <b>" + tr + "</b> fail the slope";
    if (st === 7) h += "  ·  static limit " + fL(sim.p.Ls) + ", damaged to as low as <b>" + fL(d3.min(sim.lim, function (v) { return v[1]; })) + "</b>";
    d3.select("#readout").html(h);
  }

  // ---- Animation and transport ----------------------------------------------
  function tick(now) {
    if (!playing) return;
    var dt = Math.min((now - lastNow) / 1000, 0.1); lastNow = now;
    tc = Math.min(sim.p.tmax, tc + dt * sim.p.tmax / (+state.stage <= 5 ? state.drawEarly : state.drawTime));
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
  function setStage(n) {
    n = Math.max(1, Math.min(STAGES.length, n));
    if (n === +state.stage) return;
    state.stage = String(n); controls.refresh(); onChange("stage");
  }
  function fullscreen() {
    var el = document.getElementById("stage");
    if (document.fullscreenElement) document.exitFullscreen(); else if (el.requestFullscreen) el.requestFullscreen();
  }
  function saveSvg() {
    var src = '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg.node());
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([src], { type: "image/svg+xml" }));
    a.download = "relief_limit_step" + state.stage + ".svg";
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 1000);
  }

  var tp = d3.select("#transport"), stageBtns = [], playBtn;
  STAGES.forEach(function (s, i) {
    stageBtns.push(tp.append("button").attr("class", "stg").attr("title", s[0]).text(i + 1).on("click", function () { setStage(i + 1); this.blur(); }));
  });
  tp.append("span").attr("class", "gap");
  playBtn = tp.append("button").on("click", function () { if (playing) pause(); else play(); this.blur(); });
  tp.append("button").text("Replay").on("click", function () { replay(); this.blur(); });
  tp.append("button").text("Full screen").on("click", function () { fullscreen(); this.blur(); });
  tp.append("button").text("Save SVG").attr("title", "The frame as shown, 1920 × 1080, self-contained").on("click", function () { saveSvg(); this.blur(); });
  tp.append("span").attr("class", "keys").html("<kbd>←</kbd> <kbd>→</kbd> step  ·  <kbd>space</kbd> play / pause  ·  <kbd>R</kbd> replay  ·  <kbd>F</kbd> full screen  ·  <kbd>1</kbd>–<kbd>7</kbd>");
  function buttons() {
    playBtn.text(playing ? "Pause" : "Play");
    stageBtns.forEach(function (b, i) { b.classed("current", i + 1 === +state.stage); });
  }

  d3.select(window).on("keydown", function (ev) {
    var tag = ev.target && ev.target.tagName;
    if (tag === "INPUT" || tag === "SELECT" || tag === "TEXTAREA" || ev.metaKey || ev.ctrlKey || ev.altKey) return;
    var k = ev.key;
    if (EMBED && (k === "ArrowRight" || k === "ArrowLeft" || k === "PageDown" || k === "PageUp" || k === " ")) {
      window.parent.postMessage({ type: "deck-key", key: k }, "*");
      ev.preventDefault();
      return;
    }
    if (k === "ArrowRight" || k === "PageDown") setStage(+state.stage + 1);
    else if (k === "ArrowLeft" || k === "PageUp") setStage(+state.stage - 1);
    else if (k === " ") { if (playing) pause(); else play(); }
    else if (k === "r" || k === "R") replay();
    else if (k === "f" || k === "F") fullscreen();
    else if (k >= "1" && k <= "7") setStage(+k);
    else return;
    ev.preventDefault();
  });

  // ---- Controls -------------------------------------------------------------
  function onChange(id) {
    var old = sim ? sim.path : null;
    rebuild();
    if (id === "stage") {
      if (!samePath(old, sim.path)) { tc = 0; draw(); playing = false; play(); } // a new history: draw it
      else draw();                                                              // same history: keep the cursor, new elements fade in
    } else if (id !== "theme" && id !== "drawTime" && id !== "drawEarly" && id !== "frame") {
      playing = false; tc = sim.p.tmax; draw();                                 // tuning a parameter: show the whole record
    } else draw();
    syncFrame(); readout(); buttons(); url.write(state);
  }
  var controls = GS.ui.buildControls(d3.select("#controls"), specs, state, onChange);

  rebuild(); draw(); syncFrame(); readout(); buttons(); url.write(state);
  play();

  window.addEventListener("message", function (ev) {
    var d = ev.data;
    if (!EMBED || !d || d.type !== "relief-step") return;
    if (+d.step !== +state.stage) setStage(+d.step); else if (d.replay) replay();
  });
})();
