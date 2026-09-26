/*
 * d1/app.js — controls + D3 panels for the column tier. Depends on core/*.js,
 * ui/*.js (window.GS) and d3 v7. Display units are mm/yr, km, kyr, Myr; the model
 * runs in m and yr (core/units.js converts at this boundary only).
 *
 * Two scenarios share one page: "reservoir" (the legacy relief reservoir) and
 * "columns" (two lockstep columns A/B: lithology persistence + ELA buzzsaw).
 * Every control and panel carries a `method` tooltip with its equations.
 */
(function () {
  "use strict";
  var GS = window.GS, U = GS.units;

  var C = { uplift: "#2c7fb8", relief: "#238b45", erosion: "#d94801", flux: "#888", eq: "#999", net: "#6a51a3",
            A: "#238b45", B: "#8c2d04", ela: "#4292c6", body: "rgba(140,45,4,0.18)", diff: "#6a51a3" };

  // ---- Methodology text ----------------------------------------------------
  // M(title, [equations], text, source, tag) → tooltip HTML. Equations are HTML.
  function M(title, eqs, text, src, tag) {
    var h = "<h4>" + (tag ? '<span class="tag ' + tag + '">' + tag + "</span>" : "") + title + "</h4>";
    (eqs || []).forEach(function (e) { h += '<span class="eq">' + e + "</span>"; });
    if (text) h += "<div>" + text + "</div>";
    if (src) h += '<div class="src">' + src + "</div>";
    return h;
  }
  var SRC_RES = "Relief-reservoir model (legacy uplift_sandbox); cf. Whipple 2004 style steady-state relief scaling.";
  var m = {
    shape: M("Uplift forcing U(t)", ["U(t) = U<sub>peak</sub> · s(t)"],
      "s(t) ∈ [0, 1] is the event shape on [0, T]: bell (Gaussian, σ = T/6), plateau (25 % ramp, 50 % hold, 25 % ramp), pulse (10 % rise, exponential decay τ = 0.3 T), arc, triangle or square pulse. The run continues to 2 T so the relaxation tail is visible.",
      "core/forcing.js", "asserted"),
    peakUplift: M("Peak rock-uplift rate", ["U<sub>peak</sub> (mm/yr = km/Myr = 10<sup>−3</sup> m/yr)"],
      "Rock uplift relative to a fixed base level. Applied identically to every column.", "core/forcing.js", "asserted"),
    duration: M("Event duration T", ["s(t) = 0 for t > T (except the square pulse: on until T)"],
      "Length of the uplift event; the dashed vertical line marks its end.", null, "asserted"),
    reliefLimit: M("Relief limit R<sub>lim</sub> → erodibility k", ["k = U<sub>peak</sub> / R<sub>lim</sub><sup>n</sup>", "R<sub>eq</sub>(U) = (U / k)<sup>1/n</sup>"],
      "The equilibrium relief the landscape would reach under sustained peak uplift. It sets k, the constant of the erosion law, in an intuitive unit. The column does not resolve the drainage network that actually sets k, so this is asserted.",
      SRC_RES, "asserted"),
    exponent: M("Erosion nonlinearity n", ["E(R) = k R<sup>n</sup>"],
      "How steeply erosion grows with relief. Large n gives a hard relief limit (R<sub>eq</sub> ∝ U<sup>1/n</sup>).", SRC_RES, "asserted"),
    initialRelief: M("Initial relief / elevation", ["R(0) = R<sub>0</sub>  (columns: z<sub>A</sub>(0) = z<sub>B</sub>(0) = z<sub>0</sub>)"],
      "Starting state above base level. In the columns scenario the resistant body is positioned relative to this initial surface.", null, "asserted"),
    bodyDepth: M("Resistant body: depth of its top below the initial surface", ["top: z<sub>m</sub> = z<sub>0</sub> − d", "z<sub>m</sub> = z − U<sub>cum</sub>(t)  (material frame)"],
      "The body is defined in material coordinates, so it rides up with uplift and is exposed once cumulative erosion of column B reaches d, whatever the uplift rate (checked in py/test_d1_columns.py).",
      "core/lithology.js layer generator", "asserted"),
    bodyThick: M("Resistant body thickness T", ["exposed while  z<sub>0</sub> − d − T ≤ z<sub>B</sub> − U<sub>cum</sub> ≤ z<sub>0</sub> − d", "lifetime at the surface ≈ T · c / E<sub>ref</sub>"],
      "Hypothesis B needs the body to survive the uplift event: its lifetime once exposed is its thickness divided by its (reduced) erosion rate. Compare that with the event duration.", null, "asserted"),
    contrast: M("Subaerial erodibility contrast c", ["E<sub>sub,B</sub> = D<sub>B</sub> · f(r, c),   f = 1 / (1 + r (c − 1))"],
      "r ∈ [0, 1] is the resistance at B's surface (1 inside the body); D<sub>B</sub> is the erosional demand set by the coupling control below. While the body is at the surface B lowers at D<sub>B</sub>/c, so its lifetime is T·c/D<sub>B</sub>.",
      "core/lithology.js erodibilityFactor", "asserted"),
    coupling: M("What sets the erosional demand on column B?", [
        "landscape:  D<sub>B</sub> = k z<sub>A</sub><sup>n</sup>", "own relief:  D<sub>B</sub> = k z<sub>B</sub><sup>n</sup>", "max:  D<sub>B</sub> = k max(z<sub>A</sub>, z<sub>B</sub>)<sup>n</sup>"],
      "The one real design choice of this tier, so it is a switch. <b>landscape</b>: B is a low-relief plateau inside A's landscape; its erosion is A's rate reduced by c, whatever height it reaches. Consequence: B can climb without limit while protected, and once stripped the upland <i>freezes</i> (E<sub>B</sub> = E<sub>A</sub>); with A glacially capped, A's small subaerial rate lets B run away. <b>own relief</b>: B is its own relief reservoir; its erosion grows with its height, so while protected it equilibrates at c<sup>1/n</sup>·R<sub>eq</sub> and then erodes at U: the body is consumed in ≈T/U and persistence needs great thickness. <b>max</b>: landscape while B is low, own relief once it stands above A, so the upland is bounded and decays after stripping. Which is right is the hypothesis-B question; the sandbox shows the three side by side.",
      "docs/d1.md", "asserted"),
    contrastG: M("Glacial erodibility contrast c<sub>g</sub>", ["E<sub>gl,B</sub> = E<sub>gl</sub>(z<sub>B</sub>) · f(r, c<sub>g</sub>)"],
      "Same body, its own contrast for glacial erosion (a quartzite resists abrasion differently than it resists fluvial incision).", null, "asserted"),
    glacierOn: M("Elevation-dependent glacial erosion (buzzsaw)", ["E<sub>gl</sub>(z) = P · g(z),  g = exp(−½ u<sup>2</sup>) below the ELA,  φ + (1 − φ) exp(−½ u<sup>2</sup>) above", "u = (z − ELA(t)) / w"],
      "The column asserts the ice: erosion peaks at the equilibrium-line altitude and fades over a half-width w. With P ≫ U the column is capped where E<sub>gl</sub> + E<sub>sub</sub> = U, on the lower flank of the window: z<sub>cap</sub> = ELA − w √(2 ln(P / (U − E<sub>sub</sub>))). Hypothesis A in one dimension.",
      "Buzzsaw concept: Brozović et al. 1997; Egholm et al. 2009; Mitchell & Montgomery 2006. Gaussian window is an assertion of this tier.", "asserted"),
    elaBase: M("ELA (interglacial)", ["ELA(t) = ELA<sub>0</sub> − A<sub>ELA</sub> · ½(1 − cos(2π t / P<sub>ELA</sub>))"], "Equilibrium-line altitude above base level in the warm state.", null, "asserted"),
    elaAmp: M("ELA cycle amplitude", ["ELA<sub>min</sub> = ELA<sub>0</sub> − A<sub>ELA</sub>"], "Glacial-cycle depression of the ELA (0 = constant ELA).", null, "asserted"),
    elaPeriod: M("ELA cycle period", ["P<sub>ELA</sub>"], "e.g. 100 kyr for late-Pleistocene cycles, 41 kyr for obliquity-paced ones.", null, "asserted"),
    buzzPeak: M("Peak glacial erosion rate P", ["E<sub>gl</sub>(ELA) = P"], "Maximum lowering rate at the ELA. Real glacial erosion rates span 0.01–10 mm/yr (Hallet et al. 1996).", null, "asserted"),
    buzzWidth: M("Buzzsaw half-width w", ["E<sub>gl</sub>(ELA ± w) = P · e<sup>−½</sup>"], "Vertical extent of the erosion window round the ELA.", null, "asserted"),
    buzzAbove: M("Erosion floor above the ELA φ", ["E<sub>gl</sub>(z ≫ ELA) → φ · P"],
      "φ = 0: the symmetric Gaussian window, so ground that rises above the window escapes glacial erosion. φ > 0: deep ice keeps eroding high ground at a fraction φ of the peak rate — the hypothesis-A ingredient (planation under deep ice). With φ·P/c<sub>g</sub> > U no column can escape.", "Egholm et al. 2009 (ice-cap erosion of high surfaces)", "asserted"),

    pUplift: M("Uplift forcing", ["U(t) = U<sub>peak</sub> s(t)"], "The shared forcing. Dashed vertical line: end of the event.", "core/forcing.js"),
    pRelief: M("Relief reservoir", ["dR/dt = U(t) − k R<sup>n</sup>", "R<sub>eq</sub>(t) = (U(t) / k)<sup>1/n</sup>  (dashed)"],
      "RK4 with 2000 steps; the stiff R<sup>n</sup> term is why RK4 rather than Euler. Relief chases the instantaneous equilibrium and lags it.", SRC_RES, "computed"),
    pErosion: M("Erosion vs uplift flux", ["E = k R<sup>n</sup>"], "When the two curves cross, relief peaks (dR/dt = 0).", null, "computed"),
    pNet: M("Net rate", ["dR/dt = U − E"], "Above zero relief grows, below zero it decays.", null, "computed"),
    pElev: M("Column elevations", ["dz<sub>A</sub>/dt = U − k z<sub>A</sub><sup>n</sup> − E<sub>gl</sub>(z<sub>A</sub>)", "dz<sub>B</sub>/dt = U − D<sub>B</sub> f(r<sub>B</sub>, c) − E<sub>gl</sub>(z<sub>B</sub>) f(r<sub>B</sub>, c<sub>g</sub>)"],
      "A (green): weak rock, the landscape's relief reservoir. B (brown): same forcing, with a resistant body (shaded band, drawn in the present frame: it rises with U<sub>cum</sub>). While B's surface lies inside the band its erosion is reduced. Blue dashed: ELA(t) when the glacier is on. RK4 on (z<sub>A</sub>, z<sub>B</sub>, U<sub>cum</sub>, E<sub>cum,A</sub>, E<sub>cum,B</sub>).",
      "py/d1_columns.py mirrors this; py/test_d1_columns.py checks z = z<sub>0</sub> + U<sub>cum</sub> − E<sub>cum</sub>, the T·c/U lifetime, exposure at d/U, the analytic buzzsaw cap.", "computed"),
    pRates: M("Erosion rates", ["E<sub>A</sub> = k z<sub>A</sub><sup>n</sup> + E<sub>gl</sub>(z<sub>A</sub>)", "E<sub>B</sub> = D<sub>B</sub> f(r<sub>B</sub>, c) + E<sub>gl</sub>(z<sub>B</sub>) f(r<sub>B</sub>, c<sub>g</sub>)"],
      "Total lowering rate of each column against the uplift rate (dashed grey). B's rate drops by the contrast while the body is at its surface.", null, "computed"),
    pDiff: M("Upland height z<sub>B</sub> − z<sub>A</sub>", ["d(z<sub>B</sub> − z<sub>A</sub>)/dt = E<sub>A</sub> − E<sub>B</sub>"],
      "The hypothesis-B diagnostic: how high the resistant column stands above the weak landscape, and whether that survives the uplift event. Under the landscape coupling it can only grow or freeze (destruction by marginal retreat needs d2_across); under own-relief or max it is bounded and decays after stripping.", null, "computed")
  };

  // ---- Controls (display units) ------------------------------------------
  var SHAPES = [
    { v: "bell", t: "Bell curve" }, { v: "plateau", t: "Ramp · hold · decay" }, { v: "pulse", t: "Pulse (sharp rise, decay)" },
    { v: "arc", t: "Circular arc" }, { v: "triangle", t: "Triangle" }, { v: "constant", t: "Constant (square pulse)" }
  ];
  var SCEN = [{ v: "columns", t: "Two columns A / B (lithology + buzzsaw)" }, { v: "reservoir", t: "Relief reservoir (legacy)" }];
  var common = [
    { id: "shape", label: "Uplift shape", type: "select", val: "bell", options: SHAPES, group: "Forcing", method: m.shape },
    { id: "peakUplift", label: "Peak uplift", type: "range", min: 0.1, max: 10, step: 0.1, val: 3, unit: "mm/yr", group: "Forcing", method: m.peakUplift },
    { id: "duration", label: "Event duration", type: "range", min: 1, max: 40, step: 0.5, val: 12, unit: "Myr", group: "Forcing", method: m.duration }
  ];
  var landscape = [
    { id: "reliefLimit", label: "Relief limit @ peak U", type: "range", min: 0.5, max: 8, step: 0.1, val: 3, unit: "km", group: "Landscape (column A)",
      asserted: true, method: m.reliefLimit },
    { id: "exponent", label: "Erosion exponent n", type: "range", min: 1, max: 6, step: 0.5, val: 4, unit: "", group: "Landscape (column A)", asserted: true, method: m.exponent },
    { id: "initialRelief", label: "Initial relief", type: "range", min: 0, max: 4, step: 0.05, val: 0, unit: "km", group: "Landscape (column A)", method: m.initialRelief }
  ];
  var body = [
    { id: "bodyDepth", label: "Body top depth below initial surface", type: "range", min: 0, max: 3, step: 0.05, val: 0, unit: "km", group: "Resistant body (column B)", asserted: true, method: m.bodyDepth },
    { id: "bodyThick", label: "Body thickness", type: "range", min: 0, max: 3, step: 0.05, val: 0.5, unit: "km", group: "Resistant body (column B)", asserted: true, method: m.bodyThick },
    { id: "contrast", label: "Subaerial contrast c", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Resistant body (column B)", asserted: true, method: m.contrast },
    { id: "contrastG", label: "Glacial contrast c_g", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Resistant body (column B)", asserted: true, method: m.contrastG },
    { id: "coupling", label: "Erosional demand on B", type: "select", val: "landscape", group: "Resistant body (column B)", asserted: true, method: m.coupling,
      options: [{ v: "landscape", t: "A's landscape rate  k·z_Aⁿ" }, { v: "own", t: "B's own relief  k·z_Bⁿ" }, { v: "max", t: "max of the two" }] }
  ];
  var glacier = [
    { id: "glacierOn", label: "Glacial buzzsaw on", type: "check", val: false, group: "Glacier (asserted ice)", method: m.glacierOn },
    { id: "elaBase", label: "ELA (interglacial)", type: "range", min: 0.5, max: 5, step: 0.1, val: 2, unit: "km", group: "Glacier (asserted ice)", asserted: true, method: m.elaBase },
    { id: "elaAmp", label: "ELA cycle amplitude", type: "range", min: 0, max: 1.5, step: 0.05, val: 0, unit: "km", group: "Glacier (asserted ice)", asserted: true, method: m.elaAmp },
    { id: "elaPeriod", label: "ELA cycle period", type: "range", min: 20, max: 400, step: 10, val: 100, unit: "kyr", group: "Glacier (asserted ice)", asserted: true, method: m.elaPeriod },
    { id: "buzzPeak", label: "Peak glacial rate", type: "range", min: 0, max: 20, step: 0.25, val: 1, unit: "mm/yr", group: "Glacier (asserted ice)", asserted: true, method: m.buzzPeak },
    { id: "buzzWidth", label: "Buzzsaw half-width", type: "range", min: 50, max: 1000, step: 25, val: 300, unit: "m", group: "Glacier (asserted ice)", asserted: true, method: m.buzzWidth },
    { id: "buzzAbove", label: "Erosion floor above ELA φ", type: "range", min: 0, max: 1, step: 0.05, val: 0, unit: "", group: "Glacier (asserted ice)", asserted: true, method: m.buzzAbove }
  ];
  var scenarioSpec = { id: "scenario", label: "Scenario", type: "select", val: "columns", options: SCEN, group: "",
    method: M("Scenarios", null, "<b>Two columns</b>: the illustration device (DESIGN.md §3) — a weak and a resistant column under one forcing, with an asserted elevation-dependent glacial rate. <b>Relief reservoir</b>: the legacy 1-D model, kept as a scenario.") };
  var allSpecs = [scenarioSpec].concat(common, landscape, body, glacier);
  var url = GS.ui.urlState(GS.ui.toUrlSpecs(allSpecs));
  var state = url.read();

  // ---- Panels ---------------------------------------------------------------
  var P = GS.ui.panels({ width: 820, height: 150 });
  var charts = d3.select("#charts");
  var tMyr = function (d) { return U.toMyr(d.t); };
  var mm = function (key) { return function (d) { return U.toMmyr(d[key]); }; };
  var km = function (key) { return function (d) { return U.toKm(d[key]); }; };
  var panels = [];

  function buildPanels() {
    charts.selectAll("*").remove();
    panels = [];
    if (state.scenario === "reservoir") {
      panels.push(P.make(charts, { title: "Uplift", yLabel: "U (mm/yr)", color: C.uplift, method: m.pUplift }));
      panels.push(P.make(charts, { title: "Relief", yLabel: "R (km)", color: C.relief, method: m.pRelief }));
      panels.push(P.make(charts, { title: "Erosion  vs  uplift flux", yLabel: "rate (mm/yr)", color: C.erosion, method: m.pErosion }));
      panels.push(P.make(charts, { title: "Net rate  (uplift − erosion = dR/dt)", yLabel: "net (mm/yr)", color: C.net, showX: true, xLabel: "time (Myr)", method: m.pNet }));
    } else {
      panels.push(P.make(charts, { title: "Uplift", yLabel: "U (mm/yr)", color: C.uplift, method: m.pUplift }));
      panels.push(P.make(charts, { title: "Elevation above base level  —  A weak · B resistant body", yLabel: "z (km)", color: C.B, method: m.pElev }));
      panels.push(P.make(charts, { title: "Erosion rates  vs  uplift", yLabel: "rate (mm/yr)", color: C.erosion, method: m.pRates }));
      panels.push(P.make(charts, { title: "Upland height  z_B − z_A", yLabel: "Δz (km)", color: C.diff, showX: true, xLabel: "time (Myr)", method: m.pDiff }));
    }
    d3.select("#legend").html(state.scenario === "reservoir"
      ? '<span style="color:#238b45"><i></i>relief R</span> <span style="color:#999"><i class="dash"></i>equilibrium R</span><br/>' +
        '<span style="color:#d94801"><i></i>erosion</span> <span style="color:#888"><i class="dash"></i>uplift flux</span>'
      : '<span style="color:' + C.A + '"><i></i>column A (weak)</span> <span style="color:' + C.B + '"><i></i>column B (resistant body)</span><br/>' +
        '<span style="color:' + C.B + '"><i style="border-top:8px solid ' + C.body + '"></i>resistant body</span> ' +
        '<span style="color:' + C.ela + '"><i class="dash"></i>ELA</span> <span style="color:#888"><i class="dash"></i>uplift</span>');
  }

  function renderReservoir() {
    var result = GS.d1.runReservoir({
      duration: U.fromMyr(state.duration), peakUplift: U.fromMmyr(state.peakUplift), exponent: state.exponent,
      reliefLimit: U.fromKm(state.reliefLimit), initialRelief: U.fromKm(state.initialRelief), shape: state.shape, steps: 2000
    });
    var data = result.series, last = data[data.length - 1], tMax = U.toMyr(last.t);
    var pU = panels[0], pR = panels[1], pE = panels[2], pN = panels[3];
    var uMax = d3.max(data, mm("uplift")) || 1;
    var rMax = d3.max(data, function (d) { return U.toKm(Math.max(d.relief, d.reliefEq)); }) || 1;
    var eMax = d3.max(data, function (d) { return U.toMmyr(Math.max(d.erosion, d.uplift)); }) || 1;
    var netAbs = Math.max(Math.abs(d3.min(data, mm("net"))), Math.abs(d3.max(data, mm("net"))), 1e-3);
    panels.forEach(function (pn) { pn.x.domain([0, tMax]); P.clear(pn); });
    pU.y.domain([0, uMax * 1.1]); pR.y.domain([0, rMax * 1.1]); pE.y.domain([0, eMax * 1.1]); pN.y.domain([-netAbs * 1.1, netAbs * 1.1]);
    P.line(pU, "uplift", data, mm("uplift"), C.uplift, false, tMyr);
    P.line(pR, "reliefEq", data, km("reliefEq"), C.eq, true, tMyr);
    P.line(pR, "relief", data, km("relief"), C.relief, false, tMyr);
    P.endLabel(pR, "relief", last, km("relief"), C.relief, "relief");
    P.endLabel(pR, "reliefEq", last, km("reliefEq"), C.eq, "equilibrium");
    P.line(pE, "flux", data, mm("uplift"), C.flux, true, tMyr);
    P.line(pE, "erosion", data, mm("erosion"), C.erosion, false, tMyr);
    P.zero(pN);
    P.line(pN, "net", data, mm("net"), C.net, false, tMyr);
    panels.forEach(function (pn) { P.vline(pn, "event-end", state.duration); P.axes(pn); });
    var peakR = d3.max(data, function (d) { return d.relief; });
    d3.select("#readout").html(
      "Peak relief reached: <b>" + U.fmtLen(peakR) + "</b>  ·  relief limit @ peak U: <b>" + U.fmtLen(U.fromKm(state.reliefLimit)) + "</b>" +
      "  ·  k = " + result.k.toExponential(2) + " m<sup>1−n</sup>/yr");
  }

  function renderColumns() {
    var result = GS.d1.runColumns({
      duration: U.fromMyr(state.duration), peakUplift: U.fromMmyr(state.peakUplift), shape: state.shape,
      exponent: state.exponent, reliefLimit: U.fromKm(state.reliefLimit), initialRelief: U.fromKm(state.initialRelief),
      bodyDepth: U.fromKm(state.bodyDepth), bodyThick: U.fromKm(state.bodyThick), contrast: state.contrast, contrastG: state.contrastG,
      glacierOn: !!state.glacierOn, elaBase: U.fromKm(state.elaBase), elaAmp: U.fromKm(state.elaAmp), elaPeriod: U.fromKyr(state.elaPeriod),
      buzzPeak: U.fromMmyr(state.buzzPeak), buzzWidth: state.buzzWidth, buzzAbove: state.buzzAbove, coupling: state.coupling, steps: 2000
    });
    var data = result.series, last = data[data.length - 1], tMax = U.toMyr(last.t);
    var pU = panels[0], pZ = panels[1], pE = panels[2], pD = panels[3];
    var uMax = d3.max(data, mm("uplift")) || 1;
    var zMax = d3.max(data, function (d) { return U.toKm(Math.max(d.zA, d.zB, state.glacierOn ? d.ela : 0)); }) || 1;
    var eMax = d3.max(data, function (d) { return U.toMmyr(Math.max(d.EA, d.EB, d.uplift)); }) || 1;
    var dMax = d3.max(data, function (d) { return Math.abs(U.toKm(d.diff)); }) || 0.1;
    panels.forEach(function (pn) { pn.x.domain([0, tMax]); P.clear(pn); });
    pU.y.domain([0, uMax * 1.1]); pZ.y.domain([0, zMax * 1.1]); pE.y.domain([0, eMax * 1.1]); pD.y.domain([Math.min(0, -dMax * 1.1), dMax * 1.1]);
    P.line(pU, "uplift", data, mm("uplift"), C.uplift, false, tMyr);
    if (state.bodyThick > 0) {
      var clip = function (v) { return Math.max(0, Math.min(zMax * 1.1, U.toKm(v))); };
      P.area(pZ, "body", data, function (d) { return clip(d.bodyBot); }, function (d) { return clip(d.bodyTop); }, C.body, tMyr);
    }
    if (state.glacierOn) P.line(pZ, "ela", data, km("ela"), C.ela, true, tMyr);
    P.line(pZ, "zA", data, km("zA"), C.A, false, tMyr);
    P.line(pZ, "zB", data, km("zB"), C.B, false, tMyr);
    P.endLabel(pZ, "zA", last, km("zA"), C.A, "A");
    P.endLabel(pZ, "zB", last, km("zB"), C.B, "B");
    P.line(pE, "flux", data, mm("uplift"), C.flux, true, tMyr);
    P.line(pE, "EA", data, mm("EA"), C.A, false, tMyr);
    P.line(pE, "EB", data, mm("EB"), C.B, false, tMyr);
    P.zero(pD);
    P.line(pD, "diff", data, km("diff"), C.diff, false, tMyr);
    panels.forEach(function (pn) { P.vline(pn, "event-end", state.duration); P.axes(pn); });

    // Readout: exposure / stripping times and the persistence comparison.
    var iExp = -1, iStrip = -1;
    for (var i = 0; i < data.length; i++) {
      if (iExp < 0 && data[i].rB > 0.5) iExp = i;
      else if (iExp >= 0 && iStrip < 0 && data[i].rB < 0.5) { iStrip = i; break; }
    }
    var maxDiff = d3.max(data, function (d) { return d.diff; });
    var html = "Upland height (max): <b>" + U.fmtLen(maxDiff) + "</b>, at end <b>" + U.fmtLen(last.diff) + "</b>";
    if (state.bodyThick > 0) {
      html += "  ·  body exposed: <b>" + (iExp >= 0 ? U.fmtTime(data[iExp].t) : "never") + "</b>";
      html += "  ·  stripped: <b>" + (iStrip >= 0 ? U.fmtTime(data[iStrip].t) : (iExp >= 0 ? "survives" : "—")) + "</b>";
      if (iStrip >= 0 && iExp >= 0) html += "  (lifetime " + U.fmtTime(data[iStrip].t - data[iExp].t) + " vs event " + U.fmtTime(U.fromMyr(state.duration)) + ")";
    }
    d3.select("#readout").html(html);
  }

  function render() {
    if (state.scenario === "reservoir") renderReservoir(); else renderColumns();
    url.write(state);
  }

  function buildSidebar() {
    var root = d3.select("#controls");
    root.selectAll("*").remove();
    var specs = state.scenario === "reservoir" ? [scenarioSpec].concat(common, landscape) : allSpecs;
    GS.ui.buildControls(root, specs, state, function (id) {
      if (id === "scenario") { buildSidebar(); buildPanels(); }
      render();
    });
  }

  buildSidebar();
  buildPanels();
  render();
})();
