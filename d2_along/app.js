/*
 * d2_along/app.js — controls, animation loop and D3 panels for the along-valley tier.
 * Depends on core/*.js, ui/*.js (window.GS), d2_along/model.js and d3 v7.
 * Display units: km, mm/yr, kyr/Myr; the model runs in m and yr.
 * Every control and panel carries a `method` tooltip with its equations.
 */
(function () {
  "use strict";
  var GS = window.GS, U = GS.units;
  var C = { bed: "#5a3d1e", ice: "#4292c6", iceFill: "rgba(66,146,198,0.25)", ela: "#4292c6", litho: "#8c2d04",
            uplift: "#888", fluv: "#2c7fb8", glac: "#d94801", Us: "#6a51a3", relief: "#238b45", vol: "#4292c6" };

  function M(title, eqs, text, src, tag) {
    var h = "<h4>" + (tag ? '<span class="tag ' + tag + '">' + tag + "</span>" : "") + title + "</h4>";
    (eqs || []).forEach(function (e) { h += '<span class="eq">' + e + "</span>"; });
    if (text) h += "<div>" + text + "</div>";
    if (src) h += '<div class="src">' + src + "</div>";
    return h;
  }
  var m = {
    L: M("Valley length L", ["s ∈ [0, L], node spacing ds = L/(N − 1)"], "s = 0 is the divide (valley head), s = L the outlet at fixed base level z = 0.", null, "asserted"),
    N: M("Nodes N", ["ds = L / (N − 1)"], "Resolution of the profile. 301 nodes on 30 km is 100 m.", null, "asserted"),
    initProfile: M("Initial profile", ["steady: S = (U<sub>peak</sub> f / K A<sup>m</sup>)<sup>1/n</sup> integrated up from the outlet", "linear: z = z<sub>head</sub>(1 − s/L)", "concave: z = z<sub>head</sub>(1 − s/L)<sup>1.5</sup>"],
      "The fluvial steady state (default) starts the river in equilibrium with the peak uplift, so anything that then happens is the glacier's doing. The other two are far from equilibrium and incise fast at first. Noise adds seeded Gaussian roughness. Head elevation applies to linear and concave only.", null, "asserted"),
    zHead: M("Head elevation", ["z(0, 0) = z<sub>head</sub>"], null, null, "asserted"),
    noise: M("Profile noise σ", ["z += σ·𝒩(0,1) (seeded)"], null, null, "asserted"),
    shape: M("Uplift forcing U(t)", ["U(t) = U<sub>peak</sub> s(t)"], "Event shapes on [0, T]: bell, plateau, pulse, arc, triangle, square. The outlet does not uplift (fixed base level).", "core/forcing.js", "asserted"),
    peakUplift: M("Peak rock-uplift rate", ["U<sub>peak</sub>"], null, null, "asserted"),
    duration: M("Event duration T", ["s(t) = 0 for t > T"], null, null, "asserted"),
    pattern: M("Spatial uplift pattern", ["U(s, t) = U(t) · f(s/L)"], "uniform 1; ramp rises from the ends; gaussian peaks mid-profile; tilt rises toward the outlet (f = s/L) — or toward the head when combined with a low head.", "core/forcing.js spatialFactor", "asserted"),
    KExp: M("Fluvial erodibility K", ["E<sub>f</sub> = K A<sup>m</sup> S<sup>n</sup>  (× lithology factor)", "steady state: S = (U / K A<sup>m</sup>)<sup>1/n</sup>"],
      "Detachment-limited stream power, solved implicitly node by node from the outlet upstream (Braun & Willett 2013; Yuan et al. 2019 with deposition). Only ice-free nodes are fluvial.", "core/laws.js yuanNode; validated: steady slope–area law exact, closure to round-off", "computed"),
    m: M("Area exponent m", ["E<sub>f</sub> ∝ A<sup>m</sup>"], "Concavity θ = m/n with Hack's law gives S ∝ s<sup>−h m/n</sup>.", null, "asserted"),
    nexp: M("Slope exponent n", ["E<sub>f</sub> ∝ S<sup>n</sup>"], "n ≠ 1 uses Newton iteration inside the implicit update.", null, "asserted"),
    G: M("Deposition coefficient G", ["dz/dt = −K A<sup>m</sup> S<sup>n</sup> + G Q<sub>s</sub>/A"], "G = 0: pure detachment-limited. G > 0: eroded sediment (volume flux Q<sub>s</sub>) is re-deposited downstream; iterated Gauss–Seidel, mass-conserving by construction.", "Yuan et al. 2019", "asserted"),
    hack: M("Hack exponent h", ["A(s) = k<sub>a</sub> (s + s<sub>0</sub>)<sup>h</sup>"], "The profile does not resolve its basin: drainage area is asserted from Hack's law L = 1.4 A<sup>0.6</sup> (h = 1/0.6, k<sub>a</sub> = 5.7 in m units).", "Hack 1957", "asserted"),
    W0: M("Valley width at the head W<sub>0</sub>", ["W(s) = W<sub>0</sub> + k<sub>w</sub> s"], "Ice flux Q = ∫ b W ds and the flux per unit width q = Q/W both need a width; the profile cannot resolve it.", null, "asserted"),
    kw: M("Valley width growth k<sub>w</sub>", ["W(s) = W<sub>0</sub> + k<sub>w</sub> s"], null, null, "asserted"),
    glacierOn: M("Glacier: steady ice discharge", [
        "b(z<sub>s</sub>) = β (z<sub>s</sub> − ELA) above, 2.5 β (z<sub>s</sub> − ELA) below, clamped",
        "Q(s) = ∫<sub>0</sub><sup>s</sup> b W ds, clamped ≥ 0 → terminus",
        "Γ′ H<sub>f</sub><sup>n+2</sup> S<sub>f</sub><sup>n</sup> = Q/W<sub>f</sub> at each face, Γ′ = Γ/(1 − f<sub>s</sub>), n = 3",
        "U<sub>s</sub> = f<sub>s</sub> q / H"],
      "Ice equilibrates instantly on geomorphic timescales, so each step solves the steady glacier: balance on the ice surface → routed flux → thickness by marching upstream from the outlet, solving the face flux law for H at each node (monotone → unique root). A few outer iterations couple surface and balance. Validated against an explicit flowline SIA run to steady state: thickness within 2 %, same terminus, on linear, concave and bumped beds.",
      "SIA: Glen n = 3, A = 7.57×10<sup>−17</sup> Pa<sup>−3</sup> yr<sup>−1</sup>; flowline models: MacGregor et al. 2000, Anderson et al. 2006", "computed"),
    elaBase: M("ELA (interglacial)", ["ELA(t) = ELA<sub>0</sub> − A<sub>ELA</sub> · ½(1 − cos 2πt/P)"], "Referenced to the ice surface z + H (surface feedback kept, legacy decision).", null, "asserted"),
    elaAmp: M("ELA cycle amplitude", ["ELA<sub>min</sub> = ELA<sub>0</sub> − A<sub>ELA</sub>"], null, null, "asserted"),
    elaPeriod: M("ELA cycle period", ["P"], null, null, "asserted"),
    balGrad: M("Balance gradient β", ["b = β (z<sub>s</sub> − ELA) above the ELA; 2.5 β below"], "Accumulation clamps at 2 m/yr, ablation at 8 m/yr.", "asymmetric gradient after legacy terrain_sandbox", "asserted"),
    fs: M("Sliding fraction f<sub>s</sub>", ["U<sub>s</sub> = f<sub>s</sub> q/H,   q<sub>deformation</sub> = (1 − f<sub>s</sub>) q"], "Share of the depth-averaged speed carried by basal sliding. Erosion follows sliding, so f<sub>s</sub> = 0 gives no glacial erosion.", null, "asserted"),
    flow: M("Flow-law enhancement", ["Γ → Γ · E"], "Multiplies the Glen rate factor (softer ice → thinner, faster).", null, "asserted"),
    KgExp: M("Glacial erodibility K<sub>g</sub>", ["E<sub>g</sub> = K<sub>g</sub> U<sub>s</sub><sup>l</sup>  (× lithology factor), capped at 2 cm/yr"],
      "Abrasion proportional to sliding speed (l = 1). Self-limiting: an overdeepening flattens the surface, thickens the ice and slows sliding. Validated: K<sub>g</sub> × 100 does not run away — it cuts the bed below the ELA and the glacier shuts off.", "Hallet 1979; Humphrey & Raymond 1994", "computed"),
    lexp: M("Sliding exponent l", ["E<sub>g</sub> ∝ U<sub>s</sub><sup>l</sup>"], "1 (Humphrey & Raymond) to 2 (Hallet).", null, "asserted"),
    KqExp: M("Quarrying on convexity K<sub>q</sub>", ["E<sub>q</sub> = K<sub>q</sub> U<sub>s</sub> · max(−∂²z/∂s², 0)"], "The hypothesis-A knob: extra erosion where the bed is convex-up under sliding ice. Off by default (slider at minimum).", "DESIGN.md §4", "asserted"),
    lithoType: M("Resistant body (material frame)", ["z<sub>m</sub> = z − U<sub>cum</sub>(s)", "layer: z<sub>m</sub> ∈ [top − T, top];  dike: |s − s<sub>0</sub>| < w/2;  slab: top + dip·s"],
      "Strength lives on material coordinates and rides up with uplift; the surface value r(s) multiplies K and K<sub>g</sub> by 1/(1 + r(c − 1)). Drawn as a brown overlay where the surface is in the body.", "core/lithology.js", "asserted"),
    lithoTop: M("Body top (material height)", ["top (m above base level in the initial frame)"], null, null, "asserted"),
    lithoThick: M("Body thickness / dike width", ["T (layer, slab) or w (dike)"], null, null, "asserted"),
    lithoPos: M("Dike position / slab dip", ["dike: s<sub>0</sub> = x·L;  slab: dip = x (m/m)"], null, null, "asserted"),
    contrastK: M("Fluvial contrast c", ["K → K / (1 + r(c − 1))"], null, null, "asserted"),
    contrastKg: M("Glacial contrast c<sub>g</sub>", ["K<sub>g</sub> → K<sub>g</sub> / (1 + r(c<sub>g</sub> − 1))"], null, null, "asserted"),
    dt: M("Time step", ["per step: lithology → uplift → steady ice + erosion → fluvial"], "The fluvial solve is implicit and the ice is steady, so dt is limited only by how fast the bed changes. Steady profiles are dt-independent to 0.01 m; a 200 kyr transient differs by ~15 m between 250 and 1000 yr.", null, "asserted"),
    speed: M("Steps per frame", null, "Animation speed only.", null, null),
    pProfile: M("Long profile", ["bed z(s), ice surface z + H, ELA(t)"], "Brown: bedrock. Blue fill: ice (steady discharge). Dashed: ELA. Thick brown overlay: surface within the resistant body. The outlet is fixed base level.", null, "computed"),
    pRates: M("Erosion rates", ["E<sub>f</sub> = K A<sup>m</sup> S<sup>n</sup> (ice-free);  E<sub>g</sub> = K<sub>g</sub> U<sub>s</sub><sup>l</sup> (+ quarrying) under ice;  U(s, t) dashed"], null, null, "computed"),
    pIce: M("Ice: sliding speed and thickness", ["U<sub>s</sub> = f<sub>s</sub> q/H (m/yr);  H (m, right axis not drawn — see profile)"], "Sliding speed drives erosion. Where H = 0 there is no ice.", null, "computed"),
    pHist: M("History", ["relief = max z − min z;  ice volume = Σ H W ds"], "Time series of the run. Mass closure (uplift − erosion − Δvolume) is in the readout.", null, "computed"),
    pLink: M("Isolated column at the selected node  (d1 ↔ d2 link)", [
        "frozen at the click: z<sub>down</sub>, H<sub>down</sub>, Q, A, W, W<sub>f</sub>, ds",
        "running: U(t)·f(s), lithology r(z − U<sub>cum</sub>), slope against z<sub>down</sub>, H from the face law against z<sub>down</sub> + H<sub>down</sub>",
        "divergence = z<sub>2D</sub>(t) − z<sub>1D</sub>(t)"],
      "Click the long profile to pick a node. Solid: the node's real elevation as the profile runs. Dashed: the same node run as an isolated d1 column with the same laws, where everything that came from its neighbours is held at the click-time value and badged asserted. The two curves start together; how far and how fast they part is the answer to “how much behaviour arises from the extra dimension”. Note the ELA acts on a node only through the ice flux, which is adjacency — so the isolated column is blind to glacial cycles. “Open in d1” shows this column standalone with its constants as sliders.",
      "py/test_d1_node.py: the column reproduces the 2D node exactly when its adjacency is refreshed every step.", "computed")
  };

  // ---- Controls (display units) ------------------------------------------
  var SHAPES = [{ v: "constant", t: "Constant (square pulse)" }, { v: "bell", t: "Bell curve" }, { v: "plateau", t: "Ramp · hold · decay" }, { v: "pulse", t: "Pulse" }, { v: "arc", t: "Circular arc" }, { v: "triangle", t: "Triangle" }];
  var specs = [
    { id: "L", label: "Valley length", type: "range", min: 10, max: 60, step: 5, val: 30, unit: "km", group: "Domain", asserted: true, method: m.L },
    { id: "N", label: "Nodes", type: "select", val: "301", options: [{ v: "151", t: "151" }, { v: "301", t: "301" }, { v: "601", t: "601" }], group: "Domain", asserted: true, method: m.N },
    { id: "initProfile", label: "Initial profile", type: "select", val: "steady", options: [{ v: "steady", t: "Fluvial steady state (for U, K)" }, { v: "concave", t: "Concave (1 − s/L)^1.5" }, { v: "linear", t: "Linear" }], group: "Domain", asserted: true, method: m.initProfile },
    { id: "zHead", label: "Head elevation", type: "range", min: 0.5, max: 5, step: 0.1, val: 2.5, unit: "km", group: "Domain", asserted: true, method: m.zHead },
    { id: "noise", label: "Profile noise", type: "range", min: 0, max: 50, step: 1, val: 0, unit: "m", group: "Domain", asserted: true, method: m.noise },
    { id: "shape", label: "Uplift shape", type: "select", val: "constant", options: SHAPES, group: "Forcing", method: m.shape },
    { id: "peakUplift", label: "Peak uplift", type: "range", min: 0, max: 10, step: 0.1, val: 1, unit: "mm/yr", group: "Forcing", method: m.peakUplift },
    { id: "duration", label: "Event duration", type: "range", min: 0.5, max: 40, step: 0.5, val: 10, unit: "Myr", group: "Forcing", method: m.duration },
    { id: "pattern", label: "Spatial pattern", type: "select", val: "uniform", options: [{ v: "uniform", t: "Uniform" }, { v: "ramp", t: "Ramp from ends" }, { v: "gaussian", t: "Gaussian mid-profile" }, { v: "tilt", t: "Tilt (rises to outlet)" }], group: "Forcing", asserted: true, method: m.pattern },
    { id: "KExp", label: "Erodibility K", type: "range", min: -7, max: -4, step: 0.1, val: -5.5, log: true, group: "Fluvial", method: m.KExp },
    { id: "m", label: "Area exponent m", type: "range", min: 0.3, max: 0.7, step: 0.05, val: 0.5, unit: "", group: "Fluvial", asserted: true, method: m.m },
    { id: "nexp", label: "Slope exponent n", type: "range", min: 1, max: 2, step: 0.25, val: 1, unit: "", group: "Fluvial", asserted: true, method: m.nexp },
    { id: "G", label: "Deposition G", type: "range", min: 0, max: 2, step: 0.1, val: 0, unit: "", group: "Fluvial", asserted: true, method: m.G },
    { id: "hack", label: "Hack exponent h", type: "range", min: 1.4, max: 2.0, step: 0.05, val: 1.667, unit: "", group: "Fluvial", asserted: true, method: m.hack },
    { id: "W0", label: "Valley width at head", type: "range", min: 50, max: 1000, step: 50, val: 300, unit: "m", group: "Valley geometry", asserted: true, method: m.W0 },
    { id: "kw", label: "Width growth", type: "range", min: 0, max: 0.2, step: 0.01, val: 0.05, unit: "m/m", group: "Valley geometry", asserted: true, method: m.kw },
    { id: "glacierOn", label: "Glacier on", type: "check", val: true, group: "Glacier", method: m.glacierOn },
    { id: "elaBase", label: "ELA (interglacial)", type: "range", min: 0.2, max: 5, step: 0.1, val: 1.5, unit: "km", group: "Glacier", asserted: true, method: m.elaBase },
    { id: "elaAmp", label: "ELA cycle amplitude", type: "range", min: 0, max: 1.5, step: 0.05, val: 0, unit: "km", group: "Glacier", asserted: true, method: m.elaAmp },
    { id: "elaPeriod", label: "ELA cycle period", type: "range", min: 20, max: 400, step: 10, val: 100, unit: "kyr", group: "Glacier", asserted: true, method: m.elaPeriod },
    { id: "balGrad", label: "Balance gradient", type: "range", min: 0.002, max: 0.015, step: 0.001, val: 0.007, unit: "/yr", group: "Glacier", asserted: true, method: m.balGrad },
    { id: "fs", label: "Sliding fraction", type: "range", min: 0, max: 0.95, step: 0.05, val: 0.5, unit: "", group: "Glacier", asserted: true, method: m.fs },
    { id: "flow", label: "Flow enhancement", type: "range", min: 0.2, max: 5, step: 0.1, val: 1, unit: "×", group: "Glacier", asserted: true, method: m.flow },
    { id: "KgExp", label: "Glacial erodibility K_g", type: "range", min: -6, max: -2, step: 0.1, val: -4, log: true, group: "Glacial erosion", method: m.KgExp },
    { id: "lexp", label: "Sliding exponent l", type: "range", min: 1, max: 2, step: 0.25, val: 1, unit: "", group: "Glacial erosion", asserted: true, method: m.lexp },
    { id: "KqExp", label: "Quarrying K_q (min = off)", type: "range", min: -6, max: -1, step: 0.25, val: -6, log: true, group: "Glacial erosion", asserted: true, method: m.KqExp },
    { id: "lithoType", label: "Body type", type: "select", val: "none", options: [{ v: "none", t: "None" }, { v: "layer", t: "Horizontal layer" }, { v: "slab", t: "Dipping slab" }, { v: "dike", t: "Vertical dike" }], group: "Lithology", asserted: true, method: m.lithoType },
    { id: "lithoTop", label: "Body top", type: "range", min: 0, max: 5, step: 0.1, val: 1.5, unit: "km", group: "Lithology", asserted: true, method: m.lithoTop },
    { id: "lithoThick", label: "Thickness / width", type: "range", min: 0.05, max: 3, step: 0.05, val: 0.3, unit: "km", group: "Lithology", asserted: true, method: m.lithoThick },
    { id: "lithoPos", label: "Dike position / slab dip", type: "range", min: 0, max: 1, step: 0.05, val: 0.5, unit: "", group: "Lithology", asserted: true, method: m.lithoPos },
    { id: "contrastK", label: "Fluvial contrast", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Lithology", asserted: true, method: m.contrastK },
    { id: "contrastKg", label: "Glacial contrast", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Lithology", asserted: true, method: m.contrastKg },
    { id: "dt", label: "Time step", type: "range", min: 100, max: 2000, step: 100, val: 500, unit: "yr", group: "Run", asserted: true, method: m.dt },
    { id: "speed", label: "Steps per frame", type: "range", min: 1, max: 20, step: 1, val: 4, unit: "", group: "Run", method: m.speed }
  ];
  var url = GS.ui.urlState(GS.ui.toUrlSpecs(specs));
  var state = url.read();

  // ---- Model ---------------------------------------------------------------
  var model = null;
  function lithoSpec() {
    return GS.d2along.lithoFromControls(state.lithoType, U.fromKm(state.lithoTop), U.fromKm(state.lithoThick), state.lithoPos, U.fromKm(state.L));
  }
  function build() {
    model = GS.d2along.createModel({
      L: U.fromKm(state.L), N: +state.N, initProfile: state.initProfile, zHead: U.fromKm(state.zHead), noise: state.noise, seed: 1,
      shape: state.shape, peakUplift: U.fromMmyr(state.peakUplift), duration: U.fromMyr(state.duration), pattern: state.pattern,
      K: Math.pow(10, state.KExp), m: state.m, nexp: state.nexp, G: state.G, hack: state.hack, W0: state.W0, kw: state.kw,
      glacierOn: !!state.glacierOn, elaBase: U.fromKm(state.elaBase), elaAmp: U.fromKm(state.elaAmp), elaPeriod: U.fromKyr(state.elaPeriod),
      balGrad: state.balGrad, fs: state.fs, flow: state.flow, Kg: Math.pow(10, state.KgExp), lexp: state.lexp,
      Kq: state.KqExp <= -6 ? 0 : Math.pow(10, state.KqExp), litho: lithoSpec(), contrastK: state.contrastK, contrastKg: state.contrastKg,
      dt: state.dt
    });
    model.primeIce();
    model.record();
  }

  // ---- Panels --------------------------------------------------------------
  var P = GS.ui.panels({ width: 820, height: 170 });
  var charts = d3.select("#charts");
  var pProfile = P.make(charts, { title: "Long profile", yLabel: "z (km)", color: C.bed, method: m.pProfile });
  var pRates = P.make(charts, { title: "Erosion rates  vs  uplift", yLabel: "rate (mm/yr)", color: C.glac, method: m.pRates });
  var pIce = P.make(charts, { title: "Sliding speed", yLabel: "U_s (m/yr)", color: C.Us, showX: true, xLabel: "distance from divide s (km)", method: m.pIce });
  var pHist = P.make(charts, { title: "History", yLabel: "relief (km)", color: C.relief, showX: true, xLabel: "time (Myr)", method: m.pHist });
  var pLink = P.make(charts, { title: "Isolated column — click the profile to pick a node", yLabel: "z (km)", color: C.Us, showX: true, xLabel: "time since click (kyr)", method: m.pLink });
  var sKm = function (d) { return U.toKm(d.s); };
  var link = null; // { i, column, series: [{t, z2d, z1d, H2d, H1d}] }
  // Axis ratchet: a panel's range may grow but never shrink during a run (cleared on Reset),
  // so a transient spike (e.g. at the glacier front) does not make the axes jump.
  var seen = {};
  function ratchetMax(key, v) { if (seen[key] == null || v > seen[key]) seen[key] = v; return seen[key]; }
  function ratchetMin(key, v) { if (seen[key] == null || v < seen[key]) seen[key] = v; return seen[key]; }

  function linkNode(i) {
    delete seen.linkLo; delete seen.linkHi;
    var spec = model.nodeSpec(i);
    link = { i: spec.i, column: GS.d1.createNodeColumn(spec), series: [] };
    recordLink();
  }
  function recordLink() {
    if (!link) return;
    var st = model.state, c = link.column.state, i = link.i;
    link.series.push({ t: st.t - link.column.spec.t0, z2d: st.z[i], z1d: c.z, H2d: st.H[i], H1d: c.H, E2d: st.Ef[i] + st.Eg[i], E1d: c.Ef + c.Eg });
  }
  function stepLink(nsteps) {
    if (!link) return;
    for (var k = 0; k < nsteps; k++) link.column.step();
    recordLink();
  }
  function d1Link() {
    if (!link) return "#";
    var sp = link.column.spec, q = new URLSearchParams();
    q.set("scenario", "node");
    ["shape", "peakUplift", "duration", "KExp", "m", "nexp", "KgExp", "lexp", "fs", "flow", "contrastK", "contrastKg", "lithoType", "lithoTop", "lithoThick", "lithoPos", "L", "elaBase", "elaAmp", "elaPeriod", "dt"].forEach(function (id) { q.set(id, state[id]); });
    q.set("s", U.toKm(sp.s).toFixed(3)); q.set("z0", U.toKm(sp.z0).toFixed(4)); q.set("zDown", U.toKm(sp.zDown).toFixed(4)); q.set("HDown", sp.HDown.toFixed(2));
    q.set("QExp", sp.Q > 0 ? Math.log10(sp.Q).toFixed(3) : 0); q.set("AExp", Math.log10(sp.A / 1e6).toFixed(3)); q.set("W", sp.W.toFixed(1)); q.set("Wf", sp.Wf.toFixed(1));
    q.set("ds", sp.ds.toFixed(1)); q.set("fac", sp.fac.toFixed(4)); q.set("ucum0", sp.ucum0.toFixed(3)); q.set("t0", sp.t0.toFixed(0)); q.set("runKyr", 500);
    return "../d1/index.html?" + q.toString();
  }

  function draw() {
    var st = model.state, N = st.N, i;
    var rows = [];
    for (i = 0; i < N; i++) rows.push({ s: st.s[i], z: st.z[i], zs: st.z[i] + st.H[i], H: st.H[i], Us: st.Us[i], Ef: st.Ef[i], Eg: st.Eg[i], r: st.r[i], ela: model.ELA(st.t), u: model.U(st.t) });
    var Lkm = U.toKm(st.p.L);
    var zMax = ratchetMax("z", d3.max(rows, function (d) { return U.toKm(Math.max(d.zs, state.glacierOn ? d.ela : 0)); }) || 1);
    var eMax = ratchetMax("e", d3.max(rows, function (d) { return U.toMmyr(Math.max(d.Ef, d.Eg, d.u)); }) || 0.1);
    var usMax = ratchetMax("us", d3.max(rows, function (d) { return d.Us; }) || 1);
    [pProfile, pRates, pIce].forEach(function (pn) { pn.x.domain([0, Lkm]); P.clear(pn); });
    pProfile.y.domain([0, zMax * 1.05]); pRates.y.domain([0, eMax * 1.1]); pIce.y.domain([0, usMax * 1.1]);
    // profile
    P.area(pProfile, "ice", rows, function (d) { return U.toKm(d.z); }, function (d) { return U.toKm(d.zs); }, C.iceFill, sKm);
    if (state.glacierOn) P.line(pProfile, "ela", rows, function (d) { return U.toKm(d.ela); }, C.ela, true, sKm);
    P.line(pProfile, "ice", rows, function (d) { return U.toKm(d.zs); }, C.ice, false, sKm);
    P.line(pProfile, "bed", rows, function (d) { return U.toKm(d.z); }, C.bed, false, sKm);
    // lithology overlay: segments where r > 0.5
    var segs = [], cur = null;
    rows.forEach(function (d) { if (d.r > 0.5) { if (!cur) { cur = []; segs.push(cur); } cur.push(d); } else cur = null; });
    var lsel = pProfile.g.selectAll("path.line.litho").data(segs);
    lsel.enter().append("path").attr("class", "line litho").attr("fill", "none").attr("stroke", C.litho).attr("stroke-width", 5).attr("stroke-opacity", 0.6)
      .merge(lsel).attr("d", d3.line().x(function (d) { return pProfile.x(sKm(d)); }).y(function (d) { return pProfile.y(U.toKm(d.z)); }));
    lsel.exit().remove();
    // rates
    P.line(pRates, "u", rows, function (d) { return U.toMmyr(d.u); }, C.uplift, true, sKm);
    P.line(pRates, "Ef", rows, function (d) { return U.toMmyr(d.Ef); }, C.fluv, false, sKm);
    P.line(pRates, "Eg", rows, function (d) { return U.toMmyr(d.Eg); }, C.glac, false, sKm);
    // ice
    P.line(pIce, "Us", rows, function (d) { return d.Us; }, C.Us, false, sKm);
    if (link) [pProfile, pRates, pIce].forEach(function (pn) { P.vline(pn, "node", U.toKm(st.s[link.i])); });
    [pProfile, pRates, pIce].forEach(function (pn) { P.axes(pn); });
    // linked column overlay
    P.clear(pLink);
    if (link && link.series.length > 1) {
      var ls = link.series, tk = function (d) { return U.toKyr(d.t); };
      pLink.x.domain([0, Math.max(U.toKyr(ls[ls.length - 1].t), 1)]);
      var lo = ratchetMin("linkLo", d3.min(ls, function (d) { return U.toKm(Math.min(d.z2d, d.z1d)); }));
      var hi = ratchetMax("linkHi", d3.max(ls, function (d) { return U.toKm(Math.max(d.z2d + d.H2d, d.z1d + d.H1d)); }));
      pLink.y.domain([lo - 0.02, hi + 0.02]);
      P.line(pLink, "ice2d", ls, function (d) { return U.toKm(d.z2d + d.H2d); }, C.ice, false, tk);
      P.line(pLink, "ice1d", ls, function (d) { return U.toKm(d.z1d + d.H1d); }, C.ice, true, tk);
      P.line(pLink, "z2d", ls, function (d) { return U.toKm(d.z2d); }, C.bed, false, tk);
      P.line(pLink, "z1d", ls, function (d) { return U.toKm(d.z1d); }, C.Us, true, tk);
      P.endLabel(pLink, "z2d", ls[ls.length - 1], function (d) { return U.toKm(d.z2d); }, C.bed, "profile");
      P.endLabel(pLink, "z1d", ls[ls.length - 1], function (d) { return U.toKm(d.z1d); }, C.Us, "isolated");
      P.axes(pLink);
      var last = ls[ls.length - 1];
      d3.select("#linkReadout").html("Node at <b>" + U.fmtLen(st.s[link.i]) + "</b> · after " + U.fmtTime(last.t) + ": profile z = <b>" + U.fmtLen(last.z2d) + "</b>, isolated z = <b>" + U.fmtLen(last.z1d) +
        "</b> · divergence <b>" + (last.z2d - last.z1d).toFixed(1) + " m</b> · erosion now " + U.fmtRate(last.E2d) + " vs " + U.fmtRate(last.E1d) +
        " · ice " + last.H2d.toFixed(0) + " vs " + last.H1d.toFixed(0) + " m · <a href='" + d1Link() + "' target='_blank'>open this column in d1 ↗</a>");
    } else {
      pLink.x.domain([0, 1]); pLink.y.domain([0, 1]); P.axes(pLink);
      d3.select("#linkReadout").html(link ? "Node linked — run to see the curves." : "Click a point on the long profile to isolate that node as a d1 column and overlay the two.");
    }
    // history
    var h = st.history, tMyr = function (d) { return U.toMyr(d.t); };
    P.clear(pHist);
    pHist.x.domain([0, Math.max(U.toMyr(h[h.length - 1].t), 0.01)]);
    var rMax = ratchetMax("relief", d3.max(h, function (d) { return U.toKm(d.relief); }) || 1), vMax = ratchetMax("vol", d3.max(h, function (d) { return d.iceVol; }) || 1);
    pHist.y.domain([0, rMax * 1.1]);
    P.line(pHist, "relief", h, function (d) { return U.toKm(d.relief); }, C.relief, false, tMyr);
    P.line(pHist, "vol", h, function (d) { return rMax * 1.1 * d.iceVol / (vMax * 1.1); }, C.vol, true, tMyr);
    P.axes(pHist);
    var dg = model.diagnostics();
    d3.select("#readout").html(
      "t = <b>" + U.fmtTime(dg.t) + "</b>  ·  relief <b>" + U.fmtLen(dg.relief) + "</b>  ·  max ice <b>" + dg.maxH.toFixed(0) + " m</b>" +
      (dg.terminus >= 0 ? "  ·  terminus <b>" + U.fmtLen(dg.terminus) + "</b>" : "  ·  no ice") +
      "  ·  ice volume " + (dg.iceVol / 1e9).toFixed(2) + " km³  ·  ice iters " + st.iceIters +
      "  ·  mass closure " + (dg.upliftVol > 0 ? (dg.closure / dg.upliftVol).toExponential(1) : "—"));
  }

  // ---- Loop --------------------------------------------------------------
  var playing = false, raf = null;
  function frame() {
    if (!playing) return;
    for (var k = 0; k < state.speed; k++) model.step();
    model.record(); stepLink(state.speed);
    draw();
    raf = requestAnimationFrame(frame);
  }
  function setPlaying(v) {
    playing = v;
    d3.select("#playPause").text(playing ? "❚❚ Pause" : "▶ Play");
    if (playing) raf = requestAnimationFrame(frame); else if (raf) cancelAnimationFrame(raf);
  }
  function reset() { setPlaying(false); link = null; seen = {}; build(); draw(); url.write(state); }
  d3.select("#playPause").on("click", function () { setPlaying(!playing); });
  d3.select("#stepBtn").on("click", function () { setPlaying(false); for (var k = 0; k < state.speed; k++) model.step(); model.record(); stepLink(state.speed); draw(); });
  // click tool on the profile: pick the node under the pointer
  pProfile.svg.style("cursor", "crosshair").on("click", function (ev) {
    var xy = d3.pointer(ev, pProfile.g.node());
    var sKmClicked = pProfile.x.invert(xy[0]);
    var i = Math.round(U.fromKm(sKmClicked) / model.state.ds);
    linkNode(i); draw();
  });
  d3.select("#resetBtn").on("click", reset);

  GS.ui.buildControls(d3.select("#controls"), specs, state, function (id) {
    // domain / initial-condition changes need a rebuild; everything else applies live via rebuild-in-place of params
    if (["L", "N", "initProfile", "zHead", "noise"].indexOf(id) !== -1 || (state.initProfile === "steady" && model.state.t === 0 && ["KExp", "m", "nexp", "peakUplift", "pattern", "hack"].indexOf(id) !== -1)) { reset(); return; }
    var wasPlaying = playing;
    link = null; // a parameter change invalidates the frozen column (click again to re-link)
    var keepT = model.state.t, keepZ = model.state.z, keepU = model.state.ucum, keepH = model.state.H, keepHist = model.state.history;
    build();
    // carry the evolving state across a parameter change (same grid)
    model.state.z.set(keepZ); model.state.ucum.set(keepU); model.state.H.set(keepH); model.state.t = keepT; model.state.history = keepHist;
    draw(); url.write(state);
    if (wasPlaying) setPlaying(true);
  });
  window.GS_d2along_model = function () { return model; }; // for validation from the console
  build(); draw();
})();
