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
      "The time history of rock uplift, the engine of every scenario here. s(t) ∈ [0, 1] is the event shape on [0, T]: bell (Gaussian, σ = T/6), plateau (25 % ramp, 50 % hold, 25 % ramp), pulse (10 % rise, exponential decay τ = 0.3 T), arc, triangle or square pulse. The run continues to 2 T so you see the relaxation after the event as well as the response during it.",
      "core/forcing.js", "asserted"),
    peakUplift: M("Peak rock-uplift rate", ["U<sub>peak</sub>  (mm/yr = km/Myr = 10<sup>−3</sup> m/yr)"],
      "How fast rock rises past base level at the height of the event, applied identically to every column. Faster uplift raises the equilibrium relief (as U<sup>1/n</sup>) and shortens the time a resistant body survives, since erosion has to keep pace. Active mountain belts sit at 0.5–10 mm/yr.", "core/forcing.js", "asserted"),
    duration: M("Event duration T", ["s(t) = 0 for t > T (except the square pulse: on until T)"],
      "How long the uplift lasts. The dashed vertical line marks its end; after it relief decays under erosion alone. Whether a resistant body outlasts T is the hypothesis-B question.", null, "asserted"),
    reliefLimit: M("Relief limit R<sub>lim</sub> → erodibility k", ["k = U<sub>peak</sub> / R<sub>lim</sub><sup>n</sup>", "R<sub>eq</sub>(U) = (U / k)<sup>1/n</sup>"],
      "The relief the landscape would settle at under sustained peak uplift: an intuitive way to set k, the constant of the erosion law. A high relief limit means an inefficient landscape (weak rivers, dry climate) that must stand tall to erode as fast as it rises; a low one means an efficient landscape that stays low. The column does not resolve the drainage network that actually sets k, so this is asserted.",
      SRC_RES, "asserted"),
    exponent: M("Erosion nonlinearity n", ["E(R) = k R<sup>n</sup>"],
      "How steeply erosion grows with relief. n = 1 makes relief respond in proportion to uplift; large n makes erosion explode as relief grows, so the equilibrium relief hardly changes when uplift changes (R<sub>eq</sub> ∝ U<sup>1/n</sup>) and the landscape behaves as if it had a hard ceiling. The legacy default of 4 gives that ceiling.", SRC_RES, "asserted"),
    initialRelief: M("Initial relief / elevation", ["R(0) = R<sub>0</sub>  (columns: z<sub>A</sub>(0) = z<sub>B</sub>(0) = z<sub>0</sub>)"],
      "Where the landscape starts above base level. Starting at zero shows relief building from nothing; starting at the relief limit shows a landscape already in balance being perturbed. In the columns scenario the resistant body is positioned relative to this initial surface.", null, "asserted"),
    bodyDepth: M("Resistant body: depth of its top below the initial surface", ["top: z<sub>m</sub> = z<sub>0</sub> − d", "z<sub>m</sub> = z − U<sub>cum</sub>(t)  (material frame)"],
      "How deep the hard rock is buried at the start. The body is defined in material coordinates, so it rides up with uplift and is exposed once cumulative erosion of column B reaches d, whatever the uplift rate. Zero puts it at the surface from the start; a kilometre means the landscape must first cut down a kilometre to find it.",
      "core/lithology.js layer generator", "asserted"),
    bodyThick: M("Resistant body thickness T", ["exposed while  z<sub>0</sub> − d − T ≤ z<sub>B</sub> − U<sub>cum</sub> ≤ z<sub>0</sub> − d", "lifetime at the surface ≈ T · c / E<sub>ref</sub>"],
      "How much hard rock there is to get through. Once exposed, the body lasts its thickness divided by its (reduced) erosion rate. Hypothesis B needs that lifetime to exceed the uplift event: a thin cap is a passing episode, a thick body makes a lasting upland. The readout compares the two times.", null, "asserted"),
    contrast: M("Subaerial erodibility contrast c", ["E<sub>sub,B</sub> = D<sub>B</sub> · f(r, c),   f = 1 / (1 + r (c − 1))"],
      "How much harder the body is than the surrounding rock for rivers and hillslopes: c = 5 means it erodes five times slower. r ∈ [0, 1] is the resistance at the surface of B (1 inside the body); D<sub>B</sub> is the erosional demand set by the coupling control. While the body is at the surface B lowers at D<sub>B</sub>/c, so its lifetime is T·c/D<sub>B</sub>: a quartzite (c ~ 10) lasts ten times longer than the same thickness of shale.",
      "core/lithology.js erodibilityFactor", "asserted"),
    coupling: M("What sets the erosional demand on column B?", [
        "landscape:  D<sub>B</sub> = k z<sub>A</sub><sup>n</sup>", "own relief:  D<sub>B</sub> = k z<sub>B</sub><sup>n</sup>", "max:  D<sub>B</sub> = k max(z<sub>A</sub>, z<sub>B</sub>)<sup>n</sup>"],
      "The one real design choice of this tier, so it is a switch. <b>landscape</b>: B is a low-relief plateau inside the landscape of A; its erosion is the rate of A reduced by c, whatever height it reaches. Consequence: B can climb without limit while protected, and once stripped the upland <i>freezes</i> (E<sub>B</sub> = E<sub>A</sub>); with A glacially capped, the small subaerial rate of A lets B run away. <b>own relief</b>: B is its own relief reservoir; its erosion grows with its height, so while protected it equilibrates at c<sup>1/n</sup>·R<sub>eq</sub> and then erodes at U: the body is consumed in ≈T/U and persistence needs great thickness. <b>max</b>: landscape while B is low, own relief once it stands above A, so the upland is bounded and decays after stripping. Which is right is the hypothesis-B question; the sandbox shows the three side by side.",
      "docs/d1.md", "asserted"),
    contrastG: M("Glacial erodibility contrast c<sub>g</sub>", ["E<sub>gl,B</sub> = E<sub>gl</sub>(z<sub>B</sub>) · f(r, c<sub>g</sub>)"],
      "The resistance of the same body to glacial erosion, set separately because resistance to abrasion and plucking is not resistance to river incision. A body that is glacially hard but fluvially ordinary can survive the buzzsaw yet be cut by rivers, and vice versa.", null, "asserted"),
    glacierOn: M("Elevation-dependent glacial erosion (buzzsaw)", ["E<sub>gl</sub>(z) = P · g(z),  g = exp(−½ u<sup>2</sup>) below the ELA,  φ + (1 − φ) exp(−½ u<sup>2</sup>) above", "u = (z − ELA(t)) / w"],
      "Turns on glacial erosion as the column can represent it: not ice, but an erosion rate that depends on elevation. It peaks at the equilibrium-line altitude, where real glaciers slide fastest, and fades over a half-width w. With P ≫ U the column is capped where E<sub>gl</sub> + E<sub>sub</sub> = U, on the lower flank of the window: z<sub>cap</sub> = ELA − w √(2 ln(P / (U − E<sub>sub</sub>))). This is hypothesis A in one dimension.",
      "Brozović et al. 1997; Egholm et al. 2009; Mitchell & Montgomery 2006. The Gaussian window is an assertion of this tier.", "asserted"),
    elaBase: M("ELA (interglacial)", ["ELA(t) = ELA<sub>0</sub> − A<sub>ELA</sub> · ½(1 − cos(2π t / P<sub>ELA</sub>))"],
      "The equilibrium-line altitude in the warm state, above base level: where the glacial erosion window is centred. Lower it and the buzzsaw bites lower in the landscape; raise it above the relief limit and the columns never feel ice.", null, "asserted"),
    elaAmp: M("ELA cycle amplitude", ["ELA<sub>min</sub> = ELA<sub>0</sub> − A<sub>ELA</sub>"],
      "How far the ELA drops at the cold end of a glacial cycle. Zero holds it constant; a few hundred metres sweeps the erosion window up and down so the column is planed over a band rather than a line. Pleistocene depressions were roughly 500–1000 m.", null, "asserted"),
    elaPeriod: M("ELA cycle period", ["P<sub>ELA</sub>"],
      "The length of a glacial–interglacial cycle: 100 kyr for the late Pleistocene, 41 kyr for the obliquity-paced early Pleistocene. Shorter cycles give the buzzsaw less time at each level.", null, "asserted"),
    buzzPeak: M("Peak glacial erosion rate P", ["E<sub>gl</sub>(ELA) = P"],
      "The maximum lowering rate, at the ELA. If it is below the uplift rate the column climbs through the window; if well above, the column is pinned near the ELA. Measured glacial erosion rates span 0.01–10 mm/yr (Hallet et al. 1996).", null, "asserted"),
    buzzWidth: M("Buzzsaw half-width w", ["E<sub>gl</sub>(ELA ± w) = P · e<sup>−½</sup>"],
      "How far above and below the ELA glacial erosion remains strong. A narrow window planes a sharp level; a wide one erodes a broad band and caps the column further below the ELA for the same P.", null, "asserted"),
    buzzAbove: M("Erosion floor above the ELA φ", ["E<sub>gl</sub>(z ≫ ELA) → φ · P"],
      "How much glacial erosion continues high above the ELA. φ = 0 is the symmetric Gaussian window, so ground that rises above the window escapes glacial erosion altogether. φ > 0 says deep ice keeps eroding high ground at a fraction φ of the peak rate, the hypothesis-A ingredient (planation under deep ice). With φ·P/c<sub>g</sub> > U no column can escape upward.", "Egholm et al. 2009 (ice-cap erosion of high surfaces)", "asserted"),

    pUplift: M("Uplift forcing", ["U(t) = U<sub>peak</sub> s(t)"], "The shared forcing, as a function of time. Dashed vertical line: end of the event.", "core/forcing.js"),
    pRelief: M("Relief reservoir", ["dR/dt = U(t) − k R<sup>n</sup>", "R<sub>eq</sub>(t) = (U(t) / k)<sup>1/n</sup>  (dashed)"],
      "Relief as a reservoir that uplift fills and erosion drains. The dashed line is the relief the current uplift would sustain if held forever; the solid line chases it and lags, because the reservoir takes time to fill and drain. RK4 with 2000 steps; the stiff R<sup>n</sup> term is why RK4 rather than Euler.", SRC_RES, "computed"),
    pErosion: M("Erosion vs uplift flux", ["E = k R<sup>n</sup>"], "The two rates that compete. Where the curves cross, relief peaks (dR/dt = 0); erosion trails uplift on the way up and exceeds it on the way down.", null, "computed"),
    pNet: M("Net rate", ["dR/dt = U − E"], "The difference of the two: above zero relief grows, below zero it decays.", null, "computed"),
    pElev: M("Column elevations", ["dz<sub>A</sub>/dt = U − k z<sub>A</sub><sup>n</sup> − E<sub>gl</sub>(z<sub>A</sub>)", "dz<sub>B</sub>/dt = U − D<sub>B</sub> f(r<sub>B</sub>, c) − E<sub>gl</sub>(z<sub>B</sub>) f(r<sub>B</sub>, c<sub>g</sub>)"],
      "A (green): weak rock, the relief reservoir of the landscape. B (brown): same forcing, with a resistant body (shaded band, drawn in the present frame: it rises with U<sub>cum</sub>). While the surface of B lies inside the band its erosion is reduced. Blue dashed: ELA(t) when the glacier is on. RK4 on (z<sub>A</sub>, z<sub>B</sub>, U<sub>cum</sub>, E<sub>cum,A</sub>, E<sub>cum,B</sub>).",
      "py/d1_columns.py mirrors this; py/test_d1_columns.py checks z = z<sub>0</sub> + U<sub>cum</sub> − E<sub>cum</sub>, the T·c/U lifetime, exposure at d/U, the analytic buzzsaw cap.", "computed"),
    pRates: M("Erosion rates", ["E<sub>A</sub> = k z<sub>A</sub><sup>n</sup> + E<sub>gl</sub>(z<sub>A</sub>)", "E<sub>B</sub> = D<sub>B</sub> f(r<sub>B</sub>, c) + E<sub>gl</sub>(z<sub>B</sub>) f(r<sub>B</sub>, c<sub>g</sub>)"],
      "Total lowering rate of each column against the uplift rate (dashed grey). The rate of B drops by the contrast while the body is at its surface, and jumps back when it is stripped.", null, "computed"),
    nodeIntro: M("Isolated node of the along-valley profile", [
        "dz/dt = U(t) f − E<sub>f</sub> − E<sub>g</sub>",
        "ice-free:  E<sub>f</sub> = K f<sub>K</sub> A<sup>m</sup> ((z − z<sub>down</sub>)/ds)<sup>n</sup>  (implicit)",
        "under ice: Γ′ (½(H + H<sub>down</sub>))<sup>5</sup> ((z + H − z<sub>down</sub> − H<sub>down</sub>)/ds)<sup>3</sup> = Q/W<sub>f</sub>,  U<sub>s</sub> = f<sub>s</sub> (Q/W)/max(H, H<sub>min</sub>),  E<sub>g</sub> = K<sub>g</sub> f<sub>Kg</sub> U<sub>s</sub><sup>l</sup>"],
      "One node of the along-valley profile run on its own, with the same laws as the profile and everything that came from its neighbours frozen and badged asserted: the neighbour's erosion rate (it keeps uplifting and eroding as it was), its ice, the ice flux from upstream. It shows what a single cell can do by itself: relax its slope against its receiver, thicken its ice against the downstream surface, expose a body as it erodes. What it cannot do is see neighbours erode or a glacier advance. Perturb the constants to ask what the neighbours would have to do to change the fate of this column. The ELA sliders are inert here: the ELA reaches a node only through the flux Q, which is frozen.",
      "d1/model.js createNodeColumn; py/d1_node.py; docs/linking.md", "asserted"),
    nS: M("Position along the valley s", ["r = body(s, z − U<sub>cum</sub>)"], "Where along the valley this node sits. It enters only through the lithology field (a dike is a band in s), since area and width are given separately.", null, "asserted"),
    nZ0: M("Initial elevation z(t<sub>0</sub>)", null, "The bed of the node at the moment it was isolated. Everything is measured from base level at the outlet.", null, "asserted"),
    nZDown: M("Downstream bed z<sub>down</sub> at isolation", ["z<sub>down</sub>(t) = z<sub>down</sub>(t<sub>0</sub>) + ΔU<sub>cum</sub> − e<sub>down</sub>(t − t<sub>0</sub>);   S = (z − z<sub>down</sub>)/ds"], "The receiver's elevation when the node was isolated. From then on it keeps doing what it was doing: uplifting with the shared forcing and eroding at its isolation-time rate. In a steady profile those cancel and the receiver holds its height, so a steady node stays exactly where the profile would keep it. Lower it and the node steepens and erodes faster.", null, "asserted"),
    nHDown: M("Downstream ice thickness H<sub>down</sub>", ["surface below = z<sub>down</sub> + H<sub>down</sub>"], "How thick the ice is on the downstream cell, frozen. It sets the ice surface the face law works against: thicker downstream ice backs up the ice of this node and slows its sliding.", null, "asserted"),
    nQ: M("Ice flux Q leaving the node", ["Q = ∫<sub>0</sub><sup>s</sup> b W ds  (in the profile)"], "How much ice passes through this cell per year, frozen. It contains everything upstream: accumulation area, balance gradient, ELA. Q = 0 means no ice and the node is fluvial; larger Q means thicker, faster ice and more erosion.", null, "asserted"),
    nA: M("Drainage area A", ["A = k<sub>a</sub>(s + s<sub>0</sub>)<sup>h</sup>"], "How much basin drains through this cell, hence the discharge and erosive power of the river. Hack's law in the profile; a constant here. The steady slope goes as A<sup>−m/n</sup>.", null, "asserted"),
    nW: M("Valley width W", ["q = Q/W;  W<sub>f</sub> = ½(W + W<sub>down</sub>)"], "How wide the valley floor is: the same ice flux spread over a wider valley is thinner and slower, so wider means less glacial erosion.", null, "asserted"),
    nDs: M("Node spacing ds", ["S = (z − z<sub>down</sub>)/ds"], "The cell length over which the slope to the receiver is measured. Smaller spacing makes the same elevation drop a steeper slope.", null, "asserted"),
    nFac: M("Spatial uplift factor f", ["U(s, t) = U(t) f"], "How much of the peak uplift this position receives (1 under uniform uplift, less near the ends of a ramp or gaussian pattern).", null, "asserted"),
    nRun: M("Run length", null, "How far past t<sub>0</sub> to integrate the isolated column.", null, null),
    pNodeZ: M("Elevation of the isolated node", ["z(t), ice surface z + H, z<sub>down</sub> (frozen)"],
      "The bed and ice surface of the isolated cell through time, against its downstream bed (dashed grey), which keeps uplifting and eroding at its isolation-time rate. Under ice the bed lowers relative to the downstream surface, so H grows and sliding slows: a self-limiting overdeepening. Ice-free, the slope relaxes toward (U f / K A<sup>m</sup>)<sup>1/n</sup> and the node then rises with the block at the uplift rate.", null, "computed"),
    pNodeE: M("Erosion rates", ["E<sub>f</sub>, E<sub>g</sub>, U f"], "Fluvial and glacial lowering rates against the local uplift rate. Only one of the two acts, depending on whether the frozen flux Q is zero.", null, "computed"),
    pDiff: M("Upland height z<sub>B</sub> − z<sub>A</sub>", ["d(z<sub>B</sub> − z<sub>A</sub>)/dt = E<sub>A</sub> − E<sub>B</sub>"],
      "The hypothesis-B diagnostic: how high the resistant column stands above the weak landscape, and whether that survives the uplift event. Under the landscape coupling it can only grow or freeze (destruction by marginal retreat needs d2_across); under own-relief or max it is bounded and decays after stripping.", null, "computed")
  };

  // ---- Controls (display units) ------------------------------------------
  var SHAPES = [
    { v: "bell", t: "Bell curve" }, { v: "plateau", t: "Ramp · hold · decay" }, { v: "pulse", t: "Pulse (sharp rise, decay)" },
    { v: "arc", t: "Circular arc" }, { v: "triangle", t: "Triangle" }, { v: "constant", t: "Constant (square pulse)" }
  ];
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
  var nodeSpecs = [
    { id: "s", label: "Position s", type: "range", min: 0, max: 60, step: 0.1, val: 10, unit: "km", group: "Node (frozen adjacency)", asserted: true, method: m.nS },
    { id: "z0", label: "Initial elevation", type: "range", min: 0, max: 5, step: 0.01, val: 1.5, unit: "km", group: "Node (frozen adjacency)", asserted: true, method: m.nZ0 },
    { id: "zDown", label: "Downstream bed", type: "range", min: 0, max: 5, step: 0.01, val: 1.48, unit: "km", group: "Node (frozen adjacency)", asserted: true, method: m.nZDown },
    { id: "HDown", label: "Downstream ice", type: "range", min: 0, max: 800, step: 1, val: 0, unit: "m", group: "Node (frozen adjacency)", asserted: true, method: m.nHDown },
    { id: "eDown", label: "Downstream erosion rate", type: "range", min: 0, max: 10, step: 0.05, val: 1, unit: "mm/yr", group: "Node (frozen adjacency)", asserted: true,
      method: M("Receiver erosion rate e<sub>down</sub>", ["z<sub>down</sub>(t) = z<sub>down</sub>(t<sub>0</sub>) + ΔU<sub>cum</sub> − e<sub>down</sub>(t − t<sub>0</sub>)"],
        "How fast the downstream neighbour keeps eroding: its rate at the moment of isolation, held constant. Equal to the uplift rate means the receiver holds its height (a steady profile); less and it rises under this node, flattening the slope and starving erosion; more and it drops away, steepening the node. In the profile this rate is set by the receiver's own neighbours, which is exactly what the column cannot see.", null, "asserted") },
    { id: "QExp", label: "Ice flux Q (min = none)", type: "range", min: 0, max: 8, step: 0.05, val: 0, log: true, group: "Node (frozen adjacency)", asserted: true, method: m.nQ,
      fmt: function (v) { return v <= 0 ? "no ice" : Math.pow(10, v).toExponential(2) + " m³/yr"; } },
    { id: "AExp", label: "Drainage area A", type: "range", min: -2, max: 3, step: 0.05, val: 1.5, log: true, group: "Node (frozen adjacency)", asserted: true, method: m.nA,
      fmt: function (v) { return Math.pow(10, v).toPrecision(3) + " km²"; } },
    { id: "W", label: "Valley width", type: "range", min: 50, max: 3000, step: 10, val: 800, unit: "m", group: "Node (frozen adjacency)", asserted: true, method: m.nW },
    { id: "ds", label: "Node spacing", type: "range", min: 25, max: 500, step: 5, val: 100, unit: "m", group: "Node (frozen adjacency)", asserted: true, method: m.nDs },
    { id: "fac", label: "Uplift factor at s", type: "range", min: 0, max: 1, step: 0.01, val: 1, unit: "", group: "Node (frozen adjacency)", asserted: true, method: m.nFac },
    { id: "Wf", type: "range", val: 800, min: 0, max: 1e5, hidden: true }, { id: "ucum0", type: "range", val: 0, min: -1e6, max: 1e6, hidden: true }, { id: "t0", type: "range", val: 0, min: 0, max: 1e9, hidden: true },
    { id: "runKyr", label: "Run length", type: "range", min: 50, max: 3000, step: 50, val: 500, unit: "kyr", group: "Run", method: m.nRun },
    { id: "dt", label: "Time step", type: "range", min: 100, max: 2000, step: 100, val: 500, unit: "yr", group: "Run", asserted: true, method: M("Time step", ["one step: lithology → uplift → ice or river"], "How much time each step covers. The river is implicit and the ice steady, so results barely depend on it.", null, "asserted") },
    { id: "KExp", label: "Erodibility K", type: "range", min: -7, max: -4, step: 0.1, val: -5.5, log: true, group: "Process laws (as in d2_along)", method: M("Fluvial erodibility K", ["E<sub>f</sub> = K A<sup>m</sup> S<sup>n</sup>"], "How easily the river cuts rock; lower K needs a steeper slope to keep pace with uplift. Same value as the profile it came from.", null, "asserted") },
    { id: "m", label: "Area exponent m", type: "range", min: 0.3, max: 0.7, step: 0.05, val: 0.5, unit: "", group: "Process laws (as in d2_along)", asserted: true, method: M("Area exponent m", ["E<sub>f</sub> ∝ A<sup>m</sup>"], "How strongly drainage area (discharge) drives incision at this cell.", null, "asserted") },
    { id: "nexp", label: "Slope exponent n", type: "range", min: 1, max: 2, step: 0.25, val: 1, unit: "", group: "Process laws (as in d2_along)", asserted: true, method: M("Slope exponent n", ["E<sub>f</sub> ∝ S<sup>n</sup>"], "How nonlinearly incision responds to steepness; larger n makes a steep cell relax to its equilibrium slope faster.", null, "asserted") },
    { id: "KgExp", label: "Glacial erodibility K_g", type: "range", min: -6, max: -2, step: 0.1, val: -4, log: true, group: "Process laws (as in d2_along)", method: M("Glacial erodibility", ["E<sub>g</sub> = K<sub>g</sub> U<sub>s</sub><sup>l</sup>, capped at 2 cm/yr"], "How much bed is removed per metre of sliding. Larger values deepen the cell faster, but its ice thickens against the frozen downstream surface and sliding slows, so it self-limits.", null, "asserted") },
    { id: "lexp", label: "Sliding exponent l", type: "range", min: 1, max: 2, step: 0.25, val: 1, unit: "", group: "Process laws (as in d2_along)", asserted: true, method: M("Sliding exponent l", ["E<sub>g</sub> ∝ U<sub>s</sub><sup>l</sup>"], "How nonlinearly erosion grows with sliding speed: 1 proportional, 2 the Hallet abrasion law.", null, "asserted") },
    { id: "fs", label: "Sliding fraction", type: "range", min: 0, max: 0.95, step: 0.05, val: 0.5, unit: "", group: "Process laws (as in d2_along)", asserted: true, method: M("Sliding fraction f<sub>s</sub>", ["U<sub>s</sub> = f<sub>s</sub> q/H"], "What share of ice motion is sliding rather than internal deformation. Only sliding erodes; more sliding also means thinner ice for the same flux.", null, "asserted") },
    { id: "flow", label: "Flow enhancement", type: "range", min: 0.2, max: 5, step: 0.1, val: 1, unit: "×", group: "Process laws (as in d2_along)", asserted: true, method: M("Flow-law enhancement", ["Γ → Γ · E"], "Softness of the ice: softer ice (E > 1) carries the frozen flux with a thinner, faster column.", null, "asserted") },
    { id: "Hmin", label: "Min sliding thickness", type: "range", min: 0, max: 50, step: 1, val: 10, unit: "m", group: "Process laws (as in d2_along)", asserted: true, method: M("Minimum sliding thickness H<sub>min</sub>", ["U<sub>s</sub> = f<sub>s</sub> q / max(H, H<sub>min</sub>)"], "Ice thinner than this is treated as too thin to slide erosively, which keeps q/H from diverging when the ice is a sliver.", null, "asserted") },
    { id: "Hf", label: "Fluvial shut-off thickness", type: "range", min: 10, max: 300, step: 10, val: 100, unit: "m", group: "Process laws (as in d2_along)", asserted: true, method: M("Fluvial shut-off thickness H<sub>f</sub>", ["K → K · (φ<sub>sub</sub> + (1 − φ<sub>sub</sub>) max(0, 1 − H/H<sub>f</sub>))"], "How thick the ice must be before it fully protects this cell from the river; thinner ice lets the river keep working at reduced efficiency.", null, "asserted") },
    { id: "phiSub", label: "Subglacial fluvial efficiency", type: "range", min: 0, max: 1, step: 0.05, val: 0, unit: "", group: "Process laws (as in d2_along)", asserted: true, method: M("Subglacial fluvial efficiency φ<sub>sub</sub>", ["w(H ≥ H<sub>f</sub>) = φ<sub>sub</sub>"], "How much river incision continues under thick ice (subglacial meltwater); 0 means full protection.", null, "asserted") },
    { id: "contrastK", label: "Fluvial contrast", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Lithology (as in d2_along)", asserted: true, method: M("Fluvial contrast c", ["K → K / (1 + r(c − 1))"], "How much slower the river cuts the resistant body than the surrounding rock.", null, "asserted") },
    { id: "contrastKg", label: "Glacial contrast", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Lithology (as in d2_along)", asserted: true, method: M("Glacial contrast c<sub>g</sub>", ["K<sub>g</sub> → K<sub>g</sub> / (1 + r(c<sub>g</sub> − 1))"], "How much slower ice abrades the resistant body; set separately from the fluvial contrast.", null, "asserted") },
    { id: "lithoType", label: "Body type", type: "select", val: "none", options: [{ v: "none", t: "None" }, { v: "layer", t: "Horizontal layer" }, { v: "slab", t: "Dipping slab" }, { v: "dike", t: "Vertical dike" }], group: "Lithology (as in d2_along)", asserted: true, method: M("Resistant body", ["r = body(s, z − U<sub>cum</sub>)"], "A body of harder rock in the material frame, as on the profile; it rides up with uplift and is exposed when this cell has eroded down to it.", "core/lithology.js", "asserted") },
    { id: "lithoTop", label: "Body top", type: "range", min: 0, max: 5, step: 0.1, val: 1.5, unit: "km", group: "Lithology (as in d2_along)", asserted: true, method: M("Body top", ["top (m, initial frame)"], "How high the top of the body sits in the initial material frame; below the current surface means it emerges later.", null, "asserted") },
    { id: "lithoThick", label: "Thickness / width", type: "range", min: 0.05, max: 3, step: 0.05, val: 0.3, unit: "km", group: "Lithology (as in d2_along)", asserted: true, method: M("Body thickness / dike width", ["T or w"], "How much resistant rock there is to get through before erosion returns to the background rate.", null, "asserted") },
    { id: "lithoPos", label: "Dike position / slab dip", type: "range", min: 0, max: 1, step: 0.05, val: 0.5, unit: "", group: "Lithology (as in d2_along)", asserted: true, method: M("Dike position / slab dip", ["dike: s<sub>0</sub> = x·L;  slab: dip = 0.2 x"], "Where the dike sits along the valley (fraction of L), or how steeply the slab dips.", null, "asserted") },
    { id: "L", label: "Valley length (for the dike position)", type: "range", min: 10, max: 60, step: 5, val: 30, unit: "km", group: "Lithology (as in d2_along)", asserted: true, method: M("Valley length L", ["s<sub>0</sub> = x·L"], "Only used to place a dike; the column itself has no length.", null, "asserted") }
  ];
  var SCEN = [{ v: "columns", t: "Two columns A / B (lithology + buzzsaw)" }, { v: "reservoir", t: "Relief reservoir (legacy)" }, { v: "node", t: "Isolated node of the along-valley profile" }];
  var scenarioSpec = { id: "scenario", label: "Scenario", type: "select", val: "columns", options: SCEN, group: "",
    method: M("Scenarios", null, "<b>Two columns</b>: the illustration device (DESIGN.md §3) — a weak and a resistant column under one forcing, with an asserted elevation-dependent glacial rate. <b>Relief reservoir</b>: the legacy 1-D model, kept as a scenario. <b>Isolated node</b>: one node of the along-valley profile with its adjacency frozen — opened from a click on the d2_along page, or configured here.") };
  var nodeIntroSpec = { id: "nodeIntro", label: "What this is", type: "check", val: true, hidden: true };
  var allSpecs = [scenarioSpec].concat(common, landscape, body, glacier, nodeSpecs);
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
    } else if (state.scenario === "node") {
      panels.push(P.make(charts, { title: "Isolated node: elevation, ice surface, frozen downstream bed", yLabel: "z (km)", color: C.B, method: m.nodeIntro }));
      panels.push(P.make(charts, { title: "Erosion rates  vs  uplift", yLabel: "rate (mm/yr)", color: C.erosion, method: m.pNodeE }));
      panels.push(P.make(charts, { title: "Ice thickness and sliding speed (scaled)", yLabel: "H (m)", color: C.ela, showX: true, xLabel: "time since t₀ (kyr)", method: m.pNodeZ }));
    } else {
      panels.push(P.make(charts, { title: "Uplift", yLabel: "U (mm/yr)", color: C.uplift, method: m.pUplift }));
      panels.push(P.make(charts, { title: "Elevation above base level  —  A weak · B resistant body", yLabel: "z (km)", color: C.B, method: m.pElev }));
      panels.push(P.make(charts, { title: "Erosion rates  vs  uplift", yLabel: "rate (mm/yr)", color: C.erosion, method: m.pRates }));
      panels.push(P.make(charts, { title: "Upland height  z_B − z_A", yLabel: "Δz (km)", color: C.diff, showX: true, xLabel: "time (Myr)", method: m.pDiff }));
    }
    d3.select("#legend").html(state.scenario === "node"
      ? '<span style="color:' + C.B + '"><i></i>bed z</span> <span style="color:' + C.ela + '"><i></i>ice surface</span> <span style="color:#999"><i class="dash"></i>downstream bed (frozen)</span><br/>' +
        '<span style="color:#2c7fb8"><i></i>fluvial</span> <span style="color:' + C.erosion + '"><i></i>glacial</span> <span style="color:#888"><i class="dash"></i>uplift</span> <span style="color:' + C.diff + '"><i class="dash"></i>sliding (scaled)</span>'
      : state.scenario === "reservoir"
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

  function renderNode() {
    var Lm = U.fromKm(state.L);
    var params = {
      shape: state.shape, peakUplift: U.fromMmyr(state.peakUplift), duration: U.fromMyr(state.duration),
      elaBase: U.fromKm(state.elaBase), elaAmp: U.fromKm(state.elaAmp), elaPeriod: U.fromKyr(state.elaPeriod),
      K: Math.pow(10, state.KExp), m: state.m, nexp: state.nexp, Kg: Math.pow(10, state.KgExp), lexp: state.lexp, fs: state.fs, flow: state.flow,
      eroCap: 0.02, Hmin: state.Hmin, Hf: state.Hf, phiSub: state.phiSub, contrastK: state.contrastK, contrastKg: state.contrastKg
    };
    var col = GS.d1.createNodeColumn({
      s: U.fromKm(state.s), z0: U.fromKm(state.z0), ucum0: +state.ucum0, zDown: U.fromKm(state.zDown), HDown: state.HDown,
      Q: state.QExp > 0 ? Math.pow(10, state.QExp) : 0, A: Math.pow(10, state.AExp) * 1e6, W: state.W, Wf: +state.Wf > 0 ? +state.Wf : state.W, eDown: U.fromMmyr(state.eDown),
      ds: state.ds, fac: state.fac, t0: +state.t0, dt: state.dt, params: params,
      litho: GS.d2along.lithoFromControls(state.lithoType, U.fromKm(state.lithoTop), U.fromKm(state.lithoThick), state.lithoPos, Lm)
    });
    var nsteps = Math.round(U.fromKyr(state.runKyr) / state.dt);
    for (var k = 0; k < nsteps; k++) { col.step(); col.record(); }
    var data = col.state.series, last = data[data.length - 1], t0 = +state.t0;
    var tk = function (d) { return U.toKyr(d.t - t0); };
    var pZ = panels[0], pE = panels[1], pI = panels[2];
    var zLo = d3.min(data, function (d) { return U.toKm(Math.min(d.z, d.zDown)); }), zHi = d3.max(data, function (d) { return U.toKm(d.z + d.H); });
    var eMax = d3.max(data, function (d) { return U.toMmyr(Math.max(d.Ef, d.Eg, d.u)); }) || 0.1;
    var hMax = d3.max(data, function (d) { return d.H; }) || 1, usMax = d3.max(data, function (d) { return d.Us; }) || 1;
    panels.forEach(function (pn) { pn.x.domain([0, U.toKyr(last.t - t0)]); P.clear(pn); });
    pZ.y.domain([zLo - 0.02, zHi + 0.02]); pE.y.domain([0, eMax * 1.1]); pI.y.domain([0, hMax * 1.1]);
    P.line(pZ, "zdown", data, function (d) { return U.toKm(d.zDown); }, "#999", true, tk);
    P.line(pZ, "ice", data, function (d) { return U.toKm(d.z + d.H); }, C.ela, false, tk);
    P.line(pZ, "z", data, function (d) { return U.toKm(d.z); }, C.B, false, tk);
    P.line(pE, "u", data, function (d) { return U.toMmyr(d.u); }, C.flux, true, tk);
    P.line(pE, "Ef", data, function (d) { return U.toMmyr(d.Ef); }, "#2c7fb8", false, tk);
    P.line(pE, "Eg", data, function (d) { return U.toMmyr(d.Eg); }, C.erosion, false, tk);
    P.line(pI, "H", data, function (d) { return d.H; }, C.ela, false, tk);
    P.line(pI, "Us", data, function (d) { return hMax * 1.1 * d.Us / (usMax * 1.1); }, C.diff, true, tk);
    panels.forEach(function (pn) { P.axes(pn); });
    d3.select("#readout").html("Isolated node at <b>" + U.fmtLen(U.fromKm(state.s)) + "</b>: z " + U.fmtLen(U.fromKm(state.z0)) + " → <b>" + U.fmtLen(last.z) + "</b> after " + U.fmtTime(last.t - t0) +
      " · ice " + last.H.toFixed(0) + " m · sliding " + last.Us.toFixed(1) + " m/yr · erosion " + U.fmtRate(last.Ef + last.Eg) + (col.frozen.Q > 0 ? " (glacial)" : " (fluvial)") +
      " · <a href='../d2_along/index.html'>back to the profile</a>");
  }

  function render() {
    if (state.scenario === "reservoir") renderReservoir(); else if (state.scenario === "node") renderNode(); else renderColumns();
    GS.ui.holdHeight(d3.select("#readout"));
    url.write(state);
  }

  function buildSidebar() {
    var root = d3.select("#controls");
    root.selectAll("*").remove();
    var specs = state.scenario === "reservoir" ? [scenarioSpec].concat(common, landscape)
              : state.scenario === "node" ? [scenarioSpec].concat(common, nodeSpecs)
              : [scenarioSpec].concat(common, landscape, body, glacier);
    GS.ui.buildControls(root, specs, state, function (id) {
      if (id === "scenario") { GS.ui.holdHeight(d3.select("#readout"), true); buildSidebar(); buildPanels(); }
      render();
    });
  }

  buildSidebar();
  buildPanels();
  render();
})();
