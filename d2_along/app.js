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
    L: M("Valley length L", ["s ∈ [0, L],  ds = L / (N − 1)"],
      "How long the modelled valley is, from the divide (left, s = 0) to the outlet (right, s = L), which sits at fixed base level. Longer valleys gather more drainage area and more ice, so both the river and the glacier get stronger downstream; the glacier's terminus and the fluvial relief scale with it.", null, "asserted"),
    N: M("Nodes N", ["ds = L / (N − 1)"],
      "The number of cells along the profile. More nodes resolve narrower features (cirque steps, thin termini) but each step costs proportionally more. 301 nodes on 30 km is a 100 m cell, about the scale at which the process laws here stop being meaningful.", null, "asserted"),
    initProfile: M("Initial profile", ["steady: S = (U<sub>peak</sub> f / K A<sup>m</sup>)<sup>1/n</sup> integrated up from the outlet", "linear: z = z<sub>head</sub>(1 − s/L)", "concave: z = z<sub>head</sub>(1 − s/L)<sup>1.5</sup>"],
      "The shape of the valley at t = 0. The fluvial steady state (default) starts the river in equilibrium with the peak uplift, so anything that then happens is the glacier's doing. Its slope is capped at the threshold slope S<sub>c</sub>: where U/(K A<sup>m</sup>) asks for more, the profile is at threshold and not in fluvial equilibrium, and the readout says so. Linear and concave are far from equilibrium and incise fast at first, which is useful for watching a river adjust. Head elevation applies to those two only.", null, "asserted"),
    zHead: M("Head elevation", ["z(0, 0) = z<sub>head</sub>"],
      "How high the divide starts (linear and concave profiles only). Higher heads put more of the valley above the ELA and grow a bigger initial glacier; they also mean steeper initial slopes and faster early incision.", null, "asserted"),
    noise: M("Profile noise σ", ["z += σ · 𝒩(0, 1), seeded"],
      "Small random roughness added to the initial bed. Zero gives a smooth profile; a few metres seeds the small bumps that quarrying and overdeepening can amplify, which is how a uniform bed learns to make steps.", null, "asserted"),
    shape: M("Uplift forcing shape U(t)", ["U(t) = U<sub>peak</sub> · s(t),  s ∈ [0, 1]"],
      "The time history of rock uplift. Constant holds the peak rate for the whole event; bell, plateau, pulse, arc and triangle rise and fall within the event duration so you can watch relief lag the forcing and decay after it. The outlet never uplifts: it is the fixed base level everything else is measured against.", "core/forcing.js", "asserted"),
    peakUplift: M("Peak rock-uplift rate", ["U<sub>peak</sub>  (mm/yr = km/Myr)"],
      "How fast rock rises past base level at the height of the event. This is the engine of the whole run: at steady state every erosion process together must remove exactly this much. Doubling it doubles the steady fluvial slope (for n = 1) and pushes more of the valley up through the ELA. Active ranges sit at 0.5–10 mm/yr.", null, "asserted"),
    duration: M("Event duration T", ["s(t) = 0 for t > T"],
      "How long the uplift event lasts. After it the profile relaxes under erosion alone, so relief and the glacier both decay; a short event against a slow river never reaches steady state.", null, "asserted"),
    pattern: M("Spatial uplift pattern", ["U(s, t) = U(t) · f(s/L)"],
      "Where along the valley the uplift is concentrated. Uniform lifts everything equally; ramp lifts the middle and not the ends; gaussian peaks mid-valley; tilt rises toward the outlet, i.e. uplift is strongest near base level. Non-uniform patterns change where the river steepens and where ice accumulates.", "core/forcing.js spatialFactor", "asserted"),
    KExp: M("Fluvial erodibility K", ["E<sub>f</sub> = K A<sup>m</sup> S<sup>n</sup>  (× lithology factor)", "steady state:  S = (U / K A<sup>m</sup>)<sup>1/n</sup>"],
      "How easily the river cuts rock. It sets the relief a given uplift can hold: for n = 1 the steady slope is inversely proportional to K, so lowering K by 10 makes a range 10× higher. Very low K asks for slopes steeper than rock can stand and the profile sits at threshold instead. Solved implicitly node by node from the outlet upstream, on ice-free nodes only. Typical values 10<sup>−6</sup>–10<sup>−5</sup> in these units.", "Braun & Willett 2013; Yuan et al. 2019; validated: steady slope–area law exact, closure to round-off", "computed"),
    m: M("Area exponent m", ["E<sub>f</sub> ∝ A<sup>m</sup>"],
      "How strongly discharge (via drainage area) drives incision. Larger m makes the lower valley erode much faster than the head and so makes the steady profile more concave; the concavity θ = m/n. Usually 0.4–0.6.", null, "asserted"),
    nexp: M("Slope exponent n", ["E<sub>f</sub> ∝ S<sup>n</sup>"],
      "How nonlinearly incision responds to steepness. n = 1 is the linear stream-power law; larger n makes steep reaches erode disproportionately faster, so knickpoints sharpen and relief responds less than proportionally to uplift. n ≠ 1 needs a Newton iteration inside the implicit update.", null, "asserted"),
    G: M("Deposition coefficient G", ["dz/dt = −K A<sup>m</sup> S<sup>n</sup> + G Q<sub>s</sub>/A"],
      "How much of the eroded sediment is laid back down downstream instead of being carried out. G = 0 is purely detachment-limited: everything eroded leaves at once. Larger G fills low-gradient reaches and lets a valley floor aggrade below a fast-eroding reach. Mass-conserving by construction.", "Yuan et al. 2019", "asserted"),
    Hf: M("Fluvial shut-off thickness H<sub>f</sub>", ["w(H) = φ<sub>sub</sub> + (1 − φ<sub>sub</sub>) · max(0, 1 − H/H<sub>f</sub>),   K → K · w"],
      "How thick the ice must be before it fully protects the bed from the river. Fluvial efficiency is 1 on bare ground and ramps down to the subglacial floor over the first H<sub>f</sub> metres of ice, instead of switching off at the first cell with ice. The switch was a cell-scale defect: each cell the terminus vacated as the ELA cycled took a burst of incision before ice returned, leaving a sawtooth in bed and ice thickness at the margin. The ramp reads as meltwater under thin marginal ice. 50–150 m; larger is smoother.", null, "asserted"),
    phiSub: M("Subglacial fluvial efficiency φ<sub>sub</sub>", ["w(H ≥ H<sub>f</sub>) = φ<sub>sub</sub>"],
      "How much river incision continues under thick ice: subglacial meltwater channels erode too. 0 means thick ice fully protects the bed; 0.2–0.5 lets the subglacial river keep pace with part of the uplift, which shrinks the step that forms at a terminus where glacial erosion is weak.", "subglacial meltwater erosion: Beaud et al. 2014; Herman et al. 2011", "asserted"),
    hack: M("Hack exponent h", ["A(s) = k<sub>a</sub> (s + s<sub>0</sub>)<sup>h</sup>"],
      "How fast drainage area grows downstream along the profile. The profile cannot see its basin, so area is asserted from Hack's law (L = 1.4 A<sup>0.6</sup> gives h = 1/0.6 = 1.67). Larger h means area, and so discharge and incision, grow faster with distance; the steady profile becomes more concave.", "Hack 1957", "asserted"),
    W0: M("Valley width at the head W<sub>0</sub>", ["W(s) = W<sub>0</sub> + k<sub>w</sub> s"],
      "How wide the valley floor is at the divide. Width matters twice for ice: the flux gathered per unit length is b·W, and the flux per unit width that sets thickness is Q/W. A wider head gathers more ice; a wider trunk spreads the same flux thinner and slower.", null, "asserted"),
    kw: M("Valley width growth k<sub>w</sub>", ["W(s) = W<sub>0</sub> + k<sub>w</sub> s"],
      "How much wider the valley gets per metre downstream. Zero is a canal of constant width; 0.05 makes a 300 m head into a 1.8 km trunk at 30 km. Faster widening thins and slows the trunk glacier for a given balance, weakening erosion there.", null, "asserted"),
    glacierOn: M("Glacier: steady ice discharge  (outlet: free outflow, Γ′H<sup>5</sup>S<sub>bed</sub><sup>3</sup> = q)", [
        "b(z<sub>s</sub>) = β (z<sub>s</sub> − ELA) above, 2.5 β (z<sub>s</sub> − ELA) below, clamped",
        "Q(s) = ∫<sub>0</sub><sup>s</sup> b W ds, clamped ≥ 0 → terminus",
        "Γ′ H<sub>f</sub><sup>n+2</sup> S<sub>f</sub><sup>n</sup> = Q/W<sub>f</sub> at each face, Γ′ = Γ/(1 − f<sub>s</sub>), n = 3",
        "U<sub>s</sub> = f<sub>s</sub> q / max(H, H<sub>min</sub>)"],
      "Turns the glacier on. Ice adjusts in centuries while the bed changes over hundreds of millennia, so every step solves the glacier that is in equilibrium with the current bed: mass balance on the ice surface, summed down-valley to give the flux at every point (where the sum hits zero the glacier ends), then the thickness that carries that flux, marched upstream from the outlet by solving the flow law face by face. A few outer passes couple surface and balance. Validated against an explicit flowline SIA: thickness within 2 %, same terminus.",
      "SIA: Glen n = 3, A = 7.57×10<sup>−17</sup> Pa<sup>−3</sup> yr<sup>−1</sup>; MacGregor et al. 2000; Anderson et al. 2006", "computed"),
    elaBase: M("ELA (interglacial)", ["ELA(t) = ELA<sub>0</sub> − A<sub>ELA</sub> · ½(1 − cos 2πt/P)"],
      "The equilibrium-line altitude in the warm state: above it snow accumulates, below it ice melts. Lower it and more of the valley is accumulation area, so the glacier grows, thickens and reaches further down; raise it above the divide and the glacier vanishes. Referenced to the ice surface, so thicker ice sees a colder surface (a feedback kept from the legacy model).", null, "asserted"),
    elaAmp: M("ELA cycle amplitude", ["ELA<sub>min</sub> = ELA<sub>0</sub> − A<sub>ELA</sub>"],
      "How far the ELA drops at the cold end of a glacial cycle. Zero holds the ELA constant; a few hundred metres makes the glacier advance and retreat each cycle, sweeping its erosion up and down the valley. Pleistocene cycles moved ELAs by roughly 500–1000 m.", null, "asserted"),
    elaPeriod: M("ELA cycle period", ["P"],
      "The length of one glacial–interglacial cycle. 100 kyr matches the late Pleistocene, 41 kyr the obliquity-paced early Pleistocene. Shorter periods give the glacier less time to reshape the bed per advance.", null, "asserted"),
    balGrad: M("Balance gradient β", ["b = β (z<sub>s</sub> − ELA) above the ELA; 2.5 β below"],
      "How quickly mass balance changes with height: how many metres of ice per year are gained per metre above the ELA (and lost 2.5× faster per metre below, the usual asymmetry). Steeper gradients belong to maritime climates and make bigger, faster glaciers for the same ELA; continental glaciers have gradients a few times smaller. Accumulation clamps at 2 m/yr, ablation at 8 m/yr.", "asymmetry after the legacy terrain_sandbox", "asserted"),
    fs: M("Sliding fraction f<sub>s</sub>", ["U<sub>s</sub> = f<sub>s</sub> q/H,   q<sub>deformation</sub> = (1 − f<sub>s</sub>) q"],
      "What share of the ice's motion is sliding over the bed rather than deforming internally. Only sliding erodes, so f<sub>s</sub> = 0 gives no glacial erosion at all, and a temperate glacier near 0.8 erodes hard. Sliding also thins the glacier, since less thickness is needed to carry the same flux.", null, "asserted"),
    flow: M("Flow-law enhancement", ["Γ → Γ · E"],
      "Multiplies the softness of the ice. Warmer or more damaged ice (E > 1) flows more easily, so the same flux is carried by thinner, faster ice; E < 1 is stiffer, thicker ice. 1 is temperate ice.", null, "asserted"),
    KgExp: M("Glacial erodibility K<sub>g</sub>", ["E<sub>g</sub> = K<sub>g</sub> U<sub>s</sub><sup>l</sup>  (× lithology factor), capped at 2 cm/yr"],
      "How much bed is removed per metre of sliding: abrasion by rock in the ice sole. Raising it deepens troughs and overdeepenings faster, but the process self-limits: an overdeepening flattens the ice surface, thickens the ice and slows sliding. At 100× the default the glacier cuts its bed below the ELA and shuts itself off. Typical 10<sup>−4</sup> for l = 1.", "Hallet 1979; Humphrey & Raymond 1994", "computed"),
    lexp: M("Sliding exponent l", ["E<sub>g</sub> ∝ U<sub>s</sub><sup>l</sup>"],
      "How nonlinearly erosion grows with sliding speed. l = 1 is proportional (Humphrey & Raymond); l = 2 (Hallet's abrasion law) concentrates erosion where ice moves fastest, sharpening trunk versus tributary contrasts.", null, "asserted"),
    Hmin: M("Minimum sliding thickness H<sub>min</sub>", ["U<sub>s</sub> = f<sub>s</sub> q / max(H, H<sub>min</sub>)"],
      "Ice thinner than this is treated as too thin to slide erosively. It is a regularisation with a physical reading: the last few metres of a terminus wedge are cold, crevassed or debris-laden and do not abrade. Without it q/H diverges as the ice thins to nothing at the last cell, and a sliding spike appeared at the glacier front and at the outlet. 5–20 m is reasonable.", null, "asserted"),
    eroSmooth: M("Erosion footprint [¼ ½ ¼]", ["E<sub>g,i</sub> ← ¼E<sub>i−1</sub> + ½E<sub>i</sub> + ¼E<sub>i+1</sub>"],
      "Spreads each cell's glacial erosion over its two neighbours, on the grounds that abrasion under ice a hundred metres thick acts over a patch that size, not a single 100 m cell. It also annihilates the two-cell checkerboard exactly: with erosion ∝ 1/H and a face law that averages neighbouring thicknesses, odd and even cells can otherwise decouple. Off shows the raw law.", null, "asserted"),
    KqExp: M("Quarrying on convexity K<sub>q</sub>", ["E<sub>q</sub> = K<sub>q</sub> U<sub>s</sub> · max(−∂²z/∂s², 0)"],
      "Extra erosion where the bed is convex-up under sliding ice, the plucking of steps and bumps. This is the knob hypothesis A needs: it planes off high convexities under deep ice. Off by default (slider at its minimum); with it on, riegels and cirque lips erode preferentially and small bumps in the bed get removed rather than amplified.", "DESIGN.md §4", "asserted"),
    lithoType: M("Resistant body (material frame)", ["z<sub>m</sub> = z − U<sub>cum</sub>(s)", "layer: z<sub>m</sub> ∈ [top − T, top];  dike: |s − s<sub>0</sub>| < w/2;  slab: top + dip·s"],
      "A body of harder rock inside the valley. It is defined in material coordinates, so it rides up with uplift and gets exposed once erosion has removed everything above it; where the surface is inside it, erosion is divided by the contrasts. A layer is horizontal, a slab dips along the valley, a dike is a vertical band at one position. Drawn as a brown overlay where the surface is in the body.", "core/lithology.js", "asserted"),
    lithoTop: M("Body top (material height)", ["top  (m above base level, in the initial frame)"],
      "How high the top of the body sits in the initial material frame. Put it below the initial surface and it emerges later as the valley cuts down and the rock rises; put it above and it is exposed from the start.", null, "asserted"),
    lithoThick: M("Body thickness / dike width", ["T (layer, slab)  or  w (dike)"],
      "How much resistant rock there is to get through. A thin layer is stripped and its effect is a passing episode; a thick one persists for the whole run. For hypothesis B the question is whether it outlasts the uplift event.", null, "asserted"),
    lithoPos: M("Dike position / slab dip", ["dike: s<sub>0</sub> = x·L;  slab: dip = 0.2 x (m/m, rising toward the divide)"],
      "For a dike, where along the valley it sits (0 = divide, 1 = outlet). For a slab, how steeply it dips, so that it is exposed at different heights in different reaches.", null, "asserted"),
    contrastK: M("Fluvial contrast c", ["K → K / (1 + r(c − 1))"],
      "How much harder the body is for the river than the surrounding rock: c = 5 means it erodes five times slower. Steady slopes across it are c<sup>1/n</sup> times steeper, so a dike makes a knickzone and a layer makes a bench.", null, "asserted"),
    contrastKg: M("Glacial contrast c<sub>g</sub>", ["K<sub>g</sub> → K<sub>g</sub> / (1 + r(c<sub>g</sub> − 1))"],
      "The same body's resistance to glacial abrasion, set separately because rocks resist plucking and abrasion differently than fluvial incision. A body that is glacially hard but fluvially soft is exactly the case that separates the two hypotheses.", null, "asserted"),
    Sc: M("Threshold slope S<sub>c</sub> (rockfall / step failure)", ["every cell:  z<sub>i</sub> ≤ z<sub>i+1</sub> + S<sub>c</sub> ds", "excess removed → rockfall export"],
      "The steepest slope a 100 m cell may stand at above its downstream neighbour. Anything steeper collapses to the threshold, with the debris assumed evacuated by the glacier or river below. It stands in for headwall retreat and step collapse, processes the profile cannot resolve. Without it the divide cell, which receives almost no ice flux and so cannot be glacially eroded, rises with uplift into a single-cell spire (680 m in 200 kyr in one scenario) that rings the cirque floor below. The same value caps the steady initial profile. 0.6–1.0 (30–45°) is a hillslope threshold; cirque headwalls are steeper, but not for 100 m at a time.",
      "threshold hillslopes: Burbank et al. 1996; DESIGN.md §4 hillslopes", "asserted"),
    dt: M("Time step", ["per step: lithology → uplift → steady ice + erosion → fluvial → rockfall"],
      "How much geomorphic time each step covers. The river is solved implicitly and the ice is steady, so the step is limited only by how much the bed may change per solve; steady profiles are dt-independent to 0.01 m and a 200 kyr transient differs by ~15 m between 250 and 1000 yr. Larger steps run faster and smear fast transients.", null, "asserted"),
    speed: M("Steps per frame", null, "How many model steps run between screen redraws when playing. Only affects animation speed, not the result.", null, null),
    pProfile: M("Long profile", ["bed z(s), ice surface z + H, ELA(t)"],
      "The valley seen from the side. Brown: bedrock. Blue fill: the steady glacier. Dashed: the ELA. Thick brown overlay: where the surface is inside the resistant body. The left end is the divide, the right end the outlet at fixed base level. Click anywhere to isolate that node as a d1 column.", null, "computed"),
    pRates: M("Erosion rates", ["E<sub>f</sub> = K A<sup>m</sup> S<sup>n</sup> (ice-free);  E<sub>g</sub> = K<sub>g</sub> U<sub>s</sub><sup>l</sup> (+ quarrying) under ice;  rockfall (brown dashed, clipped);  U(s, t) grey dashed"],
      "How fast each process lowers the bed at every point, against the uplift rate. Where a curve sits on the uplift line that reach is in steady state; above it the bed is being lowered net, below it is rising. Under thin marginal ice both act, the river at reduced efficiency.", null, "computed"),
    pIce: M("Sliding speed", ["U<sub>s</sub> = f<sub>s</sub> q / max(H, H<sub>min</sub>)  (m/yr)"],
      "How fast the ice slides over its bed at every point. This is what erodes: glacial erosion is proportional to it. It rises down the accumulation area as flux grows, peaks where the ice is fast and not too thick, and falls to zero at the terminus.", null, "computed"),
    pHist: M("History", ["relief = max z − min z;  ice volume = Σ H W ds"],
      "How the run has evolved so far: total relief (solid) and ice volume (dashed, scaled to the panel). Mass closure (uplift − erosion − Δvolume) is in the readout.", null, "computed"),
    pLink: M("Isolated column at the hovered node  (d1 ↔ d2 link)", [
        "frozen at isolation: the neighbour's erosion rate (z<sub>down</sub> keeps uplifting and eroding as it was), H<sub>down</sub>, Q, A, W, W<sub>f</sub>, ds",
        "running: U(t)·f(s), lithology r(z − U<sub>cum</sub>), slope against z<sub>down</sub>(t), H from the face law against z<sub>down</sub>(t) + H<sub>down</sub>",
        "divergence = z<sub>2D</sub>(t) − z<sub>1D</sub>(t)"],
      "Every node of the profile is also run as an isolated d1 column, in lockstep, from the moment of isolation (Reset, a parameter change, or the Re-isolate button). Move the mouse across the long profile and this panel shows, for the node under the pointer, its real elevation as the profile runs (solid) and the same node isolated with the same laws, everything from its neighbours frozen at the isolation time (dashed). Click to pin a node. How far and how fast the curves part is the answer to “how much behaviour arises from the extra dimension”. Two things to expect: for a profile in steady state every isolated column stays identical to its node (the receiver keeps uplifting and eroding as it was), so any gap is adjacency changing — neighbours speeding up or slowing down, or a terminus moving over the node; and the ELA reaches a node only through the ice flux, which is frozen, so the isolated column is blind to glacial cycles. “Open in d1” shows the column standalone with its constants as sliders.",
      "py/test_d1_node.py: the column reproduces the 2D node exactly when its adjacency is refreshed every step.", "computed")
  };

  // ---- Curve tooltips (hover a line) ------------------------------------
  var CT = GS.ui.curveTip;
  var ct = {
    bed: CT("Bedrock surface z(s)", C.bed, "solid", "The valley floor: uplifted each step, lowered by the river where ice-free, by sliding ice where glaciated, and by threshold failure where steeper than S<sub>c</sub>. The outlet (right) is held at base level.", "z(s, t)"),
    iceFill: CT("Ice", C.iceFill, "fill", "The steady glacier that fits the current bed: mass balance summed down-valley gives the flux, the SIA flux law marched upstream gives the thickness. Recomputed every step.", "H(s) from Γ′ H<sub>f</sub><sup>5</sup> S<sub>f</sub><sup>3</sup> = Q/W<sub>f</sub>"),
    iceSurf: CT("Ice surface z + H", C.ice, "solid", "Top of the glacier. Its slope, not the bed's, drives the flow; mass balance is evaluated at this elevation.", "z<sub>s</sub> = z + H"),
    ela: CT("Equilibrium-line altitude", C.ela, "dashed", "Above it the ice surface gains mass, below it loses. Moves with the glacial cycle when the amplitude is nonzero.", "ELA(t) = ELA<sub>0</sub> − A · ½(1 − cos 2πt/P)"),
    litho: CT("Surface inside the resistant body", C.litho, "solid", "Where the bed lies within the resistant body (material frame): erodibilities are divided by the contrasts here.", "r(s) > ½"),
    u: CT("Uplift rate U(s, t)", C.uplift, "dashed", "The shared forcing at this instant, times the spatial pattern. Where an erosion curve sits on this line that reach is in steady state.", "U(t) · f(s)"),
    Ef: CT("Fluvial erosion rate", C.fluv, "solid", "Stream-power incision on ice-free nodes (and at reduced efficiency under ice thinner than H<sub>f</sub>), solved implicitly from the outlet upstream.", "E<sub>f</sub> = K A<sup>m</sup> S<sup>n</sup> · w(H)"),
    Eg: CT("Glacial erosion rate", C.glac, "solid", "Abrasion proportional to sliding speed under ice, through a three-cell footprint, capped at 2 cm/yr; plus quarrying on convex bed if enabled.", "E<sub>g</sub> = K<sub>g</sub> U<sub>s</sub><sup>l</sup>"),
    Er: CT("Threshold-failure (rockfall) rate", "#8c6d31", "dashed", "Lowering applied where a cell stood steeper than S<sub>c</sub> above its downstream neighbour. Usually zero except at headwalls and steps; clipped to the axis when large.", "excess = z<sub>i</sub> − z<sub>i+1</sub> − S<sub>c</sub> ds"),
    Us: CT("Sliding speed U<sub>s</sub>", C.Us, "solid", "Basal sliding, the share f<sub>s</sub> of the depth-averaged speed q/H (with a minimum thickness H<sub>min</sub>). Erosion follows this curve.", "U<sub>s</sub> = f<sub>s</sub> q / max(H, H<sub>min</sub>)"),
    relief: CT("Relief", C.relief, "solid", "Highest minus lowest bed elevation on the profile, through time.", "max z − min z"),
    vol: CT("Ice volume (scaled)", C.vol, "dashed", "Total ice on the profile, scaled to the panel height (its own maximum reaches the top of the axis).", "Σ H W ds"),
    z2d: CT("Profile node: bed", C.bed, "solid", "The real elevation of the hovered node as the whole profile runs, since the moment of isolation.", "z<sub>2D</sub>(t)"),
    z1d: CT("Isolated column: bed", C.Us, "dashed", "The same node run alone with the same laws, everything from its neighbours frozen at isolation (receiver uplifting and eroding as it was, ice flux and downstream ice fixed). The gap to the solid curve is what the neighbours did.", "z<sub>1D</sub>(t)"),
    ice2d: CT("Profile node: ice surface", C.ice, "solid", "Bed plus ice thickness at the node in the full profile.", "z<sub>2D</sub> + H<sub>2D</sub>"),
    ice1d: CT("Isolated column: ice surface", C.ice, "dashed", "Bed plus the thickness that carries the frozen flux against the frozen downstream surface.", "z<sub>1D</sub> + H<sub>1D</sub>")
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
    { id: "Hf", label: "Fluvial shut-off thickness", type: "range", min: 10, max: 300, step: 10, val: 100, unit: "m", group: "Fluvial", asserted: true, method: m.Hf },
    { id: "phiSub", label: "Subglacial fluvial efficiency", type: "range", min: 0, max: 1, step: 0.05, val: 0, unit: "", group: "Fluvial", asserted: true, method: m.phiSub },
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
    { id: "Hmin", label: "Min sliding thickness", type: "range", min: 0, max: 50, step: 1, val: 10, unit: "m", group: "Glacial erosion", asserted: true, method: m.Hmin },
    { id: "eroSmooth", label: "Erosion footprint (3-cell)", type: "check", val: true, group: "Glacial erosion", asserted: true, method: m.eroSmooth },
    { id: "KqExp", label: "Quarrying K_q (min = off)", type: "range", min: -6, max: -1, step: 0.25, val: -6, log: true, group: "Glacial erosion", asserted: true, method: m.KqExp },
    { id: "lithoType", label: "Body type", type: "select", val: "none", options: [{ v: "none", t: "None" }, { v: "layer", t: "Horizontal layer" }, { v: "slab", t: "Dipping slab" }, { v: "dike", t: "Vertical dike" }], group: "Lithology", asserted: true, method: m.lithoType },
    { id: "lithoTop", label: "Body top", type: "range", min: 0, max: 5, step: 0.1, val: 1.5, unit: "km", group: "Lithology", asserted: true, method: m.lithoTop },
    { id: "lithoThick", label: "Thickness / width", type: "range", min: 0.05, max: 3, step: 0.05, val: 0.3, unit: "km", group: "Lithology", asserted: true, method: m.lithoThick },
    { id: "lithoPos", label: "Dike position / slab dip", type: "range", min: 0, max: 1, step: 0.05, val: 0.5, unit: "", group: "Lithology", asserted: true, method: m.lithoPos },
    { id: "contrastK", label: "Fluvial contrast", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Lithology", asserted: true, method: m.contrastK },
    { id: "contrastKg", label: "Glacial contrast", type: "range", min: 1, max: 20, step: 0.5, val: 5, unit: "×", group: "Lithology", asserted: true, method: m.contrastKg },
    { id: "Sc", label: "Threshold slope (rockfall)", type: "range", min: 0.3, max: 2, step: 0.05, val: 0.8, unit: "", group: "Hillslope", asserted: true, method: m.Sc },
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
      Kq: state.KqExp <= -6 ? 0 : Math.pow(10, state.KqExp), Hmin: state.Hmin, eroSmooth: !!state.eroSmooth, Hf: state.Hf, phiSub: state.phiSub, Sc: state.Sc, litho: lithoSpec(), contrastK: state.contrastK, contrastKg: state.contrastKg,
      dt: state.dt
    });
    model.primeIce();
    model.record();
  }

  // ---- Panels --------------------------------------------------------------
  // Three panel factories at one pixel scale: the profile (820 wide) and the side panel (420)
  // share the top row; the full-width panels below are 1246 = 820 + 6 + 420 so nothing rescales.
  var P = GS.ui.panels({ width: 820, height: 170 });
  var PS = GS.ui.panels({ width: 420, height: 170, margin: { top: 18, right: 56, bottom: 8, left: 52 } });
  var PW = GS.ui.panels({ width: 1246, height: 170 });
  var charts = d3.select("#charts");
  var topRow = charts.append("div").attr("class", "chart-row");
  var mainCol = topRow.append("div").attr("class", "main"), sideCol = topRow.append("div").attr("class", "side");
  var pProfile = P.make(mainCol, { title: "Long profile  (hover: isolate that node →)", yLabel: "z (km)", color: C.bed, method: m.pProfile });
  var pLink = PS.make(sideCol, { title: "Isolated vs profile at the hovered node", yLabel: "z (km)", color: C.Us, showX: true, xLabel: "time since isolation (kyr)", method: m.pLink });
  var pRates = PW.make(charts, { title: "Erosion rates  vs  uplift", yLabel: "rate (mm/yr)", color: C.glac, method: m.pRates });
  var pIce = PW.make(charts, { title: "Sliding speed", yLabel: "U_s (m/yr)", color: C.Us, showX: true, xLabel: "distance from divide s (km)", method: m.pIce });
  var pHist = PW.make(charts, { title: "History", yLabel: "relief (km)", color: C.relief, showX: true, xLabel: "time (Myr)", method: m.pHist });
  var sKm = function (d) { return U.toKm(d.s); };
  // ---- Isolated columns for EVERY node, stepped in lockstep with the profile ----
  // iso = { t0, cols[], rec: [{t, z2d: Float32Array, z1d, H2d, H1d}] }. Hover the profile to see
  // the hovered node's pair of curves; click to pin; "Re-isolate" re-bases all columns on the
  // current state. Memory: 4 arrays of N per record; records thin by 2 beyond 3000.
  var iso = null, hoverNode = -1, pinnedNode = -1;
  // Axis ratchet: a panel's range may grow but never shrink during a run (cleared on Reset),
  // so a transient spike (e.g. at the glacier front) does not make the axes jump.
  var seen = {};
  function ratchetMax(key, v) { if (seen[key] == null || v > seen[key]) seen[key] = v; return seen[key]; }
  function ratchetMin(key, v) { if (seen[key] == null || v < seen[key]) seen[key] = v; return seen[key]; }

  function isolateAll() {
    var st = model.state, N = st.N, cols = [];
    for (var i = 0; i < N - 1; i++) cols.push(GS.d1.createNodeColumn(model.nodeSpec(i)));
    iso = { t0: st.t, cols: cols, rec: [] };
    recordIso();
  }
  function recordIso() {
    if (!iso) return;
    var st = model.state, N = st.N;
    var r = { t: st.t - iso.t0, z2d: new Float32Array(N), z1d: new Float32Array(N), H2d: new Float32Array(N), H1d: new Float32Array(N) };
    for (var i = 0; i < N - 1; i++) {
      var c = iso.cols[i].state;
      r.z2d[i] = st.z[i]; r.z1d[i] = c.z; r.H2d[i] = st.H[i]; r.H1d[i] = c.H;
    }
    iso.rec.push(r);
    if (iso.rec.length > 3000) iso.rec = iso.rec.filter(function (_, k) { return k % 2 === 0; });
  }
  function stepIso(nsteps) {
    if (!iso) return;
    for (var k = 0; k < nsteps; k++) for (var i = 0; i < iso.cols.length; i++) iso.cols[i].step();
    recordIso();
  }
  function shownNode() { return hoverNode >= 0 ? hoverNode : pinnedNode; }
  function d1Link(i) {
    if (!iso || i < 0) return "#";
    var sp = iso.cols[i].spec, q = new URLSearchParams();
    q.set("scenario", "node");
    ["shape", "peakUplift", "duration", "KExp", "m", "nexp", "KgExp", "lexp", "fs", "flow", "Hmin", "Hf", "phiSub", "contrastK", "contrastKg", "lithoType", "lithoTop", "lithoThick", "lithoPos", "L", "elaBase", "elaAmp", "elaPeriod", "dt"].forEach(function (id) { q.set(id, state[id]); });
    q.set("s", U.toKm(sp.s).toFixed(3)); q.set("z0", U.toKm(sp.z0).toFixed(4)); q.set("zDown", U.toKm(sp.zDown).toFixed(4)); q.set("HDown", sp.HDown.toFixed(2)); q.set("eDown", U.toMmyr(sp.eDown).toFixed(3));
    q.set("QExp", sp.Q > 0 ? Math.log10(sp.Q).toFixed(3) : 0); q.set("AExp", Math.log10(sp.A / 1e6).toFixed(3)); q.set("W", sp.W.toFixed(1)); q.set("Wf", sp.Wf.toFixed(1));
    q.set("ds", sp.ds.toFixed(1)); q.set("fac", sp.fac.toFixed(4)); q.set("ucum0", sp.ucum0.toFixed(3)); q.set("t0", sp.t0.toFixed(0)); q.set("runKyr", 500);
    return "../d1/index.html?" + q.toString();
  }
  function drawLinkPanel(ni) {
    var st = model.state;
    PS.clear(pLink);
    if (!iso || ni < 0 || iso.rec.length < 2) {
      pLink.x.domain([0, 1]); pLink.y.domain([0, 1]); PS.axes(pLink);
      d3.select("#linkReadout").html(iso ? "Move the mouse across the long profile to isolate a node; click to pin it. Columns isolated at t = " + U.fmtTime(iso.t0) + "." : "");
      return;
    }
    var ls = iso.rec.map(function (r) { return { t: r.t, z2d: r.z2d[ni], z1d: r.z1d[ni], H2d: r.H2d[ni], H1d: r.H1d[ni] }; });
    var tk = function (d) { return U.toKyr(d.t); }, last = ls[ls.length - 1];
    pLink.x.domain([0, Math.max(U.toKyr(last.t), 1)]);
    var lo = ratchetMin("linkLo" + ni, d3.min(ls, function (d) { return U.toKm(Math.min(d.z2d, d.z1d)); }));
    var hi = ratchetMax("linkHi" + ni, d3.max(ls, function (d) { return U.toKm(Math.max(d.z2d + d.H2d, d.z1d + d.H1d)); }));
    pLink.y.domain([lo - 0.01, hi + 0.01]);
    PS.line(pLink, "ice2d", ls, function (d) { return U.toKm(d.z2d + d.H2d); }, C.ice, false, tk, ct.ice2d);
    PS.line(pLink, "ice1d", ls, function (d) { return U.toKm(d.z1d + d.H1d); }, C.ice, true, tk, ct.ice1d);
    PS.line(pLink, "z2d", ls, function (d) { return U.toKm(d.z2d); }, C.bed, false, tk, ct.z2d);
    PS.line(pLink, "z1d", ls, function (d) { return U.toKm(d.z1d); }, C.Us, true, tk, ct.z1d);
    PS.endLabel(pLink, "z2d", last, function (d) { return U.toKm(d.z2d); }, C.bed, "profile");
    PS.endLabel(pLink, "z1d", last, function (d) { return U.toKm(d.z1d); }, C.Us, "isolated");
    PS.axes(pLink, 5, 4);
    var c = iso.cols[ni].state;
    d3.select("#linkReadout").html((hoverNode >= 0 ? "Hovering" : "Pinned") + " node at <b>" + U.fmtLen(st.s[ni]) + "</b> · after " + U.fmtTime(last.t) + ": profile z = <b>" + U.fmtLen(last.z2d) +
      "</b>, isolated z = <b>" + U.fmtLen(last.z1d) + "</b> · divergence <b>" + (last.z2d - last.z1d).toFixed(1) + " m</b> · erosion now " + U.fmtRate(st.Ef[ni] + st.Eg[ni]) + " vs " + U.fmtRate(c.Ef + c.Eg) +
      " · ice " + last.H2d.toFixed(0) + " vs " + last.H1d.toFixed(0) + " m · <a href='" + d1Link(ni) + "' target='_blank'>open in d1 ↗</a>");
  }

  function draw() {
    var st = model.state, N = st.N, i;
    var rows = [];
    for (i = 0; i < N; i++) rows.push({ s: st.s[i], z: st.z[i], zs: st.z[i] + st.H[i], H: st.H[i], Us: st.Us[i], Ef: st.Ef[i], Eg: st.Eg[i], Er: st.Er[i], r: st.r[i], ela: model.ELA(st.t), u: model.U(st.t) });
    var Lkm = U.toKm(st.p.L);
    var zMax = ratchetMax("z", d3.max(rows, function (d) { return U.toKm(Math.max(d.zs, state.glacierOn ? d.ela : 0)); }) || 1);
    var eMax = ratchetMax("e", d3.max(rows, function (d) { return U.toMmyr(Math.max(d.Ef, d.Eg, d.u)); }) || 0.1);
    var usMax = ratchetMax("us", d3.max(rows, function (d) { return d.Us; }) || 1);
    pProfile.x.domain([0, Lkm]); P.clear(pProfile);
    [pRates, pIce].forEach(function (pn) { pn.x.domain([0, Lkm]); PW.clear(pn); });
    pProfile.y.domain([0, zMax * 1.05]); pRates.y.domain([0, eMax * 1.1]); pIce.y.domain([0, usMax * 1.1]);
    // profile
    P.area(pProfile, "ice", rows, function (d) { return U.toKm(d.z); }, function (d) { return U.toKm(d.zs); }, C.iceFill, sKm, ct.iceFill);
    if (state.glacierOn) P.line(pProfile, "ela", rows, function (d) { return U.toKm(d.ela); }, C.ela, true, sKm, ct.ela);
    P.line(pProfile, "ice", rows, function (d) { return U.toKm(d.zs); }, C.ice, false, sKm, ct.iceSurf);
    P.line(pProfile, "bed", rows, function (d) { return U.toKm(d.z); }, C.bed, false, sKm, ct.bed);
    // lithology overlay: segments where r > 0.5
    var segs = [], cur = null;
    rows.forEach(function (d) { if (d.r > 0.5) { if (!cur) { cur = []; segs.push(cur); } cur.push(d); } else cur = null; });
    var lsel = pProfile.g.selectAll("path.line.litho").data(segs);
    var lEnter = lsel.enter().append("path").attr("class", "line litho").attr("fill", "none").attr("stroke", C.litho).attr("stroke-width", 5).attr("stroke-opacity", 0.6).style("cursor", "help");
    GS.ui.tip.attach(lEnter, ct.litho);
    lEnter.merge(lsel).attr("d", d3.line().x(function (d) { return pProfile.x(sKm(d)); }).y(function (d) { return pProfile.y(U.toKm(d.z)); }));
    lsel.exit().remove();
    // rates
    PW.line(pRates, "u", rows, function (d) { return U.toMmyr(d.u); }, C.uplift, true, sKm, ct.u);
    PW.line(pRates, "Ef", rows, function (d) { return U.toMmyr(d.Ef); }, C.fluv, false, sKm, ct.Ef);
    PW.line(pRates, "Eg", rows, function (d) { return U.toMmyr(d.Eg); }, C.glac, false, sKm, ct.Eg);
    PW.line(pRates, "Er", rows, function (d) { return U.toMmyr(Math.min(d.Er, U.fromMmyr(eMax))); }, "#8c6d31", true, sKm, ct.Er);
    // ice
    PW.line(pIce, "Us", rows, function (d) { return d.Us; }, C.Us, false, sKm, ct.Us);
    var ni = shownNode();
    if (ni >= 0) { P.vline(pProfile, "node", U.toKm(st.s[ni])); PW.vline(pRates, "node", U.toKm(st.s[ni])); PW.vline(pIce, "node", U.toKm(st.s[ni])); }
    P.axes(pProfile); PW.axes(pRates); PW.axes(pIce);
    drawLinkPanel(ni);
    // history
    var h = st.history, tMyr = function (d) { return U.toMyr(d.t); };
    PW.clear(pHist);
    pHist.x.domain([0, Math.max(U.toMyr(h[h.length - 1].t), 0.01)]);
    var rMax = ratchetMax("relief", d3.max(h, function (d) { return U.toKm(d.relief); }) || 1), vMax = ratchetMax("vol", d3.max(h, function (d) { return d.iceVol; }) || 1);
    pHist.y.domain([0, rMax * 1.1]);
    PW.line(pHist, "relief", h, function (d) { return U.toKm(d.relief); }, C.relief, false, tMyr, ct.relief);
    PW.line(pHist, "vol", h, function (d) { return rMax * 1.1 * d.iceVol / (vMax * 1.1); }, C.vol, true, tMyr, ct.vol);
    PW.axes(pHist);
    var dg = model.diagnostics();
    d3.select("#readout").html(
      "t = <b>" + U.fmtTime(dg.t) + "</b>  ·  relief <b>" + U.fmtLen(dg.relief) + "</b>  ·  max ice <b>" + dg.maxH.toFixed(0) + " m</b>" +
      (dg.terminus >= 0 ? "  ·  terminus <b>" + U.fmtLen(dg.terminus) + "</b>" : "  ·  no ice") +
      "  ·  ice volume " + (dg.iceVol / 1e9).toFixed(2) + " km³  ·  ice iters " + st.iceIters +
      (st.steadyCapped > 0 ? "  ·  <span style='color:#b8860b'>steady profile at threshold slope over " + U.fmtLen(st.steadyCapped * st.ds) + " (K too low for U: not in fluvial equilibrium there)</span>" : "") +
      "  ·  mass closure " + (dg.upliftVol > 0 ? (dg.closure / dg.upliftVol).toExponential(1) : "—"));
    GS.ui.holdHeight(d3.select("#readout")); GS.ui.holdHeight(d3.select("#linkReadout"));
  }

  // ---- Loop --------------------------------------------------------------
  var playing = false, raf = null;
  function frame() {
    if (!playing) return;
    for (var k = 0; k < state.speed; k++) model.step();
    model.record(); stepIso(state.speed);
    draw();
    raf = requestAnimationFrame(frame);
  }
  function setPlaying(v) {
    playing = v;
    d3.select("#playPause").text(playing ? "❚❚ Pause" : "▶ Play");
    if (playing) raf = requestAnimationFrame(frame); else if (raf) cancelAnimationFrame(raf);
  }
  function reset() { setPlaying(false); seen = {}; pinnedNode = -1; hoverNode = -1; GS.ui.holdHeight(d3.selectAll("#readout, #linkReadout"), true); build(); isolateAll(); draw(); url.write(state); }
  d3.select("#resetBtn").on("click", reset);
  d3.select("#reisoBtn").on("click", function () { Object.keys(seen).forEach(function (k) { if (k.indexOf("link") === 0) delete seen[k]; }); isolateAll(); draw(); });
  d3.select("#playPause").on("click", function () { setPlaying(!playing); });
  d3.select("#stepBtn").on("click", function () { setPlaying(false); for (var k = 0; k < state.speed; k++) model.step(); model.record(); stepIso(state.speed); draw(); });
  // hover / click on the profile: show / pin the node under the pointer
  function nodeUnder(ev) {
    var xy = d3.pointer(ev, pProfile.g.node());
    var i = Math.round(U.fromKm(pProfile.x.invert(xy[0])) / model.state.ds);
    if (i < 0) i = 0; if (i > model.state.N - 2) i = model.state.N - 2;
    return i;
  }
  pProfile.svg.style("cursor", "crosshair")
    .on("mousemove", function (ev) { var i = nodeUnder(ev); if (i !== hoverNode) { hoverNode = i; draw(); } })
    .on("mouseleave", function () { hoverNode = -1; draw(); })
    .on("click", function (ev) { pinnedNode = nodeUnder(ev); hoverNode = -1; draw(); });

  GS.ui.buildControls(d3.select("#controls"), specs, state, function (id) {
    // domain / initial-condition changes need a rebuild; everything else applies live via rebuild-in-place of params
    if (["L", "N", "initProfile", "zHead", "noise"].indexOf(id) !== -1 || (state.initProfile === "steady" && model.state.t === 0 && ["KExp", "m", "nexp", "peakUplift", "pattern", "hack"].indexOf(id) !== -1)) { reset(); return; }
    var wasPlaying = playing;
    var keepT = model.state.t, keepZ = model.state.z, keepU = model.state.ucum, keepH = model.state.H, keepHist = model.state.history;
    build();
    // carry the evolving state across a parameter change (same grid)
    model.state.z.set(keepZ); model.state.ucum.set(keepU); model.state.H.set(keepH); model.state.t = keepT; model.state.history = keepHist;
    isolateAll(); // the columns must share the profile's laws: re-base them on the current state
    draw(); url.write(state);
    if (wasPlaying) setPlaying(true);
  });
  window.GS_d2along_model = function () { return model; }; // for validation from the console
  build(); isolateAll(); draw();
})();
