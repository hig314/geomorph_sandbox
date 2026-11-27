/*
 * relief_limit/glacier_model.js — "Glacier buttressing": relief against a threshold that a
 * glacier both buttresses and weakens, over two glacial cycles (the last 200 kyr). Pure
 * model, no DOM. A stand-alone slide beside the seven-step figure of model.js.
 *
 * Time t in kyr from 200 ka (t = 0) to the present (t = 200). Relief and threshold in units
 * of the intact threshold with no ice, T0 = 1. Ice thickness H in units of the largest
 * glacial thickness.
 *
 *   ice         H(t) = max(Hmin, envelope(t) + fluctuation(t)), a polyline on a dtH grid.
 *               envelope: glacial (MIS 6) – interglacial (MIS 5e) – glacial (to the LGM) –
 *               interglacial (Holocene), piecewise linear through the keys of the ice history
 *               p.ice: "thick" (ICE_KEYS) never zero, so ice is present throughout; "retreat"
 *               (RETREAT_KEYS) gone early in both interglacials, re-advancing in the late
 *               Holocene, gone at the end. fluctuation: internal dynamics, seeded, knots every
 *               1.5–4 kyr with Gaussian values σ·(0.4 + 0.6·envelope)·min(1, envelope/taper)
 *               (no taper for "thick"), cosine-interpolated
 *   threshold   T(t) = T0 − Dmg + B·b(t)
 *               buttressing  b = (1 − e^(−H/Hs)) / (1 − e^(−1/Hs)), Hs = 0.3, B fixed, at the
 *                            nodes of the ice polyline, linear between them: saturating, so the
 *                            response is stretched for thin ice and compressed for thick ice
 *               weakening    the stress change is the change in buttressing |d(B·b)|; it
 *                            accumulates from the last fracture, and when it reaches a quota
 *                            q = qMean·(1 + qJit·(2u − 1)) (u uniform, its own stream) the rock
 *                            fractures: Dmg += W·q, a step down in T
 *               a landslide removes the damaged rock: Dmg ← 0
 *   relief      dR/dt = E where there is ice (glacial erosion, constant), 0 on an ice-free grid
 *               cell (H = 0 at both ends); with "thick" ice is present the whole time
 *   failure     where R reaches T (erosion, or T falling as the ice thins) or a fracture steps
 *               T to or below R:  R ← max(T − d, floor), then Dmg ← 0; d from a truncated
 *               power law on [dMin, dMax] (its own stream, as model.js)
 *
 * H and b are linear on each grid cell, so T and R are linear between events and every event time
 * is exact. Mirrored in py/glacier_buttress.py; validated in py/test_glacier_buttress.py
 * (budget, threshold reconstruction, fracture quotas against ∫|d(B·b)|, failure rule, floor,
 * W = 0 and E = 0 limits, a converging time-stepped reference, fingerprints).
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});
  var EPS = 1e-9;

  // [age ka, envelope]: build-up to MIS 6, termination II, the last interglacial, the stepped
  // build-up of the last glacial to the LGM, termination I, the Holocene, and climate change:
  // the ice plunges in the last 600 years. A thick glacier: interglacials keep about 40 % of
  // the glacial maximum, so the plunge takes it lower than any interglacial.
  var ICE_KEYS = [[200, 0.65], [186, 0.73], [176, 0.69], [160, 0.80], [148, 0.86], [140, 0.97],
                  [135, 0.93], [129, 0.41], [117, 0.38], [110, 0.57], [103, 0.50], [92, 0.64],
                  [85, 0.57], [70, 0.73], [60, 0.79], [50, 0.71], [35, 0.76], [24, 1.00],
                  [19, 0.94], [11, 0.40], [0.6, 0.38], [0, 0.04]];
  // The same glacials, but the glacier retreats completely early in both interglacials: gone by
  // 128 ka, back briefly in a small stade (125–122 ka, to 0.28 at 123.5 ka), gone again until it
  // re-advances from 118 ka, gone by 10 ka until a late-Holocene re-advance from
  // 4.5 ka to a maximum at 0.3 ka (Little Ice Age), then the final retreat.
  var RETREAT_KEYS = [[200, 0.65], [186, 0.73], [176, 0.69], [160, 0.80], [148, 0.86], [140, 0.97],
                      [135, 0.93], [128, 0], [125, 0], [123.5, 0.28], [122, 0], [118, 0], [110, 0.50], [103, 0.45], [92, 0.62],
                      [85, 0.55], [70, 0.73], [60, 0.79], [50, 0.71], [35, 0.76], [24, 1.00],
                      [19, 0.94], [10, 0], [4.5, 0], [0.3, 0.45], [0, 0]];
  var ICE = { thick: { keys: ICE_KEYS, Hmin: 0.05, taper: 0 },
              retreat: { keys: RETREAT_KEYS, Hmin: 0, taper: 0.15 } };

  var DEFAULTS = {
    tmax: 200,     // kyr
    dtH: 0.25,     // kyr, the ice polyline's grid
    T0: 1,         // intact threshold with no ice
    B: 0.5,        // buttressing at the largest thickness (fixed)
    Hs: 0.3,       // saturation thickness of buttressing: stretched for thin ice, compressed for thick
    ice: "thick",  // ice history: "thick" or "retreat" (ICE)
    sigma: 0.09,   // amplitude of the internal fluctuations
    E: 1,          // erosion, as a multiple of Eref
    Eref: 0.012,   // relief per kyr at E = 1
    W: 1,          // weakening: drop of T per unit change in buttressing at a fracture
    qMean: 0.025,  // mean change in buttressing between fractures (threshold units)
    qJit: 0.25,    // quota spread: q uniform in qMean·(1 ± qJit)
    R0: 1.3,       // above the mean relief and the ice-free threshold: a slope held up by the ice
    floor: 0.7,    // frictional floor: broken rock still stands at its friction angle, so no landslide
                   // leaves less relief than this (py/xsection_experiment.py); was 0.165 until 2026-11-27
    dMin: 0.05, dMax: 0.5, dAlpha: 1,
    seed: 7
  };

  function gauss(rnd) {
    var u = 1 - rnd(), v = rnd();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function envelope(t, tmax, keys) {
    var a = tmax - t, i, a0, a1;
    keys = keys || ICE_KEYS;
    for (i = 0; i + 1 < keys.length; i++) {
      a0 = keys[i][0]; a1 = keys[i + 1][0];
      if (a <= a0 && a >= a1) return keys[i][1] + (keys[i + 1][1] - keys[i][1]) * (a0 - a) / (a0 - a1);
    }
    return a < keys[keys.length - 1][0] ? keys[keys.length - 1][1] : keys[0][1];
  }

  function makeIce(p) {
    var rnd = GS.noise.mulberry32(p.seed | 0), sc = ICE[p.ice], tmax = p.tmax, kt = [-2 * rnd()], kv = [], i, j = 0, n, t, s, w, env, amp, out = [];
    while (kt[kt.length - 1] <= tmax) kt.push(kt[kt.length - 1] + 1.5 + 2.5 * rnd());
    for (i = 0; i < kt.length; i++) kv.push(gauss(rnd));
    n = Math.round(tmax / p.dtH);
    for (i = 0; i <= n; i++) {
      t = tmax * i / n;
      while (j + 2 < kt.length && kt[j + 1] <= t) j++;
      s = (t - kt[j]) / (kt[j + 1] - kt[j]); w = 0.5 - 0.5 * Math.cos(Math.PI * s);
      env = envelope(t, tmax, sc.keys);
      amp = p.sigma * (0.4 + 0.6 * env) * (sc.taper > 0 ? Math.min(1, env / sc.taper) : 1);
      out.push([t, Math.max(sc.Hmin, env + amp * (kv[j] + (kv[j + 1] - kv[j]) * w))]);
    }
    return out;
  }

  function makeSizes(p) {
    var lo = Math.min(p.dMin, p.dMax), hi = Math.max(p.dMin, p.dMax), a = p.dAlpha;
    if (hi <= lo * (1 + 1e-9) || a <= 0) return function () { return lo; };
    var rnd = GS.noise.mulberry32((p.seed | 0) + 1000003), span = 1 - Math.pow(lo / hi, a);
    return function () { return lo * Math.pow(1 - rnd() * span, -1 / a); };
  }

  function makeQuotas(p) {
    var rnd = GS.noise.mulberry32((p.seed | 0) + 2000003), qm = p.qMean, j = p.qJit;
    return function () { return qm * (1 + j * (2 * rnd() - 1)); };
  }

  // buttressing per unit B: saturating in H, 0 at H = 0 and 1 at the largest thickness
  function butt(H, p) { return (1 - Math.exp(-H / p.Hs)) / (1 - Math.exp(-1 / p.Hs)); }

  function simulate(params) {
    var p = Object.assign({}, DEFAULTS, params || {});
    var ice = makeIce(p), nextD = makeSizes(p), nextQ = makeQuotas(p);
    var E = p.E * p.Eref, B = p.B, T0 = p.T0, W = p.W, floor = Math.max(p.floor, 0);
    var Dmg = 0, acc = 0, q = nextQ(), t = 0, H = ice[0][1], Bt = B * butt(H, p), Bt1, sb;
    var R = Math.min(p.R0, T0 + Bt);
    var path = [[t, R]], thr = [[t, T0 + Bt]], fracs = [], fails = [];
    var i, t1, H1, s, Ec, gap, closing, tf, tq, tn, dt, Told;

    function fail(cause) {
      var T = T0 - Dmg + Bt, d = nextD(), Rn = Math.max(T - d, floor);
      if (Rn < R) {
        fails.push({ t: t, from: R, to: Rn, d: d, lim: T, dmg: Dmg, cause: cause });
        R = Rn; path.push([t, R]);
      }
      if (Dmg > 0) { thr.push([t, T]); Dmg = 0; thr.push([t, T0 + Bt]); } // the landslide removes the damaged rock
    }

    for (i = 0; i + 1 < ice.length; i++) {
      t1 = ice[i + 1][0]; H1 = ice[i + 1][1];
      s = t1 > t ? (H1 - H) / (t1 - t) : 0;       // dH/dt on this cell
      Bt1 = B * butt(H1, p);
      sb = t1 > t ? (Bt1 - Bt) / (t1 - t) : 0;   // d(buttressing)/dt on this cell
      Ec = H > 0 || H1 > 0 ? E : 0;              // glacial erosion only where there is ice
      for (;;) {
        gap = T0 - Dmg + Bt - R;                // ≥ 0
        closing = Ec - sb;                         // d(R − T)/dt
        tf = closing > 0 ? t + gap / closing : Infinity;
        tq = sb !== 0 ? t + (q - acc) / Math.abs(sb) : Infinity;
        tn = Math.min(tf, tq, t1);
        if (tn >= t1 - EPS && tf > t1 - EPS && tq > t1 - EPS) {
          acc += Math.abs(sb) * (t1 - t); R += Ec * (t1 - t); t = t1; H = H1; Bt = Bt1;
          R = Math.min(R, T0 - Dmg + Bt);
          break;
        }
        dt = tn - t;
        acc += Math.abs(sb) * dt; R += Ec * dt; H += s * dt; Bt += sb * dt; t = tn;
        path.push([t, R]); thr.push([t, T0 - Dmg + Bt]);
        if (tq <= tf) {                            // fracture
          Told = T0 - Dmg + Bt;
          Dmg += W * q; acc = 0;
          fracs.push({ t: t, drop: W * q, q: q, from: Told, to: Told - W * q });
          q = nextQ();
          thr.push([t, T0 - Dmg + Bt]);
          if (R >= T0 - Dmg + Bt - EPS) fail("fracture");
        } else {
          R = T0 - Dmg + Bt;
          fail("limit");
        }
      }
      path.push([t, R]); thr.push([t, T0 - Dmg + Bt]);
    }
    return { path: path, thr: thr, ice: ice, fracs: fracs, fails: fails, p: p };
  }

  // Value at t on a time-ordered polyline (after any step at exactly t).
  function valueAt(poly, t) {
    var lo = 0, hi = poly.length - 1, mid, a, b;
    while (lo < hi) { mid = (lo + hi + 1) >> 1; if (poly[mid][0] <= t) lo = mid; else hi = mid - 1; }
    if (lo === poly.length - 1) return poly[lo][1];
    a = poly[lo]; b = poly[lo + 1];
    return b[0] === a[0] ? a[1] : a[1] + (b[1] - a[1]) * (t - a[0]) / (b[0] - a[0]);
  }

  // ---- The cross-section: a glacial trough read from the same history ---------------------
  // A half-section z(x), valley centre at x = 0, crest at x = 1 (x in units of the half-spacing,
  // z in units of T0), mirrored about the centre and the crest when drawn. Two shapes, and the
  // surface is the lower of them: a U cut by the glacier from the floor up, and a straight slope
  // left by landslides from the crest down.
  //   incision    the glacier cuts the floor down and a parabolic trough z_floor + C·x² (C the
  //               highest threshold, T0 + B: a full trough stands at it) climbs from it; above
  //               where it meets the older surface, that surface is untouched:
  //               z(x) = min(relict(x), z_floor + C·x²),  z_floor = z_crest − R(t)
  //               (the parabola as a polyline on SEC_N cells; the min is exact for polylines).
  //               The U climbs the wall as the floor goes down: the glacier oversteepens the wall
  //   landslide   R → R′: the crest comes down by R − R′ and the floor stays. The failure surface
  //               is a plane from the new crest down to the wall, at the slope SEC_S or, where a
  //               plane that steep would miss the wall, at the steepest slope that still meets it
  //               (s_tan: the tangent from the new crest to the trough):
  //               z(x) = min(old(x), z_crest′ − s·(1 − x))  for x ≥ x_cut,   s = min(SEC_S, s_tan),
  //               x_cut the highest point where the plane meets the old surface; old(x) below it.
  //               A wedge off the upper wall: the trough below x_cut is never touched, so the
  //               floor stays a U. The debris is not kept in the section (it lands on the glacier
  //               and is carried away)
  //   ice         the glacier fills the trough to min(SEC_ICE·C·H(t), SEC_CAP·R) above the floor,
  //               flat-topped: never over the crest
  //   weakening   damage D / T_intact, drawn as a lens under the upper wall (glacier_app.js)
  //   datum       the floor: drawn at −SEC_DATUM·C, so a full trough is about centred
  // Mirrored in py/glacier_buttress.py and checked there. Nothing here feeds back on R(t).
  var SEC_N = 32, SEC_ICE = 0.6, SEC_CAP = 0.9, SEC_S = 1.2, SEC_DATUM = 0.45;

  function secArea(poly, xTo) {
    var a = 0, i, p0, p1, zc;
    if (xTo === undefined) xTo = 1;
    for (i = 0; i + 1 < poly.length; i++) {
      p0 = poly[i]; p1 = poly[i + 1];
      if (xTo <= p0[0]) break;
      if (xTo < p1[0]) {
        zc = p0[1] + (p1[1] - p0[1]) * (xTo - p0[0]) / (p1[0] - p0[0]);
        return a + 0.5 * (p0[1] + zc) * (xTo - p0[0]);
      }
      a += 0.5 * (p0[1] + p1[1]) * (p1[0] - p0[0]);
    }
    return a;
  }
  function secValue(poly, x) {
    var i, p0, p1;
    for (i = 0; i + 1 < poly.length; i++) {
      p0 = poly[i]; p1 = poly[i + 1];
      if (x <= p1[0]) return p1[0] === p0[0] ? p0[1] : p0[1] + (p1[1] - p0[1]) * (x - p0[0]) / (p1[0] - p0[0]);
    }
    return poly[poly.length - 1][1];
  }
  // Pointwise minimum of two polylines on [0, 1], crossings inserted (exact).
  function secMin(a, b) {
    var seen = {}, xs = [], i, x, x0, za, zb, d0, d1, xc, out = [], clean, v;
    a.concat(b).forEach(function (p) { if (!seen[p[0]]) { seen[p[0]] = true; xs.push(p[0]); } });
    xs.sort(function (u, w) { return u - w; });
    for (i = 0; i < xs.length; i++) {
      x = xs[i]; za = secValue(a, x); zb = secValue(b, x);
      if (i) {
        x0 = xs[i - 1]; d0 = secValue(a, x0) - secValue(b, x0); d1 = za - zb;
        if (d0 * d1 < 0) { xc = x0 + (x - x0) * d0 / (d0 - d1); out.push([xc, secValue(a, xc)]); }
      }
      out.push([x, Math.min(za, zb)]);
    }
    clean = [out[0]];
    for (i = 1; i < out.length; i++) {
      v = out[i];
      if (v[0] - clean[clean.length - 1][0] > 1e-12) clean.push(v);
      else clean[clean.length - 1] = [clean[clean.length - 1][0], Math.min(clean[clean.length - 1][1], v[1])];
    }
    return clean;
  }
  function secTrough(zv, C) {
    var out = [], i;
    for (i = 0; i <= SEC_N; i++) out.push([i / SEC_N, zv + C * (i / SEC_N) * (i / SEC_N)]);
    return out;
  }
  function secSurface(rel, R, C) {
    var zv = rel[rel.length - 1][1] - R;
    if (rel[0][1] - zv <= EPS) return rel.map(function (v) { return v.slice(); });
    return secMin(rel, secTrough(zv, Math.max(C, R)));
  }
  // {rel, xcut, slope}: the wedge above a plane from the new crest (floor + Rn) is removed. The
  // plane has slope min(SEC_S, s_tan), s_tan the steepest slope at which a line from the new
  // crest still meets the old surface below it; xcut is the highest point where it does. Below
  // xcut the old surface is kept.
  function secLandslide(cur, Rn) {
    var zf = cur[0][1], zc = zf + Rn, st = -Infinity, sl, nw, i, xc, out = [], j;
    for (j = 0; j < cur.length; j++) if (cur[j][0] < 1 - 1e-9) st = Math.max(st, (zc - cur[j][1]) / (1 - cur[j][0]));
    sl = Math.min(SEC_S, st);
    nw = secMin(cur, [[0, zc - sl], [1, zc]]); i = nw.length - 2;
    while (i > 0 && nw[i][1] < secValue(cur, nw[i][0]) - 1e-12) i--;
    xc = nw[i][0];
    for (j = 0; j < cur.length; j++) if (cur[j][0] < xc - 1e-12) out.push(cur[j].slice());
    return { rel: out.concat(nw.slice(i)), xcut: xc, slope: sl };
  }
  function secC(p) { return p.T0 + p.B; }
  // [{t, rel, pre}]: the relict surface from each landslide on; at t = 0 a trough of relief R0
  // under a crest at 0, with the failure slope from the crest down to it.
  function sectionStates(sim) {
    var C = secC(sim.p), R0 = sim.path[0][1], i, f, cur, ls, g, j;
    var rel = secLandslide(secMin([[0, 0], [1, 0]], secTrough(-R0, Math.max(C, R0))), R0).rel;
    var out = [{ t: 0, rel: rel, pre: null, ghost: null }];
    for (i = 0; i < sim.fails.length; i++) {
      f = sim.fails[i];
      if (f.from - f.to <= EPS) continue;
      cur = secSurface(rel, f.from, C); ls = secLandslide(cur, f.to); rel = ls.rel;
      g = [[ls.xcut, secValue(cur, ls.xcut)]];             // ghost: the part of the old surface that failed
      for (j = 0; j < cur.length; j++) if (cur[j][0] > ls.xcut + 1e-12) g.push(cur[j].slice());
      out.push({ t: f.t, rel: rel, pre: cur, xcut: ls.xcut, slope: ls.slope, ghost: g });
    }
    return out;
  }
  // The section at time t: pts (floor at −SEC_DATUM·C), relief R, the surface the last
  // landslide failed from (ghost, same datum) and the time since (age), the damage standing in
  // the rock as a fraction of the intact threshold (weak), the ice H and its level zIce.
  function sectionAt(sim, states, t) {
    var p = sim.p, C = secC(p), k = 0, st, R, cur, m, ice = sim.ice, n = ice.length - 1, i, ba, bb, Ti, H, pts;
    while (k + 1 < states.length && states[k + 1].t <= t) k++;
    st = states[k]; R = valueAt(sim.path, t);
    cur = secSurface(st.rel, R, C); m = cur[0][1] + SEC_DATUM * C;
    function shift(v) { return [v[0], v[1] - m]; }
    pts = cur.map(shift);
    i = Math.min(Math.max(Math.floor(t / p.dtH), 0), n - 1);
    ba = p.B * butt(ice[i][1], p); bb = p.B * butt(ice[i + 1][1], p);
    Ti = p.T0 + ba + (bb - ba) * (t - ice[i][0]) / (ice[i + 1][0] - ice[i][0]);
    H = valueAt(ice, t);
    // rock uplift relative to the floor so far = glacial erosion so far (the texture's offset)
    var up = R - sim.path[0][1], j;
    for (j = 0; j < sim.fails.length && sim.fails[j].t <= t; j++) up += sim.fails[j].from - sim.fails[j].to;
    return { up: up, pts: pts, R: R, ghost: st.ghost ? st.ghost.map(shift) : null, age: t - st.t,
             weak: Math.max(0, (Ti - valueAt(sim.thr, t)) / Ti), H: H, zIce: pts[0][1] + Math.min(SEC_ICE * C * H, SEC_CAP * R) };
  }
  // Same line as py/glacier_buttress.py section_fingerprint() prints.
  function sectionFingerprint(extra) {
    var s = simulate(extra), states = sectionStates(s), out = String(states.length - 1);
    [37.3, 101.9, 171.3, 199.9].forEach(function (t) {
      var sec = sectionAt(s, states, t), pts = sec.pts, sum = 0, i;
      for (i = 0; i < pts.length; i++) sum += Math.abs(pts[i][1]) + pts[i][0];
      out += " " + sum.toFixed(9) + " " + pts.length + " " + (sec.zIce + sec.weak).toFixed(9);
      out += " " + (sec.ghost ? sec.ghost.length : 0);
    });
    return out;
  }

  // Same line as py/glacier_buttress.py fingerprint() prints.
  function fingerprint(extra) {
    var s = simulate(extra), sp = 0, st = 0, si = 0, i;
    for (i = 0; i < s.path.length; i++) sp += s.path[i][1];
    for (i = 0; i < s.thr.length; i++) st += s.thr[i][1];
    for (i = 0; i < s.ice.length; i++) si += s.ice[i][1];
    return s.path.length + " " + s.thr.length + " " + s.fracs.length + " " + s.fails.length + " " +
      sp.toFixed(9) + " " + st.toFixed(9) + " " + si.toFixed(9);
  }

  GS.glacier = { DEFAULTS: DEFAULTS, ICE_KEYS: ICE_KEYS, ICE: ICE, butt: butt, simulate: simulate, makeIce: makeIce, valueAt: valueAt, fingerprint: fingerprint,
                 sectionStates: sectionStates, sectionAt: sectionAt, sectionFingerprint: sectionFingerprint, SEC_ICE: SEC_ICE, SEC_S: SEC_S };
})(typeof window !== "undefined" ? window : this);
