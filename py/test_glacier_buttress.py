"""
test_glacier_buttress.py — the glacier buttressing slide (relief_limit/glacier_model.js mirror).
Run: /opt/anaconda3/bin/python3 py/test_glacier_buttress.py
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(__file__))
import glacier_buttress as gb

fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  [" + detail + "]" if detail else ""))
    if not ok: fails.append(name)

def ice_at(ice, t):
    n = len(ice) - 1; tmax = ice[-1][0]
    i = min(int(t / tmax * n), n - 1)
    while i > 0 and ice[i][0] > t: i -= 1
    while i < n - 1 and ice[i + 1][0] < t: i += 1
    a, b = ice[i], ice[i + 1]
    return a[1] + (b[1] - a[1]) * (t - a[0]) / (b[0] - a[0])

def butt_poly(s):
    """the buttressing B·b(H) at the ice nodes (the model interpolates it linearly)"""
    return [[v[0], s["p"]["B"] * gb.butt(v[1], s["p"])] for v in s["ice"]]

def stress(ice, ta, tb):
    """exact integral of |dy| over [ta, tb] on a polyline"""
    s = 0.0
    for a, b in zip(ice, ice[1:]):
        lo, hi = max(a[0], ta), min(b[0], tb)
        if hi > lo: s += abs(b[1] - a[1]) / (b[0] - a[0]) * (hi - lo)
    return s

def erosion_on(ice, t):
    """glacial erosion acts on a grid cell with ice at either end (the cell holding t)"""
    n = len(ice) - 1; i = min(max(int(t / ice[-1][0] * n), 0), n - 1)
    return ice[i][1] > 0 or ice[i + 1][1] > 0

def ice_time(ice):
    return sum(b[0] - a[0] for a, b in zip(ice, ice[1:]) if a[1] > 0 or b[1] > 0)

SETS = [{}, {"W": 0}, {"W": 1}, {"E": 0}, {"E": 3}, {"E": 3, "W": 1}, {"seed": 3}, {"seed": 41, "W": 0.6},
        {"dMax": 1.2, "W": 0.8}, {"floor": 0.5, "E": 2}, {"sigma": 0.25, "seed": 11},
        {"ice": "retreat"}, {"ice": "retreat", "W": 0, "E": 0.25}, {"ice": "retreat", "E": 3, "W": 1.5, "seed": 5}]

worst = dict(budget=0, thr=0, quota=0, slope=0)
for c in SETS:
    s = gb.simulate(c); p = s["p"]; path = s["path"]; thr = s["thr"]; ice = s["ice"]
    E = p["E"] * p["Eref"]; tag = str(c)
    # relief: time monotone; every segment grows at exactly E (0 where the ice is gone) or is a vertical drop
    mono = all(b[0] >= a[0] for a, b in zip(path, path[1:]))
    for a, b in zip(path, path[1:]):
        if b[0] > a[0] + 1e-12:
            worst["slope"] = max(worst["slope"], abs((b[1] - a[1]) / (b[0] - a[0]) - (E if erosion_on(ice, 0.5 * (a[0] + b[0])) else 0.0)))
        else:
            mono = mono and b[1] <= a[1] + 1e-12
    check(tag + " time monotone, drops only downward", mono)
    # budget
    drops = sum(f["frm"] - f["to"] for f in s["fails"])
    r = abs(path[0][1] + E * ice_time(ice) - drops - path[-1][1]); worst["budget"] = max(worst["budget"], r)
    # threshold reconstructed from the event lists: T0 + B·H − Dmg, Dmg = Σ fracture drops since the last failure
    ev = sorted([(f["t"], 0, f["drop"]) for f in s["fracs"]] + [(f["t"], 1, 0) for f in s["fails"]])
    ok = True
    for (t, T) in thr:
        dm = 0.0
        for (te, kind, dr) in ev:
            if te < t - 1e-9: dm = 0.0 if kind == 1 else dm + dr
        alt = dm  # just before any event at t; just after a fracture or failure at t the value differs
        dm2 = 0.0
        for (te, kind, dr) in ev:
            if te <= t + 1e-9: dm2 = 0.0 if kind == 1 else dm2 + dr
        base = p["T0"] + ice_at(butt_poly(s), t)
        e = min(abs(T - (base - alt)), abs(T - (base - dm2)))
        # several events at one instant: any intermediate partial sum is a legitimate vertex
        if e > 1e-9:
            part = alt; best = e
            for (te, kind, dr) in ev:
                if abs(te - t) <= 1e-9:
                    part = 0.0 if kind == 1 else part + dr
                    best = min(best, abs(T - (base - part)))
            e = best
        worst["thr"] = max(worst["thr"], e)
    # relief never above the threshold (at every relief vertex, the threshold after events there)
    above = 0.0
    for (t, R) in path:
        Ts = [T for (tt, T) in thr if abs(tt - t) <= 1e-9]
        if Ts: above = max(above, R - max(Ts))
    check(tag + " relief never above the threshold", above < 1e-9, "max excess %.1e" % above)
    # each fracture: the stress change since the previous one equals its quota
    prev = 0.0
    for f in s["fracs"]:
        worst["quota"] = max(worst["quota"], abs(stress(butt_poly(s), prev, f["t"]) - f["q"])); prev = f["t"]
        ok = ok and abs(f["q"] - p["qMean"]) <= p["qJit"] * p["qMean"] + 1e-12
        ok = ok and abs(f["drop"] - p["W"] * f["q"]) < 1e-12
    check(tag + " fracture drop = W × quota", ok)
    # failures: at the threshold, to max(T − d, floor), d within the cut-offs, damage cleared
    okf = all(abs(f["frm"] - f["lim"]) < 1e-9 or f["frm"] >= f["lim"] - 1e-9 for f in s["fails"])
    okf = okf and all(abs(f["to"] - max(f["lim"] - f["d"], p["floor"])) < 1e-12 for f in s["fails"])
    okf = okf and all(p["dMin"] - 1e-12 <= f["d"] <= p["dMax"] + 1e-12 for f in s["fails"])
    check(tag + " failures at the threshold, to max(T − d, floor)", okf, "%d failures, %d fractures" % (len(s["fails"]), len(s["fracs"])))
    check(tag + " relief never below the floor", min(v[1] for v in path) >= min(p["floor"], p["R0"]) - 1e-12)
    check(tag + " ice ≥ Hmin", min(v[1] for v in ice) >= gb.ICE[p["ice"]]["Hmin"] - 1e-15)

print("worst: budget %.1e, threshold reconstruction %.1e, quota vs ∫|d(B·b)| %.1e, growth slope − E %.1e" %
      (worst["budget"], worst["thr"], worst["quota"], worst["slope"]))
check("relief budget R0 + E·(time with ice) − Σdrops = R_end", worst["budget"] < 1e-9)
check("threshold = T0 + B·b(H) − damage since the last failure", worst["thr"] < 1e-9)
check("each fracture after exactly its quota of buttressing change ∫|d(B·b)|", worst["quota"] < 1e-9)
check("relief grows at exactly E with ice, not at all without", worst["slope"] < 1e-9)

# W = 0: no damage, the threshold is the buttressed one; fractures still drawn at the same times
s0 = gb.simulate({"W": 0}); s1 = gb.simulate({})
check("W = 0: threshold = T0 + B·b(H) (node values, linear between) everywhere",
      max(abs(T - (1 + ice_at(butt_poly(s0), t))) for t, T in s0["thr"]) < 1e-9)
# buttressing saturates: b(0) = 0, b(1) = 1, increasing, slope falling with thickness
pp = gb.simulate()["p"]; hs = [i / 100 for i in range(101)]; bs = [gb.butt(h, pp) for h in hs]
sl = [b - a for a, b in zip(bs, bs[1:])]
check("buttressing saturating: b(0) = 0, b(1) = 1, rising, slope falling", abs(bs[0]) < 1e-15 and abs(bs[-1] - 1) < 1e-15 and all(v > 0 for v in sl) and all(b < a for a, b in zip(sl, sl[1:])),
      "slope at H = 0 vs 1: %.2f vs %.2f" % (sl[0] * 100, sl[-1] * 100))
pl = dict(pp); pl["Hs"] = 1e4
check("large Hs is linear buttressing", max(abs(gb.butt(h, pl) - h) for h in hs) < 1e-4)
check("ice independent of W and E", s0["ice"] == s1["ice"] == gb.simulate({"E": 2.5})["ice"])
# E = 0: relief never rises
sE = gb.simulate({"E": 0})
check("E = 0: relief non-increasing", all(b[1] <= a[1] + 1e-12 for a, b in zip(sE["path"], sE["path"][1:])))
# with no weakening every failure is at the threshold reached by erosion or by thinning ice
check("W = 0: no failure caused by a fracture", all(f["cause"] == "limit" for f in s0["fails"]))
th = [f for f in s0["fails"] if f["cause"] == "limit"]
thin = 0
for f in th:
    H0 = ice_at(s0["ice"], max(f["t"] - 0.25, 0)); H1 = ice_at(s0["ice"], min(f["t"] + 0.25, 200))
    thin += H1 < H0
check("W = 0: failures fall where the ice is thinning (debuttressing), most of them", thin >= 0.6 * len(th), "%d of %d" % (thin, len(th)))

# the retreating glacier: gone early in both interglacials, back in the late Holocene, gone at the end
sr = gb.simulate({"ice": "retreat"}); ir = sr["ice"]
gone = lambda a0, a1: all(v[1] == 0 for v in ir if 200 - a0 <= v[0] <= 200 - a1)
stade = max(v[1] for v in ir if 200 - 125 <= v[0] <= 200 - 122)
check("retreat: no ice 127.5-125.5 ka, a small stade 125-122 ka, none 121.5-118.5 and 9-5 ka, ice in the late Holocene, none at the end",
      gone(127.5, 125.5) and 0.15 < stade < 0.45 and gone(121.5, 118.5) and gone(9, 5) and max(v[1] for v in ir if v[0] >= 199) > 0.3 and ir[-1][1] == 0,
      "stade max %.2f, ice-free %.2f kyr" % (stade, 200 - ice_time(ir)))
check("retreat: no fracture where the ice is gone (buttressing unchanged)",
      all(erosion_on(ir, f["t"] - 1e-9) or erosion_on(ir, f["t"] + 1e-9) for f in sr["fracs"]))
check("retreat: without ice the threshold is at most T0 (no buttressing)",
      all(T <= 1 + 1e-12 for t, T in sr["thr"] if not erosion_on(ir, t) and not erosion_on(ir, t - 0.25)))
check("thick history unchanged by the retreat option (fingerprint)", gb.fingerprint() == "906 997 70 21 1153.612080056 1396.570513716 546.738034185")

# independent time-stepped integration of the same rules (same random streams, same order)
def stepped(params, dt):
    p = dict(gb.DEFAULTS); p.update(params)
    ice = gb.make_ice(p); nd = gb.make_sizes(p); nq = gb.make_quotas(p)
    E = p["E"] * p["Eref"]; t = 0.0; bp = [[v[0], p["B"] * gb.butt(v[1], p)] for v in ice]; R = min(p["R0"], 1 + bp[0][1]); D = 0.0; acc = 0.0; q = nq(); out = []; bo = bp[0][1]
    n = int(round(p["tmax"] / dt))
    for k in range(1, n + 1):
        t = p["tmax"] * k / n; bn = ice_at(bp, t); acc += abs(bn - bo); bo = bn; R += (E if erosion_on(ice, t - 0.5 * dt) else 0.0) * dt
        if acc >= q: D += p["W"] * q; acc -= q; q = nq()
        T = 1 - D + ice_at(bp, t)
        if R >= T:
            d = nd(); Rn = max(T - d, p["floor"])
            if Rn < R: out.append(t); R = Rn
            D = 0.0
    return out, R
for c in [{}, {"W": 0}, {"E": 2}, {"ice": "retreat"}]:
    ex = gb.simulate(c); tf = [f["t"] for f in ex["fails"]]; errs = []
    for dt in [0.01, 0.0025, 0.000625, 0.00015625]:
        ts, Rend = stepped(c, dt)
        errs.append((len(ts) == len(tf), max([abs(a - b) for a, b in zip(ts, tf)] + [0]), abs(Rend - ex["path"][-1][1])))
    # first order, in steps (each error is set by where single events fall on the step grid):
    # failure-time errors fall at every refinement, R_end errors never rise, and both at an
    # order ≥ 0.5 over the 64-fold refinement (≥ 8×; the four settings give 0.6-1.0)
    ok = all(e[0] for e in errs) and all(b[1] < a[1] for a, b in zip(errs, errs[1:])) and \
        all(b[2] <= a[2] * (1 + 1e-9) for a, b in zip(errs, errs[1:])) and \
        errs[-1][1] <= errs[0][1] / 8 and errs[-1][2] <= errs[0][2] / 8
    check("time-stepped reference converges to the exact solution %s" % c, ok,
          "same failure count; worst time error " + " → ".join("%.2g" % e[1] for e in errs) + " kyr, R_end error " +
          " → ".join("%.1e" % e[2] for e in errs) + " (dt 0.01 … 0.00015625); order %.2f, %.2f" %
          (math.log(errs[0][1] / errs[-1][1]) / math.log(64), math.log(errs[0][2] / errs[-1][2]) / math.log(64)))

# ---- the cross-section (glacial trough) -------------------------------------------------
wm = wr = 0.0; ok_mono = ok_ridge = ok_above = ok_ice = ok_lens = ok_trough = True; nls = 0
for c in SETS:
    s = gb.simulate(c); states = gb.section_states(s); C = gb.sec_C(s["p"])
    for k in range(1, len(states)):
        pre, rel = states[k]["pre"], states[k]["rel"]; f = [f for f in s["fails"] if f["t"] == states[k]["t"] and f["frm"] - f["to"] > gb.EPS][0]
        nls += 1
        wm = max(wm, abs(gb.sec_area(pre) - gb.sec_area(rel)))                        # cut = fill
        wr = max(wr, abs(pre[-1][1] - pre[0][1] - f["frm"]), abs(rel[-1][1] - rel[0][1] - f["to"]))
        ok_ridge = ok_ridge and rel[-1][1] <= pre[-1][1] + 1e-12 and rel[0][1] >= pre[0][1] - 1e-12
        # outside the cut/fill line the old surface is kept: every vertex is on it, or is the new floor or ridge
        ok_above = ok_above and all(abs(gb.sec_value(pre, v[0]) - v[1]) < 1e-12 or v[0] in (0.0, 1.0) for v in rel)
    for t in [0.0, 3.7, 29.9, 64.2, 101.9, 123.4, 150.05, 171.3, 188.8, 199.9, 200.0]:
        sec = gb.section_at(s, states, t); p = sec["pts"]
        wr = max(wr, abs(p[-1][1] - p[0][1] - sec["R"]), abs(gb.sec_area(p)))         # relief = R(t); datum = mean
        ok_mono = ok_mono and all(b[0] > a[0] and b[1] >= a[1] - 1e-12 for a, b in zip(p, p[1:])) and p[0][0] == 0 and p[-1][0] == 1
        # carving only removes rock: the surface is never above the one the last landslide left
        kk = max(j for j in range(len(states)) if states[j]["t"] <= t); cur = gb.sec_surface(states[kk]["rel"], sec["R"], C)
        ok_trough = ok_trough and all(v[1] <= gb.sec_value(states[kk]["rel"], v[0]) + 1e-12 for v in cur)
        ok_ice = ok_ice and abs(sec["zIce"] - p[0][1] - gb.SEC_ICE * C * sec["H"]) < 1e-12 and (sec["H"] > 0 or sec["zIce"] == p[0][1])
        ok_lens = ok_lens and 0 <= sec["weak"] < 1 and (c.get("W", 1) != 0 or sec["weak"] < 1e-12)
check("section: every landslide conserves area (cut = fill)", wm < 1e-12, "%d landslides, worst %.1e" % (nls, wm))
check("section: floor-to-ridge height = R(t) before and after every landslide and at sampled times; mean elevation is the datum", wr < 1e-9, "worst %.1e" % wr)
check("section: surface rises from the valley centre to the ridge", ok_mono)
check("section: carving only removes rock (never above the surface the last landslide left)", ok_trough)
check("section: a landslide never raises the ridge nor lowers the floor", ok_ridge)
check("section: outside the cut/fill line the old surface is unchanged", ok_above)
check("section: ice fills the trough to 0.6·C·H above the floor; none when the glacier is gone", ok_ice)
check("section: damage lens 0 ≤ D/T_intact < 1, and 0 with no weakening", ok_lens)
print("section fingerprint (GS.glacier.sectionFingerprint() in the browser):")
for c in [{}, {"ice": "retreat"}, {"ice": "retreat", "W": 0, "E": 0.25}]:
    print("  ", c, gb.section_fingerprint(c))

print("fingerprint (GS.glacier.fingerprint() in the browser):")
for c in [{}, {"W": 0}, {"E": 2, "W": 0.6}, {"ice": "retreat"}, {"ice": "retreat", "W": 0, "E": 0.25}]:
    print("  ", c, gb.fingerprint(c))
print("\n%d failed" % len(fails) if fails else "\nall passed")
sys.exit(1 if fails else 0)
