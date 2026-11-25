"""
test_relief_limit.py — the "Limits to relief" event model (relief_limit/model.js mirror).
Run: /opt/anaconda3/bin/python3 py/test_relief_limit.py
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(__file__))
import relief_limit as rl

fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  [" + detail + "]" if detail else ""))
    if not ok: fails.append(name)

def drops(s): return sum(f["from"] - f["to"] for f in s["fails"])

# ---- 0. PRNG: first values of mulberry32(7) are fixed (the JS page prints the same) ----
r = rl.mulberry32(7); v = [r() for _ in range(3)]
check("mulberry32 in [0,1) and deterministic", all(0 <= x < 1 for x in v) and v == [rl.mulberry32(7)() for _ in range(1)] + v[1:],
      "first three: " + ", ".join("%.10f" % x for x in v))

# ---- 1. invariants that hold in every stage ----
for st in range(1, 8):
    s = rl.simulate({"stage": st}); p = s["p"]; path = s["path"]
    tmono = all(path[i + 1][0] >= path[i][0] for i in range(len(path) - 1))
    # every segment is growth at exactly g, a hold at the ceiling (stage 1), or a vertical drop
    ok = True
    for (t0, r0), (t1, r1) in zip(path, path[1:]):
        if t1 > t0 + 1e-12:
            sl = (r1 - r0) / (t1 - t0)
            ok = ok and (abs(sl - p["g"]) < 1e-6 or (st == 1 and abs(sl) < 1e-9))
        else:
            ok = ok and r1 <= r0 + 1e-12
    check("stage %d: time monotone, segments are growth at g or drops" % st, tmono and ok, "%d vertices, %d failures" % (len(path), len(s["fails"])))
    # budget: relief gained = g t - sum of landslide drops (stage 1 holds instead)
    if st > 1:
        bud = p["R0"] + p["g"] * p["tmax"] - drops(s) - path[-1][1]
        check("stage %d: relief budget R0 + g t - drops = R_end" % st, abs(bud) < 1e-9, "residual %.2e" % bud)
    # never above the ceiling in force at the time
    lim = s["lim"]
    def ceil_at(t):
        c = lim[0][1]
        for (tt, cc) in lim:
            if tt <= t: c = cc
        return c
    check("stage %d: relief never exceeds the ceiling" % st, all(r_ <= ceil_at(t_) + 1e-9 for t_, r_ in path))
    # every earthquake leaves relief at or below its depressed limit
    okq = all(rl.relief_at(path, q["t"]) <= q["y"] + 1e-9 for q in s["quakes"])
    check("stage %d: after each earthquake R <= Ls - D" % st, okq, "%d earthquakes, %d triggered" % (len(s["quakes"]), sum(q["trig"] for q in s["quakes"])))

# ---- 2. stage 1 flatlines; stages 2-3 are exact saw-teeth ----
s = rl.simulate({"stage": 1}); p = s["p"]
check("stage 1: reaches L at (L-R0)/g and holds", abs(s["path"][1][0] - 0.7) < 1e-12 and s["path"][-1] == (p["tmax"], p["L"]) and not s["fails"])
for st, key in ((2, "dSmall"), (3, "dLarge")):
    s = rl.simulate({"stage": st}); p = s["p"]; d = p[key]
    ts = [f["t"] for f in s["fails"]]
    per = [b - a for a, b in zip(ts, ts[1:])]
    n_exp = sum(1 for j in range(10000) if 0.7 + j * d / p["g"] < p["tmax"] - rl.EPS)   # one landing on the end of the record is not counted
    check("stage %d: period d/g, drop d, count" % st, max(abs(x - d / p["g"]) for x in per) < 1e-12 and
          all(abs(f["from"] - f["to"] - d) < 1e-12 for f in s["fails"]) and len(ts) == n_exp, "%d teeth of %.2f" % (len(ts), d))

# ---- 3. stages 3, 4, 5 are the same relief history when T = d/g ----
a, b, c = (rl.simulate({"stage": k})["path"] for k in (3, 4, 5))
same = len(a) == len(b) == len(c) and max(max(abs(x[0] - y[0]), abs(x[1] - y[1])) for x, y in zip(a, b)) < 1e-9 and b == c
check("stages 3 = 4 = 5 (the coincident earthquakes reproduce the saw-tooth)", same)
s4 = rl.simulate({"stage": 4})
check("stage 4: every earthquake triggers, at the functional limit", all(q["trig"] and abs(q["y"] - 1.0) < 1e-12 for q in s4["quakes"]), "%d earthquakes" % len(s4["quakes"]))

# ---- 4. stipulated coincidence in stages 4-5: every earthquake strikes exactly as relief reaches L ----
okc = True; nq = 0
for extra in ({}, {"g": 0.37}, {"g": 0.62, "dLarge": 0.41, "R0": 0.85}, {"dLarge": 0.07, "tmax": 9.0}, {"dLarge": 0.6, "floor": 0.7, "g": 0.5}, {"R0": 1.0}):
    for st in (4, 5):
        s = rl.simulate(dict(extra, stage=st)); p = s["p"]; qk = s["quakes"]; nq += len(qk)
        okc = okc and len(qk) >= 2 and all(q["trig"] for q in qk) and len(s["fails"]) == len(qk)
        okc = okc and all(abs(f["from"] - p["L"]) < 1e-9 and abs(f["t"] - q["t"]) < 1e-15 for f, q in zip(s["fails"], qk))
        okc = okc and max(r for _, r in s["path"]) <= p["L"] + 1e-9
        s3 = rl.simulate(dict(extra, stage=3))["path"]      # the same relief history as the unexplained limit of stage 3
        okc = okc and all(abs(rl.relief_at(s["path"], tt) - rl.relief_at(s3, tt)) < 1e-9 for tt in [p["tmax"] * (j + 0.37) / 211 for j in range(211)])
check("stages 4-5: each earthquake finds relief exactly at L and fails it, for any uplift rate, landslide size, start and floor", okc, "%d earthquakes over 12 runs" % nq)
check("stages 4-5 with g = 0: relief never reaches L, so no earthquake is drawn", rl.simulate({"stage": 5, "g": 0.0})["quakes"] == [])

# ---- 5. static failure: with no earthquake large enough, relief saw-tooths at Ls ----
s = rl.simulate({"stage": 6, "fmax": 1e-6, "tmax": 8})
check("stage 6, negligible earthquakes: failures at the static limit", len(s["fails"]) > 5 and all(f["cause"] == "limit" and abs(f["from"] - 1.6) < 1e-12 for f in s["fails"]))

# ---- 6. stage 7 damage bookkeeping ----
s = rl.simulate({"stage": 7}); p = s["p"]
lim = s["lim"]
# replay the bookkeeping from the event lists alone: damage accumulates over surviving
# earthquakes and is cleared by every landslide
ev = sorted([(q["t"], 0, q) for q in s["quakes"]] + [(f["t"], 1, f) for f in s["fails"]], key=lambda e: (e[0], e[1]))
acc = 0.0; ok = True
for t_, kind, e in ev:
    if kind == 0:
        ok = ok and abs(e["top"] - (p["Ls"] - acc)) < 1e-12        # the limit each earthquake pulls down from
        if not e["trig"]: acc += p["damage"] * e["D"]
    else:
        acc = 0.0
ups = [b[1] - a[1] for a, b in zip(lim, lim[1:]) if b[1] > a[1] + 1e-15]
check("stage 7: Ls = Ls0 - k * sum(D of surviving quakes since the last landslide)",
      ok and abs(lim[-1][1] - (p["Ls"] - acc)) < 1e-12 and min(c for _, c in lim) < p["Ls"] - 0.05, "lowest %.3f, final %.3f" % (min(c for _, c in lim), lim[-1][1]))
check("stage 7: every landslide restores Ls to its original value",
      len(ups) > 0 and all(abs(b[1] - p["Ls"]) < 1e-15 for a, b in zip(lim, lim[1:]) if b[1] > a[1] + 1e-15) and
      all(abs(rl.simulate({"stage": 7})["lim"][i][1] - c) < 1e-15 for i, (_, c) in enumerate(lim)) and
      all(any(abs(f["t"] - a[0]) < 1e-12 for f in s["fails"]) for a, b in zip(lim, lim[1:]) if b[1] > a[1] + 1e-15), "%d restorations, %d landslides" % (len(ups), len(s["fails"])))
s = rl.simulate({"stage": 7, "fmax": 0.03, "damage": 1.0, "rate": 40, "tmax": 8})
check("stage 7: a failure at the damaged static limit also restores it", any(f["cause"] == "limit" and f["from"] < 1.6 - 1e-6 for f in s["fails"]) and
      all(abs(b[1] - 1.6) < 1e-15 for a, b in zip(s["lim"], s["lim"][1:]) if b[1] > a[1] + 1e-15), "%d static failures" % sum(f["cause"] == "limit" for f in s["fails"]))
s = rl.simulate({"stage": 7}); p = s["p"]
s6 = rl.simulate({"stage": 6})
check("stages 6 and 7 share the earthquake sequence (same seed)", [(q["t"], q["f"]) for q in s6["quakes"]] == [(q["t"], q["f"]) for q in s["quakes"]])
check("every depression is the fraction f of the limit as it stands (D = f * top)",
      all(abs(q["D"] - q["f"] * q["top"]) < 1e-15 and abs(q["y"] - q["top"] * (1 - q["f"])) < 1e-12 for st in (4, 5, 6, 7) for q in rl.simulate({"stage": st})["quakes"]) and
      any(q["top"] < 1.6 - 1e-3 for q in s["quakes"]))
check("stage 7 with k = 0 is stage 6", rl.simulate({"stage": 7, "damage": 0})["path"] == s6["path"])

# ---- 7. Poisson / exponential statistics of the sequence ----
s = rl.simulate({"stage": 6, "tmax": 2000}); q = s["quakes"]; p = s["p"]
n = len(q); b = p["b"]; dM = p["dM"]
check("earthquake rate (long run)", abs(n / 2000 - p["rate"]) < 0.5, "rate %.2f" % (n / 2000))
# doubly truncated Gutenberg-Richter: P(M > Mmax - x) = (10^(-b (dM - x)) - 10^(-b dM)) / (1 - 10^(-b dM))
okg = True; det = []
for x_ in (1.5, 1.0, 0.5):
    exp_ = (10 ** (-b * (dM - x_)) - 10 ** (-b * dM)) / (1 - 10 ** (-b * dM)); got = sum(1 for e in q if e["m"] < x_) / n
    okg = okg and abs(got - exp_) < 4 * math.sqrt(exp_ * (1 - exp_) / n) + 1e-4; det.append("%.4f/%.4f" % (got, exp_))
check("magnitudes follow the truncated Gutenberg-Richter law (b = %g)" % b, okg, "P(within 1.5, 1, 0.5 of Mmax) sampled/exact: " + ", ".join(det))
check("f = fmax 10^(-c (Mmax - M)), within [fmax 10^(-c dM), fmax]",
      all(abs(e["f"] - p["fmax"] * 10 ** (-p["c"] * e["m"])) < 1e-15 and p["fmax"] * 10 ** (-p["c"] * dM) - 1e-12 <= e["f"] <= p["fmax"] for e in q))
# the implied size distribution is a power law of exponent b/c: one decade in f is c^-1 magnitude units
big = sum(1 for e in q if e["f"] > 0.3); small = sum(1 for e in q if e["f"] > 0.1)
check("N(>f) ~ f^(-b/c): ratio N(>0.1)/N(>0.3) near 3^(b/c) = %.1f" % (3 ** (b / p["c"])), abs(small / big / ((0.1 ** -2 - 0.6 ** -2) / (0.3 ** -2 - 0.6 ** -2)) - 1) < 0.15, "%d / %d" % (small, big))

# ---- 8. an independent time-stepped integration of the same rules agrees once dt resolves the events ----
#         (samples that land on an event time may sit on the other side of a drop, hence the mean; the
#         reference applies a ceiling failure before any earthquake in the same step, so a step holding
#         both can misorder them: 4e-3 is reported, the two finer steps are the test)
for st in (2, 3, 5, 6, 7):
    s = rl.simulate({"stage": st}); errs = []
    for dt in (4e-3, 1e-3, 2.5e-4):
        ref = rl.reference({"stage": st}, dt)
        e = sum(abs(r_ - rl.relief_at(s["path"], t_)) for t_, r_ in ref) / len(ref)
        errs.append(e)
    check("stage %d: time-stepped reference agrees (mean |dR| < 1e-3 at dt 1e-3 and 2.5e-4)" % st, max(errs[1:]) < 1e-3, "dt 4e-3/1e-3/2.5e-4 -> " + " / ".join("%.1e" % e for e in errs))

# ---- 8b. uplift rate: the same invariants away from g = 1 ----
for gg in (0.2, 0.6, 1.7):
    s = rl.simulate({"stage": 3, "g": gg, "tmax": 12.0}); ts = [f["t"] for f in s["fails"]]
    okp = max(abs((b_ - a_) - 0.25 / gg) for a_, b_ in zip(ts, ts[1:])) < 1e-12
    s = rl.simulate({"stage": 7, "g": gg}); p = s["p"]
    bud = p["R0"] + gg * p["tmax"] - drops(s) - s["path"][-1][1]
    okq = all(rl.relief_at(s["path"], q["t"]) <= q["y"] + 1e-9 for q in s["quakes"])
    okg = all(abs((r1 - r0) / (t1 - t0) - gg) < 1e-6 for (t0, r0), (t1, r1) in zip(s["path"], s["path"][1:]) if t1 > t0 + 1e-9)
    ref = rl.reference({"stage": 7, "g": gg}, 2.5e-4)
    e = sum(abs(r_ - rl.relief_at(s["path"], t_)) for t_, r_ in ref) / len(ref)
    check("g = %.1f: period d/g (stage 3); stage 7 budget, growth at g, R <= Ls - D, reference" % gg, okp and abs(bud) < 1e-9 and okq and okg and e < 1e-3,
          "%d landslides, mean relief %.2f, reference %.1e" % (len(s["fails"]), rl.mean_relief(s["path"]), e))
check("the earthquake sequence does not depend on g (stage 7)", [(q["t"], q["f"]) for q in rl.simulate({"stage": 7, "g": 2.5})["quakes"]] == [(q["t"], q["f"]) for q in rl.simulate({"stage": 7})["quakes"]])

# ---- 8c. zero uplift: relief never grows, the run terminates, failures need an earthquake ----
for st in range(1, 8):
    s = rl.simulate({"stage": st, "g": 0.0}); pth = s["path"]
    okz = all(r1 <= r0 + 1e-15 for (_, r0), (_, r1) in zip(pth, pth[1:])) and pth[-1][0] == s["p"]["tmax"] and all(f["cause"] == "quake" for f in s["fails"])
    check("g = 0, stage %d: relief non-increasing, no failure without an earthquake" % st, okz, "%d landslides, R %.2f -> %.2f" % (len(s["fails"]), pth[0][1], pth[-1][1]))
s = rl.simulate({"stage": 7, "g": 0.0, "R0": 1.0, "damage": 0.3})
check("g = 0, stage 7: damage alone brings the limit down to the relief (landslides, R0 = L)", len(s["fails"]) >= 1 and min(c for _, c in s["lim"]) < 1.6 - 0.2,
      "%d landslides, limit as low as %.2f" % (len(s["fails"]), min(c for _, c in s["lim"])))
check("g = 0 with R0 = L in stage 3: holds at the limit, no division by zero", rl.simulate({"stage": 3, "g": 0.0, "R0": 1.0})["path"] == [(0.0, 1.0), (4.0, 1.0)])

# ---- 8d. landslide sizes: constant in stages 2-5, truncated power law in 6-7 ----
for st in (3, 4, 5):
    check("stage %d: every landslide has d = dLarge" % st, all(f["d"] == 0.25 for f in rl.simulate({"stage": st})["fails"]))
for st in (6, 7):
    s = rl.simulate({"stage": st}); p = s["p"]
    check("stage %d: each failure leaves relief its own d below the limit in force, dMin <= d <= dMax" % st,
          all(abs(f["to"] - max(f["lim"] - f["d"], p["floor"])) < 1e-12 and p["dMin"] - 1e-12 <= f["d"] <= p["dMax"] + 1e-12 for f in s["fails"]) and
          len(set(round(f["d"], 9) for f in s["fails"])) == len(s["fails"]), "sizes " + " ".join("%.2f" % f["d"] for f in s["fails"]))
    check("stage %d with dDist = constant: every d = dLarge" % st, all(f["d"] == 0.25 for f in rl.simulate({"stage": st, "dDist": "constant"})["fails"]))
check("landslide sizes use their own stream: the earthquakes are unchanged by dAlpha", [(q["t"], q["f"]) for q in rl.simulate({"stage": 7, "dAlpha": 2.5})["quakes"]] == [(q["t"], q["f"]) for q in rl.simulate({"stage": 7})["quakes"]])
for al in (0.5, 1.0, 2.0):
    p = dict(rl.DEFAULTS); p.update({"stage": 7, "dAlpha": al}); nd = rl.make_sizes(p); xs = [nd() for _ in range(40000)]
    lo, hi = p["dMin"], p["dMax"]; okd = lo <= min(xs) and max(xs) <= hi; det = []
    for x_ in (0.1, 0.2, 0.4):
        exp_ = (x_ ** -al - hi ** -al) / (lo ** -al - hi ** -al); got = sum(1 for v in xs if v > x_) / len(xs)
        okd = okd and abs(got - exp_) < 4 * math.sqrt(exp_ * (1 - exp_) / len(xs)); det.append("%.3f/%.3f" % (got, exp_))
    check("truncated power law, alpha = %g: P(d > 0.1, 0.2, 0.4) sampled/exact" % al, okd, ", ".join(det))
p = dict(rl.DEFAULTS); p.update({"stage": 7, "dMin": 0.3, "dMax": 0.3}); nd = rl.make_sizes(p)
check("dMin = dMax degenerates to a constant size", all(nd() == 0.3 for _ in range(5)))

# ---- 8e. the floor: no landslide cuts relief below it ----
okf = True; low = 9.0
for st in range(2, 8):
    for extra in ({}, {"g": 0.0, "R0": 1.0, "damage": 0.4}, {"g": 0.3, "dMax": 1.2, "dAlpha": 0.2, "fmax": 0.95}, {"floor": 0.5, "dLarge": 0.6}):
        s = rl.simulate(dict(extra, stage=st)); fl = s["p"]["floor"]
        okf = okf and all(f["to"] >= fl - 1e-15 and f["to"] < f["from"] for f in s["fails"]) and min(r for _, r in s["path"]) >= min(fl, s["p"]["R0"]) - 1e-15
        okf = okf and all(q["trig"] == any(abs(f["t"] - q["t"]) < 1e-15 and f["cause"] == "quake" for f in s["fails"]) for q in s["quakes"])
        low = min(low, min(r for _, r in s["path"]))
check("no landslide leaves relief below the floor; an earthquake is marked as failing the slope only if a landslide happened", okf, "lowest relief reached %.3f" % low)
s = rl.simulate({"stage": 7, "g": 0.0, "R0": 1.0, "damage": 0.4})
check("g = 0, heavy weakening: relief ends on the floor, not at zero; the limit stops sagging once it reaches the relief",
      abs(s["path"][-1][1] - 0.176) < 1e-12 and s["lim"][-1][1] > 0.176 - 1e-9, "R_end %.3f, limit_end %.3f" % (s["path"][-1][1], s["lim"][-1][1]))
check("floor = 0 reproduces the unfloored rule", all(abs(f["to"] - max(f["lim"] - f["d"], 0.0)) < 1e-12 for f in rl.simulate({"stage": 7, "floor": 0.0, "g": 0.2, "dMax": 1.2})["fails"]))

# ---- 9. time-mean relief: a long saw-tooth averages L - d/2 ----
s = rl.simulate({"stage": 3, "R0": 1.0, "tmax": 250.0})
check("mean relief of the saw-tooth = L - d/2", abs(rl.mean_relief(s["path"]) - (1 - 0.125)) < 1e-9, "%.6f" % rl.mean_relief(s["path"]))

# ---- 10. fingerprints for the JS port (compare in the browser console: GS.relief.fingerprint(stage)) ----
print("\nfingerprints (stage: vertices, failures, quakes, triggered, sum R over vertices, final ceiling, mean relief)")
for st in range(1, 8):
    s = rl.simulate({"stage": st})
    print("  %d: %d %d %d %d %.9f %.9f %.9f" % (st, len(s["path"]), len(s["fails"]), len(s["quakes"]), sum(q["trig"] for q in s["quakes"]),
                                             sum(r_ for _, r_ in s["path"]), s["lim"][-1][1], rl.mean_relief(s["path"])))

s = rl.simulate({"stage": 7, "g": 0.6, "damage": 0.3})
print("  7 with g = 0.6, damage = 0.3 (GS.relief.fingerprint(7, {g: 0.6, damage: 0.3})):\n     %d %d %d %d %.9f %.9f %.9f" % (
    len(s["path"]), len(s["fails"]), len(s["quakes"]), sum(q["trig"] for q in s["quakes"]), sum(r_ for _, r_ in s["path"]), s["lim"][-1][1], rl.mean_relief(s["path"])))

# ---- 11. the cross-section: geometry read from the same history ----
def slopes(poly): return [(z1 - z0) / (x1 - x0) for (x0, z0), (x1, z1) in zip(poly, poly[1:])]
stress = [{}, {"g": 0.3}, {"floor": 0.0, "dMax": 1.2}, {"damage": 0.6, "rate": 60}, {"Ls": 1.1}, {"Ls": 2.5, "fmax": 0.9}, {"R0": 0.0}, {"R0": 1.0, "g": 0.0}]
worst_mass = 0.0; worst_rel = 0.0; ok_shape = True; ok_cutfill = True; ok_ridge = True; nls = 0; nB = 0
for st in range(1, 8):
    for ex in (stress if st >= 6 else stress[:3] + stress[6:]):
        s = rl.simulate(dict(ex, stage=st)); p = s["p"]; S = p["L"] if st <= 3 else p["Ls"]
        states = rl.section_states(s)
        for a_, b_ in zip(states, states[1:]):
            pre, post = b_["pre"], b_["rel"]; nls += 1
            worst_mass = max(worst_mass, abs(rl.sec_area(pre) - rl.sec_area(post)))
            # the valley is filled, the ridge is never raised, and the new surface is the old one or a single straight cut/fill line
            ok_ridge = ok_ridge and post[0][1] > pre[0][1] and post[-1][1] <= pre[-1][1] + 1e-9
            if post[-1][1] > pre[-1][1] - 1e-9 and len(post) > 2: nB += 1
            x2 = post[1][0]
            if post[-1][1] < pre[-1][1] - 1e-9:           # the line runs valley to ridge: nothing of the old surface is left
                ok_cutfill = ok_cutfill and len(post) == 2
                continue
            ok_cutfill = ok_cutfill and all(abs(rl.sec_value(pre, x) - z) < 1e-9 for x, z in post[2:]) and abs(rl.sec_value(pre, x2) - post[1][1]) < 1e-9
        ts = [p["tmax"] * i / 400.0 for i in range(401)] + [f["t"] for f in s["fails"]]
        for t in ts:
            sec = rl.section_at(s, states, t); pts = sec["pts"]; sl = slopes(pts)
            worst_rel = max(worst_rel, abs(pts[-1][1] - pts[0][1] - rl.relief_at(s["path"], t)), abs(rl.sec_area(pts)))
            ok_shape = ok_shape and pts[0][0] == 0.0 and pts[-1][0] == 1.0 and all(x1 > x0 for (x0, _), (x1, _) in zip(pts, pts[1:])) \
                and all(-1e-9 <= v <= S + 1e-9 for v in sl) and all(b2 <= a2 + 1e-7 for a2, b2 in zip(sl, sl[1:])) and 0.0 <= sec["weak"] < 1.0
check("section: every landslide conserves area (cut = fill)", worst_mass < 1e-12 and nls > 300, "%d landslides, worst residual %.1e" % (nls, worst_mass))
check("section: ridge-to-valley height equals R(t), and the mean elevation is the datum, at every sampled time in every step", worst_rel < 1e-9, "worst %.1e" % worst_rel)
check("section: surface is a concave polyline rising from valley to ridge, nowhere steeper than the limiting gradient", ok_shape)
check("section: a landslide raises the valley floor and never the ridge", ok_ridge, "%d of %d leave the ridge untouched" % (nB, nls))
check("section: above the cut/fill line the old surface is unchanged", ok_cutfill)

# step 1: a V at the limiting gradient cut into a flat surface; the flat is consumed when R = L
s = rl.simulate({"stage": 1}); states = rl.section_states(s)
okv = True
for t in (0.0, 0.2, 0.5, 0.69):
    pts = rl.section_at(s, states, t)["pts"]; R = 0.3 + t
    okv = okv and len(pts) == 3 and abs(pts[1][0] - R) < 1e-12 and abs(pts[2][1] - pts[1][1]) < 1e-12 and abs(pts[2][1] - R * R / 2) < 1e-12
pts = rl.section_at(s, states, 2.0)["pts"]
check("section, step 1: V of half-width w R/L in a flat surface (which stands R^2/2L above the mean); at the limit a straight slope from -L/2 to +L/2",
      okv and len(pts) == 2 and abs(pts[0][1] + 0.5) < 1e-12 and abs(pts[1][1] - 0.5) < 1e-12)

# steps 2-3: the landslide pivots the slope about its midpoint; recovery steepens it from the valley up
for st, d in ((2, 0.05), (3, 0.25)):
    s = rl.simulate({"stage": st}); states = rl.section_states(s)
    okp = all(len(b_["rel"]) == 2 and abs(b_["rel"][0][1] - b_["pre"][0][1] - d / 2) < 1e-12 and abs(b_["pre"][1][1] - b_["rel"][1][1] - d / 2) < 1e-12 for b_ in states[1:])
    t0 = s["fails"][0]["t"]; okk = True
    for fr in (0.25, 0.5, 0.9):
        sec = rl.section_at(s, states, t0 + fr * d); pts = sec["pts"]; R = 1 - d + fr * d
        okk = okk and len(pts) == 3 and abs(pts[1][0] - (R - (1 - d)) / d) < 1e-9 and abs(slopes(pts)[0] - 1) < 1e-9 and abs(slopes(pts)[1] - (1 - d)) < 1e-9
    check("section, step %d: ridge down d/2, valley up d/2, slope L-d over w; then the limiting gradient climbs to x = w (R - R1)/(L - R1)" % st, okp and okk and len(states) - 1 == len(s["fails"]))

# steps 3 and 4 share R(t) but not the section: in 4 the limiting gradient is the static one, never reached
s3 = rl.simulate({"stage": 3}); s4 = rl.simulate({"stage": 4}); t = s3["fails"][2]["t"] - 1e-6
a3 = rl.section_at(s3, rl.section_states(s3), t)["pts"]; a4 = rl.section_at(s4, rl.section_states(s4), t)["pts"]
check("section, step 4: just before a failure the slope is straight in step 3 and still two-part in step 4 (lower part at Ls/w)", len(a3) <= 3 and abs(slopes(a3)[0] - 1) < 1e-6 and len(a4) >= 3 and abs(slopes(a4)[0] - 1.6) < 1e-9)

# step 7: weak follows the damage and is cleared by every landslide
s = rl.simulate({"stage": 7}); states = rl.section_states(s)
wk = [rl.section_at(s, states, f["t"])["weak"] for f in s["fails"]]
wmax = max(rl.section_at(s, states, 4.0 * i / 2000)["weak"] for i in range(2001))
check("section, step 7: weak = (Ls0 - Ls)/Ls0, zero after every landslide; zero throughout step 6", max(wk) == 0.0 and abs(wmax - (1.6 - min(c for _, c in s["lim"])) / 1.6) < 1e-3 and
      all(rl.section_at(s6, rl.section_states(s6), 4.0 * i / 50)["weak"] == 0.0 for i in range(51)), "largest %.3f" % wmax)

print("\nsection fingerprints (stage: landslide surfaces, then sum of |z| and vertex count of the section at t = 0.37, 1.91, 3.3; GS.relief.sectionFingerprint(stage))")
for st in range(1, 8):
    print("  " + rl.section_fingerprint(st))

print("\n%d failed" % len(fails) if fails else "\nall passed")
sys.exit(1 if fails else 0)
