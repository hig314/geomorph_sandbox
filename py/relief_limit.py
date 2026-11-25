"""
relief_limit.py — mirror of relief_limit/model.js ("Limits to relief", R(t)).

Relief grows at a constant rate g until a limit fails it. Seven stages build the idea up:
  1  constant limit L, relief flatlines at it
  2  each arrival at L is a landslide: relief drops by a small increment d and regrows
  3  the same with a large d
  4  a static limit Ls stands above; L is where constant earthquakes (depression
     D = Ls - L) put the limit at the moment they strike. Their timing is stipulated to
     coincide exactly with relief reaching L, so the recurrence is d/g, not a parameter
  5  the same model, with each earthquake drawn (a drop-line from Ls to Ls - D)
  6  earthquakes are a Poisson process with Gutenberg-Richter magnitudes; each depresses
     the limit by a fraction f(M) of the limit as it stands; one fails the slope only if
     Ls (1 - f) reaches the relief
  7  an earthquake that does not fail the slope damages it: Ls drops by k D; the next
     landslide removes the damaged rock and restores Ls

One rule covers every stage. With ceiling C (L in stages 1-3, Ls after):
  between events   dR/dt = g, and R = C is a failure
  earthquake       D = f C (a fraction of the current limit); Leff = C - D; failure if R >= Leff
  failure          R <- max(Leff - d, 0)   (Leff = C for a failure at the ceiling)
  landslide size   d constant in stages 2-5; in stages 6-7 each failure draws its own d from
                   a truncated power law, N(>d) ~ d^(-alpha) on [dMin, dMax] (log-log linear
                   in size-frequency space), from a stream separate from the earthquakes
  floor            no landslide takes relief below the floor: R <- max(Leff - d, floor); a slope
                   already at or below it does not fail (and an earthquake that would have
                   failed it leaves no damage either)
  g = 0            relief does not grow; only earthquakes (and damage) can fail the slope
  d = 0            no landslide: relief holds at the ceiling (stage 1)
  damage (7)       earthquake without failure: C <- C - k D;  any failure: C <- Ls

Earthquake sizes (stages 6-7): magnitudes from the doubly truncated Gutenberg-Richter law
on [Mmax - dM, Mmax], N(>M) ~ 10^(-b M); shaking amplitude at a fixed distance scales as
10^(c M), and the depression fraction is taken proportional to it:
  f = fmax 10^(-c (Mmax - M))
so f is a truncated power law, N(>f) ~ f^(-b/c). Distance to the source is not sampled.
Stages 4-5 use the constant fraction f = (Ls - L) / Ls.

The solution is piecewise linear, so it is advanced event to event and is exact: there is
no time step. Dimensionless: relief in units of L, time in units of L/g.
"""
import math

EPS = 1e-9

DEFAULTS = {
    "stage": 3,
    "g": 1.0,        # growth rate of relief (uplift + incision), L per time unit
    "L": 1.0,        # functional limit
    "R0": 0.3,       # initial relief
    "tmax": 4.0,     # length of the record
    "dSmall": 0.05,  # landslide increment, stage 2
    "dLarge": 0.25,  # landslide increment, stages 3-5 (and 6-7 when dDist = "constant")
    "dDist": "power",  # landslide sizes in stages 6-7: "power" (truncated power law) or "constant"
    "dMin": 0.05,    # smallest landslide increment, stages 6-7
    "dMax": 0.6,     # largest landslide increment, stages 6-7
    "dAlpha": 1.0,   # exponent of the cumulative size distribution, N(>d) ~ d^(-alpha)
    "floor": 0.176,  # relief no landslide can cut below (the page sets it to 10 % of the graph height)
    "Ls": 1.6,       # static limit
    "rate": 30.0,    # earthquakes per time unit, stages 6-7
    "b": 1.0,        # Gutenberg-Richter b-value, stages 6-7
    "dM": 2.0,       # magnitude range drawn, Mmax - Mmin, stages 6-7
    "fmax": 0.6,     # fraction of the current limit removed by the largest earthquake (Mmax)
    "c": 0.5,        # shaking amplitude ~ 10^(c M) at a fixed distance
    "damage": 0.15,  # k: drop of Ls per unit depression, until the next landslide, stage 7
    "seed": 82,
}


def mulberry32(a):
    """Seeded PRNG in [0, 1); bit-for-bit the function in core/noise.js."""
    M = 0xFFFFFFFF
    state = [a & M]

    def imul(x, y):
        return ((x & M) * (y & M)) & M

    def rnd():
        state[0] = (state[0] + 0x6D2B79F5) & M
        a_ = state[0]
        t = imul(a_ ^ (a_ >> 15), 1 | a_)
        t = ((t + imul(t ^ (t >> 7), 61 | t)) & M) ^ t
        return ((t ^ (t >> 14)) & M) / 4294967296.0
    return rnd


def make_quakes(p):
    """[(t, f, m)] in time order, t <= tmax: f = fraction of the current limit removed,
    m = Mmax - M (magnitude units below the largest earthquake; 0 in stages 4-5)."""
    st = p["stage"]; out = []
    if st <= 3:
        return out
    if st <= 5:
        # stipulated coincidence: an earthquake strikes at each instant relief reaches L, so
        # the recurrence is the regrowth time of one landslide, (L - max(L - d, floor)) / g
        f = (p["Ls"] - p["L"]) / p["Ls"]
        if p["g"] <= 0 or p["R0"] > p["L"]:
            return out                              # relief never arrives at L: no coincidence to draw
        T = (p["L"] - max(p["L"] - p["dLarge"], max(p["floor"], 0.0))) / p["g"]
        t = (p["L"] - p["R0"]) / p["g"]; k = 0
        while t + k * T <= p["tmax"]:
            out.append((t + k * T, f, 0.0)); k += 1
            if T <= 0: break
        return out
    rnd = mulberry32(int(p["seed"])); t = 0.0
    beta = p["b"] * math.log(10.0); span = 1.0 - math.exp(-beta * p["dM"])
    while True:
        t += -math.log(1.0 - rnd()) / p["rate"]
        m = p["dM"] + math.log(1.0 - rnd() * span) / beta      # Mmax - M, in [0, dM]
        if t > p["tmax"]:
            return out
        out.append((t, p["fmax"] * math.pow(10.0, -p["c"] * m), m))


def make_sizes(p):
    """Returns next_d(): the increment of the next landslide. Constant in stages 1-5 (0 in
    stage 1) or when dDist = "constant"; otherwise a truncated power law on [dMin, dMax],
    P(d > x) = (x^-a - hi^-a) / (lo^-a - hi^-a), drawn by inversion from its own stream."""
    st = int(p["stage"])
    const = 0.0 if st == 1 else (p["dSmall"] if st == 2 else p["dLarge"])
    lo = min(p["dMin"], p["dMax"]); hi = max(p["dMin"], p["dMax"]); a = p["dAlpha"]
    if st <= 5 or p["dDist"] != "power":
        return lambda: const
    if hi <= lo * (1 + 1e-9) or a <= 0:
        return lambda: lo
    rnd = mulberry32(int(p["seed"]) + 1000003); span = 1.0 - math.pow(lo / hi, a)
    return lambda: lo * math.pow(1.0 - rnd() * span, -1.0 / a)


def simulate(params=None):
    p = dict(DEFAULTS); p.update(params or {})
    st = int(p["stage"]); g = p["g"]; tmax = p["tmax"]
    next_d = make_sizes(p); hold = st == 1; floor = max(p["floor"], 0.0)
    C = p["L"] if st <= 3 else p["Ls"]
    k = p["damage"] if st == 7 else 0.0
    C0 = C
    qs = make_quakes(p)
    t = 0.0; R = min(p["R0"], C)
    path = [(t, R)]; lim = [(0.0, C)]; fails = []; quakes = []
    qi = 0
    while True:
        tq = qs[qi][0] if qi < len(qs) else float("inf")
        tn = min(tq, tmax)
        # grow to tn; every arrival at the ceiling on the way is a failure
        while True:
            th = t + (C - R) / g if g > 0 else float("inf")
            if th >= tn - EPS:
                break
            path.append((th, C)); t = th
            if hold:                        # no landslide: hold at the ceiling
                R = C; t = tn; break
            d = next_d()
            R = min(max(C - d, floor), C)
            fails.append({"t": t, "from": C, "to": R, "cause": "limit", "d": d, "lim": C})
            path.append((t, R))
            if C < C0:                      # the landslide removes the damaged rock
                lim.append((t, C)); lim.append((t, C0)); C = C0
        R = min(R + g * (tn - t), C); t = tn
        path.append((t, R))
        if tq > tmax:
            break
        fq = qs[qi][1]; mq = qs[qi][2]; qi += 1
        D = fq * C                          # relative to the limit as it stands
        Leff = max(C - D, 0.0)
        reach = R >= Leff - EPS             # the depressed limit reaches the relief
        q = {"t": t, "D": D, "f": fq, "m": mq, "top": C, "y": Leff, "trig": False}
        if reach:
            d = next_d()
            Rn = max(Leff - d, floor)
            if Rn < R:
                q["trig"] = True
                fails.append({"t": t, "from": R, "to": Rn, "cause": "quake", "d": d, "lim": Leff})
                R = Rn; path.append((t, R))
                if C < C0:                  # the landslide removes the damaged rock
                    lim.append((t, C)); lim.append((t, C0)); C = C0
        elif k > 0:
            C2 = C - k * D                  # > R, since k <= 1 and R < C - D
            lim.append((t, C)); lim.append((t, C2)); C = C2
        quakes.append(q)
    lim.append((tmax, C))
    return {"path": path, "lim": lim, "quakes": quakes, "fails": fails, "p": p}


def relief_at(path, t):
    """R(t) on the piecewise-linear path (after any drop at exactly t)."""
    lo, hi = 0, len(path) - 1
    while lo < hi:                           # last vertex with time <= t
        mid = (lo + hi + 1) // 2
        if path[mid][0] <= t: lo = mid
        else: hi = mid - 1
    if lo == len(path) - 1:
        return path[lo][1]
    (t0, r0), (t1, r1) = path[lo], path[lo + 1]
    return r0 if t1 == t0 else r0 + (r1 - r0) * (t - t0) / (t1 - t0)


def mean_relief(path):
    """Time-averaged relief over the record (trapezoids; drops have zero width)."""
    a = 0.0
    for (t0, r0), (t1, r1) in zip(path, path[1:]):
        a += 0.5 * (r0 + r1) * (t1 - t0)
    T = path[-1][0] - path[0][0]
    return a / T if T > 0 else path[0][1]


def reference(params, dt):
    """Independent time-stepped integration of the same rules (forward Euler, events
    applied at the first step on or after their time). Converges to simulate() as dt -> 0.
    Returns [(t, R)] sampled every step."""
    p = dict(DEFAULTS); p.update(params or {})
    st = int(p["stage"]); g = p["g"]
    next_d = make_sizes(p); hold = st == 1; floor = max(p["floor"], 0.0)
    C = p["L"] if st <= 3 else p["Ls"]
    k = p["damage"] if st == 7 else 0.0
    qs = make_quakes(p); qi = 0; C0 = C
    R = min(p["R0"], C); out = [(0.0, R)]
    n = int(round(p["tmax"] / dt))
    for i in range(1, n + 1):
        t = i * dt
        R += g * dt
        if R >= C and g > 0:
            R = C if hold else min(max(C - next_d(), floor), C) + (R - C)
            if not hold: C = C0
        while qi < len(qs) and qs[qi][0] <= t:
            D = qs[qi][1] * C; qi += 1
            Leff = max(C - D, 0.0)
            # relief at the quake itself, not at the end of the step
            Rq = R - g * (t - qs[qi - 1][0])
            if Rq >= Leff - EPS:
                Rn = max(Leff - next_d(), floor)
                if Rn < Rq:
                    R = Rn + g * (t - qs[qi - 1][0])
                    C = C0
            elif k > 0:
                C -= k * D
        out.append((t, R))
    return out


# ---- The cross-section ---------------------------------------------------------------
# A ridge-valley half-section z(x), valley at x = 0, ridge at x = 1 (x in units of the
# half-spacing w, z in units of L), read from the same relief history. Three rules:
#   incision    the valley floor cuts down and a slope at the limiting gradient S = C0
#               (the intact ceiling: relief C0 over the half-spacing) rises from it; above
#               where that slope meets it, the older (relict) surface is untouched:
#               z(x) = min(relict(x), z_valley + S x),   z_valley = z_ridge - R(t)
#   landslide   a failure taking relief from R to R' replaces the surface by a straight
#               line of equal area: cut above where it crosses the old surface, fill below.
#               If that line would stand above the old ridge, the ridge is left alone, the
#               valley floor rises by R - R', and the line climbs from it to the point x2
#               on the old surface at which cut = fill; above x2 the old surface remains
#   datum       drawn about the mean elevation of the section, which a landslide conserves
# The surface is a concave polyline throughout, so every step is exact except the root x2
# (bisection to round-off). Nothing here feeds back on R(t).

def sec_area(poly, x_to=1.0):
    """Integral of the polyline from 0 to x_to (trapezoids)."""
    a = 0.0
    for (x0, z0), (x1, z1) in zip(poly, poly[1:]):
        if x_to <= x0:
            break
        if x_to < x1:
            zc = z0 + (z1 - z0) * (x_to - x0) / (x1 - x0)
            return a + 0.5 * (z0 + zc) * (x_to - x0)
        a += 0.5 * (z0 + z1) * (x1 - x0)
    return a


def sec_value(poly, x):
    for (x0, z0), (x1, z1) in zip(poly, poly[1:]):
        if x <= x1:
            return z0 if x1 == x0 else z0 + (z1 - z0) * (x - x0) / (x1 - x0)
    return poly[-1][1]


def sec_surface(rel, R, S):
    """The surface for relief R over the relict polyline rel: min(rel, z_valley + S x)."""
    n = len(rel); zv = rel[n - 1][1] - R
    if rel[0][1] - zv <= EPS:
        return list(rel)                              # no incision since the relict surface formed
    for i in range(n - 1):
        (x0, z0), (x1, z1) = rel[i], rel[i + 1]
        g0 = z0 - zv - S * x0; g1 = z1 - zv - S * x1
        if g1 <= EPS or i == n - 2:                   # the limiting slope meets the relict surface here
            xk = x1 if g0 - g1 <= 0 else x0 + (x1 - x0) * g0 / (g0 - g1)
            xk = min(max(xk, x0), x1)
            if x1 - xk < 1e-9:
                return [(0.0, zv)] + list(rel[i + 1:])
            return [(0.0, zv), (xk, z0 + (z1 - z0) * (xk - x0) / (x1 - x0))] + list(rel[i + 1:])


def sec_landslide(cur, Rn):
    """The surface after a landslide that leaves relief Rn on the surface cur (equal area)."""
    zr = cur[-1][1]
    a = sec_area(cur) - 0.5 * Rn                      # straight line, valley to ridge
    if a + Rn <= zr + EPS:
        return [(0.0, a), (1.0, a + Rn)]
    a = zr - Rn; lo = 0.0; hi = 1.0                   # ridge untouched: find x2 with cut = fill
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if mid * 0.5 * (a + sec_value(cur, mid)) - sec_area(cur, mid) > 0: lo = mid
        else: hi = mid
    x2 = 0.5 * (lo + hi)
    return [(0.0, a), (x2, sec_value(cur, x2))] + [v for v in cur if v[0] > x2 + 1e-9]


def section_states(sim):
    """[{t, rel, pre}]: the relict surface in force from each landslide on (pre = the surface
    it failed from), starting from a flat surface at t = 0."""
    p = sim["p"]; S = p["L"] if int(p["stage"]) <= 3 else p["Ls"]
    rel = [(0.0, 0.0), (1.0, 0.0)]
    out = [{"t": 0.0, "rel": rel, "pre": None}]
    for f in sim["fails"]:
        if f["from"] - f["to"] <= EPS:
            continue
        cur = sec_surface(rel, f["from"], S)
        rel = sec_landslide(cur, f["to"])
        out.append({"t": f["t"], "rel": rel, "pre": cur})
    return out


def ceiling_at(lim, t):
    c = lim[0][1]
    for (tt, cc) in lim:
        if tt <= t: c = cc
        else: break
    return c


def section_at(sim, states, t):
    """The section at time t: pts (about its mean elevation), the relief R, the surface the
    last landslide failed from (ghost, same datum) with the time since it (age), and
    weak = (C0 - C) / C0, the damage standing in the slope."""
    p = sim["p"]; S = p["L"] if int(p["stage"]) <= 3 else p["Ls"]
    k = 0
    while k + 1 < len(states) and states[k + 1]["t"] <= t:
        k += 1
    st = states[k]; R = relief_at(sim["path"], t)
    cur = sec_surface(st["rel"], R, S); m = sec_area(cur)
    return {"pts": [(x, z - m) for x, z in cur], "R": R,
            "ghost": [(x, z - m) for x, z in st["pre"]] if st["pre"] else None, "age": t - st["t"],
            "weak": max(0.0, (S - ceiling_at(sim["lim"], t)) / S)}


def section_fingerprint(stage, extra=None):
    s = simulate(dict(extra or {}, stage=stage)); states = section_states(s); out = "%d: %d" % (stage, len(states) - 1)
    for t in (0.37, 1.91, 3.3):
        pts = section_at(s, states, t)["pts"]
        out += " %.9f %d" % (sum(abs(z) + x for x, z in pts), len(pts))
    return out
