"""
glacier_buttress.py — mirror of relief_limit/glacier_model.js: relief against a threshold
that a glacier both buttresses and weakens, over two glacial cycles (the last 200 kyr).

Time t in kyr from 200 ka (t = 0) to the present (t = 200). Relief and threshold in units of
the intact threshold with no ice, T0 = 1. Ice thickness H in units of the largest glacial
thickness (0..1).

  ice         H(t) = max(Hmin, envelope(t) + fluctuation(t)), a polyline on a dtH grid.
              envelope: glacial (MIS 6) - interglacial (MIS 5e) - glacial (to the LGM) -
              interglacial (Holocene), piecewise linear through the keys of the ice history
              (p["ice"]): "thick" (ICE_KEYS) never reaches zero, so ice is present throughout;
              "retreat" (RETREAT_KEYS) is gone early in both interglacials, re-advances in the
              late Holocene and retreats completely at the end. fluctuation: internal dynamics,
              seeded, knots every 1.5-4 kyr with Gaussian values
              sigma*(0.4 + 0.6*envelope)*min(1, envelope/taper) (no taper for "thick"),
              cosine-interpolated between knots.
  threshold   T(t) = T0 - Dmg + B*b(t),  Dmg = Dc*(1 - exp(-D/Dc)),  Dc = T0 - floor:
              damage D takes the cohesion, a share of what is left at each fracture, never the
              friction (the floor); for small D, Dmg = D
              buttressing: b = (1 - exp(-H/Hs)) / (1 - exp(-1/Hs)), Hs = 0.3, at the nodes of
                           the ice polyline, linear between them: saturating, so the response
                           is stretched where the ice is thin and compressed where it is thick
              weakening:   the stress change is the change in buttressing |d(B*b)|; it
                           accumulates from the last fracture, and when it reaches a quota
                           q = qMean*(1 + qJit*(2u - 1)) (u uniform, its own stream) the rock
                           fractures: D += W*q, a step down in T
              a landslide removes the damaged rock: D <- 0
  relief      dR/dt = E where there is ice (glacial erosion, constant), 0 on an ice-free grid
              cell (both ends H = 0); with "thick" ice is present the whole time
  failure     where R reaches T (erosion, or T falling as the ice thins) or a fracture
              steps T to or below R:  R <- max(T - d, floor + (1 - dFrac)*(R - floor)), then Dmg <- 0;
              d from a truncated power law on [dMin, dMax] (its own stream)

H and b are linear on each grid cell, so T and R are linear between events and every event time
is found exactly (no time step beyond the ice polyline itself).
"""
import math
from relief_limit import mulberry32

EPS = 1e-9

# (age ka, envelope): glacial build-up to MIS 6, termination II, the last interglacial,
# the stepped build-up of the last glacial to the LGM, termination I, the Holocene, and
# climate change: the ice plunges in the last 600 years. A thick glacier: interglacials keep
# about 40 % of the glacial maximum, so the plunge takes it lower than any interglacial.
ICE_KEYS = [(200, 0.65), (186, 0.73), (176, 0.69), (160, 0.80), (148, 0.86), (140, 0.97),
            (135, 0.93), (129, 0.41), (117, 0.38), (110, 0.57), (103, 0.50), (92, 0.64),
            (85, 0.57), (70, 0.73), (60, 0.79), (50, 0.71), (35, 0.76), (24, 1.00),
            (19, 0.94), (11, 0.40), (0.6, 0.38), (0, 0.04)]
# The same glacials, but the glacier retreats completely early in both interglacials: gone by
# 128 ka, back briefly in a small stade (125-122 ka, to 0.28 at 123.5 ka), gone again until it
# re-advances from 118 ka, gone by 10 ka until a late-Holocene re-advance from
# 4.5 ka to a maximum at 0.3 ka (Little Ice Age), then the final retreat.
RETREAT_KEYS = [(200, 0.65), (186, 0.73), (176, 0.69), (160, 0.80), (148, 0.86), (140, 0.97),
                (135, 0.93), (128, 0.0), (125, 0.0), (123.5, 0.28), (122, 0.0), (118, 0.0), (110, 0.50), (103, 0.45), (92, 0.62),
                (85, 0.55), (70, 0.73), (60, 0.79), (50, 0.71), (35, 0.76), (24, 1.00),
                (19, 0.94), (10, 0.0), (4.5, 0.0), (0.3, 0.45), (0, 0.0)]
ICE = {"thick": dict(keys=ICE_KEYS, Hmin=0.05, taper=0.0),
       "retreat": dict(keys=RETREAT_KEYS, Hmin=0.0, taper=0.15)}

DEFAULTS = dict(
    tmax=200.0,   # kyr
    dtH=0.25,     # kyr, the ice polyline's grid
    T0=1.0,       # intact threshold with no ice
    B=0.5,        # buttressing at the largest thickness (fixed)
    Hs=0.3,       # saturation thickness of buttressing: stretched for thin ice, compressed for thick
    ice="thick",  # ice history: "thick" or "retreat" (ICE)
    sigma=0.09,   # amplitude of the internal fluctuations
    E=1.0,        # erosion, as a multiple of Eref
    Eref=0.012,   # relief per kyr at E = 1
    W=1.0,        # weakening: drop of T per unit change in buttressing at a fracture
    qMean=0.025,  # mean change in buttressing between fractures (threshold units)
    qJit=0.25,    # quota spread: q uniform in qMean*(1 ± qJit)
    R0=1.3,      # above the mean relief and the ice-free threshold: a slope held up by the ice
    dFrac=0.5,    # a landslide takes at most this fraction of the relief standing above the floor, so
                  # relief nears the floor by ever smaller drops and never sits on it
    floor=0.7,    # frictional floor: broken rock still stands at its friction angle, so no landslide
                  # leaves less relief than this (py/xsection_experiment.py); was 0.165 until 2026-11-27
    dMin=0.05, dMax=0.5, dAlpha=1.0,
    seed=7,
)


def gauss(rnd):
    u = 1.0 - rnd(); v = rnd()
    return math.sqrt(-2.0 * math.log(u)) * math.cos(2.0 * math.pi * v)


def envelope(t, tmax, keys=ICE_KEYS):
    a = tmax - t
    for (a0, h0), (a1, h1) in zip(keys, keys[1:]):
        if a <= a0 and a >= a1:
            return h0 + (h1 - h0) * (a0 - a) / (a0 - a1)
    return keys[-1][1] if a < keys[-1][0] else keys[0][1]


def make_ice(p):
    rnd = mulberry32(int(p["seed"])); sc = ICE[p["ice"]]
    tmax = p["tmax"]; kt = [-2.0 * rnd()]; kv = []
    while kt[-1] <= tmax:
        kt.append(kt[-1] + 1.5 + 2.5 * rnd())
    for t in kt:
        kv.append(gauss(rnd))
    n = int(round(tmax / p["dtH"])); out = []; j = 0
    for i in range(n + 1):
        t = tmax * i / n
        while j + 2 < len(kt) and kt[j + 1] <= t:
            j += 1
        s = (t - kt[j]) / (kt[j + 1] - kt[j]); w = 0.5 - 0.5 * math.cos(math.pi * s)
        env = envelope(t, tmax, sc["keys"])
        amp = p["sigma"] * (0.4 + 0.6 * env) * (min(1.0, env / sc["taper"]) if sc["taper"] > 0 else 1.0)
        out.append([t, max(sc["Hmin"], env + amp * (kv[j] + (kv[j + 1] - kv[j]) * w))])
    return out


def make_sizes(p):
    lo, hi, a = min(p["dMin"], p["dMax"]), max(p["dMin"], p["dMax"]), p["dAlpha"]
    if hi <= lo * (1 + 1e-9) or a <= 0:
        return lambda: lo
    rnd = mulberry32(int(p["seed"]) + 1000003); span = 1.0 - math.pow(lo / hi, a)
    return lambda: lo * math.pow(1.0 - rnd() * span, -1.0 / a)


def make_quotas(p):
    rnd = mulberry32(int(p["seed"]) + 2000003); qm = p["qMean"]; j = p["qJit"]
    return lambda: qm * (1.0 + j * (2.0 * rnd() - 1.0))


def butt(H, p):
    """buttressing per unit B: saturating in H, 0 at H = 0 and 1 at the largest thickness"""
    return (1.0 - math.exp(-H / p["Hs"])) / (1.0 - math.exp(-1.0 / p["Hs"]))


def simulate(params=None):
    p = dict(DEFAULTS); p.update(params or {})
    ice = make_ice(p); nextD = make_sizes(p); nextQ = make_quotas(p)
    E = p["E"] * p["Eref"]; B = p["B"]; T0 = p["T0"]; W = p["W"]; floor = max(p["floor"], 0.0)
    Dc = max(T0 - floor, 0.0)                         # the cohesive part of the threshold: all that damage can take
    Dmg = 0.0; Draw = 0.0; acc = 0.0; q = nextQ()     # Dmg the threshold lost, Draw the damage (sum of W*q) behind it
    t = 0.0; H = ice[0][1]; Bt = B * butt(H, p)
    R = min(p["R0"], T0 + Bt)
    path = [[t, R]]; thr = [[t, T0 + Bt]]; fracs = []; fails = []

    def fail(cause):
        nonlocal R, Dmg, Draw
        T = T0 - Dmg + Bt
        d = nextD(); Rn = max(T - d, floor + (1.0 - p["dFrac"]) * max(R - floor, 0.0))
        if Rn < R:
            fails.append(dict(t=t, frm=R, to=Rn, d=d, lim=T, dmg=Dmg, cause=cause))
            R = Rn; path.append([t, R])
        if Dmg > 0:
            thr.append([t, T]); Dmg = 0.0; Draw = 0.0; thr.append([t, T0 + Bt])

    for i in range(len(ice) - 1):
        t1, H1 = ice[i + 1]
        s = (H1 - H) / (t1 - t) if t1 > t else 0.0   # dH/dt on this cell
        Bt1 = B * butt(H1, p)
        sb = (Bt1 - Bt) / (t1 - t) if t1 > t else 0.0  # d(buttressing)/dt on this cell
        Ec = E if (H > 0 or H1 > 0) else 0.0          # glacial erosion only where there is ice
        while True:
            gap = T0 - Dmg + Bt - R                  # >= 0
            closing = Ec - sb                           # d(R - T)/dt
            tf = t + gap / closing if closing > 0 else math.inf
            tq = t + (q - acc) / abs(sb) if sb != 0 else math.inf
            tn = min(tf, tq, t1)
            if tn >= t1 - EPS and tf > t1 - EPS and tq > t1 - EPS:
                acc += abs(sb) * (t1 - t); R += Ec * (t1 - t); t = t1; H = H1; Bt = Bt1
                R = min(R, T0 - Dmg + Bt)
                break
            dt = tn - t
            acc += abs(sb) * dt; R += Ec * dt; H += s * dt; Bt += sb * dt; t = tn
            path.append([t, R]); thr.append([t, T0 - Dmg + Bt])
            if tq <= tf:                                # fracture
                Told = T0 - Dmg + Bt
                Draw += W * q; Dn = Dc * (1.0 - math.exp(-Draw / Dc)) if Dc > 0 else 0.0; acc = 0.0
                fracs.append(dict(t=t, drop=Dn - Dmg, q=q, frm=Told, to=Told - (Dn - Dmg))); Dmg = Dn
                q = nextQ()
                thr.append([t, T0 - Dmg + Bt])
                if R >= T0 - Dmg + Bt - EPS:
                    fail("fracture")
            else:
                R = T0 - Dmg + Bt
                fail("limit")
        path.append([t, R]); thr.append([t, T0 - Dmg + Bt])
    return dict(path=path, thr=thr, ice=ice, fracs=fracs, fails=fails, p=p)


# ---- The cross-section: a glacial trough read from the same history -------------------------
# A half-section z(x), valley centre at x = 0, crest at x = 1 (x in units of the half-spacing,
# z in units of T0), mirrored about the centre and the crest when drawn. Two shapes, and the
# surface is the lower of them: a U cut by the glacier from the floor up, and a straight slope
# left by landslides from the crest down.
#   incision    the glacier cuts the floor down and a parabolic trough z_floor + C*x^2 (C the
#               highest threshold, T0 + B: a full trough stands at it) climbs from it; above
#               where it meets the older surface, that surface is untouched:
#               z(x) = min(relict(x), z_floor + C*x^2),  z_floor = z_crest - R(t)
#               (the parabola as a polyline on SEC_N cells; the min is exact for polylines).
#               The U climbs the wall as the floor goes down: the glacier oversteepens the wall
#   landslide   R -> R': the crest comes down by R - R' and the floor stays. The failure surface
#               is a plane from the new crest down to the wall, at the slope SEC_S or, where a
#               plane that steep would miss the wall, at the steepest slope that still meets it
#               (s_tan: the tangent from the new crest to the trough):
#               z(x) = min(old(x), z_crest' - s*(1 - x))  for x >= x_cut,   s = min(SEC_S, s_tan),
#               x_cut the highest point where the plane meets the old surface; old(x) below it.
#               A wedge off the upper wall: the trough below x_cut is never touched, so the
#               floor stays a U. The debris is not kept in the section (it lands on the glacier
#               and is carried away)
#   ice         the glacier fills the trough to min(SEC_ICE*C*H(t), SEC_CAP*R) above the floor,
#               flat-topped: never over the crest
#   weakening   a lens under the upper wall (glacier_app.js), depth (D / T_intact)*R*6.75*s^2*(1 - s)
#   datum       the floor: drawn at -SEC_DATUM*C, so a full trough is about centred
# Nothing here feeds back on R(t).
SEC_N = 32
SEC_ICE = 0.6
SEC_CAP = 0.9   # the ice never fills more than this fraction of the relief
SEC_DATUM = 0.45  # the floor is drawn at -SEC_DATUM*C
SEC_S = 1.2     # slope of a landslide's failure surface (dz/dx in section units; 43 deg as drawn)

def value_at(poly, t):
    """value on a time-ordered polyline at t, after any step at exactly t (as glacier_model.js)"""
    lo, hi = 0, len(poly) - 1
    while lo < hi:
        mid = (lo + hi + 1) >> 1
        if poly[mid][0] <= t: lo = mid
        else: hi = mid - 1
    if lo == len(poly) - 1: return poly[lo][1]
    a, b = poly[lo], poly[lo + 1]
    return a[1] if b[0] == a[0] else a[1] + (b[1] - a[1]) * (t - a[0]) / (b[0] - a[0])

def sec_area(poly, x_to=1.0):
    a = 0.0
    for p0, p1 in zip(poly, poly[1:]):
        if x_to <= p0[0]: break
        if x_to < p1[0]:
            zc = p0[1] + (p1[1] - p0[1]) * (x_to - p0[0]) / (p1[0] - p0[0])
            return a + 0.5 * (p0[1] + zc) * (x_to - p0[0])
        a += 0.5 * (p0[1] + p1[1]) * (p1[0] - p0[0])
    return a

def sec_value(poly, x):
    for p0, p1 in zip(poly, poly[1:]):
        if x <= p1[0]:
            return p0[1] if p1[0] == p0[0] else p0[1] + (p1[1] - p0[1]) * (x - p0[0]) / (p1[0] - p0[0])
    return poly[-1][1]

def sec_min(a, b):
    """pointwise minimum of two polylines on [0, 1], crossings inserted (exact)"""
    xs = sorted(set([v[0] for v in a] + [v[0] for v in b]))
    out = []
    for i, x in enumerate(xs):
        za, zb = sec_value(a, x), sec_value(b, x)
        if i:
            x0 = xs[i - 1]; da0 = sec_value(a, x0) - sec_value(b, x0); da1 = za - zb
            if da0 * da1 < 0:
                xc = x0 + (x - x0) * da0 / (da0 - da1)
                out.append([xc, sec_value(a, xc)])
        out.append([x, min(za, zb)])
    clean = [out[0]]
    for v in out[1:]:
        if v[0] - clean[-1][0] > 1e-12: clean.append(v)
        else: clean[-1] = [clean[-1][0], min(clean[-1][1], v[1])]
    return clean

def sec_trough(zv, C):
    return [[i / SEC_N, zv + C * (i / SEC_N) ** 2] for i in range(SEC_N + 1)]

def sec_surface(rel, R, C):
    zr = rel[-1][1]; zv = zr - R
    if rel[0][1] - zv <= EPS: return [v[:] for v in rel]
    return sec_min(rel, sec_trough(zv, max(C, R)))

def sec_landslide(cur, Rn):
    """(surface, x_cut, slope): the wedge above a plane from the new crest (floor + Rn) is
    removed. The plane has slope min(SEC_S, s_tan), s_tan the steepest slope at which a line
    from the new crest still meets the old surface below it; x_cut is the highest point where
    it does. Below x_cut the old surface is kept."""
    zf = cur[0][1]; zc = zf + Rn
    st = max((zc - v[1]) / (1.0 - v[0]) for v in cur if v[0] < 1.0 - 1e-9)
    sl = min(SEC_S, st)
    new = sec_min(cur, [[0.0, zc - sl], [1.0, zc]]); i = len(new) - 2
    while i > 0 and new[i][1] < sec_value(cur, new[i][0]) - 1e-12: i -= 1
    xc = new[i][0]
    return [v[:] for v in cur if v[0] < xc - 1e-12] + new[i:], xc, sl

def sec_C(p):
    return p["T0"] + p["B"]

def section_states(sim):
    """[{t, rel, pre}]: the relict surface from each landslide on; at t = 0 a trough of relief R0
    under a crest at 0, with the failure slope from the crest down to it"""
    p = sim["p"]; C = sec_C(p); R0 = sim["path"][0][1]
    rel, _, _ = sec_landslide(sec_min([[0.0, 0.0], [1.0, 0.0]], sec_trough(-R0, max(C, R0))), R0)
    out = [dict(t=0.0, rel=rel, pre=None, ghost=None)]
    for f in sim["fails"]:
        if f["frm"] - f["to"] <= EPS: continue
        cur = sec_surface(rel, f["frm"], C); rel, xc, sl = sec_landslide(cur, f["to"])
        # ghost: the part of the old surface that failed (x_cut to the crest)
        out.append(dict(t=f["t"], rel=rel, pre=cur, xcut=xc, slope=sl,
                        ghost=[[xc, sec_value(cur, xc)]] + [v[:] for v in cur if v[0] > xc + 1e-12]))
    return out

def section_at(sim, states, t):
    p = sim["p"]; C = sec_C(p); k = 0
    while k + 1 < len(states) and states[k + 1]["t"] <= t: k += 1
    st = states[k]; R = value_at(sim["path"], t)
    cur = sec_surface(st["rel"], R, C); m = cur[0][1] + SEC_DATUM * C
    pts = [[v[0], v[1] - m] for v in cur]
    # intact threshold at t: buttressing at the ice nodes, linear between (as the model)
    ice = sim["ice"]; n = len(ice) - 1; i = min(max(int(t / p["dtH"]), 0), n - 1)
    ba, bb = p["B"] * butt(ice[i][1], p), p["B"] * butt(ice[i + 1][1], p)
    Ti = p["T0"] + ba + (bb - ba) * (t - ice[i][0]) / (ice[i + 1][0] - ice[i][0])
    H = value_at(ice, t)
    # rock uplift relative to the floor so far = glacial erosion so far (the texture's offset)
    up = R - sim["path"][0][1] + sum(f["frm"] - f["to"] for f in sim["fails"] if f["t"] <= t)
    return dict(up=up, pts=pts, R=R, ghost=[[v[0], v[1] - m] for v in st["ghost"]] if st["ghost"] else None,
                age=t - st["t"], weak=max(0.0, (Ti - value_at(sim["thr"], t)) / Ti),
                H=H, zIce=pts[0][1] + min(SEC_ICE * C * H, SEC_CAP * R))

def section_fingerprint(extra=None):
    s = simulate(extra); states = section_states(s); out = "%d" % (len(states) - 1)
    for t in [37.3, 101.9, 171.3, 199.9]:
        sec = section_at(s, states, t); pts = sec["pts"]
        out += " %.9f %d %.9f" % (sum(abs(v[1]) + v[0] for v in pts), len(pts), sec["zIce"] + sec["weak"])
        out += " %d" % (len(sec["ghost"]) if sec["ghost"] else 0)
    return out


def fingerprint(extra=None):
    s = simulate(extra)
    sp = sum(v[1] for v in s["path"]); st = sum(v[1] for v in s["thr"]); si = sum(v[1] for v in s["ice"])
    return "%d %d %d %d %.9f %.9f %.9f" % (len(s["path"]), len(s["thr"]), len(s["fracs"]), len(s["fails"]), sp, st, si)


if __name__ == "__main__":
    s = simulate()
    print(fingerprint())
    for f in s["fails"]:
        print("fail t=%.2f (%.1f ka) %s %.3f -> %.3f dmg %.3f" % (f["t"], 200 - f["t"], f["cause"], f["frm"], f["to"], f["dmg"]))
    print("fractures", len(s["fracs"]))
    print("section", section_fingerprint(), "|", section_fingerprint({"ice": "retreat"}))
