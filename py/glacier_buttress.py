"""
glacier_buttress.py — mirror of relief_limit/glacier_model.js: relief against a threshold
that a glacier both buttresses and weakens, over two glacial cycles (the last 200 kyr).

Time t in kyr from 200 ka (t = 0) to the present (t = 200). Relief and threshold in units of
the intact threshold with no ice, T0 = 1. Ice thickness H in units of the largest glacial
thickness (0..1).

  ice         H(t) = max(Hmin, envelope(t) + fluctuation(t)), a polyline on a dtH grid.
              envelope: glacial (MIS 6) - interglacial (MIS 5e) - glacial (to the LGM) -
              interglacial (Holocene), piecewise linear through ICE_KEYS; it never reaches
              zero, so ice is present throughout. fluctuation: internal dynamics, seeded,
              knots every 1.5-4 kyr with Gaussian values sigma*(0.4 + 0.6*envelope),
              cosine-interpolated between knots.
  threshold   T(t) = T0 - Dmg + B*b(t)
              buttressing: b = (1 - exp(-H/Hs)) / (1 - exp(-1/Hs)), Hs = 0.3, at the nodes of
                           the ice polyline, linear between them: saturating, so the response
                           is stretched where the ice is thin and compressed where it is thick
              weakening:   the stress change is the change in buttressing |d(B*b)|; it
                           accumulates from the last fracture, and when it reaches a quota
                           q = qMean*(1 + qJit*(2u - 1)) (u uniform, its own stream) the rock
                           fractures: Dmg += W*q, a step down in T
              a landslide removes the damaged rock: Dmg <- 0
  relief      dR/dt = E (glacial erosion, constant: ice is present the whole time)
  failure     where R reaches T (erosion, or T falling as the ice thins) or a fracture
              steps T to or below R:  R <- max(T - d, floor), then Dmg <- 0;
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

DEFAULTS = dict(
    tmax=200.0,   # kyr
    dtH=0.25,     # kyr, the ice polyline's grid
    T0=1.0,       # intact threshold with no ice
    B=0.5,        # buttressing at the largest thickness (fixed)
    Hs=0.3,       # saturation thickness of buttressing: stretched for thin ice, compressed for thick
    Hmin=0.05,    # thinnest ice
    sigma=0.09,   # amplitude of the internal fluctuations
    E=1.0,        # erosion, as a multiple of Eref
    Eref=0.012,   # relief per kyr at E = 1
    W=1.0,        # weakening: drop of T per unit change in buttressing at a fracture
    qMean=0.025,  # mean change in buttressing between fractures (threshold units)
    qJit=0.25,    # quota spread: q uniform in qMean*(1 ± qJit)
    R0=1.3,      # above the mean relief and the ice-free threshold: a slope held up by the ice
    floor=0.165,  # 10 % of the graph height, 1.1*(T0 + B)
    dMin=0.05, dMax=0.5, dAlpha=1.0,
    seed=7,
)


def gauss(rnd):
    u = 1.0 - rnd(); v = rnd()
    return math.sqrt(-2.0 * math.log(u)) * math.cos(2.0 * math.pi * v)


def envelope(t, tmax):
    a = tmax - t
    for (a0, h0), (a1, h1) in zip(ICE_KEYS, ICE_KEYS[1:]):
        if a <= a0 and a >= a1:
            return h0 + (h1 - h0) * (a0 - a) / (a0 - a1)
    return ICE_KEYS[-1][1] if a < ICE_KEYS[-1][0] else ICE_KEYS[0][1]


def make_ice(p):
    rnd = mulberry32(int(p["seed"]))
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
        env = envelope(t, tmax)
        fl = p["sigma"] * (0.4 + 0.6 * env) * (kv[j] + (kv[j + 1] - kv[j]) * w)
        out.append([t, max(p["Hmin"], env + fl)])
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
    Dmg = 0.0; acc = 0.0; q = nextQ()
    t = 0.0; H = ice[0][1]; Bt = B * butt(H, p)
    R = min(p["R0"], T0 + Bt)
    path = [[t, R]]; thr = [[t, T0 + Bt]]; fracs = []; fails = []

    def fail(cause):
        nonlocal R, Dmg
        T = T0 - Dmg + Bt
        d = nextD(); Rn = max(T - d, floor)
        if Rn < R:
            fails.append(dict(t=t, frm=R, to=Rn, d=d, lim=T, dmg=Dmg, cause=cause))
            R = Rn; path.append([t, R])
        if Dmg > 0:
            thr.append([t, T]); Dmg = 0.0; thr.append([t, T0 + Bt])

    for i in range(len(ice) - 1):
        t1, H1 = ice[i + 1]
        s = (H1 - H) / (t1 - t) if t1 > t else 0.0   # dH/dt on this cell
        Bt1 = B * butt(H1, p)
        sb = (Bt1 - Bt) / (t1 - t) if t1 > t else 0.0  # d(buttressing)/dt on this cell
        while True:
            gap = T0 - Dmg + Bt - R                  # >= 0
            closing = E - sb                            # d(R - T)/dt
            tf = t + gap / closing if closing > 0 else math.inf
            tq = t + (q - acc) / abs(sb) if sb != 0 else math.inf
            tn = min(tf, tq, t1)
            if tn >= t1 - EPS and tf > t1 - EPS and tq > t1 - EPS:
                acc += abs(sb) * (t1 - t); R += E * (t1 - t); t = t1; H = H1; Bt = Bt1
                R = min(R, T0 - Dmg + Bt)
                break
            dt = tn - t
            acc += abs(sb) * dt; R += E * dt; H += s * dt; Bt += sb * dt; t = tn
            path.append([t, R]); thr.append([t, T0 - Dmg + Bt])
            if tq <= tf:                                # fracture
                Told = T0 - Dmg + Bt
                Dmg += W * q; acc = 0.0
                fracs.append(dict(t=t, drop=W * q, q=q, frm=Told, to=Told - W * q))
                q = nextQ()
                thr.append([t, T0 - Dmg + Bt])
                if R >= T0 - Dmg + Bt - EPS:
                    fail("fracture")
            else:
                R = T0 - Dmg + Bt
                fail("limit")
        path.append([t, R]); thr.append([t, T0 - Dmg + Bt])
    return dict(path=path, thr=thr, ice=ice, fracs=fracs, fails=fails, p=p)


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
