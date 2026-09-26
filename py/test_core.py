"""
test_core.py — validation of core/ laws and the d1 reservoir scheme.

Run:  /opt/anaconda3/bin/python3 py/test_core.py
Each check prints PASS/FAIL with the number it measured. Exit code 1 on any FAIL.
"""
import math, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import gs_core as gc
from d1_reservoir import run_reservoir

fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  [" + detail + "]" if detail else ""))
    if not ok:
        fails.append(name)

# ---- 1. forcing shapes: bounded, zero outside the window, closed-form integrals ----
T = 12e6
closed = {"constant": 1.0, "triangle": 0.5, "plateau": 0.75, "arc": math.pi / 4,
          "bell": (math.sqrt(2 * math.pi) / 6) * math.erf(3 / math.sqrt(2))}
for name, frac in closed.items():
    f = gc.make_series({"shape": name, "peak": 1.0, "duration": T})
    I = gc.integrate(f, 0, T, 20000) / T
    ts = np.linspace(-T, 3 * T, 4001)
    vals = np.array([f(t) for t in ts])
    check("forcing %s integral" % name, abs(I - frac) < 2e-3, "%.4f vs %.4f" % (I, frac))
    check("forcing %s bounded" % name, vals.min() >= 0 and vals.max() <= 1 + 1e-12)
    if name != "constant":
        outside = np.array([f(t) for t in ts if t < 0 or t > T])
        # the bell is a Gaussian with ±3σ inside the window: its tail is exp(-4.5) ≈ 0.011 at the edge
        tol = 0.0112 if name == "bell" else 1e-6
        check("forcing %s ~zero outside" % name, outside.max() < tol, "%.2e" % outside.max())
pulse = gc.make_series({"shape": "pulse", "peak": 1.0, "duration": T})
check("forcing pulse peak at 10 %", abs(pulse(0.1 * T) - 1) < 1e-12)
saw = gc.make_series({"shape": "sawtooth", "peak": -400.0, "base": 1500.0, "period": 1e5})
check("forcing sawtooth ELA range", abs(saw(0) - 1500) < 1e-9 and abs(saw(0.9e5) - 1100) < 1e-6)
sine = gc.make_series({"shape": "sine", "period": 1e5})
check("forcing sine periodic", abs(sine(1.3e5) - sine(0.3e5)) < 1e-12)
check("spatial gaussian centre", abs(gc.spatial_factor(0.5, 0.5, "gaussian") - 1) < 1e-12)
check("spatial ramp", abs(gc.spatial_factor(0.05, 0.5, "ramp", 0.1) - 0.5) < 1e-12 and gc.spatial_factor(0.5, 0.5, "ramp", 0.1) == 1)

# ---- 2. noise: deterministic, in [0,1], 3D reduces sensibly ----
v = [gc.fbm3(x * 0.37, x * 0.11, x * 0.23, 7, 4, 0.5, 2) for x in range(200)]
check("noise fbm3 range", min(v) >= 0 and max(v) <= 1, "%.3f..%.3f" % (min(v), max(v)))
check("noise deterministic", gc.hash3(3, 4, 5, 9) == gc.hash3(3, 4, 5, 9) and gc.hash3(3, 4, 5, 9) != gc.hash3(3, 4, 6, 9))
v2 = [gc.fbm2(x * 0.37, x * 0.11, 7, 4, 0.5, 2) for x in range(200)]
check("noise fbm2 range", min(v2) >= 0 and max(v2) <= 1)
# legacy hash2 reference values (computed from the JS formula by hand-mirroring imul)
check("noise hash2 stable", abs(gc.hash2(0, 0, 0)) < 1e-12 and 0 <= gc.hash2(-5, 12345, 42) < 1)

# ---- 3. lithology in the material frame: a buried layer is reached at the right time ----
# Layer top 300 m below the initial surface (z0 = 1000), 200 m thick. Constant uplift U
# and constant erosion E: surface z = z0 + (U−E)t, U_cum = U t, so zm = z0 − E t.
# The body is exposed when E t = 300 and stripped when E t = 500, whatever U is.
z0, U, E = 1000.0, 2e-3, 1e-3
layer = gc.litho_make({"type": "layer", "top": z0 - 300, "thick": 200})
def r_at(t):
    z = z0 + (U - E) * t; ucum = U * t
    return layer(0, 0, gc.material_z(z, ucum))
check("lithology layer exposure time", r_at(299e3) == 0 and r_at(301e3) == 1 and r_at(499e3) == 1 and r_at(501e3) == 0)
soft = gc.litho_make({"type": "layer", "top": 0, "thick": 100, "soft": 20})
check("lithology soft edge", abs(soft(0, 0, -10) - 0.5) < 1e-12 and soft(0, 0, -50) == 1)
slab = gc.litho_make({"type": "slab", "top": 0, "thick": 100, "dipX": 0.1})
check("lithology slab dips", slab(1000, 0, 50) == 1 and slab(0, 0, 50) == 0)
blob = gc.litho_make({"type": "blob", "seed": 3, "scale": 1000, "threshold": 0.55, "soft": 0})
frac = np.mean([blob(x * 97.0, y * 53.0, z * 71.0) for x in range(20) for y in range(20) for z in range(10)])
check("lithology blob volume fraction sane", 0.05 < frac < 0.6, "%.2f" % frac)
patch = gc.litho_make({"type": "patch", "seed": 3, "scale": 1000, "threshold": 0.55, "soft": 0})
check("lithology patch persistent in zm", patch(100, 200, 0) == patch(100, 200, -5000))
union = gc.litho_make([{"type": "layer", "top": 0, "thick": 100}, {"type": "dike", "x0": 500, "width": 50}])
check("lithology union", union(500, 0, 999) == 1 and union(0, 0, -50) == 1 and union(0, 0, 999) == 0)
check("lithology factors", gc.erodibility_factor(0, 5) == 1 and abs(gc.erodibility_factor(1, 5) - 0.2) < 1e-12
      and gc.strength_factor(1, 3) == 3 and abs(gc.erodibility_factor(0.5, 3) - 0.5) < 1e-12)

# ---- 4. Yuan node update: residual of the implicit equation is zero; limits behave ----
h0, hr, Kp, dt, dist = 500.0, 480.0, 2e-3, 1000.0, 50.0
for nexp in (1.0, 1.5, 2.0):
    nh = gc.yuan_node(h0, h0, hr, Kp, dt, dist, 0.0, nexp)
    s = max(0.0, (nh - hr) / dist)
    resid = nh - h0 + Kp * dt * s ** nexp
    check("yuan residual n=%g" % nexp, abs(resid) < 1e-6, "resid %.2e, nh %.3f" % (resid, nh))
    check("yuan monotone n=%g" % nexp, hr <= nh <= h0)
# explicit vs implicit agree for small dt
nh_i = gc.yuan_node(h0, h0, hr, Kp, 1.0, dist, 0.0, 1.0)
nh_e = h0 - Kp * 1.0 * (h0 - hr) / dist
check("yuan small-dt matches explicit", abs(nh_i - nh_e) < 1e-6, "%.3e" % abs(nh_i - nh_e))
nh_big = gc.yuan_node(h0, h0, hr, 1.0, 1e9, dist, 0.0, 1.0)
check("yuan huge dt approaches receiver", hr <= nh_big < hr + 1e-3, "nh-hr %.2e" % (nh_big - hr))
check("yuan deposit-only branches", gc.yuan_node(h0, h0, hr, Kp, dt, dist, 3.0, 1.0, submerged=True) == h0 + 3.0
      and gc.yuan_node(470.0, 470.0, hr, Kp, dt, dist, 2.0, 1.0) == 472.0)
# deposition is mass conserving in the linear branch: Δh = dep − incision, incision ≥ 0
nh_d = gc.yuan_node(h0, h0, hr, Kp, dt, dist, 5.0, 1.0)
nh_0 = gc.yuan_node(h0, h0, hr, Kp, dt, dist, 0.0, 1.0)
check("yuan deposition raises bed", nh_d > nh_0)

# ---- 5. creep laws ----
check("creep linear sign", gc.creep_flux_linear(0.01, 0.1) < 0)
check("creep roering diverges near Sc", gc.creep_flux_roering(0.01, 0.6, 0.7) < 3 * gc.creep_flux_linear(0.01, 0.6)
      and gc.creep_flux_roering(0.01, 0.7, 0.7) == -0.01 * 0.7 / 0.02)

# ---- 6. glacier laws: SIA inversion round-trips; mass balance shape; sliding ----
for H in (50.0, 200.0, 800.0):
    for S in (0.005, 0.05, 0.3):
        q = gc.sia_flux(H, S)
        H2 = gc.sia_thickness(q, S)
        check("sia round trip H=%g S=%g" % (H, S), abs(H2 - H) < 1e-9 * H, "q=%.3g m2/yr" % q)
check("sia slope floor", gc.sia_thickness(1e4, 0.0) == gc.sia_thickness(1e4, 1e-3))
check("sia zero flux", gc.sia_thickness(0.0, 0.1) == 0)
# deformation speed check against the textbook depth-averaged Glen-law velocity
# u = 2A/(n+2) (ρg S)^n H^(n+1): 300 m / 3 % → ~5 m/yr; 400 m / 5 % → ~70 m/yr (sliding adds to this)
for H, S, lo, hi in ((300.0, 0.03, 2, 10), (400.0, 0.05, 30, 150)):
    u = gc.sia_flux(H, S) / H
    ua = 2 * gc.ICE["A"] / (gc.ICE["N"] + 2) * (gc.ICE["RHOG"] * S) ** gc.ICE["N"] * H ** (gc.ICE["N"] + 1)
    check("sia deformation speed H=%g S=%g" % (H, S), lo < u < hi and abs(u - ua) < 1e-9 * ua, "u=%.1f m/yr" % u)
check("mass balance asymmetry", abs(gc.mass_balance(1100, 1000, 0.007) - 0.7) < 1e-12 and abs(gc.mass_balance(900, 1000, 0.007) + 1.75) < 1e-12)
check("mass balance clamps", gc.mass_balance(3000, 1000, 0.007, accMax=2.0) == 2.0 and gc.mass_balance(-3000, 1000, 0.007, abMax=8.0) == -8.0)
check("sliding zero when no ice", gc.sliding_speed(100.0, 0.0, 0.5) == 0 and gc.sliding_speed(100.0, 50.0, 0.5) == 1.0)
check("abrasion power law", gc.abrasion(1e-4, 100.0, 1.0) == 1e-2 and gc.abrasion(1e-4, -5.0, 1.0) == 0)
check("quarrying only on convexity", gc.quarrying(1.0, 10.0, 0.01) == 0.1 and gc.quarrying(1.0, 10.0, -0.01) == 0)
check("buzzsaw peaks at ELA", gc.buzzsaw_rate(1000, 1000, 300, 2e-3) == 2e-3 and gc.buzzsaw_rate(1300, 1000, 300, 2e-3) < 1.3e-3)

# ---- 7. d1 reservoir: steady state, dt-independence, reference ODE, float32, legacy units ----
res = run_reservoir({"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 3000.0, "duration": 40e6, "tailFactor": 1.0, "steps": 4000})
Rend = res["series"][-1]["relief"]
check("reservoir steady state = relief limit", abs(Rend - 3000) < 1.0, "R=%.2f m" % Rend)
a = run_reservoir({"steps": 400}); b = run_reservoir({"steps": 4000})
ra = np.interp([s["t"] for s in b["series"]], [s["t"] for s in a["series"]], [s["relief"] for s in a["series"]])
rb = np.array([s["relief"] for s in b["series"]])
check("reservoir dt-independent (400 vs 4000 steps)", np.max(np.abs(ra - rb)) < 2.0, "max diff %.3f m" % np.max(np.abs(ra - rb)))
try:
    from scipy.integrate import solve_ivp
    p = b["params"]; k = b["k"]; n = p["exponent"]
    U = gc.make_series({"shape": p["shape"], "peak": p["peakUplift"], "duration": p["duration"]})
    sol = solve_ivp(lambda t, y: [U(t) - k * max(y[0], 0) ** n], [0, p["duration"] * p["tailFactor"]], [0.0],
                    t_eval=[s["t"] for s in b["series"]], rtol=1e-10, atol=1e-6, max_step=p["duration"] / 200)
    err = np.max(np.abs(sol.y[0] - rb))
    check("reservoir vs scipy reference", err < 0.5, "max diff %.3f m" % err)
except ImportError:
    print("SKIP scipy reference")
# float32: rerun the RK4 loop in float32 arithmetic and compare
p = b["params"]; k = np.float32(b["k"]); n = np.float32(p["exponent"])
U = gc.make_series({"shape": p["shape"], "peak": p["peakUplift"], "duration": p["duration"]})
dt = np.float32(p["duration"] * p["tailFactor"] / p["steps"]); R = np.float32(0)
f = lambda t, R: np.float32(U(float(t))) - k * np.float32(max(R, 0)) ** n
r32 = []
for i in range(p["steps"] + 1):
    t = np.float32(i) * dt
    r32.append(float(R))
    k1 = f(t, R); k2 = f(t + dt / 2, R + (dt / 2) * k1); k3 = f(t + dt / 2, R + (dt / 2) * k2); k4 = f(t + dt, R + dt * k3)
    R = np.float32(R + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4))
    if R < 0: R = np.float32(0)
err32 = np.max(np.abs(np.array(r32) - rb))
check("reservoir float32 vs float64", err32 < 0.5, "max diff %.3f m" % err32)
# legacy units (km, Myr) give the same answer once converted: R km→m, U km/Myr→m/yr, t Myr→yr
def legacy_run(duration=12.0, peak=3.0, rlim=3.0, n=4.0, steps=2000):
    kk = peak / rlim ** n; tot = duration * 2; dtt = tot / steps; R = 0.0; out = []
    Uf = lambda t: peak * math.exp(-0.5 * ((t - duration / 2) / (duration / 6)) ** 2)
    g = lambda t, R: Uf(t) - kk * max(R, 0) ** n
    for i in range(steps + 1):
        t = i * dtt; out.append(R)
        k1 = g(t, R); k2 = g(t + dtt / 2, R + dtt / 2 * k1); k3 = g(t + dtt / 2, R + dtt / 2 * k2); k4 = g(t + dtt, R + dtt * k3)
        R = max(0.0, R + dtt / 6 * (k1 + 2 * k2 + 2 * k3 + k4))
    return np.array(out) * 1000.0
si = np.array([s["relief"] for s in run_reservoir({"steps": 2000})["series"]])
errL = np.max(np.abs(legacy_run() - si))
check("reservoir SI == legacy km/Myr", errL < 1e-6, "max diff %.2e m" % errL)

print()
if fails:
    print("FAILED: " + ", ".join(fails)); sys.exit(1)
print("all checks passed")
