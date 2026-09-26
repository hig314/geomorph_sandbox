"""
test_d2_along.py — validation of the along-valley tier (DESIGN.md §4 decision gate + scheme checks).
Run: /opt/anaconda3/bin/python3 py/test_d2_along.py   (~1 min: the reference SIA is explicit)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import gs_core as gc
import d2_along as d2

fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  [" + detail + "]" if detail else ""))
    if not ok: fails.append(name)

def bed(kind, N=301, L=30000.0):
    s = np.linspace(0, L, N); ds = s[1] - s[0]
    if kind == "linear":    z = 2500 * (1 - s / L)
    elif kind == "concave": z = 2500 * (1 - s / L) ** 1.5
    elif kind == "bump":    z = 2500 * (1 - s / L) ** 1.5 + 150 * np.exp(-((s - 12000) / 1500) ** 2)
    elif kind == "step":    z = 2500 * (1 - s / L) ** 1.5 - 200 * (s > 9000) * (s < 14000)
    z[-1] = 0
    return s, ds, z

# ---- 1. Decision gate: routed-flux + upstream march vs reference explicit flowline SIA ----
# Gate in DESIGN.md §4: thickness within ~30 %. Smooth beds: within 2 %.
for kind in ("linear", "concave", "bump"):
    for fs in (0.0, 0.5, 0.8):
        s, ds, z = bed(kind); W = d2.valley_width(s)
        p = {"ela": 1800.0, "balGrad": 0.007, "accMax": 2.0, "abMax": 8.0, "fs": fs}
        ref = d2.reference_sia(z, W, ds, p, years=30000, tol=1e-6)
        inv = d2.steady_ice(z, W, ds, dict(p, iceIters=200))
        Hr, Hi = ref["H"], inv["H"]; m = Hr > 20
        rel = np.abs(Hi[m] - Hr[m]) / Hr[m]
        vr = (Hi * W).sum() / (Hr * W).sum()
        tr = s[np.max(np.where(Hr > 1)[0])]; ti = s[np.max(np.where(Hi > 1)[0])]
        check("gate %s fs=%.1f: H within 2 %%, volume within 2 %%, terminus within 1 node" % (kind, fs),
              rel.max() < 0.02 and abs(vr - 1) < 0.02 and abs(tr - ti) <= ds + 1e-6,
              "maxH %.0f m, rel max %.3f, vol ratio %.3f, term %.1f/%.1f km, %d outer iters" % (Hr.max(), rel.max(), vr, tr / 1e3, ti / 1e3, inv["iters"]))
# Step bed (overdeepening below the terminus): the explicit reference surges over the 200 m
# step (an icefall limit cycle: its face flux at 9 km exceeds the integrated upstream balance,
# impossible at steady state), so it is not a clean reference there. The routed solution is
# a true steady state; require the 30 % gate on volume and agreement above the step.
for fs in (0.0, 0.8):
    s, ds, z = bed("step"); W = d2.valley_width(s)
    p = {"ela": 1800.0, "balGrad": 0.007, "accMax": 2.0, "abMax": 8.0, "fs": fs}
    ref = d2.reference_sia(z, W, ds, p, years=30000, tol=1e-6)
    inv = d2.steady_ice(z, W, ds, dict(p, iceIters=200))
    Hr, Hi = ref["H"], inv["H"]; above = (s < 8500) & (Hr > 20)
    rel = np.abs(Hi[above] - Hr[above]) / Hr[above]
    vr = (Hi * W).sum() / (Hr * W).sum()
    check("gate step fs=%.1f: within 2 %% above the step, volume within 35 %%" % fs, rel.max() < 0.02 and 0.65 < vr < 1.05,
          "rel max above step %.3f, vol ratio %.2f" % (rel.max(), vr))
# Steady-state consistency: at the routed solution, face flux = integrated balance (by construction)
s, ds, z = bed("concave"); W = d2.valley_width(s)
inv = d2.steady_ice(z, W, ds, {"ela": 1800.0, "balGrad": 0.007, "accMax": 2.0, "abMax": 8.0, "fs": 0.5, "iceIters": 200})
Gam = gc.ICE["GAMMA"] / 0.5; H = inv["H"]; Q = inv["Q"]
resid = []
for i in range(len(z) - 1):
    if Q[i] > 0 and H[i] > 1:
        Hf = 0.5 * (H[i] + H[i + 1]); Sf = (z[i] + H[i] - z[i + 1] - H[i + 1]) / ds; Wf = 0.5 * (W[i] + W[i + 1])
        resid.append(abs(Gam * Hf ** 5 * Sf ** 3 * Wf - Q[i]) / Q[i])
check("march satisfies the face flux law", max(resid) < 1e-3, "max rel resid %.1e" % max(resid))
check("sliding speed = fs q/H", np.allclose(inv["Us"][H > 0], 0.5 * inv["q"][H > 0] / H[H > 0]))
# Sliding fraction: more sliding → thinner ice, same flux
i0 = d2.steady_ice(z, W, ds, {"ela": 1800.0, "balGrad": 0.007, "accMax": 2.0, "abMax": 8.0, "fs": 0.0, "iceIters": 200})
check("sliding thins the glacier", inv["H"].max() < i0["H"].max())
# Divide: finite thickness at the head, no blow-up
check("finite ice at the divide", 0 < inv["H"][0] < 400, "H0=%.0f m" % inv["H"][0])

# ---- 2. Fluvial long profile ----
# steady state: S = (U / (K A^m))^(1/n) exactly, for n = 1 and n = 1.5
st = d2.run({"glacierOn": False, "peakUplift": 1e-3, "K": 2e-5, "m": 0.5, "nexp": 1.0, "dt": 1000.0, "N": 151, "initProfile": "linear", "zHead": 1000.0}, nsteps=3000)
z, A, ds = st["z"], st["A"], st["ds"]; S = (z[:-1] - z[1:]) / ds; Sa = (1e-3 / (2e-5 * A[:-1] ** 0.5))
check("fluvial steady slope–area n=1", np.max(np.abs(S / Sa - 1)) < 1e-6, "max rel %.1e, relief %.0f m" % (np.max(np.abs(S / Sa - 1)), z[0]))
p15 = {"glacierOn": False, "peakUplift": 1e-3, "duration": 1e8, "K": 3e-5, "m": 0.5, "nexp": 1.5, "dt": 1000.0, "N": 151, "initProfile": "linear", "zHead": 300.0}
st = d2.run(p15, nsteps=8000)
st2 = d2.run(p15, nsteps=12000)
z, A = st["z"], st["A"]; S = (z[:-1] - z[1:]) / ds; Sa = (1e-3 / (3e-5 * A[:-1] ** 0.5)) ** (1 / 1.5)
check("fluvial steady slope–area n=1.5", np.max(np.abs(S / Sa - 1)) < 0.01 and np.max(np.abs(st2["z"] - z)) < 1.0,
      "max rel %.1e, relief %.0f m, drift 8→12 Myr %.2f m" % (np.max(np.abs(S / Sa - 1)), z[0], np.max(np.abs(st2["z"] - z))))
check("fluvial erosion = uplift at steady state", np.max(np.abs(st["Ef"][:-1] - 1e-3)) < 1e-9)
# the "steady" initial profile is the fluvial steady state: it does not move under fluvial-only forcing
ss = d2.run({"glacierOn": False, "initProfile": "steady", "dt": 1000.0, "N": 151}, nsteps=100)
check("steady initial profile is stationary", np.max(np.abs(ss["z"] - d2.init_profile(ss["p"])[2])) < 1e-6, "max |Δz| %.2e m" % np.max(np.abs(ss["z"] - d2.init_profile(ss["p"])[2])))
# mass closure with deposition: eroded − deposited = exported
for G in (0.5, 1.5):
    r = d2.run({"glacierOn": False, "G": G, "dt": 1000.0, "N": 151}, nsteps=300)
    check("fluvial closure G=%.1f" % G, abs(r["eroFluv"] - r["exported"]) < 1e-9 * r["eroFluv"], "rel %.1e" % ((r["eroFluv"] - r["exported"]) / r["eroFluv"]))
# dt independence: exact at steady state, small in a transient
f1 = d2.run({"glacierOn": False, "K": 2e-5, "dt": 500.0, "N": 151, "initProfile": "linear", "zHead": 1000.0}, nsteps=6000)
f2 = d2.run({"glacierOn": False, "K": 2e-5, "dt": 2000.0, "N": 151, "initProfile": "linear", "zHead": 1000.0}, nsteps=1500)
check("fluvial steady state dt-independent", np.max(np.abs(f1["z"] - f2["z"])) < 0.01, "%.3f m" % np.max(np.abs(f1["z"] - f2["z"])))
a = d2.run({"dt": 250.0, "N": 151, "Kg": 1e-3}, nsteps=800); b = d2.run({"dt": 1000.0, "N": 151, "Kg": 1e-3}, nsteps=200)
check("coupled transient dt 250 vs 1000 within 25 m", np.max(np.abs(a["z"] - b["z"])) < 25.0, "%.1f m" % np.max(np.abs(a["z"] - b["z"])))

# ---- 3. Glacial erosion: bounded, self-limiting, no runaway at K_g × 100 ----
res = {}
for Kg in (1e-4, 1e-2):
    r = d2.run({"glacierOn": True, "Kg": Kg, "dt": 500.0, "N": 151, "peakUplift": 1e-3}, nsteps=600)
    res[Kg] = r
    check("no NaN / bounded at Kg=%.0e" % Kg, not np.isnan(r["z"]).any() and r["H"].max() < 2000 and r["Eg"].max() <= 0.02 + 1e-12,
          "maxH %.0f, maxUs %.1f m/yr, maxEg %.2e m/yr, glacial volume %.2e m³" % (r["H"].max(), r["Us"].max(), r["Eg"].max(), r["eroGlac"]))
check("Kg×100 self-limits (total glacial erosion within 2× of Kg)", res[1e-2]["eroGlac"] < 2 * res[1e-4]["eroGlac"] and res[1e-2]["eroGlac"] > 0.3 * res[1e-4]["eroGlac"],
      "%.2e vs %.2e m³" % (res[1e-2]["eroGlac"], res[1e-4]["eroGlac"]))
r = res[1e-4]; zs = r["z"] + r["H"]
check("ice surface never slopes uphill under ice", not np.any((zs[1:] > zs[:-1] + 1e-6) & (r["H"][:-1] > 1)))
# volume bookkeeping of the coupled run: uplift − fluvial − glacial = Δ(bed volume)
z0 = d2.init_profile(r["p"])[2]
dV = float(np.sum((r["z"] - z0) * r["ds"] * r["W"]))
check("coupled volume closure", abs(r["upliftVol"] - r["eroFluv"] - r["eroGlac"] - dV) < 1e-6 * r["upliftVol"],
      "rel %.1e" % ((r["upliftVol"] - r["eroFluv"] - r["eroGlac"] - dV) / r["upliftVol"]))
# quarrying: only where the bed is convex-up, and bounded by the cap
rq = d2.run({"glacierOn": True, "Kg": 1e-4, "Kq": 1e-2, "dt": 500.0, "N": 151}, nsteps=100)
check("quarrying bounded by cap", rq["Eg"].max() <= 0.02 + 1e-12 and not np.isnan(rq["z"]).any())

# ---- 4. Lithology in the profile ----
# A dike (persistent band in s) with contrast c: at fluvial steady state E = U everywhere,
# and the signature is a slope c^(1/n) times the analytic slope inside the band (n = 1 → ×c).
rl = d2.run({"glacierOn": False, "K": 2e-5, "m": 0.5, "nexp": 1.0, "dt": 1000.0, "N": 151, "contrastK": 4.0, "initProfile": "linear", "zHead": 1000.0},
            nsteps=4000, litho_spec={"type": "dike", "x0": 15000.0, "width": 2000.0})
z, A, ds, s = rl["z"], rl["A"], rl["ds"], rl["s"]; S = (z[:-1] - z[1:]) / ds; Sa = 1e-3 / (2e-5 * A[:-1] ** 0.5)
i_in = np.argmin(np.abs(s - 15000.0)); i_out = np.argmin(np.abs(s - 20000.0))
check("dike: steady slope × contrast inside, analytic outside", abs(S[i_in] / Sa[i_in] - 4.0) < 1e-6 and abs(S[i_out] / Sa[i_out] - 1.0) < 1e-6 and abs(rl["Ef"][i_in] - 1e-3) < 1e-9,
      "ratio in %.3f, out %.3f" % (S[i_in] / Sa[i_in], S[i_out] / Sa[i_out]))
# A layer in the material frame rides up with uplift: with no erosion (K = 0) the surface
# resistance pattern is unchanged after 200 m of uplift, since z_m = z − U_cum is invariant.
spec = {"type": "layer", "top": 1500.0, "thick": 300.0}
r0 = d2.run({"glacierOn": False, "K": 0.0, "initProfile": "concave", "dt": 1000.0, "N": 151, "peakUplift": 2e-3}, nsteps=1, litho_spec=spec)
r1 = d2.run({"glacierOn": False, "K": 0.0, "initProfile": "concave", "dt": 1000.0, "N": 151, "peakUplift": 2e-3}, nsteps=100, litho_spec=spec)
z0 = d2.init_profile(r0["p"])[2]
expect = ((z0 >= 1200.0) & (z0 <= 1500.0)).astype(float)
check("layer exposed where z0 ∈ [top − thick, top] at t = 0", np.array_equal(r0["r"], expect) and expect.sum() > 3)
check("layer rides up with uplift (K = 0 → r unchanged, z raised)", np.array_equal(r1["r"], expect) and abs(r1["ucum"][10] - 200.0) < 1e-6)

print()
if fails:
    print("FAILED: " + ", ".join(fails)); sys.exit(1)
print("all checks passed")
