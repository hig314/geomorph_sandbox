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
thick = H > 10
check("sliding speed = fs q/H (H > Hmin)", np.allclose(inv["Us"][thick], 0.5 * inv["q"][thick] / H[thick]))
check("sliding regularised below Hmin", np.all(inv["Us"][(H > 0) & (H <= 10)] <= 0.5 * inv["q"][(H > 0) & (H <= 10)] / 10 + 1e-9))
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
ss = d2.run({"glacierOn": False, "initProfile": "steady", "K": 2e-5, "dt": 1000.0, "N": 151}, nsteps=100)   # K high enough that no node hits the slope cap
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
check("coupled volume closure", abs(r["upliftVol"] - r["eroFluv"] - r["eroGlac"] - r["eroRock"] - dV) < 1e-6 * r["upliftVol"],
      "rel %.1e" % ((r["upliftVol"] - r["eroFluv"] - r["eroGlac"] - r["eroRock"] - dV) / r["upliftVol"]))
# quarrying: only where the bed is convex-up, and bounded by the cap
rq = d2.run({"glacierOn": True, "Kg": 1e-4, "Kq": 1e-2, "dt": 500.0, "N": 151}, nsteps=100)
check("quarrying bounded by cap", rq["Eg"].max() <= 0.02 + 1e-12 and not np.isnan(rq["z"]).any())

# ---- 3b. Cell-scale regularisations ----
# (a) glacier reaching the outlet: free-outflow boundary → no thickness/sliding spike at the last cells
lo = d2.run({"elaBase": 300.0, "dt": 500.0, "N": 151}, nsteps=20)
Hl, Ul = lo["H"], lo["Us"]
check("glacier reaches the outlet", Hl[-1] > 1 and Hl[-2] > 1, "H outlet %.0f, H[-2] %.0f, H[-3] %.0f" % (Hl[-1], Hl[-2], Hl[-3]))
check("no outlet spike: H and Us smooth over the last 4 cells", np.max(Hl[-4:]) < 1.3 * np.min(Hl[-4:]) and np.max(Ul[-4:]) < 1.3 * np.min(Ul[-4:]),
      "Us last 4: %s" % ", ".join("%.1f" % v for v in Ul[-4:]))
# (b) checkerboard: the reported scenario (extreme K, low ELA, strong erosion). Interior bed roughness
#     (second difference, 3 nodes off each end — the fixed outlet makes a legitimate step) with the
#     [¼ ½ ¼] footprint vs the raw law, which grows a cell-scale sawtooth.
def checker(r):
    """Checkerboard signature: fraction of interior ice-covered nodes where the bed's second
    difference flips sign from one node to the next, and the largest |Δ²z| there."""
    z = r["z"]; d = z[:-2] - 2 * z[1:-1] + z[2:]; m = r["H"][1:-1] > 1; m[:6] = False; m[-3:] = False
    dd = d[m]
    flips = np.sum(np.sign(dd[1:]) != np.sign(dd[:-1])) / max(len(dd) - 1, 1)
    return flips, (np.max(np.abs(dd)) if len(dd) else 0.0)
cb = {"L": 30000.0, "N": 301, "initProfile": "steady", "peakUplift": 2.8e-3, "K": 10 ** -6.5, "G": 0.6, "hack": 1.4,
      "elaBase": 500.0, "elaAmp": 150.0, "elaPeriod": 1e5, "Kg": 10 ** -3.2, "dt": 500.0}
for sm in (True, False):
    fl, mx = checker(d2.run(dict(cb, eroSmooth=sm), nsteps=400))
    check("no checkerboard in the reported scenario (footprint %s)" % ("on" if sm else "off"), fl < 0.2 and mx < 100.0,
          "sign-flip fraction %.2f, max interior |Δ²z| %.0f m" % (fl, mx))
# the footprint removes the two-cell mode of the erosion field exactly
alt = np.array([1.0 if i % 2 else 0.0 for i in range(20)]); sm_ = alt.copy(); sm_[1:-1] = 0.25 * alt[:-2] + 0.5 * alt[1:-1] + 0.25 * alt[2:]
check("[¼ ½ ¼] footprint annihilates a two-cell mode", np.allclose(sm_[1:-1], 0.5))
r_fp = d2.run(dict(cb, eroSmooth=True), nsteps=100)
check("footprint conserves volume", abs(r_fp["upliftVol"] - r_fp["eroFluv"] - r_fp["eroGlac"] - r_fp["eroRock"] - float(np.sum((r_fp["z"] - d2.init_profile(r_fp["p"])[2]) * r_fp["ds"] * r_fp["W"]))) < 1e-6 * r_fp["upliftVol"])
# (d) headwall: the divide cell gets almost no ice flux and would rise into a single-cell spire;
#     threshold-slope rockfall keeps every ice-free inter-cell slope ≤ Sc, and the cirque floor below stays smooth
hw = {"N": 301, "dt": 500.0, "Kg": 1e-3}
r_sc = d2.run(dict(hw, Sc=0.8), nsteps=600); r_no = d2.run(dict(hw, Sc=1e9), nsteps=600)
S_sc = (r_sc["z"][:-1] - r_sc["z"][1:]) / r_sc["ds"]
check("threshold failure caps every inter-cell slope at Sc", np.max(S_sc) <= 0.8 + 1e-6, "max slope %.2f; head drop %.0f m (was %.0f m without)" % (np.max(S_sc), r_sc["z"][0] - r_sc["z"][1], r_no["z"][0] - r_no["z"][1]))
dz2 = lambda r: np.abs(r["z"][1:9] - 2 * r["z"][2:10] + r["z"][3:11]).max()
check("cirque floor below the headwall stays smooth (|Δ²z| < 15 m)", dz2(r_sc) < 15.0, "|Δ²z| %.0f m with rockfall, %.0f m without" % (dz2(r_sc), dz2(r_no)))
rep = d2.run(dict(cb, Sc=0.8), nsteps=400)
check("reported scenario: head drop capped (ice does not exempt the spire)", rep["z"][0] - rep["z"][1] <= 0.8 * rep["ds"] + 1e-6, "drop %.0f m, rock vol %.2e" % (rep["z"][0] - rep["z"][1], rep["eroRock"]))
check("rockfall volume closes", abs(r_sc["upliftVol"] - r_sc["eroFluv"] - r_sc["eroGlac"] - r_sc["eroRock"] - float(np.sum((r_sc["z"] - d2.init_profile(r_sc["p"])[2]) * r_sc["ds"] * r_sc["W"]))) < 1e-6 * r_sc["upliftVol"] and r_sc["eroRock"] > 0)
# (c) steady profile slope cap
sp = d2.make_state(dict(d2.DEFAULTS, K=1e-7, N=151))
check("steady profile capped at the threshold slope", np.max(-np.diff(sp["z"]) / sp["ds"]) <= 0.8 + 1e-9, "max slope %.3f, head %.0f m" % (np.max(-np.diff(sp["z"]) / sp["ds"]), sp["z"][0]))

# (e) ice-margin handoff: the reported low-Kg scenario with a cycling ELA. With a binary fluvial
#     on/off at the margin, ice thickness zig-zagged cell to cell (5 sign changes over the last
#     ~15 ice cells) and the bed carried a 31 m sawtooth; the efficiency ramp removes the zig-zag.
mg = {"L": 30000.0, "N": 301, "initProfile": "steady", "peakUplift": 1e-3, "K": 1e-7, "m": 0.6, "G": 1.8, "hack": 2.0,
      "W0": 300.0, "kw": 0.2, "elaBase": 1700.0, "elaAmp": 150.0, "elaPeriod": 1e5, "Kg": 10 ** -5.3, "dt": 500.0}
r_mg = d2.run(mg, nsteps=600); Hm = r_mg["H"]; zm = r_mg["z"]
term = np.max(np.where(Hm > 1)[0]); win = slice(max(3, term - 10), term + 6)
dHm = np.diff(Hm[win]); hz = int(np.sum(np.sign(dHm[1:]) != np.sign(dHm[:-1])))
d2m = zm[:-2] - 2 * zm[1:-1] + zm[2:]; d2w = np.abs(d2m[win.start - 1:win.stop - 1]).max()
check("ice margin: thickness has no cell-scale zig-zag", hz <= 1, "%d sign changes in ΔH over the margin (was 5)" % hz)
check("ice margin: bed |Δ²z| < 15 m", d2w < 15.0, "%.1f m (was 31 m)" % d2w)
check("fluvial spike at the terminus ≤ 1.6 U", r_mg["Ef"].max() <= 1.6e-3, "Ef max %.2e (was 2.5e-3)" % r_mg["Ef"].max())
check("margin closure", abs(r_mg["upliftVol"] - r_mg["eroFluv"] - r_mg["eroGlac"] - r_mg["eroRock"] - float(np.sum((zm - d2.init_profile(r_mg["p"])[2]) * r_mg["ds"] * r_mg["W"]))) < 1e-6 * r_mg["upliftVol"])
# the ramp is inert on bare ground and under thick ice (default phiSub = 0)
check("fluvial weight: 1 on bare ground, 0 under thick ice", d2.fluvial_weight(np.array([0.0, 50.0, 100.0, 500.0]), r_mg["p"]).tolist() == [1.0, 0.5, 0.0, 0.0])

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
