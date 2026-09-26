"""
test_d1_columns.py — validation of the two-column A/B scheme (d1 step 2).
Run: /opt/anaconda3/bin/python3 py/test_d1_columns.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import gs_core as gc
from d1_reservoir import run_reservoir
from d1_columns import run_columns

fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  [" + detail + "]" if detail else ""))
    if not ok: fails.append(name)

# 1. No lithology, no glacier: both columns reproduce the reservoir model exactly.
r = run_reservoir({"steps": 2000})
c = run_columns({"steps": 2000, "bodyThick": 0.0, "glacierOn": False})
zA = np.array([s["zA"] for s in c["series"]]); zB = np.array([s["zB"] for s in c["series"]])
R = np.array([s["relief"] for s in r["series"]])
check("columns == reservoir when uniform", np.max(np.abs(zA - R)) < 1e-9 and np.max(np.abs(zB - R)) < 1e-9, "max diff %.2e" % np.max(np.abs(zA - R)))

# 2. Bookkeeping: z = z0 + U_cum − E_cum for each column at every step.
for col in ("A", "B"):
    resid = max(abs(s["z" + col] - (c["params"]["initialRelief"] + s["ucum"] - s["ecum" + col])) for s in c["series"])
    check("bookkeeping column %s" % col, resid < 1e-6, "%.2e m" % resid)

# 3. Persistence timescale. A at steady state (z0 = relief limit, constant U so E_A = U).
#    Body at the surface, thickness T, contrast c: E_B = U/c → stripped at t* = T c / U,
#    and z_B − z_A = U (1 − 1/c) t* = T (c − 1) at that moment; frozen afterwards.
U, T, cc = 3e-3, 600.0, 4.0
p = {"shape": "constant", "peakUplift": U, "reliefLimit": 3000.0, "initialRelief": 3000.0, "duration": 4e6,
     "tailFactor": 1.0, "steps": 4000, "bodyDepth": 0.0, "bodyThick": T, "bodySoft": 1.0, "contrast": cc}
c = run_columns(p)
t = np.array([s["t"] for s in c["series"]]); diff = np.array([s["diff"] for s in c["series"]]); rB = np.array([s["rB"] for s in c["series"]])
EB = np.array([s["EB"] for s in c["series"]])
t_star = T * cc / U
i_exp = np.argmax(rB > 0.5)                       # body top sits at the initial surface: exposed within one step
i_strip = i_exp + np.argmax(rB[i_exp:] < 0.5)
check("body exposed at once", i_exp <= 2, "step %d" % i_exp)
check("body stripped at T c / U", abs(t[i_strip] - t_star) < 3 * (t[1] - t[0]), "t=%.0f kyr vs %.0f kyr" % (t[i_strip] / 1e3, t_star / 1e3))
i_mid = (i_exp + i_strip) // 2
check("E_B = U/c while resistant", abs(EB[i_mid] - U / cc) < 1e-9, "%.3e" % EB[i_mid])
check("upland height at stripping = T (c − 1)", abs(diff[i_strip] - T * (cc - 1)) < 5.0, "%.1f m vs %.1f m" % (diff[i_strip], T * (cc - 1)))
check("upland height frozen after stripping", abs(diff[-1] - diff[i_strip]) < 5.0, "%.1f m → %.1f m" % (diff[i_strip], diff[-1]))
check("A stays at relief limit", abs(c["series"][-1]["zA"] - 3000.0) < 1e-6)

# 4. Buried body: exposure needs cumulative erosion = depth, independent of U.
for U2 in (1e-3, 4e-3):
    c2 = run_columns({"shape": "constant", "peakUplift": U2, "reliefLimit": 3000.0, "initialRelief": 3000.0, "duration": 4e6,
                      "tailFactor": 1.0, "steps": 4000, "bodyDepth": 900.0, "bodyThick": 300.0, "bodySoft": 1.0, "contrast": 3.0})
    tt = np.array([s["t"] for s in c2["series"]]); rr = np.array([s["rB"] for s in c2["series"]])
    i_exp = np.argmax(rr > 0.5)
    check("exposure at depth/U (U=%g)" % U2, abs(tt[i_exp] - 900.0 / U2) < 2 * (tt[1] - tt[0]), "%.0f kyr vs %.0f kyr" % (tt[i_exp] / 1e3, 900.0 / U2 / 1e3))

# 5. Buzzsaw caps elevation near the ELA when the glacial rate dominates.
c3 = run_columns({"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 6000.0, "duration": 10e6, "tailFactor": 1.0, "steps": 4000,
                  "bodyThick": 0.0, "glacierOn": True, "elaBase": 2000.0, "buzzPeak": 2e-2, "buzzWidth": 300.0})
zmax = max(s["zA"] for s in c3["series"])
# steady cap where E_gl(z) + E_sub(z) = U on the lower flank of the Gaussian: z_cap = ELA − w·sqrt(2 ln(P/(U − E_sub)))
Esub = gc.reservoir_erosion(zmax, c3["k"], 4.0)
z_cap = 2000.0 - 300.0 * np.sqrt(2 * np.log(2e-2 / (3e-3 - Esub)))
check("buzzsaw caps column at analytic level", abs(zmax - z_cap) < 20.0, "z_max=%.0f m vs %.0f m" % (zmax, z_cap))
c3b = run_columns({"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 6000.0, "duration": 10e6, "tailFactor": 1.0, "steps": 4000,
                   "bodyThick": 0.0, "glacierOn": False})
check("without glacier column rises above ELA", max(s["zA"] for s in c3b["series"]) > 4000.0)
# Resistant body under ice resists the buzzsaw by contrastG: B ends higher than A.
c3c = run_columns({"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 6000.0, "duration": 10e6, "tailFactor": 1.0, "steps": 4000,
                   "bodyDepth": 0.0, "bodyThick": 5000.0, "contrast": 1.0, "contrastG": 20.0,
                   "glacierOn": True, "elaBase": 2000.0, "buzzPeak": 2e-2, "buzzWidth": 300.0})
check("glacially resistant column stands above ELA cap", c3c["series"][-1]["zB"] > c3c["series"][-1]["zA"] + 200.0,
      "zB−zA=%.0f m" % (c3c["series"][-1]["zB"] - c3c["series"][-1]["zA"]))

# 6. ELA cycle: ELA series dips by elaAmp with the right period.
c4 = run_columns({"steps": 4000, "glacierOn": True, "elaBase": 2000.0, "elaAmp": 500.0, "elaPeriod": 1e5, "duration": 1e6, "tailFactor": 1.0})
ela = np.array([s["ela"] for s in c4["series"]])
check("ELA cycle range", abs(ela.max() - 2000.0) < 1e-6 and abs(ela.min() - 1500.0) < 1.0, "%.0f..%.0f" % (ela.min(), ela.max()))

# 7. dt-independence across the layer discontinuity (soft edge 10 m).
a = run_columns(dict(p, steps=500, bodySoft=10.0)); b = run_columns(dict(p, steps=8000, bodySoft=10.0))
ta = [s["t"] for s in a["series"]]; tb = [s["t"] for s in b["series"]]
da = np.interp(tb, ta, [s["diff"] for s in a["series"]]); db = np.array([s["diff"] for s in b["series"]])
check("dt-independent (500 vs 8000 steps)", np.max(np.abs(da - db)) < 10.0, "max diff %.2f m" % np.max(np.abs(da - db)))

# 8. Coupling "own": while resistant, B is its own reservoir with k/c → equilibrium c^(1/n)·R_lim,
#    and the body is consumed at ≈U (lifetime ≈ T/U, the "no persistence" alternative).
p8 = {"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 3000.0, "initialRelief": 3000.0, "duration": 6e6, "tailFactor": 1.0,
      "steps": 6000, "bodyDepth": 0.0, "bodyThick": 60000.0, "bodySoft": 1.0, "contrast": 16.0, "coupling": "own"}
c8 = run_columns(p8)
check("coupling own: B equilibrates at c^(1/n) R_lim", abs(c8["series"][-1]["zB"] - 2 * 3000.0) < 5.0, "z_B=%.0f m" % c8["series"][-1]["zB"])
check("coupling own: B erodes at U once there", abs(c8["series"][-1]["EB"] - 3e-3) < 1e-6)
c8b = run_columns(dict(p8, bodyThick=600.0))
tt = np.array([s["t"] for s in c8b["series"]]); rr = np.array([s["rB"] for s in c8b["series"]])
i_e = np.argmax(rr > 0.5); i_s = i_e + np.argmax(rr[i_e:] < 0.5)
life = tt[i_s] - tt[i_e]
check("coupling own: body lifetime between T/U and T c/U", 600.0 / 3e-3 < life < 600.0 * 16 / 3e-3, "%.0f kyr (T/U = %.0f, T c/U = %.0f)" % (life / 1e3, 600.0 / 3e-3 / 1e3, 600.0 * 16 / 3e-3 / 1e3))
# "max" equals "landscape" while z_B ≤ z_A, and bounds B afterwards
c9a = run_columns(dict(p, coupling="max")); c9b = run_columns(dict(p, coupling="landscape"))
d9 = np.array([s["diff"] for s in c9a["series"]])
check("coupling max: bounded upland (< landscape freeze)", d9[-1] < c9b["series"][-1]["diff"] and d9.max() > 100.0, "max %.0f m" % d9.max())
# defaults unchanged by the new switches
c10 = run_columns({"steps": 500}); c10b = run_columns({"steps": 500, "coupling": "landscape", "buzzAbove": 0.0})
check("defaults unchanged", all(abs(a["zB"] - b["zB"]) < 1e-12 for a, b in zip(c10["series"], c10b["series"])))
# 9. Above-ELA floor: far above the ELA the glacial rate is φ·P; below, unchanged.
c11 = run_columns({"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 20000.0, "initialRelief": 8000.0, "duration": 1e3, "tailFactor": 1.0,
                   "steps": 2, "bodyThick": 0.0, "glacierOn": True, "elaBase": 2000.0, "buzzPeak": 2e-3, "buzzWidth": 300.0, "buzzAbove": 0.4})
check("buzzsaw floor above ELA = φP", abs(c11["series"][0]["EgA"] - 0.4 * 2e-3) < 1e-12, "%.2e" % c11["series"][0]["EgA"])
c12 = run_columns({"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 20000.0, "initialRelief": 1700.0, "duration": 1e3, "tailFactor": 1.0,
                   "steps": 2, "bodyThick": 0.0, "glacierOn": True, "elaBase": 2000.0, "buzzPeak": 2e-3, "buzzWidth": 300.0, "buzzAbove": 0.4})
check("buzzsaw below ELA unaffected by φ", abs(c12["series"][0]["EgA"] - gc.buzzsaw_rate(1700.0, 2000.0, 300.0, 2e-3)) < 1e-15)
# with a floor, a resistant column cannot escape: it is held where φP/c_g + E_sub = U or stays capped
c13 = run_columns({"shape": "constant", "peakUplift": 3e-3, "reliefLimit": 6000.0, "duration": 10e6, "tailFactor": 1.0, "steps": 4000,
                   "bodyDepth": 0.0, "bodyThick": 50000.0, "contrast": 1.0, "contrastG": 4.0,
                   "glacierOn": True, "elaBase": 2000.0, "buzzPeak": 2e-2, "buzzWidth": 300.0, "buzzAbove": 1.0})
check("floor φ=1 stops escape above the window", c13["series"][-1]["zB"] < 2000.0 + 3 * 300.0, "z_B=%.0f m" % c13["series"][-1]["zB"])

print()
if fails:
    print("FAILED: " + ", ".join(fails)); sys.exit(1)
print("all checks passed")
