"""
test_d1_node.py — the isolated node column (d1 ↔ d2_along link).
Run: /opt/anaconda3/bin/python3 py/test_d1_node.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import gs_core as gc, d2_along as d2
from d1_node import make_node_column

fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  [" + detail + "]" if detail else ""))
    if not ok: fails.append(name)

def spec_from(st, p, i, litho=None):
    s = st["s"]; W = st["W"]
    return {"s": float(s[i]), "z0": float(st["z"][i]), "ucum0": float(st["ucum"][i]), "zDown": float(st["z"][i + 1]),
            "HDown": float(st["H"][i + 1]), "Q": float(st["q"][i] * W[i]) if "q" in st else 0.0, "A": float(st["A"][i]),
            "W": float(W[i]), "Wf": float(0.5 * (W[i] + W[i + 1])), "ds": float(st["ds"]),
            "fac": gc.spatial_factor(float(s[i]) / p["L"], 0.5, p["pattern"], p["ramp"]), "t0": float(st["t"]), "dt": p["dt"],
            "params": p, "litho": litho}

# ---- 0. the shared face law reproduces the d2 march exactly ----
N = 151; L = 30000.0; s = np.linspace(0, L, N); ds = s[1] - s[0]
z = 2500 * (1 - s / L) ** 1.5; z[-1] = 0; W = d2.valley_width(s)
inv = d2.steady_ice(z, W, ds, {"ela": 1800.0, "balGrad": 0.007, "accMax": 2.0, "abMax": 8.0, "fs": 0.5, "iceIters": 200})
Gam = gc.ICE["GAMMA"] / 0.5; H = inv["H"]; Q = inv["Q"]
ok = True
for i in range(N - 2, -1, -1):
    if Q[i] > 0:
        Hi = gc.face_thickness(z[i], z[i + 1] + H[i + 1], H[i + 1], Q[i] / (0.5 * (W[i] + W[i + 1])), ds, Gam, 3)
        ok = ok and abs(Hi - H[i]) < 0.03   # outer iteration tolerance is 0.01 m, bisection 0.001 m
check("face_thickness reproduces the march node by node", ok)

# ---- 1. Consistency: with neighbours refreshed from the true 2D run every step, the node
#         column reproduces the 2D node exactly (same laws, same order of operations). ----
p = dict(d2.DEFAULTS); p.update({"N": 151, "dt": 500.0})
st = d2.make_state(p)
U = gc.make_series({"shape": p["shape"], "peak": p["peakUplift"], "duration": p["duration"]})
ELA = gc.make_series({"shape": "sine", "peak": -p["elaAmp"], "base": p["elaBase"], "period": p["elaPeriod"]})
st = d2.make_state(p); d2.prime_ice(st, p, ELA)
ice_nodes = np.where(st["H"] > 1)[0]
NODES = (int(ice_nodes[len(ice_nodes) // 2]), int(ice_nodes[-1]) - 1, int(ice_nodes[-1]) + 30)   # under ice, near the terminus, fluvial
for i_node in NODES:
    st = d2.make_state(p); d2.prime_ice(st, p, ELA); d2.step(st, p, U, ELA, None)
    col = make_node_column(spec_from(st, p, i_node))
    maxdiff = 0.0
    for k in range(30):
        # refresh the frozen adjacency from the 2D state BEFORE this step (what the node "sees")
        col["frozen"]["zDown"] = float(st["z"][i_node + 1]); col["frozen"]["HDown"] = float(st["H"][i_node + 1])
        # the 2D step computes ice on the uplifted bed; mirror that by using the flux the 2D step will use:
        d2.step(st, p, U, ELA, None)
        col["frozen"]["Q"] = float(st["q"][i_node] * st["W"][i_node])
        # what the 2D step used for this node's receiver: under ice, the march saw the downstream bed after
        # uplift but before erosion; the fluvial sweep saw the receiver AFTER its own update (Gauss–Seidel).
        if col["frozen"]["Q"] > 0:
            col["frozen"]["zDown"] = float(st["z"][i_node + 1] + (st["Eg"][i_node + 1] + st["Ef"][i_node + 1]) * p["dt"])
        else:
            col["frozen"]["zDown"] = float(st["z"][i_node + 1])
        col["frozen"]["HDown"] = float(st["H"][i_node + 1])
        col["step"]()
        maxdiff = max(maxdiff, abs(col["state"]["z"] - st["z"][i_node]))
    check("node column tracks 2D node %d when adjacency is refreshed" % i_node, maxdiff < 0.05, "max |Δz| %.3f m (H=%.0f, Ef=%.1e)" % (maxdiff, st["H"][i_node], st["Ef"][i_node]))

# ---- 2. Isolated fluvial column relaxes to z_down + ds (U f / (K A^m))^(1/n) ----
st = d2.run({"glacierOn": False, "N": 151, "dt": 1000.0}, nsteps=1)
p2 = st["p"]; i = 100
col = make_node_column(spec_from(st, p2, i))
for k in range(3000): col["step"]()
z_eq = col["frozen"]["zDown"] + col["frozen"]["ds"] * (p2["peakUplift"] * col["frozen"]["fac"] / (p2["K"] * col["frozen"]["A"] ** p2["m"])) ** (1 / p2["nexp"])
check("isolated fluvial column → z_down + ds·(U/KA^m)^(1/n)", abs(col["state"]["z"] - z_eq) < 1e-3, "z %.3f vs %.3f" % (col["state"]["z"], z_eq))
check("isolated fluvial column erosion = uplift at equilibrium", abs(col["state"]["Ef"] - p2["peakUplift"]) < 1e-9)

# ---- 3. Isolated glacial column: self-limiting overdeepening against a frozen downstream surface ----
st = d2.run({"N": 151, "dt": 500.0}, nsteps=2)
p3 = st["p"]; i = NODES[0]
col = make_node_column(spec_from(st, p3, i))
H0 = col["frozen"]  # noqa
zs = []
for k in range(400): col["step"](); zs.append(col["state"]["z"])
cs = col["state"]
check("isolated glacial column: ice thickens as bed lowers (H > H_2D)", cs["H"] > st["H"][i] and cs["Us"] < st["Us"][i], "H %.0f→%.0f, Us %.1f→%.1f" % (st["H"][i], cs["H"], st["Us"][i], cs["Us"]))
z_at_solve = cs["z"] + cs["Eg"] * p3["dt"]   # H was solved before this step's erosion lowered the bed
qf = col["frozen"]["Q"] / col["frozen"]["Wf"]
resid = Gam * (0.5 * (cs["H"] + col["frozen"]["HDown"])) ** 5 * ((z_at_solve + cs["H"] - col["frozen"]["zDown"] - col["frozen"]["HDown"]) / col["frozen"]["ds"]) ** 3 - qf
check("isolated glacial column: face law holds", abs(resid) < 1e-3 * qf, "rel resid %.1e" % (resid / qf))
check("isolated glacial column bookkeeping", abs(cs["z"] - (col["frozen"] and (spec_from(st, p3, i)["z0"]) + (cs["ucum"] - spec_from(st, p3, i)["ucum0"]) - cs["ecum"])) < 1e-6)

# ---- 4. Divergence from the 2D node grows in time (the emergence measure is non-trivial) ----
st = d2.run({"N": 151, "dt": 500.0}, nsteps=2); p4 = st["p"]; i = NODES[1]
col = make_node_column(spec_from(st, p4, i))
U = gc.make_series({"shape": p4["shape"], "peak": p4["peakUplift"], "duration": p4["duration"]})
div = []
for k in range(200):
    d2.step(st, p4, U, ELA, None); col["step"]()
    div.append(abs(col["state"]["z"] - st["z"][i]))
check("isolated vs 2D diverge monotonically-ish (final > early)", div[-1] > 5 * div[10] and div[-1] > 1.0, "divergence %.1f m after %.0f kyr" % (div[-1], 200 * 500 / 1e3))

print()
if fails:
    print("FAILED: " + ", ".join(fails)); sys.exit(1)
print("all checks passed")
