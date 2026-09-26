"""
d1_columns.py — Python mirror of d1/model.js runColumns(): two lockstep columns
A (weak) and B (with a resistant body) under one uplift U(t) and one ELA(t).

Per column, elevation above base level z evolves as
    dz/dt = U(t) − E_sub − E_gl
  E_sub  subaerial lowering. Column A is the landscape's relief reservoir,
         E_ref = k z_A^n (dR/dt = U − k Rⁿ, legacy uplift_sandbox). Column B's demand
         D_B is an asserted coupling (p["coupling"]):
           "landscape"  D_B = k z_A^n          B sits in A's landscape (low-relief plateau)
           "own"        D_B = k z_B^n          B is its own relief reservoir
           "max"        D_B = k max(z_A, z_B)^n
         and its rock resists it: E_sub,B = D_B · f(r_B, c), f = 1/(1 + r(c−1)).
  E_gl   elevation-dependent glacial rate ("buzzsaw"): a Gaussian window of
         half-width w round the ELA, with an optional floor φ above the ELA (deep ice
         still erodes high ground): E_gl = P·g(z)·f(r, c_g),
           g = exp(−½u²) for z ≤ ELA;  g = φ + (1−φ)·exp(−½u²) for z > ELA;  u = (z−ELA)/w.
  r_B    resistance at B's surface, from a horizontal body in the material frame
         zm = z − U_cum: top at z0 − depth, thickness T, soft edges (core/lithology).
RK4 on the vector (z_A, z_B, U_cum, E_cum,A, E_cum,B). Units m, yr.
"""
import math
import gs_core as gc

DEFAULTS = {
    "duration": 12e6, "tailFactor": 2.0, "peakUplift": 3e-3, "shape": "bell",
    "exponent": 4.0, "reliefLimit": 3000.0, "initialRelief": 0.0,
    "bodyDepth": 0.0, "bodyThick": 500.0, "bodySoft": 10.0, "contrast": 5.0, "contrastG": 5.0,
    "glacierOn": False, "elaBase": 2000.0, "elaAmp": 0.0, "elaPeriod": 1e5,
    "buzzPeak": 1e-3, "buzzWidth": 300.0, "buzzAbove": 0.0,
    "coupling": "landscape",
    "steps": 2000,
}

def run_columns(params=None):
    p = dict(DEFAULTS); p.update(params or {})
    U = gc.make_series({"shape": p["shape"], "peak": p["peakUplift"], "duration": p["duration"]})
    ela_cycle = gc.make_series({"shape": "sine", "peak": -p["elaAmp"], "base": p["elaBase"], "period": p["elaPeriod"]})
    k = gc.reservoir_k(p["peakUplift"], p["reliefLimit"], p["exponent"])
    n = p["exponent"]; z0 = p["initialRelief"]
    body = gc.litho_make({"type": "layer", "top": z0 - p["bodyDepth"], "thick": p["bodyThick"], "soft": p["bodySoft"]})
    gl = p["glacierOn"]; P = p["buzzPeak"]; w = p["buzzWidth"]; phi = p["buzzAbove"]; c = p["contrast"]; cg = p["contrastG"]
    coupling = p["coupling"]

    def buzz(z, ela):
        e = gc.buzzsaw_rate(z, ela, w, P)
        return e if z <= ela else phi * P + (1.0 - phi) * e

    def rates(t, y):
        zA, zB, ucum = y[0], y[1], y[2]
        u = U(t); ela = ela_cycle(t)
        Eref = gc.reservoir_erosion(zA, k, n)
        if coupling == "own":
            DB = gc.reservoir_erosion(zB, k, n)
        elif coupling == "max":
            DB = gc.reservoir_erosion(max(zA, zB), k, n)
        else:
            DB = Eref
        rB = body(0.0, 0.0, gc.material_z(zB, ucum))
        EgA = buzz(zA, ela) if gl else 0.0
        EgB = (buzz(zB, ela) * gc.erodibility_factor(rB, cg)) if gl else 0.0
        EsubB = DB * gc.erodibility_factor(rB, c)
        EA = Eref + EgA; EB = EsubB + EgB
        return [u - EA, u - EB, u, EA, EB], (u, ela, Eref, EsubB, EgA, EgB, rB)

    sim = p["duration"] * p["tailFactor"]; dt = sim / p["steps"]
    y = [z0, z0, 0.0, 0.0, 0.0]
    series = []
    for i in range(p["steps"] + 1):
        t = i * dt
        d, aux = rates(t, y)
        u, ela, Eref, EsubB, EgA, EgB, rB = aux
        top = z0 - p["bodyDepth"] + y[2]
        series.append({"t": t, "uplift": u, "ela": ela, "zA": y[0], "zB": y[1], "ucum": y[2],
                       "ecumA": y[3], "ecumB": y[4], "EA": Eref + EgA, "EB": EsubB + EgB,
                       "EsubA": Eref, "EsubB": EsubB, "EgA": EgA, "EgB": EgB, "rB": rB,
                       "bodyTop": top, "bodyBot": top - p["bodyThick"], "diff": y[1] - y[0]})
        # RK4 on the vector
        k1 = d
        k2, _ = rates(t + dt / 2, [y[j] + dt / 2 * k1[j] for j in range(5)])
        k3, _ = rates(t + dt / 2, [y[j] + dt / 2 * k2[j] for j in range(5)])
        k4, _ = rates(t + dt, [y[j] + dt * k3[j] for j in range(5)])
        y = [y[j] + dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(5)]
        if y[0] < 0: y[0] = 0.0
        if y[1] < 0: y[1] = 0.0
    return {"params": p, "k": k, "series": series}
