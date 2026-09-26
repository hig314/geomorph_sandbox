"""
d1_reservoir.py — Python mirror of d1/model.js runReservoir().

Relief reservoir dR/dt = U(t) − k Rⁿ integrated with RK4, in metres and years.
"""
import gs_core as gc

def run_reservoir(params=None):
    p = {"duration": 12e6, "tailFactor": 2.0, "peakUplift": 3e-3, "exponent": 4.0,
         "reliefLimit": 3000.0, "initialRelief": 0.0, "steps": 2000, "shape": "bell"}
    p.update(params or {})
    U = gc.make_series({"shape": p["shape"], "peak": p["peakUplift"], "duration": p["duration"]})
    k = gc.reservoir_k(p["peakUplift"], p["reliefLimit"], p["exponent"])
    n = p["exponent"]
    sim_time = p["duration"] * p["tailFactor"]
    dt = sim_time / p["steps"]
    dRdt = lambda t, R: U(t) - gc.reservoir_erosion(R, k, n)
    series = []
    R = p["initialRelief"]
    for i in range(p["steps"] + 1):
        t = i * dt
        u = U(t); E = gc.reservoir_erosion(R, k, n)
        series.append({"t": t, "uplift": u, "relief": R, "erosion": E, "net": u - E,
                       "reliefEq": gc.reservoir_equilibrium(u, k, n)})
        R = gc.rk4(dRdt, t, R, dt)
        if R < 0:
            R = 0.0
    return {"params": p, "k": k, "series": series}
