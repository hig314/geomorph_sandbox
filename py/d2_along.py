"""
d2_along.py — Python mirror of d2_along/model.js: the along-valley (long profile)
tier, z(s, t), plus the reference flowline SIA it is validated against.

Coordinates: s along the valley, node i at s = i·ds, i = 0 (divide / valley head)
… N−1 (outlet, fixed base level z = 0). Units m, yr. numpy-vectorised, but every
operation mirrors a JS loop one-for-one (same stencils, same clamps).

Asserted (DESIGN.md §2 ladder): drainage area A(s) by Hack's law, valley width W(s).
Computed: ice flux q (∫ balance·W), thickness H (flow-law inversion), sliding U_s,
surface slope, erosion.

Glacier — steady discharge (DESIGN.md §4), not time-stepped SIA:
  1. b(z_s) mass balance on the ice SURFACE z_s = z + H (legacy decision kept)
  2. Q(s) = ∫₀ˢ b W ds, clamped ≥ 0 (ablation exhausts the flux → terminus)
  3. thickness from the flux law, marched UPSTREAM from the outlet: at each face
     Γ' Hf^(n+2) Sf^n = q_f with Γ' = Γ/(1 − f_s), solved for H_i given H_{i+1}
     (monotone → unique root). Pointwise inversion H = (q/ΓSⁿ)^(1/(n+2)) was tried
     first and diverged (slope floor → H explodes); the march is the BVP fallback of
     DESIGN.md §6. Outer iteration: surface → balance → flux → march, under-relaxed.
     f_s = sliding fraction of the depth-averaged speed q/H
  4. U_s = f_s q / H;  E = K_g U_s^l · litho, capped;  optional quarrying K_q U_s · convexity

Reference: reference_sia() — explicit conservative flux-form flowline SIA (the legacy
terrain_sandbox scheme in 1-D: CFL-subcycled, upwind guard, surface-referenced
balance, no diffusivity cap) run to steady state on a fixed bed.
"""
import numpy as np
import gs_core as gc

ICE = gc.ICE

# ---------------------------------------------------------------------------
# asserted geometry
# ---------------------------------------------------------------------------
def hack_area(s, ka=5.7, h=1.667, s0=200.0):
    """Drainage area (m²) at distance s (m) from the divide: A = ka (s + s0)^h.
    ka = 5.7, h = 1.667 reproduce Hack (1957) L = 1.4 A^0.6 (km, km²)."""
    return ka * (s + s0) ** h

def valley_width(s, W0=300.0, kw=0.05):
    """Valley width (m): W = W0 + kw s."""
    return W0 + kw * s

# ---------------------------------------------------------------------------
# glacier: steady-discharge thickness inversion
# ---------------------------------------------------------------------------
def surface_slope(zs, ds, smin):
    """Downstream-positive surface slope, centred; one-sided at the ends; floored."""
    S = np.empty_like(zs)
    S[1:-1] = (zs[:-2] - zs[2:]) / (2 * ds)
    S[0] = (zs[0] - zs[1]) / ds
    S[-1] = (zs[-2] - zs[-1]) / ds
    return np.maximum(S, smin)

def route_flux(b, W, ds):
    """Q(s) = ∫ b W ds downstream, clamped at zero (terminus)."""
    src = b * W * ds
    Q = np.empty_like(b)
    acc = 0.0
    for i in range(len(b)):
        acc = acc + src[i]
        if acc < 0.0:
            acc = 0.0
        Q[i] = acc
    return Q

def march_thickness(z, W, Q, ds, Gam, n, N):
    """Given the volume flux Q_i leaving node i (m³/yr), march upstream from the outlet
    solving the face flux law for H_i given H_{i+1} (core face_thickness: monotone → unique
    root). Nodes whose outgoing flux is zero carry no ice (terminus and ice-free ground)."""
    H = np.zeros(N)
    for i in range(N - 2, -1, -1):
        if Q[i] <= 0.0:
            H[i] = 0.0
            continue
        Wf = 0.5 * (W[i] + W[i + 1])
        H[i] = gc.face_thickness(z[i], z[i + 1] + H[i + 1], H[i + 1], Q[i] / Wf, ds, Gam, n)
    return H

def steady_ice(z, W, ds, p, H0=None):
    """Steady ice on bed z by routed flux + upstream march. Returns dict(H, q, Q, Us, b, S, iters)."""
    n = ICE["N"]; Gam = ICE["GAMMA"] * p.get("flow", 1.0) / (1.0 - p.get("fs", 0.5))
    fs = p.get("fs", 0.5); smin = p.get("smin", ICE["SLOPE_MIN"])
    relax = p.get("relax", 0.5); iters = p.get("iceIters", 30); tol = p.get("iceTol", 0.01)
    N = len(z)
    H = np.zeros_like(z) if H0 is None else H0.copy()
    it = 0
    for it in range(1, iters + 1):
        zs = z + H
        b = np.array([gc.mass_balance(v, p["ela"], p["balGrad"], p.get("accMax"), p.get("abMax")) for v in zs])
        Q = route_flux(b, W, ds)
        Hn = march_thickness(z, W, Q, ds, Gam, n, N)
        dH = Hn - H
        H = H + relax * dH
        if np.max(np.abs(dH)) < tol:
            break
    zs = z + H
    b = np.array([gc.mass_balance(v, p["ela"], p["balGrad"], p.get("accMax"), p.get("abMax")) for v in zs])
    Q = route_flux(b, W, ds); q = Q / W
    S = surface_slope(zs, ds, smin)
    Us = np.where(H > 0, fs * q / np.maximum(H, 1e-9), 0.0)
    return {"H": H, "q": q, "Q": Q, "Us": Us, "b": b, "S": S, "iters": it}

# ---------------------------------------------------------------------------
# reference: explicit flowline SIA to steady state
# ---------------------------------------------------------------------------
def reference_sia(z, W, ds, p, years=20000.0, dt_max=5.0, tol=1e-5, check=200, H0=None):
    """Explicit conservative flux-form SIA on a fixed bed; returns H at (windowed) steady state.
    Sliding enters as q_total = q_def / (1 − f_s), consistent with the inversion."""
    n = ICE["N"]; Gam = ICE["GAMMA"] * p.get("flow", 1.0) / (1.0 - p.get("fs", 0.5))
    N = len(z); H = np.zeros(N) if H0 is None else H0.copy()
    Wf = 0.5 * (W[:-1] + W[1:])
    t = 0.0; volPrev = H.sum(); sub = 0
    while t < years:
        zs = z + H
        Hf = 0.5 * (H[:-1] + H[1:])
        S = (zs[1:] - zs[:-1]) / ds                    # +ve = surface rising downstream
        D = Gam * Hf ** (n + 2) * np.abs(S) ** (n - 1)
        # upwind guard: no flux if the higher-surface cell is ice-free
        up = np.where(S > 0, H[1:], H[:-1])
        D[up <= 1e-3] = 0.0
        qf = -D * S                                     # flux per unit width across face i+1/2
        Dmax = D.max()
        dt = min(dt_max, 0.2 * ds * ds / Dmax) if Dmax > 0 else dt_max
        b = np.array([gc.mass_balance(v, p["ela"], p["balGrad"], p.get("accMax"), p.get("abMax")) for v in zs])
        div = np.zeros(N)
        F = Wf * qf                                     # volume flux across faces
        div[1:-1] = (F[1:] - F[:-1]) / (W[1:-1] * ds)
        div[0] = F[0] / (W[0] * ds)                     # divide: no inflow face
        div[-1] = (0.0 - F[-1]) / (W[-1] * ds)          # outlet: free outflow (ice leaves)
        H = H + dt * (b - div)
        H[H < 0] = 0.0
        H[-1] = 0.0                                     # outlet boundary: ice exits the domain
        t += dt; sub += 1
        if sub % check == 0:
            vol = H.sum()
            if abs(vol - volPrev) <= tol * (vol + 1.0):
                break
            volPrev = vol
    zs = z + H
    Q = np.zeros(N); Q[1:-1] = 0.5 * (F[:-1] + F[1:]); Q[0] = 0.5 * F[0]; Q[-1] = F[-1]
    Us = np.where(H > 0, p.get("fs", 0.5) * (Q / W) / np.maximum(H, 1e-9), 0.0)
    return {"H": H, "Q": Q, "Us": Us, "t": t, "sub": sub}

# ---------------------------------------------------------------------------
# the geomorphic model
# ---------------------------------------------------------------------------
DEFAULTS = {
    "L": 30000.0, "N": 301,
    "initProfile": "steady", "zHead": 2500.0, "noise": 0.0, "seed": 1,
    "shape": "constant", "peakUplift": 1e-3, "duration": 1e7, "pattern": "uniform", "ramp": 0.1,
    "K": 3.2e-6, "m": 0.5, "nexp": 1.0, "G": 0.0, "ka": 5.7, "hack": 1.667, "W0": 300.0, "kw": 0.05,
    "glacierOn": True, "elaBase": 1500.0, "elaAmp": 0.0, "elaPeriod": 1e5, "balGrad": 0.007,
    "accMax": 2.0, "abMax": 8.0, "fs": 0.5, "flow": 1.0, "Kg": 1e-4, "lexp": 1.0, "eroCap": 0.02,
    "Kq": 0.0, "iceIters": 30, "relax": 0.5, "smin": 1e-3, "iceTol": 0.01,
    "litho": None, "contrastK": 5.0, "contrastKg": 5.0,
    "dt": 500.0,
}

def init_profile(p):
    N = p["N"]; L = p["L"]; ds = L / (N - 1)
    s = np.arange(N) * ds
    if p["initProfile"] == "linear":
        z = p["zHead"] * (1 - s / L)
    elif p["initProfile"] == "steady" and p["K"] > 0:
        # fluvial steady state for the peak uplift: S_i = (U f_i / (K A_i^m))^(1/n), integrated from the outlet
        A = hack_area(s, p["ka"], p["hack"])
        z = np.zeros(N)
        for i in range(N - 2, -1, -1):
            f = gc.spatial_factor(s[i] / L, 0.5, p["pattern"], p["ramp"])
            z[i] = z[i + 1] + ds * (p["peakUplift"] * f / (p["K"] * A[i] ** p["m"])) ** (1.0 / p["nexp"])
    else:  # concave: z ∝ (1 − s/L)^1.5 … a fluvial-looking profile
        z = p["zHead"] * (1 - s / L) ** 1.5
    if p["noise"] > 0:
        rng = np.random.default_rng(p["seed"])
        z = z + p["noise"] * rng.standard_normal(N)
    z[-1] = 0.0
    return s, ds, z

def make_state(p):
    s, ds, z = init_profile(p)
    return {"t": 0.0, "s": s, "ds": ds, "z": z, "ucum": np.zeros_like(z), "H": np.zeros_like(z),
            "Us": np.zeros_like(z), "Ef": np.zeros_like(z), "Eg": np.zeros_like(z),
            "A": hack_area(s, p["ka"], p["hack"]), "W": valley_width(s, p["W0"], p["kw"]),
            "exported": 0.0, "eroFluv": 0.0, "eroGlac": 0.0, "upliftVol": 0.0}

def fluvial_step(st, p, iceMask):
    """Implicit Yuan erosion–deposition along the line (receiver = downstream node)."""
    z = st["z"]; A = st["A"]; W = st["W"]; ds = st["ds"]; dt = p["dt"]; N = len(z)
    Kf = p["K"] * st["fK"]
    h0 = z.copy(); G = p["G"]
    # Under-relaxation only while iterating the deposition feedback. With G = 0 the
    # downstream-first sweep IS the exact implicit solve, and relaxing it (as the legacy
    # DEM code did, OMEGA = 0.6 on a single pass) applies only 60 % of the incision:
    # steady slopes came out ~2× the analytic U/(K A^m) and depended on dt.
    OMEGA = 0.6 if G > 0 else 1.0
    maxIter = 1 if G <= 0 else min(25, int(np.ceil(2 + G * 20)))
    Qs = np.zeros(N)
    for it in range(maxIter):
        # (a) sediment flux downstream: Qs_out = Qs_in + (eroded − deposited) volume rate
        Qs[:] = 0.0
        for i in range(N - 1):
            Qs[i + 1] += Qs[i] + (h0[i] - z[i]) * ds * W[i] / dt
        # (b) implicit update from the outlet upstream (receiver already updated)
        maxd = 0.0
        for i in range(N - 2, -1, -1):
            if iceMask[i]:
                continue
            hr = z[i + 1]
            dep = dt * G * Qs[i] / A[i]
            Kp = Kf[i] * A[i] ** p["m"]
            nh = gc.yuan_node(h0[i], z[i], hr, Kp, dt, ds, dep, p["nexp"], False)
            nh = z[i] + OMEGA * (nh - z[i])
            d = abs(nh - z[i]); maxd = max(maxd, d)
            z[i] = nh
        if G <= 0 or maxd < 0.02:
            break
    # export at the outlet: everything that arrives at node N−1 leaves
    Qs[:] = 0.0
    for i in range(N - 1):
        Qs[i + 1] += Qs[i] + (h0[i] - z[i]) * ds * W[i] / dt
    st["exported"] += Qs[N - 1] * dt
    st["Ef"] = (h0 - z) / dt
    st["eroFluv"] += float(np.sum((h0 - z) * ds * W))

def step(st, p, U_of_t, ELA_of_t, litho_fn):
    z = st["z"]; s = st["s"]; ds = st["ds"]; dt = p["dt"]; N = len(z); L = p["L"]
    # lithology at the surface (material frame)
    if litho_fn is not None:
        r = np.array([litho_fn(s[i], 0.0, gc.material_z(z[i], st["ucum"][i])) for i in range(N)])
    else:
        r = np.zeros(N)
    st["r"] = r
    st["fK"] = 1.0 / (1.0 + r * (p["contrastK"] - 1.0))
    st["fKg"] = 1.0 / (1.0 + r * (p["contrastKg"] - 1.0))
    # 1. uplift (outlet fixed)
    u = U_of_t(st["t"])
    fac = np.array([gc.spatial_factor(s[i] / L, 0.5, p["pattern"], p["ramp"]) for i in range(N)])
    du = u * fac * dt; du[-1] = 0.0
    z += du; st["ucum"] += du
    st["upliftVol"] += float(np.sum(du * ds * st["W"]))
    # 2. glacier (steady) → mask + erosion
    if p["glacierOn"]:
        pi = dict(p); pi["ela"] = ELA_of_t(st["t"])
        ice = steady_ice(z, st["W"], ds, pi, st["H"])
        H = ice["H"]; Us = ice["Us"]
        iceMask = H > 1.0
        E = np.array([gc.cap_rate(gc.abrasion(p["Kg"] * st["fKg"][i], Us[i], p["lexp"]), p["eroCap"]) for i in range(N)])
        if p["Kq"] > 0:
            curv = np.zeros(N); curv[1:-1] = (z[:-2] - 2 * z[1:-1] + z[2:]) / (ds * ds)   # >0 convex-up? sign below
            conv = -curv                                                               # convex up: z'' < 0
            E = E + np.array([gc.quarrying(p["Kq"] * st["fKg"][i], Us[i], conv[i]) for i in range(N)])
            E = np.minimum(E, p["eroCap"])
        E[~iceMask] = 0.0; E[-1] = 0.0
        z -= E * dt
        st["H"] = H; st["Us"] = Us; st["Eg"] = E; st["q"] = ice["q"]; st["b"] = ice["b"]
        st["eroGlac"] += float(np.sum(E * dt * ds * st["W"]))
    else:
        iceMask = np.zeros(N, dtype=bool); st["H"][:] = 0.0; st["Us"][:] = 0.0; st["Eg"][:] = 0.0
    # 3. fluvial on ice-free nodes
    fluvial_step(st, p, iceMask)
    st["t"] += dt

def prime_ice(st, p, ELA_of_t):
    """Steady ice on the initial bed, no erosion (mirrors model.primeIce)."""
    if not p["glacierOn"]:
        return
    pi = dict(p); pi["ela"] = ELA_of_t(st["t"])
    ice = steady_ice(st["z"], st["W"], st["ds"], pi, st["H"])
    st["H"] = ice["H"]; st["Us"] = ice["Us"]; st["q"] = ice["q"]; st["b"] = ice["b"]

def run(params=None, nsteps=100, litho_spec=None):
    p = dict(DEFAULTS); p.update(params or {})
    st = make_state(p)
    U = gc.make_series({"shape": p["shape"], "peak": p["peakUplift"], "duration": p["duration"]})
    ELA = gc.make_series({"shape": "sine", "peak": -p["elaAmp"], "base": p["elaBase"], "period": p["elaPeriod"]})
    litho_fn = gc.litho_make(litho_spec) if litho_spec else None
    prime_ice(st, p, ELA)
    for _ in range(nsteps):
        step(st, p, U, ELA, litho_fn)
    st["p"] = p
    return st
