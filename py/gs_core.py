"""
gs_core.py — Python mirror of core/*.js (forcing, noise, lithology, laws).

Scalar, plain-Python mirrors of the JS functions, written to match the JS line
for line so a scheme validated here is the scheme that ships. Units: metres,
years. Run py/test_core.py with /opt/anaconda3/bin/python3.
"""
import math

# ---------------------------------------------------------------------------
# forcing
# ---------------------------------------------------------------------------
def _bell(t, p):
    T = p["duration"]; mu = T / 2; sigma = T / 6; z = (t - mu) / sigma
    return math.exp(-0.5 * z * z)

def _arc(t, p):
    x = (t - p["duration"] / 2) / (p["duration"] / 2)
    if x <= -1 or x >= 1:
        return 0.0
    return math.sqrt(1 - x * x)

def _constant(t, p):
    return 0.0 if t > p["duration"] else 1.0

def _triangle(t, p):
    x = abs((t - p["duration"] / 2) / (p["duration"] / 2))
    return max(0.0, 1 - x)

def _plateau(t, p):
    T = p["duration"]; r = T * 0.25
    if t <= 0 or t >= T:
        return 0.0
    if t < r:
        return t / r
    if t > T - r:
        return (T - t) / r
    return 1.0

def _pulse(t, p):
    T = p["duration"]
    if t <= 0:
        return 0.0
    rise = T * 0.1
    if t < rise:
        return t / rise
    return math.exp(-(t - rise) / (T * 0.3))

def _sine(t, p):
    return 0.5 * (1 - math.cos(2 * math.pi * t / p["period"]))

def _sawtooth(t, p):
    ph = t / p["period"]; ph -= math.floor(ph)
    return ph / 0.9 if ph < 0.9 else (1 - ph) / 0.1

SHAPES = {"bell": _bell, "arc": _arc, "constant": _constant, "triangle": _triangle,
          "plateau": _plateau, "pulse": _pulse, "sine": _sine, "sawtooth": _sawtooth}

def make_series(spec):
    p = {"shape": spec.get("shape", "bell"), "peak": spec.get("peak", 1.0),
         "base": spec.get("base", 0.0), "duration": spec.get("duration", 1e6),
         "period": spec.get("period", 1e5)}
    s = SHAPES.get(p["shape"], _bell)
    def f(t):
        return p["base"] + p["peak"] * s(t, p)
    f.spec = p
    return f

def integrate(f, t0, t1, steps=1000):
    dt = (t1 - t0) / steps
    s = 0.5 * (f(t0) + f(t1))
    for i in range(1, steps):
        s += f(t0 + i * dt)
    return s * dt

def spatial_factor(u, v, pattern, ramp=0.0):
    if pattern == "gaussian":
        ru = u - 0.5; rv = v - 0.5; sig = 0.18
        return math.exp(-(ru * ru + rv * rv) / (2 * sig * sig))
    if pattern == "tilt":
        return u
    if pattern == "ramp":
        dE = min(u, v, 1 - u, 1 - v)
        return (dE / ramp if dE < ramp else 1.0) if ramp > 0 else 1.0
    return 1.0

# ---------------------------------------------------------------------------
# noise (32-bit integer arithmetic mirrored with masks)
# ---------------------------------------------------------------------------
M32 = 0xFFFFFFFF

def _imul(a, b):
    r = (a & M32) * (b & M32) & M32
    return r - (1 << 32) if r >= (1 << 31) else r  # signed 32-bit like Math.imul

def _u32(x):
    return x & M32

def _mix(h):
    h = _imul(h ^ (_u32(h) >> 13), 1274126177)
    h = h ^ (_u32(h) >> 16)
    return _u32(h) / 4294967296.0

def hash2(ix, iy, seed):
    h = (_imul(ix, 374761393) + _imul(iy, 668265263) + _imul(seed, 982451653)) & M32
    return _mix(h)

def hash3(ix, iy, iz, seed):
    h = (_imul(ix, 374761393) + _imul(iy, 668265263) + _imul(iz, 1103515245) + _imul(seed, 982451653)) & M32
    return _mix(h)

def smoothstep(t):
    return t * t * (3 - 2 * t)

def value_noise2(x, y, seed):
    x0 = math.floor(x); y0 = math.floor(y)
    fx = smoothstep(x - x0); fy = smoothstep(y - y0)
    v00 = hash2(x0, y0, seed); v10 = hash2(x0 + 1, y0, seed)
    v01 = hash2(x0, y0 + 1, seed); v11 = hash2(x0 + 1, y0 + 1, seed)
    a = v00 + (v10 - v00) * fx
    b = v01 + (v11 - v01) * fx
    return a + (b - a) * fy

def value_noise3(x, y, z, seed):
    x0 = math.floor(x); y0 = math.floor(y); z0 = math.floor(z)
    fx = smoothstep(x - x0); fy = smoothstep(y - y0); fz = smoothstep(z - z0)
    v000 = hash3(x0, y0, z0, seed); v100 = hash3(x0 + 1, y0, z0, seed)
    v010 = hash3(x0, y0 + 1, z0, seed); v110 = hash3(x0 + 1, y0 + 1, z0, seed)
    v001 = hash3(x0, y0, z0 + 1, seed); v101 = hash3(x0 + 1, y0, z0 + 1, seed)
    v011 = hash3(x0, y0 + 1, z0 + 1, seed); v111 = hash3(x0 + 1, y0 + 1, z0 + 1, seed)
    a0 = v000 + (v100 - v000) * fx; b0 = v010 + (v110 - v010) * fx
    a1 = v001 + (v101 - v001) * fx; b1 = v011 + (v111 - v011) * fx
    c0 = a0 + (b0 - a0) * fy; c1 = a1 + (b1 - a1) * fy
    return c0 + (c1 - c0) * fz

def fbm2(x, y, seed, octaves, persistence, lacunarity):
    amp = 1.0; freq = 1.0; s = 0.0; norm = 0.0
    for o in range(octaves):
        s += amp * value_noise2(x * freq, y * freq, seed + o * 1013)
        norm += amp; amp *= persistence; freq *= lacunarity
    return s / norm

def fbm3(x, y, z, seed, octaves, persistence, lacunarity):
    amp = 1.0; freq = 1.0; s = 0.0; norm = 0.0
    for o in range(octaves):
        s += amp * value_noise3(x * freq, y * freq, z * freq, seed + o * 1013)
        norm += amp; amp *= persistence; freq *= lacunarity
    return s / norm

# ---------------------------------------------------------------------------
# lithology
# ---------------------------------------------------------------------------
def material_z(z, ucum):
    return z - ucum

def _window1(v, lo, hi, soft):
    if soft > 0:
        a = min(1.0, max(0.0, (v - lo) / soft)); b = min(1.0, max(0.0, (hi - v) / soft))
        return min(a, b)
    return 1.0 if lo <= v <= hi else 0.0

def _ramp01(v):
    return 0.0 if v < 0 else (1.0 if v > 1 else v)

def litho_make(spec):
    if isinstance(spec, list):
        fns = [litho_make(s) for s in spec]
        return lambda x, y, zm: max(f(x, y, zm) for f in fns)
    spec = spec or {}
    t = spec.get("type", "none")
    if t == "layer":
        top = spec.get("top", 0.0); thick = spec.get("thick", 200.0); soft = spec.get("soft", 0.0)
        return lambda x, y, zm: _window1(zm, top - thick, top, soft)
    if t == "slab":
        top = spec.get("top", 0.0); thick = spec.get("thick", 200.0); soft = spec.get("soft", 0.0)
        dipX = spec.get("dipX", 0.0); dipY = spec.get("dipY", 0.0)
        def slab(x, y, zm):
            tt = top + dipX * x + dipY * y
            return _window1(zm, tt - thick, tt, soft)
        return slab
    if t == "dike":
        x0 = spec.get("x0", 0.0); w = spec.get("width", 100.0); soft = spec.get("soft", 0.0)
        return lambda x, y, zm: _window1(x, x0 - w / 2, x0 + w / 2, soft)
    if t == "blob":
        seed = int(spec.get("seed", 0)); scale = spec.get("scale", 2000.0); thr = spec.get("threshold", 0.6)
        oct_ = spec.get("octaves", 3); pers = spec.get("persistence", 0.5); soft = spec.get("soft", 0.05)
        def blob(x, y, zm):
            v = fbm3(x / scale, y / scale, zm / scale, seed, oct_, pers, 2)
            return _ramp01((v - thr) / soft) if soft > 0 else (1.0 if v >= thr else 0.0)
        return blob
    if t == "patch":
        seed = int(spec.get("seed", 0)); scale = spec.get("scale", 2000.0); thr = spec.get("threshold", 0.6)
        oct_ = spec.get("octaves", 3); pers = spec.get("persistence", 0.5); soft = spec.get("soft", 0.05)
        def patch(x, y, zm):
            v = fbm2(x / scale, y / scale, seed, oct_, pers, 2)
            return _ramp01((v - thr) / soft) if soft > 0 else (1.0 if v >= thr else 0.0)
        return patch
    return lambda x, y, zm: 0.0

def erodibility_factor(r, c):
    return 1.0 / (1.0 + r * (c - 1.0))

def strength_factor(r, c):
    return 1.0 + r * (c - 1.0)

# ---------------------------------------------------------------------------
# laws
# ---------------------------------------------------------------------------
ICE = {"RHOG": 8987.0, "N": 3, "A": 7.57e-17, "ABL_RATIO": 2.5, "SLOPE_MIN": 1e-3}
ICE["GAMMA"] = (2 * ICE["A"] / (ICE["N"] + 2)) * ICE["RHOG"] ** ICE["N"]

def reservoir_k(Upeak, Rlim, n):
    return Upeak / Rlim ** n

def reservoir_erosion(R, k, n):
    return k * (R if R > 0 else 0.0) ** n

def reservoir_equilibrium(U, k, n):
    return (U / k) ** (1.0 / n)

def stream_power(K, A, S, m, n):
    return K * A ** m * (S if S > 0 else 0.0) ** n

def yuan_node(h0, hcur, hr, Kp, dt, dist, dep, nexp, submerged=False):
    if dep < 0:
        dep = 0.0
    if h0 > hr and not submerged:
        if abs(nexp - 1) < 1e-9:
            f = Kp * dt / dist
            nh = (h0 + f * hr + dep) / (1 + f)
        else:
            z = hcur
            for _ in range(20):
                s = (z - hr) / dist
                if s < 0:
                    s = 0.0
                g = z - h0 + Kp * dt * s ** nexp - dep
                dg = 1 + (Kp * dt * nexp * s ** (nexp - 1)) / dist
                dz = g / dg
                z -= dz
                if abs(dz) < 1e-4:
                    break
            nh = z
        if nh < hr:
            nh = hr
    else:
        nh = h0 + dep
    return nh

def creep_flux_linear(D, S):
    return -D * S

def creep_flux_roering(D, S, Sc, fmin=0.02):
    f = 1 - (S * S) / (Sc * Sc)
    if f < fmin:
        f = fmin
    return -D * S / f

def mass_balance(zs, ela, bg, accMax=None, abMax=None, ablRatio=None):
    ds = zs - ela
    b = ds * bg if ds > 0 else ds * bg * (ablRatio if ablRatio is not None else ICE["ABL_RATIO"])
    if accMax is not None and b > accMax:
        b = accMax
    if abMax is not None and b < -abMax:
        b = -abMax
    return b

def sia_flux(H, S, Gamma=None, n=None):
    Gamma = ICE["GAMMA"] if Gamma is None else Gamma; n = ICE["N"] if n is None else n
    return Gamma * H ** (n + 2) * abs(S) ** n

def sia_thickness(q, S, Gamma=None, n=None, Smin=None):
    Gamma = ICE["GAMMA"] if Gamma is None else Gamma; n = ICE["N"] if n is None else n
    Smin = ICE["SLOPE_MIN"] if Smin is None else Smin
    s = abs(S)
    if s < Smin:
        s = Smin
    if q <= 0:
        return 0.0
    return (q / (Gamma * s ** n)) ** (1.0 / (n + 2))

def face_thickness(zi, zs_down, H_down, qf, ds, Gamma, n):
    """Thickness H_i at a node whose downstream neighbour has surface zs_down and
    thickness H_down, such that the SIA face flux carries qf (m²/yr):
        Gamma · (½(H_i + H_down))^(n+2) · ((zi + H_i − zs_down)/ds)^n = qf
    Monotone in H_i → unique root; bisection to 1e-3 m. qf ≤ 0 → 0."""
    if qf <= 0.0:
        return 0.0
    lo = max(0.0, zs_down - zi)
    hi = lo + 10.0
    def f(x):
        Hf = 0.5 * (x + H_down); Sf = (zi + x - zs_down) / ds
        return Gamma * Hf ** (n + 2) * Sf ** n - qf
    while f(hi) < 0.0 and hi < 1e5:
        hi *= 2.0
    x = 0.5 * (lo + hi)
    for _ in range(60):
        if f(x) > 0.0:
            hi = x
        else:
            lo = x
        x = 0.5 * (lo + hi)
        if hi - lo < 1e-3:
            break
    return x

def sliding_speed(q, H, fs):
    return fs * q / H if H > 0 else 0.0

def abrasion(Kg, Us, l=1.0):
    return Kg * (Us if Us > 0 else 0.0) ** l

def quarrying(Kq, Us, convexity):
    return Kq * Us * convexity if convexity > 0 else 0.0

def buzzsaw_rate(z, ela, width, peak):
    u = (z - ela) / width
    return peak * math.exp(-0.5 * u * u)

def cap_rate(E, cap):
    return cap if E > cap else E

def rk4(f, t, y, dt):
    k1 = f(t, y)
    k2 = f(t + dt / 2, y + (dt / 2) * k1)
    k3 = f(t + dt / 2, y + (dt / 2) * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
