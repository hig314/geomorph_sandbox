"""
xsection_experiment.py — an exploratory cross-valley model, to see what a glacier, uplift and
rock-slope failure do to a valley section, and so what the schematic section of
relief_limit/glacier.html should be simplified to. NOT a validated tier (d2_across is not
built); the rules are the simplest that carry the mechanisms, and the numbers are read for
shape, not magnitude.

Half-section z(x), valley centre x = 0, crest x = L, symmetric about both.
  uplift     U everywhere
  ice        flat surface at z(0) + Hmax*H(t), H(t) the figure's ice history (glacier_buttress)
  erosion    the floor is cut down at Kg and a parabolic trough z_floor + Cp*x^2 climbs from it,
             under the ice only (the trough is asserted, as in the figure: a first try with
             erosion ~ (ice depth)^2 cut a V-notch at the centre, not a trough; the real pattern
             needs the ice-flow solution of Harbor 1992, which this is not)
  failure    Culmann planar wedges, searched over toe position and plane dip. For the wedge
             above a plane of dip th through a toe at x_i (to where the plane leaves the rock,
             or the crest):
               FoS = (c*l + (W cos th + P sin th) tan phi + P cos th) / (W sin th)
             W the wedge weight, l the plane length, P = 0.5*gamma_i*d_toe^2 the ice thrust on
             the wedge's face (d_toe the ice depth at the toe). The lowest FoS below 1 fails:
             the surface drops to the plane. Debris is removed if there is ice in the valley
             (carried away), else spread level on the floor (area conserved).
  damage     optional: cohesion falls by a fraction for every unit change in ice thickness
             (as the figure's weakening) and is restored in the rock a landslide exposes; here
             applied to the whole section as one number, reset by any failure larger than Amin.

Run:  /opt/anaconda3/bin/python3 xsection_experiment.py   (writes xsection_*.png to $TMPDIR)
"""
import math, os, sys
import numpy as np
import glacier_buttress as gb

DEG = math.pi / 180


def run(ice="thick", L=3000.0, N=150, tmax=200.0, dt=0.05, U=5.0, Kg=10.0, Hmax=900.0,
        c=1.3e6, Cp=5e-4, phi=35.0, gamma=26e3, gamma_i=9e3, slope0=0.42, damage=0.0, cmin_frac=0.0, seed=7, check=0.2):
    """times in kyr, lengths in m, rates in m/kyr, c in Pa. Returns dict of histories."""
    x = np.linspace(0, L, N + 1); dx = x[1] - x[0]
    z = slope0 * x.copy()
    sim = gb.simulate(dict(ice=ice, seed=seed)); icep = sim["ice"]
    tanphi = math.tan(phi * DEG)
    ths = np.arange(26, 72, 2.0) * DEG
    out = dict(x=x, t=[], relief=[], H=[], fails=[], snaps=[], c=[])
    dmg = 0.0; Hprev = gb.value_at(icep, 0.0)
    nsteps = int(round(tmax / dt)); next_check = 0.0; snap_t = set(np.round(np.arange(0, tmax + 1e-9, 2.0), 6))

    def search(z, zs, ceff):
        best = None
        for i in range(0, N - 2, 2):
            dtoe = max(0.0, zs - z[i]); P = 0.5 * gamma_i * dtoe * dtoe
            xr = x[i:] - x[i]
            for th in ths:
                zp = z[i] + math.tan(th) * xr
                above = z[i:] - zp
                if above[1] <= 0: continue
                neg = np.nonzero(above[1:] <= 0)[0]
                j = (neg[0] + 1) if len(neg) else len(above) - 1       # exit index (relative), or the crest
                A = np.trapezoid(np.maximum(above[:j + 1], 0), dx=dx)
                if A < 2000.0: continue
                W = gamma * A; l = xr[j] / math.cos(th)
                fos = (ceff * l + (W * math.cos(th) + P * math.sin(th)) * tanphi + P * math.cos(th)) / (W * math.sin(th))
                if fos < 1 and (best is None or fos < best[0]): best = (fos, i, th, j, A)
        return best

    for n in range(nsteps + 1):
        t = n * dt; H = gb.value_at(icep, t); zs = z[0] + Hmax * H if H > 0 else -1e9
        if damage > 0:
            dmg = min(1.0 - cmin_frac, dmg + damage * abs(H - Hprev)); Hprev = H
        ceff = c * (1 - dmg)
        if t >= next_check - 1e-9:
            next_check += check
            for _ in range(20):
                b = search(z, zs, ceff)
                if b is None: break
                fos, i, th, j, A = b
                zp = z[i] + math.tan(th) * (x[i:] - x[i]); zold = z.copy()
                z[i:i + j + 1] = np.minimum(z[i:i + j + 1], zp[:j + 1])
                removed = np.trapezoid(zold - z, dx=dx)
                if H <= 0:                                              # no ice: debris level on the floor
                    lo, hi = z[0], z[0] + 2000.0
                    for _k in range(50):
                        f = 0.5 * (lo + hi)
                        if np.trapezoid(np.maximum(0, f - z[:i + 1]), dx=dx) < removed: lo = f
                        else: hi = f
                    z[:i + 1] = np.maximum(z[:i + 1], f)
                out["fails"].append(dict(t=t, toe_h=zold[i] - zold[0], toe_x=x[i], th=th / DEG, A=A, H=H,
                                         ice_h=(Hmax * H if H > 0 else 0.0), crest_drop=zold[-1] - z[-1],
                                         relief0=zold[-1] - zold[0], relief1=z[-1] - z[0], exit_x=x[i + j]))
                if damage > 0 and A > 2e4: dmg = 0.0; ceff = c
        if round(t, 6) in snap_t: out["snaps"].append((t, z - z[0], Hmax * H))
        out["t"].append(t); out["relief"].append(z[-1] - z[0]); out["H"].append(H); out["c"].append(ceff)
        # uplift and glacial erosion
        if H > 0:
            z0 = z[0] + (U - Kg) * dt                       # the floor: uplift less glacial erosion
            zu = z + U * dt; par = z0 + Cp * x * x
            z = np.where(par < zs, np.minimum(zu, par), zu)  # the U is cut only under the ice
        else:
            z = z + U * dt
    return out


def summarize(name, o):
    f = o["fails"]; r = np.array(o["relief"])
    print("== %s: %d failures; relief min %.0f mean %.0f max %.0f m" % (name, len(f), r.min(), r.mean(), r.max()))
    if not f: return
    big = [q for q in f if q["A"] > 2e4]
    th = np.array([q["th"] for q in f]); A = np.array([q["A"] for q in f])
    print("   plane dip: area-weighted mean %.0f deg (range %.0f-%.0f)" % ((th * A).sum() / A.sum(), th.min(), th.max()))
    wi = [q for q in f if q["H"] > 0]
    if wi:
        rel = np.array([(q["toe_h"] - q["ice_h"]) / q["relief0"] for q in wi]); Aw = np.array([q["A"] for q in wi])
        print("   with ice: toe height minus ice surface, as a fraction of relief: area-weighted mean %+.2f (10-90%%: %+.2f to %+.2f)"
              % ((rel * Aw).sum() / Aw.sum(), np.percentile(rel, 10), np.percentile(rel, 90)))
    toe = np.array([q["toe_h"] / q["relief0"] for q in f])
    print("   toe height / relief: area-weighted mean %.2f" % ((toe * A).sum() / A.sum()))
    cd = np.array([q["crest_drop"] / q["relief0"] for q in f])
    print("   crest drop / relief per failure: mean %.3f, max %.2f; failures that lower the crest: %d of %d" % (cd.mean(), cd.max(), (cd > 1e-6).sum(), len(f)))
    drop = np.array([(q["relief0"] - q["relief1"]) / q["relief0"] for q in f])
    print("   relief lost per failure: mean %.3f, max %.2f of relief" % (drop.mean(), drop.max()))
    print("   lowest relief / highest relief: %.2f" % (r[int(len(r) * 0.1):].min() / r.max()))


def plot(name, o, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(11, 9), gridspec_kw=dict(height_ratios=[1, 2]))
    t = np.array(o["t"]); ax[0].plot(t, o["relief"], "k"); ax[0].set_ylabel("relief (m)")
    a2 = ax[0].twinx(); a2.fill_between(t, 0, o["H"], color="#a9cbe8", alpha=0.5); a2.set_ylim(0, 4)
    for q in o["fails"]: ax[0].axvline(q["t"], color="r", alpha=min(1, q["A"] / 1e5), lw=0.6)
    ax[0].set_title(name)
    x = o["x"]; snaps = o["snaps"]; cm = plt.get_cmap("viridis")
    for k, (ts, zz, hi) in enumerate(snaps):
        if k % 5: continue
        ax[1].plot(x, zz - zz[-1] + 0, color=cm(ts / t[-1]), lw=1)        # crest-anchored
    ax[1].set_aspect("equal"); ax[1].set_xlabel("x (m), valley centre to crest; colour = time; crest held at 0")
    fig.tight_layout(); fig.savefig(path, dpi=80); plt.close(fig)


if __name__ == "__main__":
    tmp = os.environ.get("XS_OUT", os.environ.get("TMPDIR", "."))
    cases = [("thick", dict(ice="thick")), ("retreat", dict(ice="retreat")),
             ("retreat_damage", dict(ice="retreat", damage=1.0, cmin_frac=0.1)),
             ("thick_damage", dict(ice="thick", damage=1.0, cmin_frac=0.1))]
    for name, kw in cases:
        if len(sys.argv) > 1 and name not in sys.argv[1:]: continue
        o = run(**kw); summarize(name, o); plot(name, o, os.path.join(tmp, "xsection_%s.png" % name))
