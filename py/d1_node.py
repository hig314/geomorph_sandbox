"""
d1_node.py — Python mirror of d1/model.js createNodeColumn(): one node of the
along-valley profile run as an isolated column.

Everything the profile computed from ADJACENCY is frozen at its isolation-time value and
asserted: the downstream neighbour keeps doing what it was doing — it uplifts with the shared
forcing and erodes at its isolation-time rate e_down, so z_down(t) = z_down(t0) + ΔU_cum −
e_down·(t − t0) (flag extrapolate; False supplies z_down directly) — its ice thickness H_down,
the ice volume flux Q leaving the node, the drainage area A and widths W, W_f. Shared FORCING
keeps running: U(t) (× the node's spatial factor) and ELA(t) — though ELA acts on a node only
through the flux, which is frozen, so the isolated column is blind to it. For a profile in
steady state every isolated column stays identical to its node; divergence is adjacency
changing (neighbours speeding up or slowing down, a terminus moving over the node).
Two rejected variants: freezing z_down itself (correct in steady state, but a receiver that
speeds up or slows down is missed); letting z_down uplift without eroding (an infinitely
hard dam: the node erodes its slope away and then rides uplift with no erosion).
LOCAL laws and feedbacks remain: lithology in the material frame, slope against the
frozen receiver, thickness from the face flux law against the frozen downstream surface.

Per step (mirrors d2_along.step for one node):
  r = litho(s, z − U_cum);  z += U f dt;  U_cum += U f dt
  if Q > 0:  H = face_thickness(z, z_down + H_down, H_down, Q/W_f, ds, Γ′, n)
             U_s = f_s (Q/W)/H;  E_g = min(K_g f_Kg U_s^l, cap);  z −= E_g dt
  else:      z = yuan_node(z, z, z_down, K f_K A^m, dt, ds, 0, n_exp)   (E_f = Δz/dt)
"""
import gs_core as gc

def make_node_column(spec):
    """spec: s, z0, ucum0, zDown, HDown, Q, A, W, Wf, ds, fac, t0, dt, params (dict of the
    d2_along parameters), litho (spec or None)."""
    p = spec["params"]
    U = gc.make_series({"shape": p["shape"], "peak": p["peakUplift"], "duration": p["duration"]})
    ELA = gc.make_series({"shape": "sine", "peak": -p["elaAmp"], "base": p["elaBase"], "period": p["elaPeriod"]})
    litho = gc.litho_make(spec.get("litho")) if spec.get("litho") else None
    Gam = gc.ICE["GAMMA"] * p["flow"] / (1.0 - p["fs"])
    n = gc.ICE["N"]
    st = {"t": spec["t0"], "z": spec["z0"], "ucum": spec["ucum0"], "H": 0.0, "Us": 0.0, "Ef": 0.0, "Eg": 0.0, "r": 0.0,
          "ecum": 0.0, "series": []}
    frozen = {k: spec[k] for k in ("zDown", "HDown", "Q", "A", "W", "Wf", "ds", "fac", "s")}
    frozen["extrapolate"] = spec.get("extrapolate", True); frozen["eDown"] = spec.get("eDown", 0.0)

    def record():
        st["series"].append({"t": st["t"], "z": st["z"], "H": st["H"], "Us": st["Us"], "Ef": st["Ef"], "Eg": st["Eg"], "r": st["r"],
                             "u": U(st["t"]) * frozen["fac"], "ela": ELA(st["t"]),
                             "zDown": st.get("zDown", frozen["zDown"])})

    def step():
        dt = spec["dt"]
        r = litho(frozen["s"], 0.0, gc.material_z(st["z"], st["ucum"])) if litho else 0.0
        st["r"] = r
        fK = gc.erodibility_factor(r, p["contrastK"]); fKg = gc.erodibility_factor(r, p["contrastKg"])
        du = U(st["t"]) * frozen["fac"] * dt
        st["z"] += du; st["ucum"] += du
        zDown = frozen["zDown"] + (((st["ucum"] - spec["ucum0"]) - frozen["eDown"] * (st["t"] + dt - spec["t0"])) if frozen["extrapolate"] else 0.0)
        zsDown = zDown + frozen["HDown"]
        H = 0.0; Us = 0.0; Eg = 0.0
        if frozen["Q"] > 0:
            H = gc.face_thickness(st["z"], zsDown, frozen["HDown"], frozen["Q"] / frozen["Wf"], frozen["ds"], Gam, n)
            Us = gc.sliding_speed(frozen["Q"] / frozen["W"], H, p["fs"], p.get("Hmin", 10.0)) if H > 0 else 0.0
            Eg = gc.cap_rate(gc.abrasion(p["Kg"] * fKg, Us, p["lexp"]), p["eroCap"]) if H > 1.0 else 0.0
            st["z"] -= Eg * dt; st["ecum"] += Eg * dt
        # fluvial with the same under-ice efficiency ramp as the profile
        w = p.get("phiSub", 0.0) + (1.0 - p.get("phiSub", 0.0)) * min(1.0, max(0.0, 1.0 - H / p.get("Hf", 100.0)))
        Ef = 0.0
        if w > 0.0:
            h0 = st["z"]
            Kp = p["K"] * fK * w * frozen["A"] ** p["m"]
            nh = gc.yuan_node(h0, h0, zDown, Kp, dt, frozen["ds"], 0.0, p["nexp"], False)
            st["z"] = nh; Ef = (h0 - nh) / dt; st["ecum"] += h0 - nh
        st["H"] = H; st["Us"] = Us; st["Eg"] = Eg; st["Ef"] = Ef; st["zDown"] = zDown
        st["t"] += dt

    record()
    return {"state": st, "frozen": frozen, "step": step, "record": record, "U": U, "ELA": ELA}
