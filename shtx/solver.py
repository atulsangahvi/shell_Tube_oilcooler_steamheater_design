import math
from dataclasses import dataclass
from shtx.geometry import tube_outside_area, tube_flow_area, shell_crossflow_area, tube_pitch, equiv_diameter
from shtx.htc import tube_singlephase_htc, shell_htc_zukauskas, bell_delaware_apply, nusselt_film_condensation_horizontal_bank
from shtx.dp import tube_dp_darcy, shell_dp_kern_segmental
from shtx.effectiveness import eps_crossflow_unmixed, eps_1shell
from shtx.bell_estimates import estimate_leakage_fractions

@dataclass
class SegRes:
    flow_model: str
    segments: int
    duty_W: float
    U_avg: float
    Th_out: float
    Tc_out: float
    dp_shell_kpa: float
    dp_tube_kpa: float
    area_m2: float
    J_used: dict
    leak: dict
    debug: dict

def overall_Uo(hs, ht, Rfs, Rft, do, di, tw, kw):
    Ao_Ai = do / max(di, 1e-12)
    Rw = tw / max(kw, 1e-12) * Ao_Ai
    R = 1 / max(hs, 1e-12) + Rfs + Rw + Ao_Ai / max(ht, 1e-12) + Rft
    return 1 / max(R, 1e-12)

def segmented_singlephase(p_shell, p_tube, Th_in, Tc_in, mh, mc, geom, foul, J, flow_model="1-shell", N=20, kw=16.0):
    do = geom["do"]; di = geom["di"]; tw = geom["tw"]; L = geom["L"]; n = geom["n"]; passes = geom["passes"]
    Ds = geom["Ds"]; bs = geom["bs"]; bc = geom["bc"]

    layout = geom.get("layout","tri30")
    p_over_do = geom.get("p_over_do",1.25)
    Pt = tube_pitch(do, p_over_do)
    De = equiv_diameter(Pt, do, layout)
    As = shell_crossflow_area(Ds, bs, bc, do, Pt, layout)

    A = tube_outside_area(n, do, L)
    dA = A / max(N, 1)
    At = tube_flow_area(n, di, passes)

    leak = estimate_leakage_fractions(n, do, Ds, geom.get("delta_tb_m", 0.0), geom.get("delta_sb_m", 0.0))
    leak["r_leak"] = leak["A_leak"] / max(As, 1e-12)

    Th, Tc = Th_in, Tc_in
    Qt = 0.0
    Us = 0.0
    dpt = 0.0
    dps = 0.0

    # baffle compartments for dp (approx)
    Nb = max(int(round(L / max(bs, 1e-6))), 1)

    debug_mid=None

    for i in range(N):
        ph = p_shell(Th)
        pc = p_tube(Tc)

        ht, Re_t, Pr_t, v_t = tube_singlephase_htc(mc, pc["rho"], pc["mu"], pc["k"], pc["cp"], di, At)

        v_max = mh / (ph["rho"] * max(As, 1e-12))
        hid, Re_s, Pr_s = shell_htc_zukauskas(ph["rho"], ph["mu"], ph["k"], ph["cp"], do, v_max, layout)
        hs = bell_delaware_apply(hid, J["Jc"], J["Jl"], J["Jb"], J["Jr"])

        U = overall_Uo(hs, ht, foul["Rf_shell"], foul["Rf_tube"], do, di, tw, kw)
        Us += U

        Ch = mh * ph["cp"]
        Cc = mc * pc["cp"]
        Cmin = min(Ch, Cc)
        Cr = Cmin / max(max(Ch, Cc), 1e-12)
        NTU = U * dA / max(Cmin, 1e-12)

        eps = eps_crossflow_unmixed(NTU, Cr) if flow_model == "crossflow" else eps_1shell(NTU, Cr)
        dQ = eps * Cmin * (Th - Tc)
        Qt += dQ

        Th -= dQ / max(Ch, 1e-12)
        Tc += dQ / max(Cc, 1e-12)

        dL = L / max(N, 1)
        dpt += tube_dp_darcy(mc, pc["rho"], pc["mu"], di, dL, At, passes, 0.0)

        # distribute shell dp uniformly over segments (total computed using midpoint props)
        if i == int(N/2):
            debug_mid=dict(v_max=v_max, As=As, De=De, Re_shell=Re_s, Pr_shell=Pr_s, ht=ht, hs=hs, hid=hid, Re_tube=Re_t, v_tube=v_t)
            dp_total = shell_dp_kern_segmental(mh, ph["rho"], ph["mu"], v_max, De, Ds, Nb, K_window=1.6) / max(J["Jl"]*J["Jb"],1e-12)
            dps = dp_total

    return SegRes(flow_model, N, Qt, Us / max(N, 1), Th, Tc, dps / 1000.0, dpt / 1000.0, A, J, leak, debug_mid or {})

def steam_heater_simple(p_tube, tci, mc, psteam, geom, foul, kw=16.0):
    do = geom["do"]; di = geom["di"]; tw = geom["tw"]; L = geom["L"]; n = geom["n"]; passes = geom["passes"]
    A = tube_outside_area(n, do, L)
    At = tube_flow_area(n, di, passes)

    pc = p_tube(tci)
    ht, Re, _, vt = tube_singlephase_htc(mc, pc["rho"], pc["mu"], pc["k"], pc["cp"], di, At)

    sat = geom["steam_sat"](psteam)
    Ts = sat["t_sat_c"]
    hs = nusselt_film_condensation_horizontal_bank(Ts, tci + 0.35 * (Ts - tci), sat, do, max(int(n / 10), 1))

    U = overall_Uo(hs, ht, foul["Rf_shell"], foul["Rf_tube"], do, di, tw, kw)
    Cc = mc * pc["cp"]
    NTU = U * A / max(Cc, 1e-12)
    eps = 1 - math.exp(-NTU)

    tco = tci + eps * (Ts - tci)
    Q = Cc * (tco - tci)

    return dict(Q=Q, U=U, tco=tco, Ts=Ts, A=A, Re=Re, vt=vt, hs=hs, ht=ht)
