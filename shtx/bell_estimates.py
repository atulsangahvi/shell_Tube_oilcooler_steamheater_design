"""Practical (approx) estimators for Bell–Delaware leakage/bypass corrections.

These are NOT a full Bell–Delaware implementation; they are a controlled, transparent approximation
to get you started without HTRI geometry export.

Inputs are typical manufacturing/assembly parameters:
- tube-to-baffle diametral clearance
- shell-to-baffle diametral clearance
- tube layout + pitch ratio
- sealing strips (rough bypass reduction)

You can always override J-factors manually in the app.
"""

import math

def estimate_leakage_fractions(n_tubes:int, do_m:float, Ds_m:float, delta_tb_m:float, delta_sb_m:float):
    # crude leakage areas
    # tube-to-baffle leakage area ~ sum of annular gaps around tube holes
    A_tb = n_tubes * math.pi * do_m * max(delta_tb_m, 0.0)
    # shell-to-baffle leakage area ~ circumference * clearance
    A_sb = math.pi * Ds_m * max(delta_sb_m, 0.0)

    A_leak = A_tb + A_sb
    return dict(A_tb=A_tb, A_sb=A_sb, A_leak=A_leak)

def estimate_Jl_from_clearances(A_leak:float, A_cross:float):
    # map leakage ratio to a reasonable Jl range [0.4 .. 1.0]
    r = A_leak / max(A_cross, 1e-12)
    # smooth decay: small leakage -> near 1; big leakage -> down to ~0.4
    Jl = 1.0 / (1.0 + 3.0*r)
    return max(0.4, min(1.0, Jl))

def estimate_Jb_from_bypass(p_over_do:float, sealing_strips:int):
    # bypass tendency grows with pitch ratio and bundle/shell gap; sealing strips reduce it.
    # Start from a typical baseline and apply modifiers.
    base = 0.90
    pitch_penalty = max(0.0, (p_over_do - 1.25))  # 0 at 1.25
    Jb = base - 0.15*pitch_penalty
    # sealing strips: each pair provides diminishing returns
    ss = max(int(sealing_strips), 0)
    Jb += 0.03 * (1.0 - math.exp(-0.6*ss))
    return max(0.5, min(1.0, Jb))

def estimate_Jc_from_baffle_cut(baffle_cut_frac:float):
    # Jc mainly accounts for baffle cut/window effect on crossflow heat transfer.
    # Keep mild dependence and stay in [0.7..1.0]
    bc = max(0.15, min(0.45, baffle_cut_frac))
    Jc = 1.0 - 0.6*(bc - 0.25)
    return max(0.7, min(1.0, Jc))
