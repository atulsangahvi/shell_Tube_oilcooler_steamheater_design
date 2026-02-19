import math

def tube_outside_area(n, do, L):
    return n * math.pi * do * L

def tube_flow_area(n, di, passes):
    tpp = max(n / max(passes, 1), 1.0)
    return tpp * math.pi * (di**2) / 4

def approx_shell_crossflow_area(Ds, bs, bc):
    # very simple crossflow area estimate, used for velocity scaling
    return Ds * bs * max(1 - bc, 0.2)

def tube_pitch(do_m: float, p_over_do: float):
    return do_m * p_over_do
def shell_crossflow_area(Ds, bs, bc):
    # Compatibility alias for newer solver.py
    return approx_shell_crossflow_area(Ds, bs, bc)
