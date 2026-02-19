import math
from shtx.htc import friction_factor_blasius, reynolds

def tube_dp_darcy(m, rho, mu, D, L, A, passes=1, K=0.0):
    v = m / (rho * max(A, 1e-12))
    Re = reynolds(rho, v, D, mu)
    f = friction_factor_blasius(Re)
    return f * (L * passes / max(D, 1e-12)) * 0.5 * rho * v * v + K * 0.5 * rho * v * v

def shell_dp_kern_segmental(m, rho, mu, v_max, De, Ds, Nb, K_window=1.2):
    # Practical Kern-style segmental baffle pressure drop:
    # Crossflow friction + window/turning allowance.
    Re = rho * v_max * De / max(mu, 1e-12)
    # Kern friction factor (smooth tube bank, engineering approximation)
    f = 0.14 * (max(Re, 1e-6) ** -0.2)
    # crossflow length scale per baffle compartment ~ Ds/De
    dp_cross = Nb * f * (Ds / max(De, 1e-12)) * 0.5 * rho * v_max * v_max
    dp_window = Nb * K_window * 0.5 * rho * v_max * v_max
    return dp_cross + dp_window
