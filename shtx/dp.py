from shtx.htc import friction_factor_blasius, reynolds

def tube_dp_darcy(m, rho, mu, D, L, A, passes=1, K=0.0):
    v = m / (rho * max(A, 1e-12))
    Re = reynolds(rho, v, D, mu)
    f = friction_factor_blasius(Re)
    return f * (L * passes / max(D, 1e-12)) * 0.5 * rho * v * v + K * 0.5 * rho * v * v

def shell_dp_ideal_kern(m, rho, A, K=1.2):
    v = m / (rho * max(A, 1e-12))
    return K * 0.5 * rho * v * v
