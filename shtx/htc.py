import math

def reynolds(rho, v, D, mu):
    return rho * v * D / max(mu, 1e-12)

def friction_factor_blasius(Re):
    return 0.3164 / (Re**0.25) if Re > 3000 else 64.0 / max(Re, 1e-12)

def tube_singlephase_htc(m, rho, mu, k, cp, D, A):
    v = m / (rho * max(A, 1e-12))
    Re = reynolds(rho, v, D, mu)
    Pr = cp * mu / max(k, 1e-12)
    if Re < 2300:
        Nu = 3.66
    else:
        f = friction_factor_blasius(Re)
        Nu = (f / 8 * (Re - 1000) * Pr) / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2/3) - 1))
    return Nu * k / max(D, 1e-12), Re, Pr, v

def shell_htc_zukauskas(rho, mu, k, cp, do, v_max, layout="tri30"):
    # Tube bank correlation (Zukauskas-style) with staggered vs inline constants.
    # Uses do as characteristic; v_max should be maximum crossflow velocity through bundle.
    Re = rho * v_max * do / max(mu, 1e-12)
    Pr = cp * mu / max(k, 1e-12)

    # constants for tube banks (engineering-grade)
    if layout == "tri30":  # staggered
        C, m = (0.35, 0.60) if Re < 1e3 else (0.27, 0.63)
        n = 0.36
    else:  # inline
        C, m = (0.30, 0.60) if Re < 1e3 else (0.21, 0.64)
        n = 0.36

    Nu = C * (max(Re, 1e-6) ** m) * (max(Pr, 1e-6) ** n)
    h = Nu * k / max(do, 1e-12)
    return h, Re, Pr

def bell_delaware_apply(h, Jc, Jl, Jb, Jr):
    return h * Jc * Jl * Jb * Jr

def nusselt_film_condensation_horizontal_bank(t_sat, t_wall, p, do, N_rows=1):
    dt = max(t_sat - t_wall, 0.2)
    g = 9.81
    h1 = 0.725 * ((p["rho_l"] * p["rho_l"] * g * p["h_fg"] * (p["k_l"]**3)) / (p["mu_l"] * do * dt))**0.25
    return h1 * (max(N_rows, 1.0)**(-0.25))
