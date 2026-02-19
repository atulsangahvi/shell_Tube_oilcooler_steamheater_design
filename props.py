from CoolProp.CoolProp import PropsSI
import math

def water_props(T_c: float, P_kpa: float = 300.0):
    T = T_c + 273.15
    P = P_kpa * 1000.0
    rho = PropsSI("D","T",T,"P",P,"Water")
    mu  = PropsSI("V","T",T,"P",P,"Water")
    k   = PropsSI("L","T",T,"P",P,"Water")
    cp  = PropsSI("C","T",T,"P",P,"Water")
    pr  = cp*mu/max(k,1e-12)
    return dict(rho=rho, mu=mu, k=k, cp=cp, pr=pr)

def steam_sat_props(P_kpa: float):
    P = P_kpa * 1000.0
    T_sat = PropsSI("T","P",P,"Q",0,"Water")
    h_g = PropsSI("H","P",P,"Q",1,"Water")
    h_f = PropsSI("H","P",P,"Q",0,"Water")
    rho_l = PropsSI("D","P",P,"Q",0,"Water")
    mu_l  = PropsSI("V","P",P,"Q",0,"Water")
    k_l   = PropsSI("L","P",P,"Q",0,"Water")
    cp_l  = PropsSI("C","P",P,"Q",0,"Water")
    sigma = PropsSI("I","P",P,"Q",0,"Water")
    return dict(t_sat_c=T_sat-273.15, h_fg=h_g-h_f, rho_l=rho_l, mu_l=mu_l, k_l=k_l, cp_l=cp_l, sigma=sigma)

def brine_props(base: str, x_mass: float, T_c: float):
    T = T_c + 273.15
    f1 = f"INCOMP::{base}[{x_mass:.6f}]"
    try:
        rho = PropsSI("D","T",T,"P",101325.0,f1)
        mu  = PropsSI("V","T",T,"P",101325.0,f1)
        k   = PropsSI("L","T",T,"P",101325.0,f1)
        cp  = PropsSI("C","T",T,"P",101325.0,f1)
    except Exception:
        f2 = f"INCOMP::{base}-{x_mass*100:.3f}%"
        rho = PropsSI("D","T",T,"P",101325.0,f2)
        mu  = PropsSI("V","T",T,"P",101325.0,f2)
        k   = PropsSI("L","T",T,"P",101325.0,f2)
        cp  = PropsSI("C","T",T,"P",101325.0,f2)
    pr = cp*mu/max(k,1e-12)
    return dict(rho=rho, mu=mu, k=k, cp=cp, pr=pr)

def astm_d341_fit(nu40, nu100):
    def Y(nu): return math.log10(math.log10(nu + 0.7))
    T1, T2 = 313.15, 373.15
    X1, X2 = math.log10(T1), math.log10(T2)
    Y1, Y2 = Y(nu40), Y(nu100)
    B = (Y1 - Y2) / max((X2 - X1), 1e-12)
    A = Y1 + B * X1
    return A, B

def astm_d341_nu(A, B, T_c):
    T = T_c + 273.15
    Y = A - B * math.log10(T)
    return max((10 ** (10 ** Y)) - 0.7, 0.2)

def oil_props_from_d341(T_c, nu40, nu100, rho15):
    A, B = astm_d341_fit(nu40, nu100)
    nu_cSt = astm_d341_nu(A, B, T_c)
    beta = 0.0007
    rho = max(rho15 * (1 - beta * (T_c - 15.0)), 500.0)
    mu = nu_cSt * 1e-6 * rho
    cp, k = 2000.0, 0.13
    pr = cp * mu / max(k, 1e-12)
    return dict(rho=rho, mu=mu, k=k, cp=cp, pr=pr, nu_cSt=nu_cSt)
