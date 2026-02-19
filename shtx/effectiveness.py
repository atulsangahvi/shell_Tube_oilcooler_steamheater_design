import math

def eps_crossflow_unmixed(NTU, Cr):
    NTU = max(NTU, 0.0)
    Cr = max(Cr, 0.0)
    if Cr < 1e-12:
        return 1 - math.exp(-NTU)
    return 1 - math.exp((math.exp(-Cr * NTU) - 1) / Cr)

def eps_1shell(NTU, Cr):
    NTU = max(NTU, 0.0)
    Cr = min(max(Cr, 0.0), 0.999999999)
    if Cr < 1e-12:
        return 1 - math.exp(-NTU)
    s = math.sqrt(1 + Cr * Cr)
    x = NTU * s
    ex = math.exp(x)
    coth = (ex + 1) / max(ex - 1, 1e-12)
    eps = 2 / (1 + Cr + s * coth)
    return min(max(eps, 0.0), 1.0)
