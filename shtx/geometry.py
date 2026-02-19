import math

def tube_outside_area(n, do, L):
    return n * math.pi * do * L

def tube_flow_area(n, di, passes):
    tpp = max(n / max(passes, 1), 1.0)
    return tpp * math.pi * (di**2) / 4

def tube_pitch(do_m: float, p_over_do: float):
    return do_m * p_over_do

# Backward-compatible name used by older code
def approx_shell_crossflow_area(Ds, bs, bc):
    return Ds * bs * max(1 - bc, 0.2)

def equiv_diameter(*args, **kwargs):
    """
    Compatible with both call styles:
      - equiv_diameter(Pt, do, layout)
      - equiv_diameter(do, pitch, layout="30° triangular" / "90° square" / "tri30" / "sq90")
    Returns De in meters.
    """
    if len(args) >= 3:
        a0, a1, a2 = args[0], args[1], args[2]
        # If user passes layout strings like "30° triangular", normalize
        layout = str(a2)
        Pt = float(a0)
        do = float(a1)
    else:
        do = float(kwargs.get("do"))
        Pt = float(kwargs.get("pitch"))
        layout = str(kwargs.get("layout", "tri30"))

    # normalize layout
    L = layout.lower()
    tri = ("tri" in L) or L.startswith("30")

    # Kern-style equivalent diameter for tube banks
    if tri:
        A_flow = 0.866 * Pt * Pt - math.pi * do * do / 4
    else:
        A_flow = Pt * Pt - math.pi * do * do / 4

    De = 4 * A_flow / max(math.pi * do, 1e-12)
    return max(De, 1e-9)

def shell_crossflow_area(*args, **kwargs):
    """
    Compatible with both call styles:
      - shell_crossflow_area(Ds, B, bc)              [older simplified solver]
      - shell_crossflow_area(Ds, B, bc, do, Pt, layout)  [newer solver]
    Returns A_cross in m².
    """
    # Parse positional arguments
    if len(args) >= 3:
        Ds = float(args[0])
        B  = float(args[1])
        bc = float(args[2])
    else:
        Ds = float(kwargs["Ds"])
        B  = float(kwargs["B"])
        bc = float(kwargs["bc"])

    # Optional extra geometry
    do = None
    Pt = None
    layout = "tri30"

    if len(args) >= 6:
        do = float(args[3])
        Pt = float(args[4])
        layout = str(args[5])
    else:
        do = kwargs.get("do", None)
        Pt = kwargs.get("Pt", None)
        layout = str(kwargs.get("layout", layout))

    # Base strip area (window reduction)
    strip = Ds * B * max(1 - bc, 0.2)

    # If we know tube/pitch, apply a porosity factor (helps match HTRI-like vmax)
    if (do is not None) and (Pt is not None):
        do = float(do); Pt = float(Pt)
        ratio = do / max(Pt, 1e-12)
        L = layout.lower()
        tri = ("tri" in L) or L.startswith("30")
        k = 0.55 if tri else 0.50
        phi = max(0.25, min(0.90, 1 - k * ratio))
        return strip * phi

    return strip
