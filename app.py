import streamlit as st
from shtx.units import mm_to_m
from shtx.props import water_props, brine_props, oil_props_from_d341, steam_sat_props
from shtx.solver import segmented_singlephase, steam_heater_simple
from shtx.geometry import tube_pitch
from shtx.bell_estimates import (
    estimate_leakage_fractions,
    estimate_Jl_from_clearances,
    estimate_Jb_from_bypass,
    estimate_Jc_from_baffle_cut
)

st.set_page_config(page_title="S&T HX (BD inputs)", layout="wide")
st.title("Shell-and-Tube HX Designer — Segmented + Shell-side detailed inputs")

mode = st.selectbox("Mode", [
    "Steam condensing → Water heating (simple)",
    "Single-phase segmented (water/glycol/oil)"
])


st.sidebar.subheader("HTRI benchmark (optional)")
htri_file = st.sidebar.file_uploader("Upload HTRI/Xist .xls for % error check", type=["xls"])
htri_vals = None
if htri_file is not None:
    # save temp
    tmp_path = "/tmp/htri_upload.xls"
    with open(tmp_path, "wb") as f:
        f.write(htri_file.getbuffer())
    try:
        htri_vals = parse_htri_key_values(tmp_path)
        st.sidebar.success("HTRI file read OK.")
        st.sidebar.write({k: htri_vals[k] for k in ['Q_MW','U_actual','Th_in','Th_out','Tc_in','Tc_out','dp_shell_kpa','dp_tube_kpa'] if k in htri_vals})
    except Exception as e:
        st.sidebar.error(f"Could not parse HTRI file: {e}")
        htri_vals = None


st.subheader("Geometry")
g1, g2, g3, g4 = st.columns(4)
with g1:
    Ds = st.number_input("Shell ID (mm)", 150.00, 4000.00, 334.00, step=0.01, format="%.2f")
    do = st.number_input("Tube OD (mm)", 6.00, 60.00, 15.875, step=0.01, format="%.3f")
with g2:
    tw = st.number_input("Tube wall thickness (mm)", 0.20, 5.00, 1.245, step=0.01, format="%.3f")
    n  = st.number_input("Number of tubes", 10, 50000, 134, step=1)
with g3:
    L  = st.number_input("Tube effective length (m)", 0.50, 12.00, 2.34, step=0.01, format="%.2f")
    passes = st.number_input("Tube passes", 1, 12, 2, step=1)
with g4:
    bc = st.number_input("Baffle cut (%)", 15.00, 45.00, 25.00, step=0.01, format="%.2f")
    bs = st.number_input("Baffle spacing (m)", 0.05, 2.00, 0.23, step=0.01, format="%.2f")

st.subheader("Shell-side detailed geometry (optional; defaults are sensible starts)")
d1, d2, d3, d4 = st.columns(4)
with d1:
    layout_ui = st.selectbox("Tube layout", ["30° triangular", "45° rotated square", "90° square"])
layout = "tri30" if layout_ui.startswith("30") else "sq90"
with d2:
    p_over_do = st.number_input("Tube pitch ratio p/do", 1.10, 1.80, 1.25, step=0.01, format="%.2f")
with d3:
    delta_tb_mm = st.number_input("Tube-to-baffle dia. clearance δ_tb (mm)", 0.10, 3.00, 0.80, step=0.01, format="%.2f")
with d4:
    delta_sb_mm = st.number_input("Shell-to-baffle dia. clearance δ_sb (mm)", 0.50, 12.00, 4.00, step=0.01, format="%.2f")

ss1, ss2 = st.columns(2)
with ss1:
    sealing_strips = st.number_input("Sealing strip pairs (integer)", 0, 20, 0, step=1)
with ss2:
    st.caption("More sealing strips ↓ bypass, ↑ Jb (closer to 1).")

pitch_m = tube_pitch(mm_to_m(do), float(p_over_do))
st.caption(f"Computed tube pitch ≈ {pitch_m*1000:.2f} mm")

st.subheader("Fouling + Material")
f1, f2, f3 = st.columns(3)
with f1:
    Rfs = st.number_input("Shell fouling (m²·K/W)", 0.0, 0.01, 0.0, step=0.000001, format="%.6f")
with f2:
    Rft = st.number_input("Tube fouling (m²·K/W)", 0.0, 0.01, 0.000179, step=0.000001, format="%.6f")
with f3:
    kw  = st.number_input("Tube wall k (W/m·K)", 8.0, 35.0, 16.0, step=0.1, format="%.1f")

st.subheader("Bell–Delaware correction factors")
auto = st.checkbox("Auto-estimate Jc/Jl/Jb from clearances & layout (recommended to start)", value=True)

if auto:
    # Estimate Jc from baffle cut, Jb from pitch+sealing strips, Jl from leakage areas
    Ds_m = mm_to_m(Ds)
    do_m = mm_to_m(do)
    As = Ds_m * float(bs) * max(1 - (float(bc)/100.0), 0.2)  # same as approx_shell_crossflow_area
    leak = estimate_leakage_fractions(int(n), do_m, Ds_m, mm_to_m(delta_tb_mm), mm_to_m(delta_sb_mm))
    Jc_val = estimate_Jc_from_baffle_cut(float(bc)/100.0)
    Jl_val = estimate_Jl_from_clearances(leak["A_leak"], As)
    Jb_val = estimate_Jb_from_bypass(float(p_over_do), int(sealing_strips))
    Jr_val = 1.00
    st.info(f"Auto J: Jc={Jc_val:.2f}, Jl={Jl_val:.2f}, Jb={Jb_val:.2f}, Jr={Jr_val:.2f}  |  leakage ratio A_leak/A_cross={leak['A_leak']/max(As,1e-12):.3f}")
else:
    c1, c2, c3, c4 = st.columns(4)
    with c1: Jc_val = st.number_input("Jc", 0.2, 1.2, 0.90, step=0.01, format="%.2f")
    with c2: Jl_val = st.number_input("Jl", 0.2, 1.2, 0.80, step=0.01, format="%.2f")
    with c3: Jb_val = st.number_input("Jb", 0.2, 1.2, 0.90, step=0.01, format="%.2f")
    with c4: Jr_val = st.number_input("Jr", 0.2, 1.2, 1.00, step=0.01, format="%.2f")

geom = dict(
    Ds=mm_to_m(Ds),
    do=mm_to_m(do),
    di=mm_to_m(do - 2 * tw),
    tw=mm_to_m(tw),
    L=float(L),
    n=int(n),
    passes=int(passes),
    bs=float(bs),
    bc=float(bc) / 100.0,
    steam_sat=steam_sat_props,
    # store clearances for reporting
    delta_tb_m=mm_to_m(delta_tb_mm),
    delta_sb_m=mm_to_m(delta_sb_mm),
    layout=layout,
    p_over_do=float(p_over_do),
)

foul = dict(Rf_shell=float(Rfs), Rf_tube=float(Rft))
J = dict(Jc=float(Jc_val), Jl=float(Jl_val), Jb=float(Jb_val), Jr=float(Jr_val))

if mode.startswith("Single-phase"):
    st.subheader("Flow model selector (HTC correlations unchanged)")
    flow_choice = st.selectbox("Flow model", [
        "1-shell-pass / 2+ tube passes (ε–NTU 1-shell)",
        "Crossflow approximation (ε–NTU crossflow)",
    ])
    flow_model = "1-shell" if flow_choice.startswith("1-shell") else "crossflow"

    N = st.number_input("Segments", 5, 200, 20, step=1)

    a, b = st.columns(2)
    with a:
        hot = st.selectbox("Shell-side (hot) fluid", ["Water", "MEG", "MPG", "Oil (custom)"])
        Th = st.number_input("Hot inlet (°C)", -30.0, 350.0, 80.0, step=0.01, format="%.2f")
        mh = st.number_input("Hot mass flow (kg/s)", 0.01, 20000.0, 5.0, step=0.01, format="%.2f")
        if hot in ["MEG", "MPG"]:
            hot_pct = st.number_input("Hot glycol (mass %)", 0.00, 80.00, 30.00, step=0.01, format="%.2f")
        if hot.startswith("Oil"):
            nu40 = st.number_input("Oil ν40 (cSt)", 0.5, 5000.0, 46.0, step=0.01, format="%.2f")
            nu100 = st.number_input("Oil ν100 (cSt)", 0.2, 2000.0, 6.8, step=0.01, format="%.2f")
            rho15 = st.number_input("Oil ρ15 (kg/m³)", 600.0, 1100.0, 870.0, step=0.1, format="%.1f")
    with b:
        cold = st.selectbox("Tube-side (cold) fluid", ["Water", "MEG", "MPG"])
        Tc = st.number_input("Cold inlet (°C)", -30.0, 200.0, 7.0, step=0.01, format="%.2f")
        mc = st.number_input("Cold mass flow (kg/s)", 0.01, 20000.0, 8.0, step=0.01, format="%.2f")
        if cold in ["MEG", "MPG"]:
            cold_pct = st.number_input("Cold glycol (mass %)", 0.00, 80.00, 30.00, step=0.01, format="%.2f")

    def pf(ftype, pct=None):
        if ftype == "Water":
            return lambda T: water_props(T, 300.0)
        if ftype in ["MEG", "MPG"]:
            return lambda T: brine_props(ftype, (pct or 0.0) / 100.0, T)
        return lambda T: oil_props_from_d341(T, nu40, nu100, rho15)

    p_shell = pf(hot, locals().get("hot_pct", None))
    p_tube  = pf(cold, locals().get("cold_pct", None))

    if st.button("Run segmented rating"):
        res = segmented_singlephase(
            p_shell, p_tube,
            float(Th), float(Tc),
            float(mh), float(mc),
            geom, foul, J,
            flow_model, int(N), float(kw),
        )

        r1, r2, r3, r4 = st.columns(4)
        r1.metric("Duty (MW)", f"{res.duty_W/1e6:.3f}")
        r2.metric("U_avg (W/m²·K)", f"{res.U_avg:.1f}")
        r3.metric("Hot out (°C)", f"{res.Th_out:.2f}")
        r4.metric("Cold out (°C)", f"{res.Tc_out:.2f}")

        d1, d2 = st.columns(2)
        d1.metric("Shell ΔP (kPa)", f"{res.dp_shell_kpa:.2f}")
        d2.metric("Tube ΔP (kPa)", f"{res.dp_tube_kpa:.2f}")

        
if htri_vals is not None:
    st.subheader("Comparison vs HTRI (percent error)")
    comp = {
        "Duty (MW)": (res.duty_W/1e6, htri_vals.get("Q_MW")),
        "U (W/m²·K)": (res.U_avg, htri_vals.get("U_actual")),
        "Hot out (°C)": (res.Th_out, htri_vals.get("Th_out")),
        "Cold out (°C)": (res.Tc_out, htri_vals.get("Tc_out")),
        "Shell ΔP (kPa)": (res.dp_shell_kpa, htri_vals.get("dp_shell_kpa")),
        "Tube ΔP (kPa)": (res.dp_tube_kpa, htri_vals.get("dp_tube_kpa")),
    }
    rows=[]
    for k,(m,ref) in comp.items():
        pe = pct_err(m, ref)
        rows.append([k, m, ref, pe])
    dfc = pd.DataFrame(rows, columns=["Metric","Model","HTRI","% error"])
    st.dataframe(dfc, use_container_width=True)

        st.subheader("J-factors used + leakage estimate")
        st.json(res.J_used)
        st.json(res.leak)
        st.subheader("Debug (mid-segment)")
        st.json(res.debug)
else:
    m = st.number_input("Water mass flow (kg/s)", 0.1, 20000.0, 40.0, step=0.01, format="%.2f")
    t = st.number_input("Water inlet (°C)", -10.0, 200.0, 60.0, step=0.01, format="%.2f")
    p = st.number_input("Steam pressure (kPa abs)", 50.0, 10000.0, 2591.0, step=0.1, format="%.1f")

    if st.button("Run steam heater rating"):
        out = steam_heater_simple(lambda T: water_props(T, 300.0), float(t), float(m), float(p), geom, foul, float(kw))
        r1, r2, r3, r4 = st.columns(4)
        r1.metric("Duty (MW)", f"{out['Q']/1e6:.3f}")
        r2.metric("U (W/m²·K)", f"{out['U']:.1f}")
        r3.metric("Water out (°C)", f"{out['tco']:.2f}")
        r4.metric("Steam Tsat (°C)", f"{out['Ts']:.2f}")
