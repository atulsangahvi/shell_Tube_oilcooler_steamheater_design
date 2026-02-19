import streamlit as st
import hmac
from io import BytesIO
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
from datetime import datetime

from shtx.units import mm_to_m
from shtx.props import water_props, brine_props, oil_props_from_d341, steam_sat_props
from shtx.solver import segmented_singlephase, steam_heater_simple
from shtx.geometry import tube_pitch
from shtx.bell_estimates import (
    estimate_leakage_fractions,
    estimate_Jl_from_clearances,
    estimate_Jb_from_bypass,
    estimate_Jc_from_baffle_cut,
)

# -----------------------------
# Password gate (no hint shown)
# -----------------------------
def require_password():
    # Put this in Streamlit Secrets (local .streamlit/secrets.toml OR Streamlit Cloud secrets):
    # APP_PASSWORD = "Snehaanju"
    pw = st.secrets.get("APP_PASSWORD", None)
    if not pw:
        st.error("App password is not configured. Add APP_PASSWORD in Streamlit Secrets.")
        st.stop()

    if "auth_ok" not in st.session_state:
        st.session_state.auth_ok = False

    if st.session_state.auth_ok:
        return

    st.subheader("Login")
    entered = st.text_input("Password", type="password")

    if not entered:
        st.stop()

    if hmac.compare_digest(entered, pw):
        st.session_state.auth_ok = True
        st.rerun()
    else:
        st.error("Incorrect password.")
        st.stop()


# -----------------------------
# PDF report builder
# -----------------------------
def build_pdf_report(title: str, sections: list[tuple[str, dict]]) -> bytes:
    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=A4)
    w, h = A4

    x0 = 18 * mm
    y = h - 18 * mm

    def line(txt, dy=6.0):
        nonlocal y
        c.drawString(x0, y, str(txt))
        y -= dy * mm
        if y < 18 * mm:
            c.showPage()
            y = h - 18 * mm

    c.setFont("Helvetica-Bold", 14)
    line(title, dy=8.0)

    c.setFont("Helvetica", 9)
    line(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", dy=6.0)
    line("-" * 90, dy=6.0)

    for sec_title, data in sections:
        c.setFont("Helvetica-Bold", 11)
        line(sec_title, dy=7.0)
        c.setFont("Helvetica", 9)

        if not data:
            line("(none)", dy=5.5)
        else:
            for k, v in data.items():
                line(f"{k}: {v}", dy=5.5)

        line("-" * 90, dy=6.0)

    c.save()
    return buf.getvalue()


# -----------------------------
# App
# -----------------------------
st.set_page_config(page_title="S&T HX (BD inputs)", layout="wide")
require_password()

st.title("Shell-and-Tube HX Designer — Segmented + Shell-side detailed inputs")

mode = st.selectbox(
    "Mode",
    [
        "Steam condensing → Water heating (simple)",
        "Single-phase segmented (water/glycol/oil)",
    ],
)

# -----------------------------
# Geometry inputs
# -----------------------------
st.subheader("Geometry")
g1, g2, g3, g4 = st.columns(4)
with g1:
    Ds = st.number_input("Shell ID (mm)", 150.00, 4000.00, 334.00, step=0.01, format="%.2f")
    do = st.number_input("Tube OD (mm)", 6.00, 60.00, 15.875, step=0.01, format="%.3f")
with g2:
    tw = st.number_input("Tube wall thickness (mm)", 0.20, 5.00, 1.245, step=0.01, format="%.3f")
    n = st.number_input("Number of tubes", 10, 50000, 134, step=1)
with g3:
    L = st.number_input("Tube effective length (m)", 0.50, 12.00, 2.34, step=0.01, format="%.2f")
    passes = st.number_input("Tube passes", 1, 12, 2, step=1)
with g4:
    bc = st.number_input("Baffle cut (%)", 15.00, 45.00, 25.00, step=0.01, format="%.2f")
    bs = st.number_input("Baffle spacing (m)", 0.05, 2.00, 0.23, step=0.01, format="%.2f")

st.subheader("Shell-side detailed geometry (optional; defaults are sensible starts)")
d1, d2, d3, d4 = st.columns(4)
with d1:
    layout = st.selectbox("Tube layout", ["30° triangular", "45° rotated square", "90° square"])
with d2:
    p_over_do = st.number_input("Tube pitch ratio p/do", 1.10, 1.80, 1.25, step=0.01, format="%.2f")
with d3:
    delta_tb_mm = st.number_input(
        "Tube-to-baffle dia. clearance δ_tb (mm)", 0.10, 3.00, 0.80, step=0.01, format="%.2f"
    )
with d4:
    delta_sb_mm = st.number_input(
        "Shell-to-baffle dia. clearance δ_sb (mm)", 0.50, 12.00, 4.00, step=0.01, format="%.2f"
    )

ss1, ss2 = st.columns(2)
with ss1:
    sealing_strips = st.number_input("Sealing strip pairs (integer)", 0, 20, 0, step=1)
with ss2:
    st.caption("More sealing strips ↓ bypass, ↑ Jb (closer to 1).")

pitch_m = tube_pitch(mm_to_m(do), float(p_over_do))
st.caption(f"Computed tube pitch ≈ {pitch_m * 1000:.2f} mm")

# -----------------------------
# Fouling + material
# -----------------------------
st.subheader("Fouling + Material")
f1, f2, f3 = st.columns(3)
with f1:
    Rfs = st.number_input("Shell fouling (m²·K/W)", 0.0, 0.01, 0.0, step=0.000001, format="%.6f")
with f2:
    Rft = st.number_input("Tube fouling (m²·K/W)", 0.0, 0.01, 0.000179, step=0.000001, format="%.6f")
with f3:
    kw = st.number_input("Tube wall k (W/m·K)", 8.0, 35.0, 16.0, step=0.1, format="%.1f")

# -----------------------------
# Bell–Delaware factors
# -----------------------------
st.subheader("Bell–Delaware correction factors")
auto = st.checkbox("Auto-estimate Jc/Jl/Jb from clearances & layout (recommended to start)", value=True)

if auto:
    Ds_m = mm_to_m(Ds)
    do_m = mm_to_m(do)
    As = Ds_m * float(bs) * max(1 - (float(bc) / 100.0), 0.2)

    leak_tmp = estimate_leakage_fractions(int(n), do_m, Ds_m, mm_to_m(delta_tb_mm), mm_to_m(delta_sb_mm))
    Jc_val = estimate_Jc_from_baffle_cut(float(bc) / 100.0)
    Jl_val = estimate_Jl_from_clearances(leak_tmp["A_leak"], As)
    Jb_val = estimate_Jb_from_bypass(float(p_over_do), int(sealing_strips))
    Jr_val = 1.00

    st.info(
        f"Auto J: Jc={Jc_val:.2f}, Jl={Jl_val:.2f}, Jb={Jb_val:.2f}, Jr={Jr_val:.2f}  |  "
        f"leakage ratio A_leak/A_cross={leak_tmp['A_leak'] / max(As, 1e-12):.3f}"
    )
else:
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        Jc_val = st.number_input("Jc", 0.2, 1.2, 0.90, step=0.01, format="%.2f")
    with c2:
        Jl_val = st.number_input("Jl", 0.2, 1.2, 0.80, step=0.01, format="%.2f")
    with c3:
        Jb_val = st.number_input("Jb", 0.2, 1.2, 0.90, step=0.01, format="%.2f")
    with c4:
        Jr_val = st.number_input("Jr", 0.2, 1.2, 1.00, step=0.01, format="%.2f")

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
    delta_tb_m=mm_to_m(delta_tb_mm),
    delta_sb_m=mm_to_m(delta_sb_mm),
)

foul = dict(Rf_shell=float(Rfs), Rf_tube=float(Rft))
J = dict(Jc=float(Jc_val), Jl=float(Jl_val), Jb=float(Jb_val), Jr=float(Jr_val))

# -----------------------------
# Modes
# -----------------------------
if mode.startswith("Single-phase"):
    st.subheader("Flow model selector (HTC correlations unchanged)")
    flow_choice = st.selectbox(
        "Flow model",
        [
            "1-shell-pass / 2+ tube passes (ε–NTU 1-shell)",
            "Crossflow approximation (ε–NTU crossflow)",
        ],
    )
    flow_model = "1-shell" if flow_choice.startswith("1-shell") else "crossflow"

    N = st.number_input("Segments", 5, 200, 20, step=1)

    a, b = st.columns(2)
    with a:
        hot = st.selectbox("Shell-side (hot) fluid", ["Water", "MEG", "MPG", "Oil (custom)"])
        Th = st.number_input("Hot inlet (°C)", -30.0, 350.0, 80.0, step=0.01, format="%.2f")
        mh = st.number_input("Hot mass flow (kg/s)", 0.01, 20000.0, 5.0, step=0.01, format="%.2f")
        hot_pct = None
        nu40 = nu100 = rho15 = None
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
        cold_pct = None
        if cold in ["MEG", "MPG"]:
            cold_pct = st.number_input("Cold glycol (mass %)", 0.00, 80.00, 30.00, step=0.01, format="%.2f")

    def pf(ftype, pct=None):
        if ftype == "Water":
            return lambda T: water_props(T, 300.0)
        if ftype in ["MEG", "MPG"]:
            return lambda T: brine_props(ftype, (pct or 0.0) / 100.0, T)
        return lambda T: oil_props_from_d341(T, float(nu40), float(nu100), float(rho15))

    p_shell = pf(hot, hot_pct)
    p_tube = pf(cold, cold_pct)

    if st.button("Run segmented rating"):
        res = segmented_singlephase(
            p_shell,
            p_tube,
            float(Th),
            float(Tc),
            float(mh),
            float(mc),
            geom,
            foul,
            J,
            flow_model,
            int(N),
            float(kw),
        )

        r1, r2, r3, r4 = st.columns(4)
        r1.metric("Duty (MW)", f"{res.duty_W / 1e6:.3f}")
        r2.metric("U_avg (W/m²·K)", f"{res.U_avg:.1f}")
        r3.metric("Hot out (°C)", f"{res.Th_out:.2f}")
        r4.metric("Cold out (°C)", f"{res.Tc_out:.2f}")

        d1c, d2c = st.columns(2)
        d1c.metric("Shell ΔP (kPa)", f"{res.dp_shell_kpa:.2f}")
        d2c.metric("Tube ΔP (kPa)", f"{res.dp_tube_kpa:.2f}")

        st.subheader("J-factors used + leakage estimate")
        st.json(res.J_used)
        st.json(res.leak)

        # PDF
        inputs_dict = {
            "Mode": "Single-phase segmented",
            "Flow model": flow_model,
            "Segments": int(N),
            "Shell ID (mm)": float(Ds),
            "Tube OD (mm)": float(do),
            "Tube wall (mm)": float(tw),
            "Tube length (m)": float(L),
            "Tube passes": int(passes),
            "No. tubes": int(n),
            "Baffle cut (%)": float(bc),
            "Baffle spacing (m)": float(bs),
            "Tube layout": layout,
            "Pitch ratio p/do": float(p_over_do),
            "δ_tb (mm)": float(delta_tb_mm),
            "δ_sb (mm)": float(delta_sb_mm),
            "Sealing strip pairs": int(sealing_strips),
            "Rf_shell (m²K/W)": float(Rfs),
            "Rf_tube (m²K/W)": float(Rft),
            "Tube wall k (W/mK)": float(kw),
            "Shell hot fluid": hot,
            "Hot inlet (°C)": float(Th),
            "Hot m_dot (kg/s)": float(mh),
            "Tube cold fluid": cold,
            "Cold inlet (°C)": float(Tc),
            "Cold m_dot (kg/s)": float(mc),
            "Jc": float(J["Jc"]),
            "Jl": float(J["Jl"]),
            "Jb": float(J["Jb"]),
            "Jr": float(J["Jr"]),
        }
        if hot in ["MEG", "MPG"]:
            inputs_dict["Hot glycol mass %"] = float(hot_pct)
        if cold in ["MEG", "MPG"]:
            inputs_dict["Cold glycol mass %"] = float(cold_pct)
        if hot.startswith("Oil"):
            inputs_dict["Oil ν40 (cSt)"] = float(nu40)
            inputs_dict["Oil ν100 (cSt)"] = float(nu100)
            inputs_dict["Oil ρ15 (kg/m³)"] = float(rho15)

        results_dict = {
            "Duty (MW)": f"{res.duty_W / 1e6:.6f}",
            "U_avg (W/m²·K)": f"{res.U_avg:.3f}",
            "Hot outlet (°C)": f"{res.Th_out:.3f}",
            "Cold outlet (°C)": f"{res.Tc_out:.3f}",
            "Shell ΔP (kPa)": f"{res.dp_shell_kpa:.3f}",
            "Tube ΔP (kPa)": f"{res.dp_tube_kpa:.3f}",
            "Area (m²)": f"{res.area_m2:.3f}",
        }

        leak_dict = res.leak if isinstance(res.leak, dict) else {}

        pdf_bytes = build_pdf_report(
            title="Shell-and-Tube HX Report",
            sections=[
                ("Inputs", inputs_dict),
                ("Results", results_dict),
                ("Leakage / Geometry Estimate", leak_dict),
            ],
        )

        st.download_button(
            label="Download PDF report",
            data=pdf_bytes,
            file_name="shtx_report.pdf",
            mime="application/pdf",
        )

else:
    st.subheader("Steam heater inputs")
    c1, c2 = st.columns(2)
    with c1:
        m = st.number_input("Water mass flow (kg/s)", 0.1, 20000.0, 40.0, step=0.01, format="%.2f")
        t = st.number_input("Water inlet (°C)", -10.0, 200.0, 60.0, step=0.01, format="%.2f")
    with c2:
        p = st.number_input("Steam pressure (kPa abs)", 50.0, 10000.0, 2591.0, step=0.1, format="%.1f")

    if st.button("Run steam heater rating"):
        out = steam_heater_simple(
            lambda T: water_props(T, 300.0),
            float(t),
            float(m),
            float(p),
            geom,
            foul,
            float(kw),
        )

        r1, r2, r3, r4 = st.columns(4)
        r1.metric("Duty (MW)", f"{out['Q'] / 1e6:.3f}")
        r2.metric("U (W/m²·K)", f"{out['U']:.1f}")
        r3.metric("Water out (°C)", f"{out['tco']:.2f}")
        r4.metric("Steam Tsat (°C)", f"{out['Ts']:.2f}")

        inputs_dict = {
            "Mode": "Steam condensing → Water heating (simple)",
            "Shell ID (mm)": float(Ds),
            "Tube OD (mm)": float(do),
            "Tube wall (mm)": float(tw),
            "Tube length (m)": float(L),
            "Tube passes": int(passes),
            "No. tubes": int(n),
            "Baffle cut (%)": float(bc),
            "Baffle spacing (m)": float(bs),
            "Rf_shell (m²K/W)": float(Rfs),
            "Rf_tube (m²K/W)": float(Rft),
            "Tube wall k (W/mK)": float(kw),
            "Water inlet (°C)": float(t),
            "Water m_dot (kg/s)": float(m),
            "Steam pressure (kPa abs)": float(p),
        }

        results_dict = {
            "Duty (MW)": f"{out['Q'] / 1e6:.6f}",
            "U (W/m²·K)": f"{out['U']:.3f}",
            "Water outlet (°C)": f"{out['tco']:.3f}",
            "Steam Tsat (°C)": f"{out['Ts']:.3f}",
            "Area (m²)": f"{out['A']:.3f}",
            "h_shell (W/m²·K)": f"{out['hs']:.1f}",
            "h_tube (W/m²·K)": f"{out['ht']:.1f}",
            "Re_tube": f"{out['Re']:.0f}",
            "v_tube (m/s)": f"{out['vt']:.3f}",
        }

        pdf_bytes = build_pdf_report(
            title="Shell-and-Tube HX Report",
            sections=[
                ("Inputs", inputs_dict),
                ("Results", results_dict),
            ],
        )

        st.download_button(
            label="Download PDF report",
            data=pdf_bytes,
            file_name="shtx_steam_heater_report.pdf",
            mime="application/pdf",
        )
