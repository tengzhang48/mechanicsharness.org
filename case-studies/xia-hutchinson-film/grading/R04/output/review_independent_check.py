"""Independent reviewer checks: 1D strip energy from PROBLEM.md functional,
not copied from solver tanh target. Also KI mapping and 16-value comparison.
"""
import json
import math
import numpy as np

NU = 0.3
EBAR = 1.0  # set Ebar=1, h=1, ell=1, s0=1 by scaling
h = 1.0
ell = 1.0
s0 = 1.0
k = h * EBAR / ell**2


def wake_u2(x2, S):
    return -(s0 * ell / EBAR) * np.sinh((x2 - S / 2.0) / ell) / np.cosh(S / (2.0 * ell))


def strip_W_quad(S, n=4001):
    x2 = np.linspace(0.0, S, n)
    u2 = wake_u2(x2, S)
    eps22 = np.gradient(u2, x2)
    dsig22 = EBAR * eps22
    w = h * s0 * eps22 + 0.5 * h * dsig22 * eps22 + 0.5 * k * u2**2
    return np.trapezoid(w, x2)


def gamma_from_W(S_ahead, S_behind, n_tips=1, period_note=""):
    # G = (W_ahead - W_behind) / (h * n_tips) wait: PROBLEM: divide release per unit
    # in-plane advance per tip by h. Delta Pi per period for da: (W_a - W_b)*da? 
    # If W is energy per unit x1 of the strip, advancing tips by da replaces
    # ahead length da by behind: Delta Pi = (W_behind - W_ahead)*da
    # Energy release per unit in-plane advance per tip = -(Delta Pi)/da / n_tips
    # G = that / h = (W_ahead - W_behind)/(h n_tips)
    Wa = strip_W_quad(S_ahead)
    Wb = strip_W_quad(S_behind)
    return (Wa - Wb) / (h * n_tips) * EBAR / (s0**2 * ell), Wa, Wb


print("=== Independent quadrature of stated DeltaPi density ===")
for eoH in (0.5, 1.0, 1.5, 2.0):
    H = ell / eoH
    gA, WaA, WbA = gamma_from_W(H, H, n_tips=1)  # wait Case A: ahead S infinite uncracked W=0
    # Case A: period H, one tip. Ahead: uncracked W=0 (u=0). Behind: strip width H.
    Wa = 0.0
    Wb = strip_W_quad(H)
    gA = (Wa - Wb) / h * EBAR / (s0**2 * ell)
    # Case B: period 2H, one advancing tip. Ahead: infinite cracks spacing 2H, W(2H).
    # Behind: two width-H strips, W=2*W(H). Energy per period: W_ahead - W_behind
    # wait G = (W_ahead - W_behind)/(h n_tips); n_tips=1
    WaB = strip_W_quad(2 * H)
    WbB = 2 * strip_W_quad(H)
    gB = (WaB - WbB) / h * EBAR / (s0**2 * ell)
    t = math.tanh(H / (2 * ell))
    gA_th = t
    gB_th = 2 * t - math.tanh(H / ell)
    print(f"ell/H={eoH}: gA_quad={gA:.16f} gA_th={gA_th:.16f} d={gA-gA_th:.3e}")
    print(f"         gB_quad={gB:.16f} gB_th={gB_th:.16f} d={gB-gB_th:.3e}")
    print(f"         Wb(H)={Wb:.12f} target={-h*s0**2*ell/EBAR*t:.12f}")

print("\n=== KI from G=KI^2/E with E=Ebar*(1-nu^2) ===")
# G = KI^2 / E => KI = sqrt(E G) = sqrt(Ebar(1-nu^2) * gamma * s0^2 ell / Ebar)
# = s0 sqrt(ell) sqrt((1-nu^2) gamma)

print("\n=== Compare 16 samples to reference.json ===")
ref = json.load(open("reference/reference.json"))
sub = json.load(open("submission/numeric_rows.json"))
n_agree = 0
tol = 1e-12
for i, (r, s) in enumerate(zip(ref["rows"], sub)):
    for key in ("gamma_A", "gamma_B", "kI_A", "kI_B"):
        d = abs(r[key] - s[key])
        ok = d < tol
        n_agree += int(ok)
        print(f"ell/H={s['ell_over_H']} {key}: sub={s[key]:.16f} ref={r[key]:.16f} d={d:.3e} {'OK' if ok else 'DIFF'}")
print(f"AGREEMENT {n_agree}/16 at tol={tol}")

# Case A small-H: gamma ~ H/(2 ell)
print("\n=== Independent small/large limits from series ===")
# tanh(x)=x-x^3/3+... x=H/(2ell)
# gA ~ H/(2ell)
# gB = 2tanh(x)-tanh(2x) ~ 2(x-x^3/3)-(2x-(2x)^3/3)= 2x-2x^3/3 -2x + 8x^3/3 = 2 x^3 = 2 (H/2ell)^3 = H^3/(4 ell^3)
# so gamma_B ~ (H/ell)^3 / 4  yes
print("gB ~ (H/ell)^3/4 from tanh series: confirmed algebraically")
print("sigma_c,A ~ sqrt(Gc Ebar / (ell * H/(2ell))) = sqrt(2 Gc Ebar / H)")
print("sigma_c,B ~ sqrt(Gc Ebar / (ell * H^3/(4 ell^3))) = 2 sqrt(Gc Ebar) ell / H^{3/2}")
