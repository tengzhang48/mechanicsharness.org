"""Re-run previous breaking attempts and probe new v0.2 checks."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "crack_tip_asymptotics" / "checks"))
from crack_tip_checks import check_field, williams_exponent, format_report


def failed(rep):
    return {r["check"] for r in rep if r["passed"] is False}


def by_name(rep):
    return {r["check"]: r for r in rep}


def dump(title, rep):
    print("\n====", title, "====")
    print("FAIL:", failed(rep) or "(none)")
    for r in rep:
        flag = {True: "PASS", False: "FAIL", None: "n/a "}[r["passed"]]
        v = r["value"]
        extra = f" {v:.6g}" if isinstance(v, float) and v == v else ""
        print(f"  {flag} {r['check']}{extra}")
        if r.get("note"):
            print(f"       note: {r['note']}")


def williams_crack(n_theta=4001, sign=1.0):
    th = np.linspace(-np.pi, np.pi, n_theta)
    k = sign / np.sqrt(2 * np.pi)
    return (th, k * (1.25 * np.cos(th / 2) - 0.25 * np.cos(1.5 * th)),
            k * (0.75 * np.cos(th / 2) + 0.25 * np.cos(1.5 * th)),
            k * (0.25 * np.sin(th / 2) + 0.25 * np.sin(1.5 * th)))


def williams_mode_II_crack(n_theta=4001):
    th = np.linspace(-np.pi, np.pi, n_theta)
    k = 1 / np.sqrt(2 * np.pi)
    srr = k * (-1.25 * np.sin(th / 2) + 0.75 * np.sin(1.5 * th))
    stt = k * (-0.75 * np.sin(th / 2) - 0.75 * np.sin(1.5 * th))
    srt = k * (0.25 * np.cos(th / 2) + 0.75 * np.cos(1.5 * th))
    return th, srr, stt, srt


def t_stress(n_theta=4001):
    th = np.linspace(-np.pi, np.pi, n_theta)
    return th, np.cos(th) ** 2, np.sin(th) ** 2, -np.sin(th) * np.cos(th)


def hydrostatic(n_theta=4001):
    th = np.linspace(-np.pi, np.pi, n_theta)
    return th, np.ones_like(th), np.ones_like(th), np.zeros_like(th)


def williams_notch_I(beta, n_theta=4001):
    p = williams_exponent(beta, "I")
    lam = 1 - p
    s = lam + 1
    th = np.linspace(0, beta, n_theta)
    B = -np.cos((lam + 1) * beta) / np.cos((lam - 1) * beta)
    F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
    Fp = -(lam + 1) * np.sin((lam + 1) * th) - B * (lam - 1) * np.sin((lam - 1) * th)
    Fpp = -(lam + 1) ** 2 * np.cos((lam + 1) * th) - B * (lam - 1) ** 2 * np.cos((lam - 1) * th)
    return th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p


def williams_notch_II(beta, n_theta=4001):
    p = williams_exponent(beta, "II")
    lam = 1 - p
    s = lam + 1
    th = np.linspace(0, beta, n_theta)
    B = -np.sin((lam + 1) * beta) / np.sin((lam - 1) * beta)
    F = np.sin((lam + 1) * th) + B * np.sin((lam - 1) * th)
    Fp = (lam + 1) * np.cos((lam + 1) * th) + B * (lam - 1) * np.cos((lam - 1) * th)
    Fpp = -(lam + 1) ** 2 * np.sin((lam + 1) * th) - B * (lam - 1) ** 2 * np.sin((lam - 1) * th)
    return th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p


def fake_airy_cos5(n_theta=4001):
    th = np.linspace(-np.pi, np.pi, n_theta)
    F = np.cos(th / 2) ** 5
    p = 0.5
    s = 2 - p
    Fp = CubicSpline_safe(th, F)
    return th, Fp[2] + s * F, s * (s - 1) * F, -(s - 1) * Fp[1]


def CubicSpline_safe(th, F):
    from scipy.interpolate import CubicSpline
    sp = CubicSpline(th, F)
    return None, sp(th, 1), sp(th, 2)


print("williams_exponent pi/2 I", williams_exponent(np.pi / 2, "I"))
print("williams_exponent 2pi/3 I", williams_exponent(2 * np.pi / 3, "I"), "p", williams_exponent(2 * np.pi / 3, "I"))
print("williams_exponent 2pi/3 II", williams_exponent(2 * np.pi / 3, "II"))
print("williams_exponent 5pi/6 II", williams_exponent(5 * np.pi / 6, "II"))
print("williams_exponent 3pi/4 I", williams_exponent(3 * np.pi / 4, "I"))
print("williams_exponent 5pi/6 I", williams_exponent(5 * np.pi / 6, "I"))
print("williams_exponent 270deg I", williams_exponent(np.radians(270) / 2, "I"))
print("williams_exponent 270deg II", williams_exponent(np.radians(270) / 2, "II"))
print("williams_exponent 240 I", williams_exponent(np.radians(240) / 2, "I"))
print("williams_exponent 240 II", williams_exponent(np.radians(240) / 2, "II"))
print("williams_exponent 180 I", williams_exponent(np.pi / 2, "I"))
print("williams_exponent 300 I", williams_exponent(np.radians(300) / 2, "I"))
print("williams_exponent 300 II", williams_exponent(np.radians(300) / 2, "II"))

from scipy.optimize import brentq
crit = brentq(lambda x: np.tan(x) - x, 4.0, 4.6)
print("crit 2beta rad", crit, "deg", np.degrees(crit), "beta deg", np.degrees(crit) / 2)
try:
    print("p II at crit", williams_exponent(crit / 2, "II"))
except Exception as e:
    print("p II at crit ERROR", type(e).__name__, e)
print("p I 257.45", williams_exponent(np.radians(257.45) / 2, "I"))
print("p II 257.45", williams_exponent(np.radians(257.45) / 2, "II"))

# T-stress leading vs not
dump("T-stress leading=False", check_field(*t_stress(), p=0.0, law="linear", leading=False))
dump("T-stress leading=True", check_field(*t_stress(), p=0.0, law="linear", leading=True))

dump("mode II crack", check_field(*williams_mode_II_crack(), p=0.5, law="linear", mode="II"))
dump("Williams KI", check_field(*williams_crack(), p=0.5, law="linear", nu=0.3, plane="strain"))
dump("closing crack", check_field(*williams_crack(sign=-1), p=0.5, law="linear"))

# 90 deg wedge
th = np.linspace(-np.pi / 2, np.pi / 2, 4001)
dump("90 wedge uniform-along-edge field",
     check_field(th, 4 * np.sin(th) ** 2, 4 * np.cos(th) ** 2, 4 * np.sin(th) * np.cos(th),
                 p=0.0, law="linear", beta=np.pi / 2))

# p=1 crash
dump("p=1 Williams stresses", check_field(*williams_crack(), p=1.0, law="linear"))

# hydrostatic
dump("hydrostatic p=0 leading=False", check_field(*hydrostatic(), p=0.0, law="linear", leading=False))

# coarse grid
dump("Williams 21 pts full", check_field(*williams_crack(21), p=0.5, law="linear"))
th, a, b, c = williams_crack(21)
keep = th >= 0
dump("Williams 21 pts half-range", check_field(th[keep], a[keep], b[keep], c[keep], p=0.5, law="linear"))

# higher-order Williams lambda=1.5 p=-0.5
th = np.linspace(-np.pi, np.pi, 4001)
lam = 1.5
# mode I Airy ~ r^{lam+1} [cos((lam+1)th)+B cos((lam-1)th)] with beta=pi
# at beta=pi, B = -cos((lam+1)pi)/cos((lam-1)pi)
B = -np.cos((lam + 1) * np.pi) / np.cos((lam - 1) * np.pi)
s = lam + 1
F = np.cos((lam + 1) * th) + B * np.cos((lam - 1) * th)
Fp = -(lam + 1) * np.sin((lam + 1) * th) - B * (lam - 1) * np.sin((lam - 1) * th)
Fpp = -(lam + 1) ** 2 * np.cos((lam + 1) * th) - B * (lam - 1) ** 2 * np.cos((lam - 1) * th)
p = 1 - lam
dump("higher-order Williams p=-0.5 leading=True",
     check_field(th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p=p, law="linear", leading=True))
dump("higher-order Williams p=-0.5 leading=False",
     check_field(th, Fpp + s * F, s * (s - 1) * F, -(s - 1) * Fp, p=p, law="linear", leading=False))

# 120 deg notch Delta p
th, a, b, c, p = williams_notch_I(2 * np.pi / 3)
dump("120 notch exact", check_field(th, a, b, c, p=p, law="linear", beta=2 * np.pi / 3))
dump("120 notch dp=1e-4", check_field(th, a, b, c, p=p + 1e-4, law="linear", beta=2 * np.pi / 3))
dump("120 notch dp=1e-5", check_field(th, a, b, c, p=p + 1e-5, law="linear", beta=2 * np.pi / 3))

# 150 mode II
th, a, b, c, p = williams_notch_II(5 * np.pi / 6)
dump("150 mode II exact", check_field(th, a, b, c, p=p, law="linear", beta=5 * np.pi / 6, mode="II"))

# T mixed into Williams
th, a, b, c = williams_crack()
tt, ta, tb, tc = t_stress()
for eps in (1e-4, 1e-3, 1e-2):
    dump(f"Williams+eps T eps={eps}",
         check_field(th, a + eps * ta, b + eps * tb, c + eps * tc, p=0.5, law="linear"))

# F=cos^5
th, a, b, c = fake_airy_cos5()
dump("fake Airy cos^5 p=0.5", check_field(th, a, b, c, p=0.5, law="linear"))

# F=(1+cos th)^2
th = np.linspace(-np.pi, np.pi, 4001)
F = (1 + np.cos(th)) ** 2
p = 0.5
s = 2 - p
from scipy.interpolate import CubicSpline
sp = CubicSpline(th, F)
dump("fake Airy (1+cos)^2",
     check_field(th, sp(th, 2) + s * F, s * (s - 1) * F, -(s - 1) * sp(th, 1), p=0.5, law="linear"))

# leading switch: T with leading True should report non-singular
# already dumped

# wrong Williams mode II exponent on 150
th, a, b, c, p = williams_notch_II(5 * np.pi / 6)
dump("150 II with p=0.5", check_field(th, a, b, c, p=0.5, law="linear", beta=5 * np.pi / 6, mode="II"))
dump("150 II with p=p+1e-5", check_field(th, a, b, c, p=p + 1e-5, law="linear", beta=5 * np.pi / 6, mode="II"))

# 120 II should be non-singular
print("120 mode II p", williams_exponent(2 * np.pi / 3, "II"))

# 90 mode II
try:
    print("90 mode II p", williams_exponent(np.pi / 2, "II"))
except Exception as e:
    print("90 mode II error", type(e), e)

# mode II 240
print("240 mode II p", williams_exponent(np.radians(240) / 2, "II"))

# reported_J on notch
th, a, b, c, p = williams_notch_I(2 * np.pi / 3)
dump("notch reported_J=1e9", check_field(th, a, b, c, p=p, law="linear", beta=2 * np.pi / 3, reported_J=1e9))

# normalization
rep = check_field(*williams_crack(), p=0.5, law="linear")
print("\nnormalization note:", by_name(rep)["normalization (report only)"]["note"])
print("max Se value", by_name(rep)["normalization (report only)"]["value"])

# table values
print("\n=== existence table check ===")
angles = [180, 240, 257.45, 270, 300, 360]
for ang in angles:
    beta = np.radians(ang) / 2
    pi_ = williams_exponent(beta, "I") if ang != 360 else 0.5
    try:
        pii = williams_exponent(beta, "II") if ang != 360 else 0.5
    except Exception as e:
        pii = f"ERR {e}"
    print(f"2beta={ang} pI={pi_:.6f} pII={pii}")

# opening sign on notch
th, a, b, c, p = williams_notch_I(2 * np.pi / 3)
dump("120 notch opening report", check_field(th, a, b, c, p=p, law="linear", beta=2 * np.pi / 3))
dump("120 notch closing stresses", check_field(th, -a, -b, -c, p=p, law="linear", beta=2 * np.pi / 3))

# Cartesian residual floor: T-stress values
rep = check_field(*t_stress(), p=0.0, law="linear", leading=False)
print("T cartesian", by_name(rep)["equilibrium (Cartesian, finite differences)"])
print("T saint-venant", by_name(rep)["compatibility (Saint-Venant, strains from the declared law)"])
print("T angular", by_name(rep)["equilibrium (angular ODEs, declared p)"])

rep = check_field(*williams_mode_II_crack(), p=0.5, law="linear", mode="II")
print("II cartesian", by_name(rep)["equilibrium (Cartesian, finite differences)"])
print("II saint-venant", by_name(rep)["compatibility (Saint-Venant, strains from the declared law)"])
