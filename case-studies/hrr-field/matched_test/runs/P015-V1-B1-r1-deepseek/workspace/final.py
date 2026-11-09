import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

mm, N, MPa = 1.0, 1.0, 1.0
t0, g0, n = 80.0, 0.002, 3.5
q = 1.0 / (n + 1.0)
Jin = 0.120


def rhs(th, y):
    F, Fp = y
    num = -q * (1 - q) * F * (q * q * F * F + Fp * Fp) \
          + (n - 1.0) / n * q * q * F * Fp * Fp
    den = q * q * F * F + Fp * Fp / n
    return [Fp, num / den]


sol = solve_ivp(rhs, [0, np.pi], [0, 1], rtol=1e-13, atol=1e-15,
                dense_output=True, max_step=1e-4)
th = np.linspace(0, np.pi, 800001)
F, Fp = sol.sol(th)
Fs, Fps = CubicSpline(th, F), CubicSpline(th, Fp)
g = np.sqrt(q * q * F * F + Fp * Fp)

thf = np.concatenate([-th[::-1], th[1:]])
Ff = np.concatenate([-F[::-1], F[1:]])
Fpf = np.concatenate([Fp[::-1], Fp[1:]])
gf = np.sqrt(q * q * Ff * Ff + Fpf * Fpf)
integ = (n / (n + 1.0) * gf ** (1 + 1.0 / n) * np.cos(thf)
         - q * Ff * gf ** (1.0 / n - 1) * (q * Ff * np.cos(thf) - Fpf * np.sin(thf)))
Int = np.trapezoid(integ, thf)
C = (Jin / (t0 * g0 ** (-1.0 / n) * Int)) ** (n / (n + 1.0))
A = t0 * g0 ** (-1.0 / n) * C ** (1.0 / n)

print("n=%.1f q=%.12f  Int=%.12f  C=%.12e  A=%.10f MPa*mm^q" % (n, q, Int, C, A))
print("F(pi)=%.10f  F'(pi)=%.2e  F'(0)=%.10f" % (F[-1], Fp[-1], Fp[0]))
print("g(0)=%.10f  g(pi/2)=%.10f  g(pi)=%.10f" % (g[0], np.interp(np.pi/2, th, g), g[-1]))
print()
# 1
print("1) Delta w = 2 C r^q F(pi):")
for r in (0.0200, 0.0800):
    print("   r=%.4f mm : Delta w = %.10f mm" % (r, 2 * C * r ** q * F[-1]))
print()
# 2
r2, t2 = 0.0200, np.radians(60.0)
g2 = float(np.interp(t2, th, g)); F2 = float(Fs(t2)); Fp2 = float(Fps(t2))
trz = A * r2 ** (-q) * q * F2 * g2 ** (1.0 / n - 1.0)
tthz = A * r2 ** (-q) * Fp2 * g2 ** (1.0 / n - 1.0)
print("2) r=0.0200 mm, theta=60deg: F=%.8f F'=%.8f g=%.8f" % (F2, Fp2, g2))
print("   tau_rz  = %.8f MPa" % trz)
print("   tau_thz = %.8f MPa" % tthz)
print()
# 3
print("3) |gamma| = 0.0100 => r = (C g/0.01)^(1/(1-q)):")
for deg in (0.0, 90.0, 180.0):
    gg = float(np.interp(np.radians(deg), th, g))
    print("   theta=%3.0f deg : g=%.8f  r=%.10f mm" % (deg, gg, (C * gg / 0.01) ** (1 / (1 - q))))
