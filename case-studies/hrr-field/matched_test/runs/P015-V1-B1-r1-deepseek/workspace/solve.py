import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# Constitutive (antiplane shear, mode III):
#   gamma = (g0/t0) (tau/t0)^(n-1) tau,  t0=80 MPa, g0=0.002, n=3.5
# Separable crack-tip field:
#   w = C r^q F(theta),  q = 1/(n+1)  (fixed by finite nonzero J)
#   gamma_r = C q r^{q-1} F, gamma_theta = C r^{q-1} F'
#   tau_r = t0 g0^{-1/n} C^{1/n} r^{-q} q F g^{1/n-1}
#   tau_theta = t0 g0^{-1/n} C^{1/n} r^{-q} F' g^{1/n-1}
#   g = sqrt(q^2 F^2 + F'^2)
# Equilibrium ODE: (F' H)' + q(1-q) F H = 0, H = g^{1/n-1}

n = 3.5
t0 = 80.0
g0 = 0.002
Jin = 0.120
q = 1.0 / (n + 1.0)
print("q = 1/(n+1) =", q)


def ode_rhs(theta, y):
    F, Fp = y
    g = np.sqrt(q * q * F * F + Fp * Fp)
    # F'' = [ -q(1-q) F (q^2 F^2 + Fp^2) + (n-1)/n q^2 F Fp^2 ] / (q^2 F^2 + Fp^2/n)
    num = -q * (1 - q) * F * (q * q * F * F + Fp * Fp) \
          + (n - 1.0) / n * q * q * F * Fp * Fp
    den = q * q * F * F + Fp * Fp / n
    return [Fp, num / den]


def integrate(Fp0, qq):
    global q
    qsave = q
    q = qq
    sol = solve_ivp(ode_rhs, [0.0, np.pi], [0.0, Fp0],
                    rtol=1e-12, atol=1e-14, dense_output=True, max_step=1e-3)
    q = qsave
    return sol


def face_slope(qq):
    sol = integrate(1.0, qq)
    return sol.y[1, -1]


# Check the crack exponent is a root of the eigenproblem (root near 1/(n+1))
r = brentq(face_slope, 0.15, 0.35, xtol=1e-15, rtol=1e-15)
print("shooting root q* = %.12f   1/(n+1) = %.12f" % (r, q))
print("F'(pi) at q=1/(n+1):", face_slope(q))

# Use q = 1/(n+1) (the physically fixed exponent)
sol = integrate(1.0, q)
th = np.linspace(0.0, np.pi, 400001)
Y = sol.sol(th)
F, Fp = Y[0], Y[1]
g = np.sqrt(q * q * F * F + Fp * Fp)

# Extend to full circle (odd symmetry: F(-th) = -F(th))
thf = np.concatenate([-th[::-1], th[1:]])
Ff = np.concatenate([-F[::-1], F[1:]])
Fpf = np.concatenate([Fp[::-1], Fp[1:]])
gf = np.sqrt(q * q * Ff * Ff + Fpf * Fpf)

# J angular integral
# J = t0 g0^{-1/n} C^{(n+1)/n} * I
#   I = int [ n/(n+1) g^{1+1/n} cos th - q F g^{1/n-1} (q F cos th - F' sin th) ] dth
integrand = (n / (n + 1.0) * gf ** (1.0 + 1.0 / n) * np.cos(thf)
             - q * Ff * gf ** (1.0 / n - 1.0) * (q * Ff * np.cos(thf) - Fpf * np.sin(thf)))
Int = np.trapezoid(integrand, thf)
print("I =", Int)

C = (Jin / (t0 * g0 ** (-1.0 / n) * Int)) ** (n / (n + 1.0))
print("C =", C)

# ---- Requested quantities ----
# 1. relative face slip Delta w(r) = w(r,pi)-w(r,-pi) = 2 C r^q F(pi)
Fpi = F[-1]
print("F(pi) =", Fpi, "F'(pi) =", Fp[-1])
for rr in (0.0200, 0.0800):
    dw = 2.0 * C * rr ** q * Fpi
    print(f"Delta w(r={rr}) = {dw:.8f} mm")


def g_at(theta):
    return np.interp(theta, th, g)


# theta in degrees
def gdeg(deg):
    return float(np.interp(np.radians(deg), th, g))


# 2. signed stresses tau_rz, tau_thz at r=0.02, theta=60deg
r2 = 0.0200
th2 = np.radians(60.0)
g2 = gdeg(60.0)
A = t0 * g0 ** (-1.0 / n) * C ** (1.0 / n)  # MPa*mm^q
F2 = float(np.interp(th2, th, F))
Fp2 = float(np.interp(th2, th, Fp))
trz = A * r2 ** (-q) * q * F2 * g2 ** (1.0 / n - 1.0)
tthz = A * r2 ** (-q) * Fp2 * g2 ** (1.0 / n - 1.0)
print(f"tau_rz(0.02,60) = {trz:.8f} MPa")
print(f"tau_thz(0.02,60) = {tthz:.8f} MPa")

# 3. radii where |gamma| = 0.01 on theta = 0, 90, 180
gtar = 0.0100
for deg in (0.0, 90.0, 180.0):
    gg = gdeg(deg)
    rr = (C * gg / gtar) ** (1.0 / (1.0 - q))
    print(f"r(|gamma|=0.01) at theta={deg}deg : {rr:.8f} mm  (g={gg:.6f})")

# verification: J from a few contours (already r-independent by construction) and
# compare with a J integral computed directly from the Cartesian field via finite numbers.
np.save('/tmp/none.npy', np.array([0]))
