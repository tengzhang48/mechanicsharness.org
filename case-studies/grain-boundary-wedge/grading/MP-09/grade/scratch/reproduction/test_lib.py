import numpy as np
from scipy.integrate import quad
from gblib import *

print("=== TEST 1: C0 = PV int_0^inf S(1,t) dt should be -4 (derived analytically) ===")
def S1(t):
    if t == 1.0:
        return np.inf
    return 1.0/(1.0-t) - 1.0/(1.0+t) - 2.0*t*(1.0-t)/(1.0+t)**3
def S1_reg(t):
    # S1 minus its PV singularity 1/(1-t); finite limit -1/2 at t=1
    if t == 1.0:
        return -0.5
    return S1(t) - 1.0/(1.0-t)
# PV at t=1: integrate [0,2] subtracting singular part analytically
val1, _ = quad(S1_reg, 0, 2, points=[1.0], epsabs=1e-13, limit=200)
# PV int_0^2 dt/(1-t) = 0 (symmetric); then [2,inf)
val2, _ = quad(lambda t: S1(t), 2, np.inf, epsabs=1e-13, limit=200)
print("C0 numeric = %.15f  (expect -4)" % (val1+val2))

print("=== TEST 2: T_moment vs direct numerical PV integration ===")
def T_moment_num(x, m):
    # PV int_0^1 xi^m S(x,xi) dxi via singularity subtraction
    def reg(xi):
        return (xi**m - x**m)/(x-xi) - xi**m/(x+xi) - 2.0*xi**(m+1)*(x-xi)/(x+xi)**3
    val, _ = quad(reg, 0, 1, points=[x], epsabs=1e-12, limit=200)
    val += x**m*np.log(x/(1.0-x))
    return val
for m in [0, 1, 3, 7]:
    for x in [0.05, 0.3, 0.65, 0.95]:
        a = T_moment(x, m); b = T_moment_num(x, m)
        print("m=%d x=%.2f  recurrence=% .12e  quad=% .12e  diff=% .2e" % (m, x, a, b, a-b))

print("=== TEST 3: T_0 via antiderivative S_int ===")
for x in [0.1, 0.4, 0.8]:
    a = T_moment(x, 0)
    b = S_int(x, 1.0) - S_int(x, 1e-300)  # S_int(x,0)=4 exactly; use limit
    b_exact = S_int(x, 1.0) - 4.0
    print("x=%.1f  T0=% .12e  via S_int(-4 at 0)=% .12e  diff=% .2e" % (x, a, b_exact, a-b_exact))

print("=== TEST 4: Theta_k vs analytic Af for f = x^2(1-x)^3 ===")
# f = x^2 - 3x^3 + 3x^4 - x^5 ; f''' = -18 + 72x - 60x^2 = c0 T0* + c1 T1* + c2 T2*
g = np.array([-4.5, 6.0, -7.5])
def Af_exact(x):
    return -18.0*T_moment(x, 0) + 72.0*T_moment(x, 1) - 60.0*T_moment(x, 2)
for x in [0.02, 0.13, 0.37, 0.5, 0.71, 0.9, 0.99]:
    th = Theta_k(x, 2)
    val = np.dot(th, g)
    print("x=%.2f  Theta.g=% .12f  exact=% .12f  diff=% .2e" % (x, val, Af_exact(x), val-Af_exact(x)))

print("=== TEST 5: F_k properties (third derivative = T_k*, BCs) ===")
for k in [0, 1, 4, 9]:
    h = 1e-4
    xs = np.linspace(0.02, 0.98, 9)
    num3 = np.array([(F_k(x+1.5*h, k) - 3*F_k(x+0.5*h, k) + 3*F_k(x-0.5*h, k) - F_k(x-1.5*h, k))/h**3 for x in xs])
    exact = chebT(k, 2*xs-1)
    bc0 = F_k(0.0, k)
    f1 = (F_k(1.0, k)-F_k(1.0-h, k))/h
    f2 = (F_k(1.0, k)-2*F_k(1.0-h, k)+F_k(1.0-2*h, k))/h**2
    print("k=%d  max|F'''-T_k*|=%.2e   F(0)=%.1e  F'(1)~%.1e  F''(1)~%.1e"
          % (k, np.max(np.abs(num3-exact)), bc0, f1, f2))
print("ALL LIB TESTS DONE")
