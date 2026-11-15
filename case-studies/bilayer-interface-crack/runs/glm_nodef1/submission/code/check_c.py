"""Refine c = J/|K|^2 and compare with candidate closed forms."""
import numpy as np
import sys
sys.path.insert(0, '/work/code')
from eigenpair import Eig, jint

E1, nu1, E2, nu2 = 9.0, 1/3, 1.0, 1/3
eig = Eig(E1, nu1, E2, nu2)
a = 0.7+0.5j
K = np.sqrt(2*np.pi)*np.conj(a)*(1-eig.cstar)
Eb1, Eb2 = E1/(1-nu1**2), E2/(1-nu2**2)
S = 1/Eb1 + 1/Eb2
eps_std = -eig.eps_m

print("convergence of J with ntheta:")
cvals = []
for n in [20000, 200000, 2000000]:
    th = np.linspace(0, 2*np.pi, n)[:-1] + 1e-9
    z = 0.1*np.exp(1j*th)
    f = eig.fields(z, a)
    iE = np.where(f["up"], 0, 1)
    E, nu = eig.E[iE], eig.nu[iE]
    sxx, syy, sxy = f["sxx"], f["syy"], f["sxy"]
    ee11 = (1+nu)/E*((1-nu)*sxx - nu*syy)
    ee22 = (1+nu)/E*((1-nu)*syy - nu*sxx)
    ee12 = (1+nu)/E*sxy
    W = 0.5*(sxx*ee11 + syy*ee22 + 2*sxy*ee12)
    nx, ny = np.cos(th), np.sin(th)
    tx, ty = sxx*nx + sxy*ny, sxy*nx + syy*ny
    integrand = (W*nx - (tx*f["ux"] + ty*f["vx"]))*0.1
    J = np.trapz(integrand, th)
    cvals.append(J/abs(K)**2)
    print(f"  n={n}: c = {J/abs(K)**2:.10f}")

c = cvals[-1]
print(f"\nc                  = {c:.10f}")
print(f"S/2*sech^2(pi eps) = {S/2/np.cosh(np.pi*eps_std)**2:.10f}")
print(f"S/2                = {S/2:.10f}")
print(f"S/2*cosh^2(pi eps) = {S/2*np.cosh(np.pi*eps_std)**2:.10f}")
print(f"ratio to sech^2 form: {c/(S/2/np.cosh(np.pi*eps_std)**2):.10f}")

# homogeneous limits of the sech^2 form are exact; verify another nu:
for (n1_, n2_) in [(0.3, 0.3), (1/3, 1/3)]:
    eh = Eig(3.0, n1_, 5.0, n2_)
    ah = 0.3+0.9j
    Kh = np.sqrt(2*np.pi)*np.conj(ah)*(1-eh.cstar)
    th = np.linspace(0, 2*np.pi, 400000)[:-1] + 1e-9
    z = 0.07*np.exp(1j*th)
    f = eh.fields(z, ah)
    iE = np.where(f["up"], 0, 1)
    E, nu = eh.E[iE], eh.nu[iE]
    sxx, syy, sxy = f["sxx"], f["syy"], f["sxy"]
    ee11 = (1+nu)/E*((1-nu)*sxx - nu*syy)
    ee22 = (1+nu)/E*((1-nu)*syy - nu*sxx)
    ee12 = (1+nu)/E*sxy
    W = 0.5*(sxx*ee11 + syy*ee22 + 2*sxy*ee12)
    nx, ny = np.cos(th), np.sin(th)
    tx, ty = sxx*nx + sxy*ny, sxy*nx + syy*ny
    J = np.trapz((W*nx - (tx*f["ux"] + ty*f["vx"]))*0.07, th)
    Sh = 1/(3.0/(1-n1_**2)) + 1/(5.0/(1-n2_**2))
    epsh = -eh.eps_m
    print(f"nu=({n1_},{n2_}) E=(3,5): c={J/abs(Kh)**2:.8f}  "
          f"S/2*sech^2={Sh/2/np.cosh(np.pi*epsh)**2:.8f}")
