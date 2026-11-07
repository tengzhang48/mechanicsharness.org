import sys
from pathlib import Path
import numpy as np
from scipy.linalg import eig
from scipy.special import roots_legendre
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).parent / 'reproduction'))
from solver10 import Res10
from solver5 import cheb_T
from crack import solve_crack

for n in [16, 24, 32, 40, 44]:
    r = Res10(n)
    vals, vl, vr = eig(-r.H, r.F, left=True, right=True)
    ix = [i for i, v in enumerate(vals) if abs(v.imag) < 1e-7 and 0 < v.real < 20000]
    ix.sort(key=lambda i: vals[i].real)
    print('N', n, 'eigenvalues', [float(vals[i].real) for i in ix[:7]])
    residues = []
    for i in ix[:7]:
        a = vr[:, i].real
        w = vl[:, i].real
        residue_a = a * ((w @ np.ones(n+1)) / (w @ r.F @ a))
        residues.append(float(r._frows([0.5])[0] @ residue_a))
    print('residues_at_half', residues, 'partial_sums', np.cumsum(residues).tolist())
    i = ix[0]
    a = vr[:, i].real
    half, end = r._frows([0.5, 1.0]) @ a
    print('shape', half/end, 'mode1_amplitude_at_1', residues[0]/(half/end))
    ac, _, _ = solve_crack(n)
    print('p0_plus_crack', np.max(np.abs(r.solve(0)+ac)))

r = Res10(32)
lam = 11.74266781911283
u, s, vt = np.linalg.svd(-lam*r.F-r.H)
a = vt[-1]
print('main_pole_smin', s[-1], 'f1_svd', r._frows([1.0])[0] @ a)
for m in [1000, 2000, 4000]:
    t, w = roots_legendre(m)
    xi, ww = (t+1)/2, w/2
    rho = (1-xi)**-0.5 * (cheb_T(xi,32) @ a)
    errs = []
    for x in [0.3,0.5,0.7]:
        rhox = (1-x)**-0.5 * (cheb_T(np.array([x]),32) @ a)[0]
        kp = np.sum(ww*(rho-rhox)/(x-xi)) + rhox*np.log(x/(1-x))
        ki = np.sum(ww*rho*(-1/(x+xi)-2*xi*(x-xi)/(x+xi)**3))
        fx = r._frows([x])[0] @ a
        errs.append(abs(kp+ki+lam*fx)/abs(kp+ki))
    print('unweighted_quad_m', m, 'errors', errs)

for m in [120, 240, 480]:
    t, w = roots_legendre(m)
    t, w = (t+1)/2, w/2
    xi = 1-t*t
    q = cheb_T(xi,32) @ a
    errs = []
    for x in [0.3,0.5,0.7]:
        d = np.sqrt(1-x)
        qx = (cheb_T(np.array([x]),32) @ a)[0]
        kp = np.sum(2*w*(q-qx)/(t*t-d*d)) + qx*np.log((1-d)/(1+d))/d
        ki = np.sum(2*w*q*(-1/(x+xi)-2*xi*(x-xi)/(x+xi)**3))
        fx = r._frows([x])[0] @ a
        errs.append(abs(kp+ki+lam*fx)/abs(kp+ki))
    print('endpoint_regularized_m', m, 'errors', errs)

print('recorded_eigen5_partial_sum', sum([1.2266,0.0622,-0.3586,-0.0667,0.0380]))
