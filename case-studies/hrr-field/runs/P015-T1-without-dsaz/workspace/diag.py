import numpy as np
import hrr2 as H
from scipy.integrate import solve_ivp


def urp(theta, Y):
    phi, d1, d2, d3 = Y
    A, B, C, se = H.quantities(phi, d1, d2, d3)
    Ap = H.k*d1 + d3
    Bp = H.k*(H.k-1.0)*d1
    Cp = -(H.k-1.0)*d2
    sep = (0.75*2.0*(A-B)*(Ap-Bp) + 3.0*2.0*C*Cp)/(2.0*se)
    q = se**(H.n-1.0)
    qp = (H.n-1.0)*se**(H.n-2.0)*sep
    Urp = -(1.0/(2.0*H.p))*(qp*(A-B) + q*(Ap-Bp))
    return Urp, q, C


print(" t   U_rp(0)   U_t(0)   A(pi)   C(pi)")
for t in [-0.5, -0.25, -0.1, 0.0, 0.25, 0.5, 1.0]:
    Y0 = [1.0, 0.0, t, 0.0]
    s = solve_ivp(H.rhs, [0.0, np.pi], Y0, rtol=1e-11, atol=1e-13,
                  method='Radau', max_step=0.01)
    if not s.success:
        print(" t=%.2f FAIL" % t); continue
    Urp0, q0, C0 = urp(0.0, Y0)
    Ut0 = (Urp0 - 2.0*q0*C0)/(H.p+1.0)
    phi, d1, d2, d3 = s.y[:, -1]
    A = H.k*phi + d2
    C = -(H.k-1.0)*d1
    print(" %+.2f  %+.4f   %+.4f   %+.4f  %+.4f" % (t, Urp0, Ut0, A, C))