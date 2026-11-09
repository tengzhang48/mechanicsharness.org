import sys, numpy as np
sys.path.insert(0, "knowledge")
import hrr

n = 6.5
p_sing = 1.0 / (n + 1.0)
print("s =", hrr.singular_exponent(n), "  p = 1/(n+1) =", p_sing)

sol = hrr.solve_angular(n, plane_strain=True, n_theta=16001)
I = hrr.integral_I(n, plane_strain=True, sol=sol)
En = hrr.energy_coefficient_En(n, plane_strain=True, sol=sol, I=I)

th = sol["theta"]
se, st, sr = sol["se"], sol["st"], sol["sr"]
shape_full = 2.0 * np.trapezoid(se ** (n + 1.0), th)   # integral over -pi..pi

H = st[0] / I ** p_sing
C = (n / (n + 1.0)) * shape_full / I

print(f"I(6.5, ps)      = {I:.10f}")
print(f"sigma_tilde_th(0) = {st[0]:.10f}")
print(f"sigma_tilde_r(0)  = {sr[0]:.10f}")
print(f"E_n              = {En:.10f}   (check E_n/pi = {En/np.pi:.10f})")
print(f"H = st(0)/I^p    = {H:.8f}")
print(f"C = (n/(n+1))*int(se^(n+1))/I = {C:.8f}")
print(f"p                = {p_sing:.10f}")
print(f"shooting f''(0) = {sol['p']:.10f},  c = {sol['c']:.6e}")
print(f"max se = {se.max():.12f} (should be 1)")
