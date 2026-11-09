"""Verify the main quasistatic branch by continuation (nearest-root walk),
confirm it reaches (M*, L*) with delta_t = delta_c, and print final answers."""
import numpy as np
from scipy.optimize import brentq
import solve_lib as sl

def roots_at(L):
    d0 = sl.delta_end(L, 0.0); d1 = sl.delta_end(L, 1.0); q = d1 - d0
    def g(M):
        return M*M*sl.CM/(2*sl.b) - (sl.G0 + sl.bridging_work(d0 + q*M))
    r = np.roots([sl.CM/(2*sl.b) + (sl.sig/(2*sl.dc))*q**2,
                  -q*sl.sig*(1 - d0/sl.dc), -(sl.G0 + sl.bridging_work(d0))])
    return sorted(x.real for x in r if abs(x.imag) < 1e-10 and x.real > 0), d0, q

# continuation from initiation in fine steps, nearest-root selection
Ls = np.concatenate([np.linspace(1e-4, 12e-3, 1200)])
Mprev = sl.Minit
path = []
ok = True
for L in Ls:
    rr, d0, q = roots_at(L)
    if not rr:
        ok = False; break
    M = min(rr, key=lambda x: abs(x - Mprev))
    dt = d0 + q*M
    path.append((L, M, dt))
    Mprev = M
print("continuation reached L =", path[-1][0]*1e3, "mm; M =", f"{path[-1][1]:.6f}",
      "; delta_t =", f"{path[-1][2]*1e6:.3f} um; monotone M:",
      bool(np.all(np.diff([p[1] for p in path]) > 0)),
      "; monotone dt:", bool(np.all(np.diff([p[2] for p in path]) > 0)))

Mpath = np.array([p[1] for p in path]); Lpath = np.array([p[0] for p in path])
dtpath = np.array([p[2] for p in path])
i = np.argmax(Mpath)
print(f"peak of path:  M = {Mpath[i]:.8f} N m at L = {Lpath[i]*1e3:.5f} mm, "
      f"delta_t = {dtpath[i]*1e6:.4f} um")
print(f"M* from J-identity: {sl.Mstar:.8f}  -> agreement: "
      f"{abs(Mpath[i]-sl.Mstar)/sl.Mstar:.2e}")

# stability (G_M > 0 -> advance stable under moment control) along the path
for L, M, dt in path[::150] + [path[i]]:
    d0_, q_ = roots_at(L)[1:]
    GM = sl.CM*M/sl.b - q_*sl.sig*(1 - dt/sl.dc)
    print(f"   L={L*1e3:7.3f} mm  M={M:.6f}  dt={dt*1e6:8.3f} um  G_M={GM:8.2f}")

# final event quantities (exact BVP at M*)
f = lambda L: sl.delta_end(L, sl.Mstar) - sl.dc
Le = brentq(f, 10e-3, 11e-3, xtol=1e-16, rtol=8.9e-16)
dmid = sl.delta_end_exact = float(np.asarray(sl.solution(Le, sl.Mstar, [Le/2])[0]).reshape(-1)[0])
print("\n================ FINAL ================")
print(f"M*   = {sl.Mstar:.6f} N m   (on each strip)")
print(f"L*   = {Le*1e3:.4f} mm")
print(f"d(L*/2) = {dmid*1e6:.3f} um")
print(f"J* = {sl.Mstar**2*sl.CM/(2*sl.b):.4f} J/m^2 = Gamma0 + sigma_b delta_c/2 "
      f"= {sl.G0 + sl.sig*sl.dc/2:.4f}")
