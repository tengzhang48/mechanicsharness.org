"""Trace the quasistatic path M(L): J_far(M) = Gamma0 + int_0^{dt} t dd (lower
branch by continuation), and find where dt = delta_c (first connector inactive)."""
import numpy as np
from scipy.optimize import brentq
import solve_lib as sl

Ls = np.concatenate([np.linspace(0.5e-3, 20e-3, 196), np.linspace(20.5e-3, 40e-3, 80)])
Mprev = sl.Minit
print(f"{'L(mm)':>9} {'M_lo':>10} {'dt(um)':>9} {'J(M)':>8} {'G_M':>10} {'2nd root':>10}")
rows = []
for L in Ls:
    d0 = sl.delta_end(L, 0.0); d1 = sl.delta_end(L, 1.0); q = d1 - d0
    def g(M):
        return M*M*sl.CM/(2*sl.b) - (sl.G0 + sl.bridging_work(d0 + q*M))
    try:
        Mlo = brentq(g, Mprev*0.7, Mprev*1.5 + 1e-4, xtol=1e-15, rtol=8.9e-16)
    except ValueError:
        Mlo = np.nan
    dt = d0 + q*Mlo
    GM = sl.CM*Mlo/sl.b - q*sl.sig*(1 - dt/sl.dc) if Mlo == Mlo else np.nan
    r2 = np.roots([sl.CM/(2*sl.b) + (sl.sig/(2*sl.dc))*q**2,
                   -q*sl.sig*(1 - d0/sl.dc), -(sl.G0 + sl.bridging_work(d0))])
    r2p = [x.real for x in r2 if abs(x.imag) < 1e-9 and x.real > 0]
    rows.append((L, Mlo, dt, GM, sorted(r2p)))
    Mprev = Mlo if Mlo == Mlo else Mprev
    print(f"{L*1e3:9.3f} {Mlo:10.6f} {dt*1e6:9.3f} {Mlo*Mlo*sl.CM/(2*sl.b):8.4f} "
          f"{GM:10.2f} {str(sorted([round(x,5) for x in r2p])):>22}")
