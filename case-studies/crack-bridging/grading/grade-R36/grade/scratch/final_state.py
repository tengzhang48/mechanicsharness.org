"""Pin down L* where delta_t = dc along the R-curve, by continuation."""
import numpy as np
import solve as S
from scipy.optimize import brentq

def M_of_L(L, M_seed):
    def g(M):
        d, _ = S.solve_delta(L, M)
        return M**2 / (2 * S.b * S.Dtil) - S.bridging_work(d(L)[0]) - S.Gam0
    Ms = np.linspace(M_seed * 0.9, M_seed * 1.3, 800)
    vals = np.array([g(m) for m in Ms])
    for i in range(len(Ms) - 1):
        if vals[i] * vals[i + 1] < 0:
            return brentq(g, Ms[i], Ms[i + 1])
    return np.nan

if __name__ == "__main__":
    # continue from L=10.4 mm (M~0.2119) to find the dt=dc crossing
    M_seed = 0.2119
    Ls = np.linspace(10.4e-3, 11.2e-3, 41)
    prev = None
    for L in Ls:
        M = M_of_L(L, M_seed if prev is None else prev)
        if np.isnan(M):
            continue
        prev = M
        dt = S.solve_delta(L, M)[0](L)[0]
        print("L=%8.4f mm  M=%.6f  dt=%9.4f um" % (L * 1e3, M, dt * 1e6))

    def h(L, M_seed):
        M = M_of_L(L, M_seed)
        return S.solve_delta(L, M)[0](L)[0] - S.dc
    Lstar = brentq(lambda L: h(L, 0.2126), 10.5e-3, 11.0e-3, xtol=1e-13)
    Mstar = M_of_L(Lstar, 0.2126)
    d, _ = S.solve_delta(Lstar, Mstar)
    print("\nL* = %.8f mm" % (Lstar * 1e3))
    print("M* = %.8f N m" % Mstar)
    print("delta(L*) = %.6f um" % (d(Lstar)[0] * 1e6))
    print("delta(L*/2) = %.6f um" % (d(Lstar / 2)[0] * 1e6))
    print("delta(end) = %.4f um" % (d(Lstar + S.a0)[0] * 1e6))
    print("G_far = %.6f J/m^2" % (Mstar**2 / (2 * S.b * S.Dtil)))