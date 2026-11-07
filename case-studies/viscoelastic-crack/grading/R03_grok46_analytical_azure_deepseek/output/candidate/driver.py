import numpy as np
import solve
from scipy.optimize import fsolve


def Cres(V, a, eps, np_, ny, dmax=0.2):
    f, K, T0, W, S2, X, Y, K0 = solve.residuals(V, a, eps, npanel=np_, ny=ny)
    U = W[:, 0]
    d = -X
    m = (X < 0) & (d < dmax) & (d > 0)
    dd = d[m]; uu = U[m]
    A = np.vstack([dd ** 0.5, dd ** 1.5, dd ** 2.5]).T
    c, *_ = np.linalg.lstsq(A, uu, rcond=None)
    return f, c[0], K, W, S2, X, Y, K0


def solve_CE(V, np_, ny, guess=(0.53, 0.87), dmax=0.2):
    def res(z):
        a, eps = z
        f, C, K, *_ = Cres(V, a, eps, np_, ny, dmax)
        return [f - 0.1, C]
    sol, info, ier, msg = fsolve(res, guess, full_output=True, xtol=1e-11)
    return sol, info, ier


if __name__ == "__main__":
    for V in [0.1]:
        for np_ in [(100, 60, 100), (160, 100, 160), (240, 150, 240), (360, 230, 360)]:
            ny = np_[1]
            sol, info = solve_CE(V, np_, ny)
            print(np_, 'a=%.6f eps=%.6f' % (sol[0], sol[1]), 'fvec', np.round(info['fvec'], 9))