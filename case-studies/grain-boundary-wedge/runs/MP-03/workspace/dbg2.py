import numpy as np
from crack import crack_solve
from numpy.polynomial.legendre import leggauss
from gbcore import Rfun
from numpy.polynomial.chebyshev import chebval

kapC, kapf, v_of = crack_solve(N=32, verbose=False)
N, Mt = 32, 800
t, wt = leggauss(Mt); t = 0.5*(t+1); wt = 0.5*wt
def A_ahead(x, v):
    # A[w](x) for x>1, w = v(t)/t with xi=1-t^2
    vq = v(1-t**2)
    I1 = 2*np.sum(wt * vq/(x-1+t**2))
    I3 = 2*np.sum(wt * Rfun(x, 1-t**2) * vq)
    return I1 + I3
d = np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05])
xs = 1 + d
y = np.array([ (1 - A_ahead(x, v_of))*np.sqrt(2*(x-1)) for x in xs])
# fit y = k(1 + c d)
Ab = np.stack([np.ones_like(d), d], axis=1)
coef, *_ = np.linalg.lstsq(Ab, y, rcond=None)
print("s(x)*sqrt(2(x-1)) values:", y)
print(f"extrapolated kappa = {coef[0]:.6f} (cf. kappa(C)={kapC:.6f}, handbook 1.1215)")
