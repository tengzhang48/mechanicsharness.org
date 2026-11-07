"""Test whether lambda1 -> 0 or a finite limit, by inspecting the slowest
eigenvector: is it (a) the gauge (constant) direction, i.e. a discretisation
null mode, or (b) a genuine smooth eigenfunction?"""
import numpy as np
from numpy.polynomial.legendre import leggauss
import solve as S

xq, wq = leggauss(2000); xq, wq = 0.5*(xq+1), 0.5*wq
for n in (10, 16, 22, 28, 34, 40):
    s = np.arange(n+1)/n
    E = S.elasticity_matrix(s, xq, wq)
    A = -(E @ S.d2_uniform(n))
    idx = np.arange(1, n+1)
    ev, V = np.linalg.eig(A[np.ix_(idx, idx)])
    keep = np.abs(ev.imag) < 1e-3*np.abs(ev)
    pairs = sorted([(ev.real[j], j) for j in np.where(keep & (ev.real>1e-3))[0]])
    l1, j1 = pairs[0]
    v = V[:, j1].real; v = v/np.linalg.norm(v)
    # overlap with constant direction
    const = np.ones(n)/np.sqrt(n)
    ov = abs(v @ const)
    # smoothness: norm of (v - shifted)*n  (large -> oscillatory/singular)
    rough = n*np.linalg.norm(np.diff(v))
    print('n=%2d  lam1=%.5f  |<v,1>|=%.4f  rough=%.3f  first few v: %s'
          % (n, l1, ov, rough, np.round(v[:4],4)))