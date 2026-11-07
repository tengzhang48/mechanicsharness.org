import numpy as np
from crack_solver import Params, assemble
P = Params(0.2, 0.4, 0.89)
s = 100.0; N = 60
Ma, n = assemble(s, N, P, kmax=800.0)
Mb, n = assemble(s, N, P, kmax=3000.0)
Mc, n = assemble(s, N, P, kmax=3000.0, width=0.25, npts=16, ntheta=1500)
for nm, (i0, i1, j0, j1) in {'M11': (0,N,0,N), 'M12': (0,N,N,2*N), 'M21': (N,2*N,0,N), 'M22': (N,2*N,N,2*N)}.items():
    A = Ma[i0:i1, j0:j1]; B = Mb[i0:i1, j0:j1]; C = Mc[i0:i1, j0:j1]
    print(nm, 'max|a-b|', np.max(np.abs(A-B)), ' max|b-c|', np.max(np.abs(B-C)), ' max|b|', np.max(np.abs(B)))
    idx = np.unravel_index(np.argmax(np.abs(A-B)), A.shape); print('    worst a-b at', idx)
rhs = np.zeros(2*N); rhs[0] = -np.pi/2/s
for M in [Ma, Mb, Mc]:
    x = np.linalg.solve(M, rhs)*s
    a, b = x[:N], x[N:]
    print('K', np.sum(a*(n+1))/(1-P.nu), ' a[16], a[32], b[16], b[32], b[59]:', a[8], a[16], b[8], b[16], b[59])
