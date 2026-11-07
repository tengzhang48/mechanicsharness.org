import numpy as np
from crack_solver import Params, assemble
P = Params(0.2, 0.4, 0.89)
s = 100.0
for N in [40, 80]:
    M, n = assemble(s, N, P)
    print('N', N, 'cond', np.linalg.cond(M))
    # diagonal scaling
    D = np.sqrt(np.abs(np.diag(M)))
    Ms = M/D[:, None]/D[None, :]
    print('  scaled cond', np.linalg.cond(Ms))
    # check symmetry structure: M21^T should equal -s*M12 * factor?
    M11 = M[:N,:N]; M12 = M[:N,N:]; M21 = M[N:,:N]; M22 = M[N:,N:]
    print('  M11 sym err', np.max(np.abs(M11-M11.T))/np.max(np.abs(M11)), ' M22 sym err', np.max(np.abs(M22-M22.T))/np.max(np.abs(M22)))
    print('  M21^T + s M12 err', np.max(np.abs(M21.T + s*M12)))
    print('  diag M22 first/last', M22[0,0], M22[1,1], M22[-1,-1], ' diag M11 last', M11[-1,-1])
    # off-diagonal decay of M22 rows
    print('  M22 row 30 sample', M22[min(30,N-1), [0,1,2,5,10,20,29]] if N>30 else '')
