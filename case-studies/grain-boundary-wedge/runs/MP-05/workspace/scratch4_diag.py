import numpy as np
from core import *
np.set_printoptions(precision=4, suppress=False, linewidth=200)
N = 32
lam, V, res, A, B, xc = solve_eig(N)
print("cond(A)=", np.linalg.cond(A), "cond(B)=", np.linalg.cond(B))
idx = np.argsort(np.abs(lam))
for k in idx[:20]:
    print(f"lam={lam[k]!r:>44}  res={res[k]:.1e}")
print("---- all |Im| small & Re>0 sorted by Re:")
sel = [l for l in lam if abs(l.imag) < 1e-6 and l.real > 0]
sel.sort()
print(np.array(sel[:15]))
print("---- count by region: Re>0:", sum(1 for l in lam if l.real>0), " Re<0:", sum(1 for l in lam if l.real<0), " |Im|>1e-6:", sum(1 for l in lam if abs(l.imag)>1e-6))
