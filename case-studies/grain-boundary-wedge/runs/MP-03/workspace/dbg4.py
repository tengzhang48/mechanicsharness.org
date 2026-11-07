import numpy as np
from eig import eigen_s
xg, lam, V, lams = eigen_s(N=24, verbose=False)
fin = lam[np.isfinite(lam)]
print("num finite:", len(fin))
print("max real part:", fin.real.max(), " min real:", fin.real.min())
order = np.argsort(-fin.real)
for l in fin[order][:12]:
    print(f"  {l.real: .6g} {l.imag:+.3g}j")
