import numpy as np
from eig import eigen_s
xg, lam, V, lams = eigen_s(N=40, verbose=False)
fin = lam[np.isfinite(lam)]
print("num eigen:", len(lam), "finite:", len(fin))
order = np.argsort(fin.real)
for l in fin[order][:15]:
    print(f"  {l.real: .6g} {l.imag:+.3g}j")
