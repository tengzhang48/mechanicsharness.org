import numpy as np
import sys
sys.path.insert(0, '/work/code')
from central_crack import CentralCrack

mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
X = mesh.X
d = np.diff(X)
print("nx:", len(X), "min dx:", d.min(), "max dx:", d.max())
# check for duplicate/too-close breakpoints
print("dx below 1e-12:", (d < 1e-12).sum())
# element areas via corner cross products
e = mesh.elems
c = mesh.coords
p = c[e[:, [0, 1, 2, 3]]]
area2 = ((p[:, 1, 0]-p[:, 0, 0])*(p[:, 2, 1]-p[:, 0, 1])
         - (p[:, 2, 0]-p[:, 0, 0])*(p[:, 1, 1]-p[:, 0, 1]))
print("min/max quad double-area:", area2.min(), area2.max())
bad = np.where(area2 <= 0)[0]
print("non-positive quads:", len(bad), mesh.elems[bad][:3] if len(bad) else "")
# where are the smallest elements?
i = np.argmin(d)
print("smallest interval at x=", X[i], "..", X[i+1])
for x0 in [-2.0, -1.0, 0.0]:
    j = np.argsort(np.abs(X-x0))[:6]
    print(f"breakpoints near {x0}:", X[np.sort(j)])
