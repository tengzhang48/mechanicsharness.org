import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from central_crack import CentralCrack

mesh = CentralCrack(1.0, 2.0, 40.0, 3e-4, 2e-4, 1.0)
nid, hnid, j0 = mesh.nid, mesh.hnid, mesh.j0
# element containing cell at x ~ -1 (between the tips)
X = mesh.X
i = np.searchsorted(X, -1.0)
print("cell", i, "spans", X[i], X[i+1])
up = mesh.elem_of_cell[j0][i]
lo = mesh.elem_of_cell[j0-1][i]
print("upper elem nodes:", sorted(up), "-> x:", sorted(mesh.coords[up][:, 0]))
print("lower elem nodes:", sorted(lo), "-> x:", sorted(mesh.coords[lo][:, 0]))
print("shared nodes between them:", set(up) & set(lo))
# find the mid-column cell
i2 = np.searchsorted(X, -0.5)
up2 = mesh.elem_of_cell[j0][i2]
lo2 = mesh.elem_of_cell[j0-1][i2]
print("mid-crack cell", i2, X[i2], X[i2+1], "shared:", set(up2) & set(lo2))
# count shared nodes along the crack for upper/lower interface-adjacent elements
shared_count = 0
for i in range(mesh.nx-1):
    if -2.0 < X[i] < 0:
        s = set(mesh.elem_of_cell[j0][i]) & set(mesh.elem_of_cell[j0-1][i])
        if s:
            shared_count += 1
            print("  SHARED at cell", i, X[i], X[i+1], s)
print("cells with shared nodes inside crack:", shared_count)
