import numpy as np
import sys
sys.path.insert(0, '/work/code')
import fem
from fem import BilayerCrack
from central_crack import CentralCrack

mB = BilayerCrack(1.0, 2.0, 30.0, 3e-4, 2e-4)
mC = CentralCrack(1.0, 2.0, 30.0, 3e-4, 2e-4, 15.0)
print("B: nx", mB.nx, "ny", mB.ny, "nodes", len(mB.coords), "elems", len(mB.elems),
      "pairs", len(mB.pairs))
print("C: nx", mC.nx, "ny", mC.ny, "nodes", len(mC.coords), "elems", len(mC.elems),
      "pairs", len(mC.pairs))
print("B X[:5]", mB.X[:5], " X[-3:]", mB.X[-3:])
print("C X[:5]", mC.X[:5], " X[-3:]", mC.X[-3:])
print("B Y[:3]", mB.Y[:3], "Y[-2:]", mB.Y[-2:])
print("C Y[:3]", mC.Y[:3], "Y[-2:]", mC.Y[-2:])
# compare left_edges triples
for j in [0, mB.j0-1, mB.j0, mB.ny-2]:
    eb = mB.left_edges[j]; ec = mC.left_edges[j]
    print(f"j={j}: B n={eb[0]} y=({eb[1]:.4f},{eb[2]:.4f}) | "
          f"C n={ec[0]} y=({ec[1]:.4f},{ec[2]:.4f})")
# compare an interface cell element (upper and lower) in the crack region
i = 5
print("B upper cell5:", sorted(mB.elem_of_cell[mB.j0][i]),
      "  lower:", sorted(mB.elem_of_cell[mB.j0-1][i]))
print("C upper cell5:", sorted(mC.elem_of_cell[mC.j0][i]),
      "  lower:", sorted(mC.elem_of_cell[mC.j0-1][i]))
# node coordinate check for those elements
for tag, m in (("B", mB), ("C", mC)):
    e = m.elem_of_cell[m.j0][i]
    print(tag, "upper elem x,y:", np.round(m.coords[e], 5).tolist())
    e = m.elem_of_cell[m.j0-1][i]
    print(tag, "lower elem x,y:", np.round(m.coords[e], 5).tolist())
