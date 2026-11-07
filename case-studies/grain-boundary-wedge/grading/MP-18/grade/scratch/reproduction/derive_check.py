import numpy as np, march
from march import Smat, D1, D2, D3

# Derivation check: which of these operators, when used in an implicit march of
# the FULL (u,s) system, reproduces the observed decay?  Build the exact
# discrete system and its matrix exponential generator.

N = 60
M = Smat(N)                     # already /pi  (march.py convention)
D1m, D2m, D3m = D1(N), D2(N), D3(N)

# Full discrete system in u (interior unknowns):
#   s = 1 - M (D1 u) ; u_tau = D2 s   -> u_tau = D2(1 - M D1 u)
# Generator of the transient (u):  A_u = -D2 M D1
A_u = -(D2m @ M @ D1m)
ev_u = np.linalg.eigvals(A_u)
print("A_u=-D2 M D1 : slowest decay lambda =", -np.max(ev_u.real))

# Now the operator in s directly, using s = -M D1 u (transient part):
#   s_tau = -M D1 u_tau = -M D1 D2 s        -> A_s1 = -M D1 D2 (acts on s-unknowns)
A_s1 = -(M @ D1m @ D2m)
print("A_s1=-M D1 D2: slowest decay lambda =", -np.max(np.linalg.eigvals(A_s1).real))

# operator() in march.py = -M D1 D2  -> SAME as A_s1. So march.operator(N) IS A_s1.
M2, _, _, L = march.operator(N)
print("march.operator matches A_s1:", np.allclose(L, A_s1))

# The stated reduced form: s_tau = -(1/pi) Mraw D3.  Smat already includes /pi,
# so "-M D3" is the "-（1/pi) Mraw D3" with Mraw = Smat*pi.  i.e. L2 = -M D3.
L2 = -(M @ D3m)
print("L2=-M D3 : max real part (should be negative for decay):", np.max(np.linalg.eigvals(L2).real))