import numpy as np, march
from march import Smat, D1, D2, D3

# Compare L1 = -M D1 D2  with  L2 = -(1/pi) M D3.
# The correct reduced eq, from s = 1 - M w, w = D1 u, w_tau = D1 u_tau = D1 D2 s:
#   s_tau = -M w_tau = -M D1 D2 s            (L1)
# And also s = 1 - M w  =>  s_tau = -M w_tau ;  w = d_x u ; s''' = d_x^3 s... but w != s''
# So why would M D3 appear? Only if w relates to D3. Let's check M D1 D2 vs M D3 numerically.
N = 40
M = Smat(N)
L1 = M @ D1(N) @ D2(N)
L2 = (1.0 / np.pi) * M @ D3(N)
print("||L1 - L2||/||L1|| =", np.linalg.norm(L1 - L2) / np.linalg.norm(L1))
# Maybe D1@D2 == D3/h? no. check
print("||D1@D2 - D3||/||D3|| =", np.linalg.norm(D1(N) @ D2(N) - D3(N)) / np.linalg.norm(D3(N)))
# ratio of norms on a smooth test vector
x = np.arange(1, N + 1) / (N + 1)
v = np.sin(np.pi * x / 2)
print("D1 D2 v vs D3 v ratio sample:", (D1(N) @ D2(N) @ v)[:3], (D3(N) @ v)[:3])