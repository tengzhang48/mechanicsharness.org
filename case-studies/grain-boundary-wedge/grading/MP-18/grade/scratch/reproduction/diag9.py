import numpy as np, march
from march import Smat, D1, D2

# Physical closed system in u:
#   s = 1 - M (D1 u)          (M already = Smat with 1/pi)
#   u_tau = D2 s              (d_tau u = d_x^2 s)
# => u_tau = D2 (1 - M D1 u) = -D2 M D1 u + D2 1
# eigenvalue problem for the transient part: u_tau = -D2 M D1 u
# Decay operator K = -D2 M D1 ; its eigenvalues -> decay of u.
# s decays via s ~ -M D1 u, so same decay rates.
N = 40
M = Smat(N); D1m = D1(N); D2m = D2(N)
K = -(D2m @ M @ D1m)
ev = np.linalg.eigvals(K)
print("K=-D2 M D1 eigs: max real (slowest decay):", np.max(ev.real))
print("top real parts:", np.sort(ev.real)[::-1][:6])

# Also the exactly stated reduced operator L = -(1/pi) M D3, but with 1/pi separate:
from march import D3
L = -(1.0 / np.pi) * (Smat(N) * np.pi) @ D3(N)   # = -(1/pi)(pi*Mraw) D3 = -Mraw D3
# Mraw = Smat*pi
print("check: does K have negative real parts (decay)?", np.all(ev.real < 1e-9))