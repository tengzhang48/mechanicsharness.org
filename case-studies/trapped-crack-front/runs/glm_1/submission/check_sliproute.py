"""Slip-route extraction and cross-check against the traction-route R matrix.
Per input channel j (unit K_j^0, delta a = 2/D_j, a' = i delta a = 2i/D_j):
  c+^x = D2 dK2^loc + a' D3 K3^0 d_{jz}      c+^z = D3 dK3^loc - a' D2 K2^0 d_{jx}
  dK2^g = dK2^loc + a' K3^0 d_{jz}           dK3^g = dK3^loc - a' K2^0 d_{jx}
  R_slip[i,j] = (D_j/2) dK_i^g
"""
import numpy as np
import first_order as fo

D1, D3 = fo.D1, fo.D3
Dvec = np.array([D1, D3])

for ML in [5, 6, 8]:
    sh = fo.solve_shear(ML=ML, verbose=False)
    cpm = sh['cpm']               # cpm[i,j] = sqrt-edge coeff per unit input, incl. b1 -1
    R_tr = sh['R']
    R_sl = np.zeros((2, 2), complex)
    for j in range(2):
        Dj = Dvec[j]
        dK2loc = (cpm[0, j] - 2j*1.0*(j == 1))/D1
        dK3loc = (cpm[1, j] + 2j*1.0*(j == 0))/D3
        dK2g = dK2loc + (2j/D3)*(j == 1)
        dK3g = dK3loc - (2j/D1)*(j == 0)
        R_sl[0, j] = (Dj/2)*dK2g
        R_sl[1, j] = (Dj/2)*dK3g
    print(f"ML={ML}")
    print("  R traction-route:"); print(np.array_str(R_tr, precision=6, suppress_small=True))
    print("  R slip-route:");     print(np.array_str(R_sl, precision=6, suppress_small=True))
    print("  diff:"); print(np.array_str(R_tr-R_sl, precision=2, suppress_small=False))
