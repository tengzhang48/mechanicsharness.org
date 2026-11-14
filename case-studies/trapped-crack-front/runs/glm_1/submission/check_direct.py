"""Cleanest extraction: delta-K per input channel from the TOTAL ahead traction:
  T_tot(x) = T_b1(x) + T_reg(x) = A32 x^-3/2 + A12 x^-1/2 + A0 + ...
  A32 known exactly = 1/(sqrt(2pi) D_j)  (input channel j, own component).
  A12 = lim_{u->0} [ (T_tot(x) - A32 x^-3/2) x^{1/2} ]  (Richardson in u = sqrt(x)).
Then dK_i = sqrt(2pi) A12_i;  R[i,j] = (D_j/2) dK_i  (K_j=1, da = 2/D_j).
"""
import numpy as np
import first_order as fo
from kernels_x import apply_ahead_matrix

for ML in [6, 8]:
    sh = fo.solve_shear(ML=ML, verbose=False)
    basis = sh['basis']; nb = len(basis)
    xs = np.array([1e-6, 4e-6, 1.6e-5, 6.4e-5, 2.56e-4, 1.024e-3, 4.1e-3])
    u = np.sqrt(xs)
    R = np.zeros((2, 2), complex)
    for jc in range(2):
        s, res = sh['sols'][jc]
        ax, az = s[:nb], s[nb:]
        A32ex = 1.0/(fo.S2PI*fo.DVEC[jc])
        for i, ent in enumerate(['xx', 'zz']):
            Tb1 = fo.b1_ahead(ent, xs)*(jc == i)      # b1 only in own channel via diagonal kernels
            Treg = np.array([apply_ahead_matrix(basis, -x, [ent])[ent] @ (ax if i == 0 else az) for x in xs])
            Ttot = Tb1 + Treg
            # cross-kernel contribution from the OTHER component (reg only; b1-xz closed form too)
            other_ax = az if i == 0 else ax
            Treg_x = np.array([apply_ahead_matrix(basis, -x, ['xz'])[ 'xz'] @ other_ax for x in xs])
            Tb1_xz = fo.b1_ahead('xz', xs)*(1 - (jc == i))   # b1 off-diagonal via xz kernel
            Ttot = Ttot + Treg_x + Tb1_xz
            f = (Ttot - A32ex*xs**-1.5)*np.sqrt(xs)
            c = np.polyfit(u, f, 3)
            A12 = np.polyval(c, 0.0)
            R[i, jc] = fo.S2PI*(fo.DVEC[jc]/2.0)*A12
    Rsl = fo.slip_route_R(sh['cpm'])
    print(f"ML={ML}")
    print("  R direct-extraction:"); print(np.array_str(R, precision=8))
    print("  R slip:"); print(np.array_str(Rsl, precision=8))
    print("  diff:"); print(np.array_str(R-Rsl, precision=3))
