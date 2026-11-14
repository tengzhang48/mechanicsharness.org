"""With the fixed Hadamard machinery: compare numerically-computed on-crack
traction of the b1 shear columns against the closed-form residue prediction
e^{-X} * sqrt(pi) cc nu /sqrt(2) * [[-1,i],[i,1]]."""
import numpy as np
from kernels_x import BasisFun, apply_kernels_matrix
from first_order import b1_oncrack_shear

bf = BasisFun(-0.5, 0)
CC = -1/(2*0.7); NU = 0.3
print("on-crack T of [b1 e_x, b1 e_z]: numeric (Hadamard) vs closed form")
for X in [0.02, 0.05, 0.2, 1.0, 3.0, 10.0]:
    d = apply_kernels_matrix([bf], X, ['xx', 'xz', 'zz'])
    num_x = np.array([d['xx'][0], d['zz'][0]])      # T from input e_x: (T_x, T_z)
    num_z = np.array([d['xz'][0], d['zz'][0]])      # input e_z contributes via xz to T_x
    # careful: T_x(b1 e_z) uses kernel xz ; T_z(b1 e_x) also uses kernel xz (symmetric)
    Tx_from_x = d['xx'][0]; Tz_from_x = d['xz'][0]
    Tx_from_z = d['xz'][0]; Tz_from_z = d['zz'][0]
    M0 = b1_oncrack_shear(X)
    print(f"X={X:7.3f}: from e_x: num=({Tx_from_x:+.5f},{Tz_from_x:+.5f}j) "
          f"closed=({M0[0,0]:+.5f},{M0[1,0]:+.5f})")
    print(f"          from e_z: num=({Tx_from_z:+.5f}j,{Tz_from_z:+.5f}) "
          f"closed=({M0[0,1]:+.5f}j,{M0[1,1]:+.5f})")
