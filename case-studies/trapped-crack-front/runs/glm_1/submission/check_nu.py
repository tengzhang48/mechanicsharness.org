"""Verify the exact closed form R = -1/2 I + (nu/(2-nu)) [[1, i],[i, -1]]
by re-running the whole pipeline at nu = 0.2, 0.25, 0.35."""
import numpy as np
import kernels_x, m_stiffness, first_order as fo

for NU in [0.2, 0.25, 0.35]:
    kernels_x.NU = NU
    fo.NU = NU
    fo.CC = -1.0/(2*(1-NU))
    fo.D1 = 2*(1-NU)*np.sqrt(2/np.pi)
    fo.D3 = 2*np.sqrt(2/np.pi)
    fo.DVEC = np.array([fo.D1, fo.D3])
    sh = fo.solve_shear(ML=6, verbose=False)
    R = sh['R']
    gam = NU/(2-NU)
    Rex = np.array([[-0.5+gam, 1j*gam], [1j*gam, -0.5-gam]])
    print(f"nu={NU}: gamma_exact = {gam:.7f}")
    print("  R numeric:"); print(np.array_str(R, precision=6, suppress_small=True))
    print("  R - R_exact:"); print(np.array_str(R - Rex, precision=6, suppress_small=True))
