"""Eigenvalue structure of the interface-crack wedge problem (justifies the form of the traction fit):
roots of det M(lam) with 0 < Re lam < 3.2 and the interface traction s_tt + i s_rt they carry."""
import numpy as np
from asym import eig_matrix, find_root, mu_kappa, bimat_eps, _rows_stress

for (E1, nu1, E2, nu2) in [(5.0, 0.22, 1.0, 0.35), (9.0, 1/3, 1.0, 1/3)]:
    mu1, k1 = mu_kappa(E1, nu1); mu2, k2 = mu_kappa(E2, nu2)
    eps = bimat_eps(E1, nu1, E2, nu2)
    roots = []
    for re in np.arange(0.1, 3.2, 0.1):
        for im in np.arange(-0.3, 0.31, 0.1):
            try:
                l = find_root(re + 1j*im, mu1, k1, mu2, k2)
            except Exception:
                continue
            if 0.01 < l.real < 3.2 and abs(l.imag) < 1 and np.isfinite(l):
                if abs(np.linalg.det(eig_matrix(l, mu1, k1, mu2, k2))) < 1e-10:
                    if all(abs(l - r) > 1e-6 for r in roots):
                        roots.append(l)
    roots.sort(key=lambda z: (z.real, z.imag))
    print(f"E1/E2={E1/E2:g}  eps={eps:.6f}")
    for l in roots:
        M = eig_matrix(l, mu1, k1, mu2, k2)
        U, S, Vh = np.linalg.svd(M)
        nnull = int(np.sum(S < 1e-9 * S[0]))
        _, stt0, srt0 = _rows_stress(l, 0.0)
        info = []
        for kk in range(1, nnull + 1):
            v = Vh[-kk].conj()
            a, b = stt0 @ v[0:4], srt0 @ v[0:4]
            info.append(f"|s_tt|={abs(a):.2e} s_rt/s_tt={b/a if abs(a)>1e-12 else float('nan'):.4f}")
        print(f"   lam={l.real:.6f}{l.imag:+.6f}i  null-dim={nnull}  ", "; ".join(info))
