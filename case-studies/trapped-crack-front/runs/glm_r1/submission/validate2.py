import numpy as np, time, sys
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

def run(tag, **kw):
    t0=time.time()
    p = WavyFrontProblem(**kw)
    p.solve()
    th, K1, K2, K3 = p.K_along_front()
    e1 = (K1.mean()-1.0); e2=(K2.mean()-0.6); e3=(K3.mean()-0.8)
    print(f"{tag}: relres {p.fitinfo['relres']:.2e} "
          f"K1err {e1:+.4f} K2err {e2:+.4f} K3err {e3:+.4f} "
          f"(dev {np.abs(K1-K1.mean()).max():.1e}) [{time.time()-t0:.0f}s]")

print("--- V1 straight, convergence in M/Nxf/margin ---")
run("M=12 Nx=4096  ", A=0.0, harm=1, Nxf=4096, Nz=64, M=12, Nzmode=4, Lcr=2.0, Lx=3.5)
run("M=16 Nx=4096  ", A=0.0, harm=1, Nxf=4096, Nz=64, M=16, Nzmode=4, Lcr=2.0, Lx=3.5)
run("M=20 Nx=4096  ", A=0.0, harm=1, Nxf=4096, Nz=64, M=20, Nzmode=4, Lcr=2.0, Lx=3.5)
run("M=24 Nx=4096  ", A=0.0, harm=1, Nxf=4096, Nz=64, M=24, Nzmode=4, Lcr=2.0, Lx=3.5)
run("M=20 Nx=8192  ", A=0.0, harm=1, Nxf=8192, Nz=64, M=20, Nzmode=4, Lcr=2.0, Lx=3.5)
print("--- V1 with larger margin ---")
import crack3d
for mdx in [3.0, 6.0, 10.0]:
    p = WavyFrontProblem(A=0.0, harm=1, Nxf=4096, Nz=64, M=20, Nzmode=4, Lcr=2.0, Lx=3.5)
    p.solve() if mdx==3.0 else None
# patch margin via monkey method
for mdx in [6.0, 10.0]:
    p = WavyFrontProblem(A=0.0, harm=1, Nxf=4096, Nz=64, M=20, Nzmode=4, Lcr=2.0, Lx=3.5)
    import types
    def solve_m(self, margin_dx=mdx, verbose=False, rcond=1e-11):
        ii, jj = self.collocation(margin_dx=margin_dx)
        Nc = len(ii); xc = self.X[ii, jj]
        rhs = np.zeros((3, Nc))
        for c in range(3):
            rhs[c] = np.interp(xc, self.xg, self.t_win[c])
        nb = 3 * self.M * (2 * self.Nzm + 1)
        Amat = np.zeros((3 * Nc, nb))
        H, u = self._window_u()
        col = 0
        for j in range(3):
            for m in range(self.M):
                for n in range(-self.Nzm, self.Nzm + 1):
                    g = self.basis_grid(j, m, n, H, u)
                    tx, ty, tz = self.op.apply(*g)
                    Amat[:, col] = np.concatenate(
                        [tx[ii, jj].real, ty[ii, jj].real, tz[ii, jj].real])
                    col += 1
        b = np.concatenate([rhs[0], rhs[1], rhs[2]])
        scale = max(np.abs(b).max(), 1e-12) * self.cw
        mm = np.arange(self.M); sgn = (-1.0) ** mm
        ns = np.arange(-self.Nzm, self.Nzm + 1); nns = len(ns)
        rows = []
        for j in range(3):
            for i_n in range(nns):
                row = np.zeros(nb)
                for m in range(self.M):
                    row[(j * self.M + m) * nns + i_n] = sgn[m]
                rows.append(row * scale)
        Amat = np.vstack([Amat, np.array(rows)])
        b = np.concatenate([b, np.zeros(len(rows))])
        c, *_ = np.linalg.lstsq(Amat, b, rcond=rcond)
        self.c = c.reshape(3, self.M, 2 * self.Nzm + 1)
        self.fitinfo = dict(shape=Amat.shape,
                            relres=np.linalg.norm(Amat @ c - b) / np.linalg.norm(b))
        return self.fitinfo
    p.solve = types.MethodType(solve_m, p)
    p.solve()
    th, K1, K2, K3 = p.K_along_front()
    print(f"margin {mdx}: relres {p.fitinfo['relres']:.2e} K1err {K1.mean()-1:+.4f} "
          f"K2err {K2.mean()-0.6:+.4f} K3err {K3.mean()-0.8:+.4f}")
print("--- V2 shift invariance (shift=0.05, A=0, harm=1) ---")
run("shift 0.05    ", A=0.0, harm=1, shift=0.05, Nxf=4096, Nz=64, M=20, Nzmode=4, Lcr=2.0, Lx=3.5)
