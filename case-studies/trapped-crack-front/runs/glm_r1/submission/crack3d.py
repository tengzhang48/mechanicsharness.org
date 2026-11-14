"""
Spectral dislocation solver v6 for a half-plane crack x < a(z), y = 0 in an
infinite isotropic solid, loaded from infinity by the pure straight-front
K-field (uniform K1^0, K2^0, K3^0; zero T-stress).

Exact quasistatic traction-dislocation operator, Fourier mode k=(kx,kz),
kp=(kx,0,kz)/|k|, tp=(kz,0,-kx)/|k|, E'=E/(1-nu^2):
  t = -(|k|/2)[ (E'/2)(b.y)y + (E'/2)(b.kp)kp + mu (b.tp)tp ]
(kz=0 limit = classical 2D DDT: t_x=-(E'/4)|kx|b_x, t_y=-(E'/4)|kx|b_y,
 t_z=-(mu/2)|kx|b_z; DDT gives T[b0] = {0 on x<0, +t_inf on x>0}).

Density:  b = W(sigma,z) * sigma * S(sigma,z) + b_tail(x)
  sigma = sqrt(a(z)-x), u = 2 sigma/sqrt(Lcr)-1,
  W = H(a-x)*env(s)  (env: 1 for s<0.8Lcr, 0 at s=Lcr),
  S = sum_{m, n>=0} [alpha_mn cos(n theta) + beta_mn sin(n theta)] T_m(u),
  theta = 2 pi harm z / lam,
  b_tail = C0 sqrt(-x) H(-x) (1-chi(x)), chi = smootherstep, 1 at x>-c1, 0
  below -c2. Exact traction of the far piece:
  T[b_tail] = t_inf 1_{x>0} - T_grid[b0 chi].

Collocation of total traction = 0 on x < a(z):
  T_grid[W sigma S](xc) = T_grid[b0 chi](xc) - t_inf(xc) 1_{xc>0}
with error-model weights and a ridge toward the straight-crack solution
(penalizing high modes) for stability. SVD-truncated solve.

K extraction: fit the total density near the front on grid points:
  b_tot/sqrt(s) = C + d sqrt(s) + e s,  s = a(z)-x in (s_lo, s1);
K1=(E' sqrt(2pi)/8) C_y ; local frame n=(1,0,-a')/sqrt(1+a'^2), t=(a',0,1)/...:
K2=(E' sqrt(2pi)/8)(C.n), K3=(mu sqrt(2pi)/4)(C.t).
"""

import numpy as np


def material(E=1.0, nu=0.3):
    mu = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    Ep = E / (1.0 - nu ** 2)
    return dict(E=E, nu=nu, mu=mu, lam=lam, Ep=Ep)


class Operator2D:
    def __init__(self, shape, dx, dz, Ep, mu):
        Nx, Nz = shape
        kx = 2.0 * np.pi * np.fft.fftfreq(Nx, d=dx)
        kz = 2.0 * np.pi * np.fft.fftfreq(Nz, d=dz)
        KX = kx[:, None] + 0.0 * kz[None, :]
        KZ = 0.0 * kx[:, None] + kz[None, :]
        K = np.hypot(KX, KZ)
        self.zero00 = (K == 0)
        K = np.where(self.zero00, 1.0, K)
        self.kpx, self.kpz = KX / K, KZ / K
        self.tpx, self.tpz = KZ / K, -KX / K
        self.fL = -0.5 * K * (Ep / 2.0)
        self.fT = -0.5 * K * mu

    def apply(self, bx, by, bz):
        Bx = np.fft.fft2(bx, axes=(-2, -1))
        By = np.fft.fft2(by, axes=(-2, -1))
        Bz = np.fft.fft2(bz, axes=(-2, -1))
        bL = Bx * self.kpx + Bz * self.kpz
        bT = Bx * self.tpx + Bz * self.tpz
        tx = self.fL * bL * self.kpx + self.fT * bT * self.tpx
        tz = self.fL * bL * self.kpz + self.fT * bT * self.tpz
        ty = self.fL * By
        for arr in (tx, ty, tz):
            arr[self.zero00] = 0.0
        return (np.fft.ifft2(tx, axes=(-2, -1)),
                np.fft.ifft2(ty, axes=(-2, -1)),
                np.fft.ifft2(tz, axes=(-2, -1)))


def smootherstep(x, x1, x2):
    """1 for x>=x1, 0 for x<=x2 (x2<x1), C2-smooth in between."""
    t = np.clip((x - x2) / (x1 - x2), 0.0, 1.0)
    return t ** 3 * (t * (6.0 * t - 15.0) + 10.0)


class WavyFrontProblem:
    def __init__(self, A, lam=1.0, harm=1, K0=(1.0, 0.6, 0.8), shift=0.0,
                 mat=None, Lcr=2.0, Lx=3.5, Nxf=4096, Nz=64,
                 M=18, Nzmode=4, cw=30.0, s1=0.12, taper=(0.30, 0.95)):
        self.A, self.lam, self.harm, self.shift = A, lam, harm, shift
        self.K10, self.K20, self.K30 = K0
        m = mat or material()
        self.Ep, self.mu = m['Ep'], m['mu']
        s2p = np.sqrt(2.0 * np.pi)
        self.C0 = np.array([8.0 * self.K20 / (self.Ep * s2p),
                            8.0 * self.K10 / (self.Ep * s2p),
                            4.0 * self.K30 / (self.mu * s2p)])
        self.Lcr, self.Lx = Lcr, Lx
        self.Nx, self.Nz = Nxf, Nz
        self.dx = 2.0 * Lx / Nxf
        self.dz = lam / Nz
        self.xg = (np.arange(Nxf) - Nxf // 2) * self.dx
        self.zg = np.arange(Nz) * self.dz
        self.X, self.Z = np.meshgrid(self.xg, self.zg, indexing='ij')
        self.M, self.Nzm = M, Nzmode
        self.cw = cw
        self.s1 = s1
        self.theta = 2.0 * np.pi * harm * self.Z / lam
        self.a = shift + A * np.cos(self.theta)
        self.c1 = -taper[0] * Lcr
        self.c2 = -taper[1] * Lcr
        self.op = Operator2D((Nxf, Nz), self.dx, self.dz, self.Ep, self.mu)
        self._precompute()

    def _precompute(self):
        chi = smootherstep(self.xg, self.c1, self.c2)
        self.b0chi_x = np.where(self.xg < 0,
                                np.sqrt(np.maximum(-self.xg, 0.0)) * chi, 0.0)
        grid = [self.b0chi_x[:, None] * np.ones((1, self.Nz)) * self.C0[j]
                for j in range(3)]
        t_win = self.op.apply(*grid)
        self.t_win = [t_win[j][:, 0].real for j in range(3)]

    def t_inf(self, x):
        w = 1.0 / np.sqrt(2.0 * np.pi * np.abs(x))
        return (self.K20 * w, self.K10 * w, self.K30 * w)

    # ---- geometry / basis ------------------------------------------------
    def _window_u(self):
        s = self.a - self.X
        H = (s > 0).astype(float)
        sig = np.sqrt(np.maximum(s, 0.0))
        u = np.clip(2.0 * sig / np.sqrt(self.Lcr) - 1.0, -1.0, 1.0)
        env = 1.0 - smootherstep(s, self.Lcr, 0.8 * self.Lcr)
        return H * env, u, sig

    def n_unknowns(self):
        return 3 * self.M * (2 * self.Nzm + 1)

    def basis_columns(self, ii, jj):
        """Real basis columns: for each (j,m,n>=0): cos-column and (n>0)
        sin-column. Column order matches self.unpack/ridge weights."""
        W, u, sig = self._window_u()
        cols = []
        arcc = np.arccos(np.clip(u, -1, 1))
        for j in range(3):
            for m in range(self.M):
                Tm = np.cos(m * arcc) * W * sig
                for n in range(0, self.Nzm + 1):
                    f = Tm * np.exp(1j * n * self.theta)
                    g = [np.zeros_like(f) for _ in range(3)]
                    g[j] = f
                    tx, ty, tz = self.op.apply(*g)
                    v = np.concatenate([tx[ii, jj], ty[ii, jj], tz[ii, jj]])
                    cols.append(v.real)
                    if n > 0:
                        cols.append(v.imag)
        return np.array(cols).T

    def mode_weights(self):
        """Ridge weight per unknown (order matches basis_columns)."""
        mm = np.arange(self.M)
        Wm = (1.0 + mm) ** 2
        w = []
        for j in range(3):
            for m in range(self.M):
                for n in range(0, self.Nzm + 1):
                    w.append(Wm[m] * (1.0 + 2.0 * n) ** 2)
                    if n > 0:
                        w.append(Wm[m] * (1.0 + 2.0 * n) ** 2)
        return np.array(w)

    # ---- collocation -------------------------------------------------------
    def collocation(self, margin_dx=15.0, nx_per_row=50):
        margin = margin_dx * self.dx
        ii_list, jj_list = [], []
        for jr in range(self.Nz):
            arow = self.a[0, jr]
            s = arow - self.xg
            ok = np.where((s > margin) & (s < 0.999 * self.Lcr))[0]
            if len(ok) < 4:
                continue
            smax, smin = 0.999 * self.Lcr, margin
            t = np.cos(np.pi * (np.arange(nx_per_row) + 0.5) / nx_per_row)
            s_target = 0.5 * (smax + smin) + 0.5 * (smax - smin) * t
            x_target = arow - s_target
            idx = np.clip(np.searchsorted(self.xg, x_target), ok[0] + 1, ok[-1])
            good = []
            for i0 in idx:
                i = i0
                while i > ok[0] and not (s[i] > margin):
                    i -= 1
                if s[i] > margin:
                    good.append(i)
            good = np.unique(good)
            ii_list.extend(good.tolist())
            jj_list.extend([jr] * len(good))
        return np.array(ii_list, dtype=int), np.array(jj_list, dtype=int)

    # ---- solve ---------------------------------------------------------------
    def solve(self, verbose=False, cut=1e-9, margin_dx=15.0, ridge=3e-3,
              prior=None):
        ii, jj = self.collocation(margin_dx=margin_dx)
        Nc = len(ii)
        xc = self.X[ii, jj]
        sc = self.a[0, jj] - xc
        sw = margin_dx * self.dx
        w = (np.maximum(sc, sw) / sw) ** 1.5
        rhs = np.zeros((3, Nc))
        pos = xc > 0
        tinf_pos = self.t_inf(np.where(pos, xc, 1.0))
        for c in range(3):
            rhs[c] = np.interp(xc, self.xg, self.t_win[c])
            rhs[c][pos] -= tinf_pos[c][pos]
        A = self.basis_columns(ii, jj)                # (3Nc, nb)
        w3 = np.concatenate([w, w, w])
        A = A * w3[:, None]
        b = np.concatenate([rhs[0], rhs[1], rhs[2]]) * w3
        nb = A.shape[1]
        if ridge is not None and ridge > 0:
            Wmode = self.mode_weights()
            cnorm = np.median(np.linalg.norm(A, axis=0))
            rho = ridge * cnorm
            r0 = np.zeros(nb) if prior is None else np.asarray(prior).ravel()
            A = np.vstack([A, np.diag(np.sqrt(Wmode)) * rho])
            b = np.concatenate([b, np.sqrt(Wmode) * rho * r0])
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        keep = S > cut * S[0]
        c = Vt[keep].T @ ((U[:, keep].T @ b) / S[keep])
        self.c = c
        self.fitinfo = dict(shape=A.shape, nsv=int(keep.sum()),
                            relres=np.linalg.norm(A @ c - b) / np.linalg.norm(b))
        if verbose:
            print(f"  solve {A.shape}: sv kept {keep.sum()}/{len(S)}, "
                  f"relres {self.fitinfo['relres']:.2e}")
        return self.fitinfo

    def solve_prior(self, ridge=1e-4, **kw):
        p0 = WavyFrontProblem(A=0.0, lam=self.lam, harm=self.harm,
                              K0=(self.K10, self.K20, self.K30),
                              shift=self.shift, Lcr=self.Lcr, Lx=self.Lx,
                              Nxf=self.Nx, Nz=self.Nz, M=self.M,
                              Nzmode=self.Nzm, s1=self.s1,
                              taper=(-self.c1 / self.Lcr, -self.c2 / self.Lcr))
        p0.solve(ridge=ridge, **kw)
        return p0.c

    # ---- reconstruction / extraction ------------------------------------------
    def S_coeff(self):
        """Return S(u=.) coefficient array as (3, M, Nzm+1, 2): cos, sin."""
        return self.c.reshape(3, self.M, 2 * self.Nzm + 1)

    def total_density(self):
        W, u, sig = self._window_u()
        arcc = np.arccos(np.clip(u, -1, 1))
        Cf = self.S_coeff()
        b_tot = [np.zeros(self.X.shape) for _ in range(3)]
        chi = smootherstep(self.xg, self.c1, self.c2)
        far = np.where(self.xg < 0,
                       np.sqrt(np.maximum(-self.xg, 0.0)) * (1.0 - chi), 0.0)
        for j in range(3):
            acc = np.zeros(self.X.shape)
            for m in range(self.M):
                Tm = np.cos(m * arcc)
                base = Tm * W * sig
                a0 = Cf[j, m, 0]
                acc += base * a0
                for n in range(1, self.Nzm + 1):
                    ac = Cf[j, m, 2 * n - 1]
                    bs = Cf[j, m, 2 * n]
                    if ac != 0.0 or bs != 0.0:
                        acc += base * (ac * np.cos(n * self.theta)
                                       + bs * np.sin(n * self.theta))
            b_tot[j] = acc + (far * self.C0[j])[:, None]
        return b_tot

    def K_along_front(self, s_lo=None):
        b_tot = self.total_density()
        s2p = np.sqrt(2.0 * np.pi)
        C = np.zeros((3, self.Nz))
        s1 = self.s1
        if s_lo is None:
            s_lo = max(self.s1 - 0.08, 8.0 * self.dx)
        for jr in range(self.Nz):
            aj = self.a[0, jr]
            s = aj - self.xg
            sel = np.where((s > s_lo) & (s < s1))[0]
            if len(sel) < 8:
                sel = np.where(s > 4.0 * self.dx)[0][:40]
            sq = np.sqrt(s[sel])
            Am = np.column_stack([np.ones_like(sq), sq, sq ** 2])
            for j in range(3):
                yv = b_tot[j][sel, jr] / sq
                cf, *_ = np.linalg.lstsq(Am, yv, rcond=None)
                C[j, jr] = cf[0]
        aph = -self.A * (2.0 * np.pi * self.harm / self.lam) * np.sin(
            2.0 * np.pi * self.harm * self.zg / self.lam)
        rot = 1.0 / np.sqrt(1.0 + aph ** 2)
        K1 = (self.Ep * s2p / 8.0) * C[1]
        K2 = (self.Ep * s2p / 8.0) * (C[0] - aph * C[2]) * rot
        K3 = (self.mu * s2p / 4.0) * (aph * C[0] + C[2]) * rot
        th = 2.0 * np.pi * self.harm * self.zg / self.lam
        return th, K1, K2, K3

    def fit_partA(self):
        th, K1, K2, K3 = self.K_along_front()
        c, s = np.cos(th), np.sin(th)
        out = {}
        for name, K in (('1', K1), ('2', K2), ('3', K3)):
            M = np.column_stack([np.ones_like(th), c, s])
            coef, *_ = np.linalg.lstsq(M, K, rcond=None)
            out['K0_' + name] = coef[0]
            if self.A > 0:
                out['C' + name] = coef[1] * self.lam / (self.A * self.K10)
                out['S' + name] = coef[2] * self.lam / (self.A * self.K10)
        return out
