"""
Semi-infinite interface crack between two elastic layers (plane strain),
loaded far behind the tip by an axial force P and bending moment M in layer 1.

Formulation
-----------
Each layer's state = its free-beam far-field stress (an exact elasticity
solution carrying the prescribed resultants) + response to the interface
tractions.  The traction response A*tau therefore carries the whole crack-tip
problem; the far-field conditions become conditions on D := (A tau):

  D(x) = -beta                 on x > 0   (bonded: delta_tot = 0),
  D(x) -> B(x) = alpha1 x e_x + c2 x^2 e_y   as x -> -inf (beam COD far field,
                                with alpha1, c2 exact from the given stresses),
  COD(0-) = 0, tau ~ x^{-1/2+i eps} at 0+, COD ~ x^{1/2+i eps} at 0-.

Represented on x<0 as D = -beta g(x) + B(x) ghat(x) + eta(x) with cutoffs
g = e^{x/lg} (g(0-)=1, g(-inf)=0), ghat = 1 - g (ghat(0-)=0, ghat(-inf)=1),
and eta decaying (basis with |x|^{1/2+i eps} tip behavior).  The Fourier
identity is enforced at every xi != 0 through the regular parts of one-sided
transforms (constants/polynomials on a half line contribute PV 1/xi^n terms):

  A tau_hat - eta_hat + beta (g_hat - i/xi)
      = -alpha1 [x g]^ e_x - c2 [x^2 g]^ e_y + alpha1/xi^2 e_x - 2 i c2/xi^3 e_y.

Bases: complex generalized Laguerre functions
  Phi_n(x)  = x^{a} e^{-x/lt} L_n^{(a)}(x/lt),  a = -1/2 + i eps   (tau, x>0)
  Psi_n(x)  = s^{b} e^{-s/ld} L_n^{(b)}(s/ld),  b =  1/2 + i eps   (eta, x<0)
with analytic transforms (s_+- = 1 +- i xi l):
  Phi_hat_n = l^{a+1} Gamma(a+n+1) (i xi l)^n / (n! (1+i xi l)^{n+a+1})
  Psi_hat_n = l^{b+1} Gamma(b+n+1) (-i xi l)^n / (n! (1-i xi l)^{n+b+1})
The complex traction p+iq is expanded in Phi_n with complex coefficients
(no Re/Im degeneracy; orthogonal basis => well conditioned).

K = lim_{x->0+} sqrt(2 pi x) (sigma22 + i sigma12) e^{-i eps ln x}.

Units: h = 1 (so H = 2), E2 = 1 (E1 = 9), nu = 1/3 for the target problem.
"""
import numpy as np
from compliance import layer_compliance

# ----------------------------------------------------------------------------
# complex gamma (Lanczos, g=7); cross-checked against scipy at runtime
# ----------------------------------------------------------------------------
_LC = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
       771.32342877765313, -176.61502916214059, 12.507343278686905,
       -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]


def cgamma(z):
    z = complex(z)
    if z.real < 0.5:
        return np.pi / (np.sin(np.pi * z) * cgamma(1.0 - z))
    x = z - 1.0
    s = complex(_LC[0])
    for i in range(1, 9):
        s += _LC[i] / (x + i)
    t = x + 7.5
    return np.sqrt(2 * np.pi) * t ** (x + 0.5) * np.exp(-t) * s


# ----------------------------------------------------------------------------
# complex generalized Laguerre polynomials L_n^{(alpha)}(t) by recurrence
# ----------------------------------------------------------------------------
def laguerre_complex(n, alpha, t):
    """Evaluate L_n^{(alpha)}(t) for complex alpha and real t>=0."""
    if n == 0:
        return np.ones_like(t) * 1.0
    if n == 1:
        return 1.0 + alpha - t
    Lm2 = np.ones_like(t) * 1.0
    Lm1 = 1.0 + alpha - t
    for k in range(1, n):
        Lk = ((2 * k + 1 + alpha - t) * Lm1 - (k + alpha) * Lm2) / (k + 1)
        Lm2, Lm1 = Lm1, Lk
    return Lm1


class Problem:
    """Bilayer with an interface crack: geometry/material data and exact
    far-field / energy-release quantities."""

    def __init__(self, h, H, E1, E2, nu):
        self.h, self.H, self.nu = h, H, nu
        self.mu1 = E1 / (2 * (1 + nu))
        self.mu2 = E2 / (2 * (1 + nu))
        self.kap = 3 - 4 * nu
        self.eps = np.log((self.kap / self.mu1 + 1 / self.mu2) /
                          (self.kap / self.mu2 + 1 / self.mu1)) / (2 * np.pi)
        self.Ep1 = E1 / (1 - nu ** 2)      # plane-strain modulus
        self.Ep2 = E2 / (1 - nu ** 2)
        self.I1 = h ** 3 / 12
        self.I2 = H ** 3 / 12

    # ---- combined interfacial compliance in the Fourier domain -----------
    def A(self, xi):
        M1 = layer_compliance(xi, self.h, self.mu1, self.nu)
        M2 = layer_compliance(xi, self.H, self.mu2, self.nu)
        R = np.diag([1.0, -1.0])
        S = np.diag([-1.0, 1.0])
        return M1 - R @ M2 @ S

    # ---- exact far-field (beam) quantities -------------------------------
    def farfield(self, P, M):
        """alpha1 (x-slope of delta_x), c2 (x^2 coeff of delta_y) as x-> -inf."""
        k1 = M / (self.Ep1 * self.I1)                        # curvature layer 1
        k2 = -(M + P * (self.h + self.H) / 2) / (self.Ep2 * self.I2)
        e1 = -P / (self.Ep1 * self.h)                        # axial strains
        e2 = P / (self.Ep2 * self.H)
        # interface-fiber axial strains incl. rotations of the sections
        a1 = (e1 + k1 * self.h / 2) - (e2 - k2 * self.H / 2)
        c2 = (k1 - k2) / 2.0
        return a1, c2

    def G_exact(self, P, M):
        """Exact energy release rate (steady-state argument, dead loads):
        G = strain energy per unit length far behind (ahead is stress free)."""
        N1, N2 = -P, P
        M1s = M
        M2s = -(M + P * (self.h + self.H) / 2)
        return (N1 ** 2 / (2 * self.Ep1 * self.h) + M1s ** 2 / (2 * self.Ep1 * self.I1)
                + N2 ** 2 / (2 * self.Ep2 * self.H) + M2s ** 2 / (2 * self.Ep2 * self.I2))

    def Lam_sum(self):
        return (self.kap + 1) / self.mu1 + (self.kap + 1) / self.mu2

    def Kabs_exact(self, P, M):
        """|K| from exact G via G = (Lam1+Lam2)|K|^2/(16 cosh^2(pi eps))."""
        return np.sqrt(16 * np.cosh(np.pi * self.eps) ** 2
                       * self.G_exact(P, M) / self.Lam_sum())

    def resultants_exact(self, P, M):
        """Interface-traction resultants implied by equilibrium of each layer
        (moment about the tip): int q = P, int p = 0, int x p = -(M + P h/2)."""
        return P, 0.0, -(M + P * self.h / 2)


# ----------------------------------------------------------------------------
# Laguerre-type basis: values, transforms, moments
# ----------------------------------------------------------------------------

class LagBasis:
    """Phi_n(x) = x^w e^{-x/l} L_n^{(w)}(x/l) on x>0  (w complex)."""

    def __init__(self, l, w):
        self.l = float(l)
        self.w = complex(w)

    def val(self, n, x):
        t = x / self.l
        return x ** self.w * np.exp(-t) * laguerre_complex(n, self.w, t)

    def that(self, n, xi):
        """Fourier transform (convention e^{-i xi x}) of Phi_n."""
        w = self.w
        l = self.l
        s = 1.0 + 1j * xi * l
        return (l ** (w + 1) * cgamma(w + n + 1) / np.math.factorial(n)
                * (1j * xi * l) ** n / s ** (n + w + 1))

    def that_conj(self, n, xi):
        """conj( Phi_hat_n(-xi) ) = same formula with conjugated order."""
        w = self.w.conjugate()
        l = self.l
        s = 1.0 + 1j * xi * l
        return (l ** (w + 1) * cgamma(w + n + 1) / np.math.factorial(n)
                * (1j * xi * l) ** n / s ** (n + w + 1))

    def moments(self, nmax):
        """int Phi_n dx and int x Phi_n dx (complex), by adaptive quadrature
        (Gauss-Laguerre nodes are inaccurate for the singular complex weight)."""
        from scipy.integrate import quad
        m0 = np.zeros(nmax + 1, complex)
        m1 = np.zeros(nmax + 1, complex)
        for n in range(nmax + 1):
            f = lambda x: x ** self.w * np.exp(-x / self.l) * laguerre_complex(n, self.w, x / self.l)
            m0[n] = (quad(lambda x: f(x).real, 0, 400 * self.l, limit=500)[0]
                     + 1j * quad(lambda x: f(x).imag, 0, 400 * self.l, limit=500)[0])
            m1[n] = (quad(lambda x: (x * f(x)).real, 0, 400 * self.l, limit=500)[0]
                     + 1j * quad(lambda x: (x * f(x)).imag, 0, 400 * self.l, limit=500)[0])
        return m0, m1


class ReflLagBasis(LagBasis):
    """Psi_n(x) = |x|^w e^{-|x|/l} L_n^{(w)}(|x|/l) on x<0 (w = 1/2 + i eps)."""

    def val(self, n, x):
        s = -x
        t = s / self.l
        return s ** self.w * np.exp(-t) * laguerre_complex(n, self.w, t)

    def that(self, n, xi):
        w = self.w
        l = self.l
        s = 1.0 - 1j * xi * l
        return (l ** (w + 1) * cgamma(w + n + 1) / np.math.factorial(n)
                * (-1j * xi * l) ** n / s ** (n + w + 1))

    def that_conj(self, n, xi):
        w = self.w.conjugate()
        l = self.l
        s = 1.0 - 1j * xi * l
        return (l ** (w + 1) * cgamma(w + n + 1) / np.math.factorial(n)
                * (-1j * xi * l) ** n / s ** (n + w + 1))


# ----------------------------------------------------------------------------
# main solve
# ----------------------------------------------------------------------------

def default_grid():
    return np.concatenate([np.geomspace(2.5e-3, 1.0, 150),
                           np.geomspace(1.0, 600.0, 170)])


def solve_case(pr, P, M, Nt=12, Nd=12, lt=0.9, ld=1.1, lg=1.0,
               xigrid=None, weight='A', rcond=1e-10, impose_res=True,
               wres=1e4):
    eps = pr.eps
    a1, c2 = pr.farfield(P, M)
    if xigrid is None:
        xigrid = default_grid()
    a_g = 1.0 / lg

    tauB = LagBasis(lt, -0.5 + 1j * eps)     # traction basis, x > 0
    etaB = ReflLagBasis(ld, 0.5 + 1j * eps)  # COD correction basis, x < 0

    # unknown layout (real vector): tau: a_n, b_n (c_n = a_n + i b_n), n<Nt
    #   eta: d_n, e_n ; beta: 2  => 2Nt + 2Nd + 2
    ncols = 2 * Nt + 2 * Nd + 2
    i_ta, i_tb = 0, Nt
    i_ed, i_ee = 2 * Nt, 2 * Nt + Nd
    i_bx, i_by = ncols - 2, ncols - 1

    nxi = len(xigrid)
    A_sys = np.zeros((4 * nxi, ncols))
    b_sys = np.zeros(4 * nxi)

    # precompute basis transforms per xi
    for m, xi in enumerate(xigrid):
        A = pr.A(xi)
        w = np.linalg.norm(A) if weight == 'A' else 1.0
        r0 = 4 * m

        def add(col, cvec2):
            A_sys[r0, col] = np.real(cvec2[0]) / w
            A_sys[r0 + 1, col] = np.imag(cvec2[0]) / w
            A_sys[r0 + 2, col] = np.real(cvec2[1]) / w
            A_sys[r0 + 3, col] = np.imag(cvec2[1]) / w

        # tau columns.  tau_c = p + iq = sum c_n Phi_n (c_n = a_n + i b_n),
        # complex-valued function; the transforms of the REAL functions are
        #   p_hat = [T(xi) + conj(T(-xi))]/2,  q_hat = [T(xi) - conj(T(-xi))]/(2i)
        # with T(xi) = sum c_n Phi_hat_n(xi).  (NOT Re/Im of T!  Phi_hat(-xi)
        # != conj(Phi_hat(xi)) for complex order alpha.)
        for n in range(Nt):
            F = tauB.that(n, xi)
            G = tauB.that_conj(n, xi)
            Sm = 0.5 * (F + G)          # (F+G)/2
            Df = -0.5j * (F - G)        # -i(F-G)/2
            add(i_ta + n, A @ np.array([Df, Sm]))   # coefficient a_n
            add(i_tb + n, A @ np.array([Sm, -Df]))  # coefficient b_n
        # eta columns: eta = eta_x + i eta_y = sum (d_n + i e_n) Psi_n,
        # same structure as the tau side:
        #   eta_x_hat = d*Sm + e*(-Df),  eta_y_hat = d*Df + e*Sm
        for n in range(Nd):
            H = etaB.that(n, xi)
            Kc = etaB.that_conj(n, xi)
            Sm = 0.5 * (H + Kc)
            Df = -0.5j * (H - Kc)
            add(i_ed + n, -np.array([Sm, Df]))     # coefficient d_n (eta_x)
            add(i_ee + n, np.array([Df, -Sm]))     # coefficient e_n (eta_y)
        # beta columns
        gcol = 1.0 / (a_g - 1j * xi) - 1j / xi
        add(i_bx, np.array([gcol, 0j]))
        add(i_by, np.array([0j, gcol]))
        # RHS
        part = np.array([-a1 * (-(a_g - 1j * xi) ** -2) + a1 / xi ** 2,
                         -c2 * (2.0 * (a_g - 1j * xi) ** -3) - 2j * c2 / xi ** 3])
        b_sys[r0] = np.real(part[0]) / w
        b_sys[r0 + 1] = np.imag(part[0]) / w
        b_sys[r0 + 2] = np.real(part[1]) / w
        b_sys[r0 + 3] = np.imag(part[1]) / w

    # ---- exact resultant constraints (equilibrium; they are equivalent to
    # the pole-matching conditions at xi -> 0 but far better conditioned as
    # explicit constraints):  int p = 0, int q = P, int x p = -(M + P h/2).
    if impose_res:
        m0, m1 = tauB.moments(Nt - 1)
        C = np.zeros((3, ncols))
        # int p = Re sum c_n m0_n ; int q = Im sum c_n m0_n
        # int x p = Re sum c_n m1_n
        C[0, i_ta:i_ta + Nt] = m0.real
        C[0, i_tb:i_tb + Nt] = -m0.imag
        C[1, i_ta:i_ta + Nt] = m0.imag
        C[1, i_tb:i_tb + Nt] = m0.real
        C[2, i_ta:i_ta + Nt] = m1.real
        C[2, i_tb:i_tb + Nt] = -m1.imag
        d = np.array([0.0, P, -(M + P * pr.h / 2.0)])
        A_sys = np.vstack([A_sys, wres * C])
        b_sys = np.concatenate([b_sys, wres * d])

    sol, res, rank, sv = np.linalg.lstsq(A_sys, b_sys, rcond=rcond)
    resid_rms = float(np.sqrt(np.mean((A_sys @ sol - b_sys) ** 2)))
    cond = float(sv[0] / sv[rank - 1])

    # ---- post-processing -------------------------------------------------
    cn = sol[i_ta:i_ta + Nt] + 1j * sol[i_tb:i_tb + Nt]      # tau coeffs
    dn = sol[i_ed:i_ed + Nd] + 1j * sol[i_ee:i_ee + Nd]      # eta coeffs

    def tauc(x):
        tot = 0j
        for n in range(Nt):
            tot += cn[n] * tauB.val(n, x)
        return tot

    def K_of(x):
        return np.sqrt(2 * np.pi * x) * tauc(x) * np.exp(-1j * eps * np.log(x))

    K = K_of(1e-3)
    K_drift = max(abs(K_of(1e-4) - K), abs(K_of(3.2e-3) - K)) / abs(K)

    # resultants (checks; not imposed in the solve)
    m0, m1 = tauB.moments(Nt - 1)
    msum0 = np.sum(cn * m0)
    msum1 = np.sum(cn * m1)
    rq = msum0.imag            # int q
    rp = msum0.real            # int p
    rxp = msum1.real           # int x p

    # COD reconstruction on x<0:  COD = beta(1-g) + B*ghat + eta
    bx, by = sol[i_bx], sol[i_by]

    def cod(x):
        g = np.exp(x / lg)
        dx = bx * (1 - g) + a1 * x * (1 - g)
        dy = by * (1 - g) + c2 * x * x * (1 - g)
        for n in range(Nd):
            z = dn[n] * etaB.val(n, x)
            dx += z.real
            dy += z.imag
        return dx, dy

    return dict(
        K=K, Kabs=abs(K), psi_deg=float(np.degrees(np.angle(K))),
        K_drift=float(K_drift),
        res_q=float(rq), res_p=float(rp), res_xp=float(rxp),
        resid_rms=resid_rms, cond=cond,
        beta=(float(bx), float(by)),
        alpha1=float(a1), c2=float(c2),
        cod_minus20=cod(-20.0), cod_minus10=cod(-10.0),
        tauc_1e3=tauc(1e-3), tauc_1e2=tauc(1e-2), tauc_1e1=tauc(1e-1),
    )
