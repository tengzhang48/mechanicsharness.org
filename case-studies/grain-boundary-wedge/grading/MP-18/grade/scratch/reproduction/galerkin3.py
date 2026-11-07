import numpy as np
from phys_oper import M_phys
from spectral_galerkin import cheb_basis_with_BCs
from numpy.polynomial import chebyshev as C2


def op_poly_apply(ns, xfine, nodes):
    """Given BC-satisfying Chebyshev basis coeffs ns (nb x (Nc+1)), compute the
    operator  (A psi_i)(x) = -(1/pi) Int S(x,y) psi_i'''(y) dy  on the fine grid
    xfine, using M_phys (validated exact hats).  psi_i''' is a polynomial; we
    represent psi_i''' as a combination of hats on `nodes` (fine grid) EXACTLY by
    interpolation of the polynomial at the nodes (high degree poly on fine grid).
    Returns Apsi (nb, len(xfine)) and Psi (nb, len(xfine))."""
    nb = ns.shape[0]
    xcols = nodes[1:-1]
    M = M_phys(xcols, nodes)          # (M, M): rows xcols, cols = hat coeffs at nodes
    Psi = np.zeros((nb, len(xfine)))
    Apsi = np.zeros((nb, len(xfine)))
    # we need M applied to values of psi''' sampled on the FULL node hat basis.
    # M maps hat-coefficients (at nodes[1..M]) to values at xcols.  Sampling
    # psi''' at nodes[1..M] gives its hat representation (interpolation).
    for i in range(nb):
        c3 = C2.chebder(C2.chebder(C2.chebder(ns[i])))
        p3_vals = C2.chebval(2 * nodes[1:-1] - 1, c3)     # psi_i''' at interior nodes
        Aw = M @ p3_vals                                   # (Int S psi_i''') at xcols
        # evaluate Aw at xfine by interpolation (fine grid, dense)
        Apsi[i] = np.interp(xfine, xcols, Aw)
        Psi[i] = C2.chebval(2 * xfine - 1, ns[i])
    return Psi, Apsi


if __name__ == "__main__":
    import scipy.linalg as sla
    from numpy.polynomial.legendre import leggauss

    def run(Nc, Mfine=400, nx=400):
        ns = cheb_basis_with_BCs(Nc)
        nb = ns.shape[0]
        xf = np.linspace(0, 1, Mfine + 2)
        nodes = xf
        Psi, Apsi = op_poly_apply(ns, xf[1:-1], nodes)
        xfine = xf[1:-1]
        # Galerkin: G_ij = sum_x w * Psi_j * Apsi_i ; Mass_ij = <Psi_j,Psi_i>
        w = np.ones(len(xfine))
        w[0] = w[-1] = 0.5
        w = w / len(w) * 1.0
        # trapezoid weights on uniform grid length Mfine
        dx = 1.0 / (Mfine + 1)
        w = np.full(len(xfine), dx)
        G = (Psi * w) @ Apsi.T
        Mass = (Psi * w) @ Psi.T
        return ns, G, Mass

    for Nc in [8, 12, 16, 20]:
        ns, G, Mass = run(Nc)
        ev = sla.eig(G, Mass, right=False).real
        ev = np.sort(ev)[::-1]
        print("Nc=%2d nb=%2d top eigs: %s" % (Nc, ns.shape[0], np.round(ev[:5], 5)))