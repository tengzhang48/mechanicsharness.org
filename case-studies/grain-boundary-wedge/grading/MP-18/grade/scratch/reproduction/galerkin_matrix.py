import numpy as np
from numpy.polynomial import chebyshev as C
from numpy.polynomial.legendre import leggauss
from spectral_galerkin import cheb_basis_with_BCs


def cheb_coeff_to_d3_coeff(c):
    """Third derivative (d^3/dx^3) coefficient vector of the Chebyshev series c,
    in the SAME t=2x-1 basis (so values can be evaluated at x)."""
    return C.chebder(C.chebder(C.chebder(c)))


def S_phys(x, y):
    return 1.0 / (y - x) + 1.0 / (x + y) + 2.0 * y * (x - y) / (x + y) ** 3


def galerkin_matrix(Nc, nx=400, ngauss=400):
    """Return G_ij = -<psi_j, (1/pi) Int S_phys(x,y) psi_i'''(y) dy>
    evaluated with high-order Gauss quadrature in x and y (PV handled by shifting
    the x nodes off the singular diagonal since we integrate x LAST against a
    smooth weight; we use the product rule with the exact 'outer' integral).

    Approach: G_ij = -(1/pi) * sum over y-quad of  psi_i'''(y_w) *
                     [ Int_0^1 S_phys(x, y_w) psi_j(x) dx ].
    The bracket is a smooth function of y_w computed by the exact antiderivative
    in x (roles: now x is the integration variable, kernel S_phys(x,y_w)).
    """
    ns = cheb_basis_with_BCs(Nc)
    nb = ns.shape[0]
    # y quadrature nodes (avoid endpoints, include PV handled inside bracket)
    gy, wy = leggauss(ngauss)
    yq = 0.5 * (gy + 1.0)          # in (0,1)
    wy = 0.5 * wy

    # x-quadrature for the bracket Int S_phys(x, yq) psi_j(x) dx
    gx, wx = leggauss(ngauss)
    xq = 0.5 * (gx + 1.0)
    wx = 0.5 * wx

    # precompute psi_j at xq and psi_i''' at yq
    PsiX = np.array([C.chebval(2 * xq - 1, ns[j]) for j in range(nb)])       # (nb,nx)
    Psi3Y = np.zeros((nb, ngauss))
    for i in range(nb):
        c3 = cheb_coeff_to_d3_coeff(ns[i])
        Psi3Y[i] = C.chebval(2 * yq - 1, c3)

    # bracket B[j, k] = Int_0^1 S_phys(x, yq[k]) psi_j(x) dx
    B = np.zeros((nb, ngauss))
    for k in range(ngauss):
        yk = yq[k]
        Sk = S_phys(xq, yk)        # (nx,)   -- PV at x=yq[k] not on grid generally
        for j in range(nb):
            B[j, k] = np.sum(wx * Sk * PsiX[j])

    G = np.zeros((nb, nb))
    Mass = np.zeros((nb, nb))
    for i in range(nb):
        for j in range(nb):
            G[i, j] = -np.sum(wy * B[j, :] * Psi3Y[i, :]) / np.pi
            Mass[i, j] = np.sum(wx * PsiX[j] * PsiX[i])
    return ns, G, Mass


if __name__ == "__main__":
    import scipy.linalg as sla
    for Nc in [12, 16, 20, 24]:
        ns, G, Mass = galerkin_matrix(Nc=Nc, nx=1200, ngauss=1200)
        ev = sla.eig(G, Mass, right=False)
        ev = np.sort(ev.real)[::-1]
        print("Nc=%2d nb=%2d  top real eigs: %s" % (Nc, ns.shape[0], np.round(ev[:5], 5)))