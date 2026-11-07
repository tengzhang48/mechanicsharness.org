"""
Time-dependent solver for the grain-boundary stress relaxation, using the
*two-field* form of the model.  This avoids ever forming the doubly
differentiated kernel.

  g(s,tau) = 1 - (2/pi) int_0^1 S(s,x) u_s(x,tau) dx        (elasticity)
  u_tau    = -(1/2) g_ss(s,tau)                              (diffusion)

Unknown: u on Chebyshev-Gauss-Lobatto nodes (interior values; u itself never
needs boundary data because only u_s enters).  The elasticity integral is
evaluated with the p.v. kernel S by singularity subtraction, the diffusion
operator is the Chebyshev D2 matrix.
"""

import numpy as np
from gb import cheb, S_kernel


def elasticity_matrix(nodes, xq, wq):
    """E such that g(node_i) = 1 - sum_j E[i,j] w(node_j).
    Uses  S(a,x) = 1/(a-x) - 1/(a+x) - 2x(a-x)/(a+x)^3  with a p.v. at x=a.
    Applying it to w_x = sum_j phi_j'(x) w_j, the p.v. of the 1/(a-x) term is
    handled by subtracting phi_j'(a):

      fp int_0^1 phi_j'(x)/(a-x) dx
        = int_0^1 [phi_j'(x) - phi_j'(a)]/(a-x) dx + phi_j'(a) * log(a/(1-a))

    The two exceptional nodes are handled explicitly:
      a = 0 : S(0,x) = -1/x, regular (no p.v., the 1/(a-x) term has no pole),
      a = 1 : p.v. is at the *endpoint* x=1, so the crossing log term is zero
              and the subtraction alone gives the finite part.
    """
    from gb2 import lagrange_basis, nodal_derivatives
    _, dphi, _ = lagrange_basis(nodes)
    DPq = dphi(xq)                     # (Nq, n) phi_j'(x)
    DPn, _ = nodal_derivatives(nodes)
    n = len(nodes)
    E = np.zeros((n, n))
    for i, a in enumerate(nodes):
        at = np.abs(xq - a) > 1e-13
        xa, wa = xq[at], wq[at]
        dpa = DPn[i, :]
        if a == 0.0:
            # The s = 0 row is never used: the surface condition g(0) = 0 is
            # imposed directly, so set it to zero.
            E[i, :] = 0.0
            continue
        pv = wa @ ((DPq[at, :] - dpa[None, :]) / (a - xa)[:, None])
        cross = dpa * (np.log(a / (1.0 - a)) if a < 1.0 else 0.0)
        # regular tail -1/(a+x) - 2x(a-x)/(a+x)^3
        reg = (-DPq[at, :] / (a + xa)[:, None]
               - 2.0 * xa[:, None] * (a - xa)[:, None] * DPq[at, :]
               / (a + xa)[:, None] ** 3)
        E[i, :] = (2.0 / np.pi) * (pv + cross + wa @ reg)
    return E