import numpy as np
from fd_nonuniform import build_D1D2
from statement_oper import Smat_statement
import oper


def operator_graded(nodes):
    """Reduced decay operator  A_u = -D2 M D1  on the interior nodes nodes[1..M].
    Equivalently in s: A_s = -M D1 D2.  Both give the same nonzero eigenvalues
    (they are similar via M^{-1} up to the constant mode); we use A_u which is
    the physically clean chain  u_tau = D2(1 - M D1 u).
    M is built from the PROBLEM-STATEMENT kernel S(x,y) (validated against
    direct PV-aware Gauss quadrature)."""
    xcols = np.asarray(nodes[1:-1], float)
    D1, D2 = build_D1D2(nodes)
    M = Smat_statement(xcols, nodes)
    A_u = -(D2 @ M @ D1)
    return M, D1, D2, A_u


if __name__ == "__main__":
    import grids
    # validation: on a graded grid, check M matches direct PV quadrature
    x = grids.graded_nodes(12, 3.0)
    xcols = x[1:-1]
    M = Smat_statement(xcols, x)
    Mn = oper.M_hat_numeric(xcols, x, ngauss=1000)
    print("graded M exact vs numeric PV: max diff", np.max(np.abs(M - Mn)))
    # eigenvalues
    M, D1, D2, A = operator_graded(x)
    ev = np.linalg.eigvals(A)
    print("N=12 graded p=3 slowest decay:", -np.max(ev.real))
    print("all real parts sorted desc:", np.round(np.sort(ev.real)[::-1][:6], 6))