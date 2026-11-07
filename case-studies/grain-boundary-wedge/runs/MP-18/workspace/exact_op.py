import numpy as np
from phys_oper import F1, F2

# ============ EXACT application of the integral operator to a trial s =========
# (Ms)(x) = (1/pi) Int_0^1 S_phys(x,y) s'''(y) dy.
# For a trial s given at nodes, we want Ms at the SAME nodes, using the exact
# y-antiderivative.  Write s''' = g.  Then Int_0^1 S(x,y) g(y) dy.  g is only
# known at nodes; represent g by its piecewise-linear hat interpolant g(y) =
# sum_j g_j phi_j(y). Then Int S g dy = sum_j g_j * (Int S phi_j dy) = M_phys @ g.
# So:  (Ms)(x_i) = [M_phys @ (s''' at nodes)]_i .   EXACT if s''' sampled exactly.

def s_double_prime_hat(xcols, nodes):
    """Return D2 (values of s'') at interior nodes; we use the nonuniform D2."""
    from fd_nonuniform import build_D1D2
    _, D2 = build_D1D2(nodes)
    return D2


if __name__ == "__main__":
    # Build the operator matrix  A = -M_phys @ D3nonuniform  and solve eigenvalues
    from fd_nonuniform import build_D1D2
    from phys_oper import M_phys
    from d3_nonuniform import build_D3_nonuniform
    import numpy as np

    N = 30
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D3 = build_D3_nonuniform(nodes)
    f = x ** 3
    print("D3 x^3 (interior should be 6):", np.round(D3 @ f, 3))