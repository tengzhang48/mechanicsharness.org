import numpy as np

# Exact antiderivatives of the statement kernel S(x,y) w.r.t. y (sympy-verified):
#   F1(x,y) = Int S dy
#   F2(x,y) = Int F1 dy
#   Int (alpha + beta y) S dy = alpha*F1 + beta*(y*F1 - F2)

def F1(x, y):
    Lm = np.log(np.abs(y - x) + 1e-300)
    Lp = np.log(np.abs(y + x) + 1e-300)
    return 2.0 * x * (3.0 * y + 2.0 * x) / (y + x) ** 2 + (Lp - Lm)


def F2(x, y):
    Lm = np.log(np.abs(y - x) + 1e-300)
    Lp = np.log(np.abs(y + x) + 1e-300)
    return (2.0 * x ** 2 + (y + x) * (-y * Lm + x * Lm + (y + 7.0 * x) * Lp)) / (y + x)


def int_linear_S(x, lo, hi, alpha, beta):
    """Exact Int_lo^hi (alpha + beta*y) S(x,y) dy."""
    def G(y):
        return alpha * F1(x, y) + beta * (y * F1(x, y) - F2(x, y))
    return G(hi) - G(lo)


def Smat_statement(xcols, nodes):
    """(1/pi) PV Int_0^1 S(x_i,y) phi_j(y) dy for hat basis on arbitrary nodes."""
    N = len(xcols)
    M = np.zeros((N, N))
    for i in range(N):
        z = xcols[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            # left piece [a,b]: phi = (y-a)/(b-a) = (-a/(b-a)) + (1/(b-a)) y
            val = int_linear_S(z, a, b, -a / (b - a), 1.0 / (b - a))
            # right piece [b,c]: phi = (c-y)/(c-b) = (c/(c-b)) + (-1/(c-b)) y
            val += int_linear_S(z, b, c, c / (c - b), -1.0 / (c - b))
            M[i, j] = val / np.pi
    return M


if __name__ == "__main__":
    import oper
    N = 10
    h = 1.0 / (N + 1)
    xcols = np.arange(1, N + 1) * h
    nodes = np.concatenate(([0.0], xcols, [1.0]))
    Ms = Smat_statement(xcols, nodes)
    Mn = oper.M_hat_numeric(xcols, nodes, ngauss=800)
    print("exact-antideriv vs numeric-PV-quad: max abs diff =", np.max(np.abs(Ms - Mn)),
          " rel =", np.max(np.abs(Ms - Mn)) / np.max(np.abs(Ms)))