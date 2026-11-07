import numpy as np


def M_phys_fast(xcols, nodes):
    """Vectorised version of phys_oper.M_phys using the exact antiderivatives."""
    xcols = np.asarray(xcols, float)
    nodes = np.asarray(nodes, float)
    N = len(xcols)
    a = nodes[:-2]
    b = nodes[1:-1]
    c = nodes[2:]
    # left piece [a,b]: alpha_l=-a/(b-a), beta_l=1/(b-a); right [b,c]: alpha_r=c/(c-b), beta_r=-1/(c-b)
    alpha_l = -a / (b - a)
    beta_l = 1.0 / (b - a)
    alpha_r = c / (c - b)
    beta_r = -1.0 / (c - b)

    def G(x, y, alpha, beta):
        Lm = np.log(np.abs(y - x) + 1e-300)
        Lp = np.log(np.abs(y + x) + 1e-300)
        F1 = -2.0 * x * (3.0 * y + 2.0 * x) / (y + x) ** 2 - (Lp - Lm)
        F2 = (-2.0 * x ** 2 + (y + x) * (y * Lm - x * Lm - (y + 7.0 * x) * Lp)) / (y + x)
        return alpha * F1 + beta * (y * F1 - F2)

    X = xcols[:, None]     # (N,1) collocation
    Al = alpha_l[None, :]
    Bl = beta_l[None, :]
    Ar = alpha_r[None, :]
    Br = beta_r[None, :]
    A_ = a[None, :]
    B_ = b[None, :]
    C_ = c[None, :]
    val = (G(X, B_, Al, Bl) - G(X, A_, Al, Bl)
           + G(X, C_, Ar, Br) - G(X, B_, Ar, Br))
    return val / np.pi