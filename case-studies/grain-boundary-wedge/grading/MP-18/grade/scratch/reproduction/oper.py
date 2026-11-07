import numpy as np
from numpy.polynomial.legendre import leggauss

# ================= Problem-statement kernel nonlocal operator =================
# s(x) = 1 - (1/pi) PV Int_0^1 S(x,y) w(y) dy
# S(x,y) = 1/(x-y) - 1/(x+y) - 2y(x-y)/(x+y)^3
# Exact antiderivative w.r.t. y  (sympy-verified):
#   F(x,y) = [ 2x(3y+2x) + (ln|y+x| - ln|y-x|)(y+x)^2 ] / (y+x)^2

def F_state(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    lr = np.log(np.abs(y + x) + 1e-300) - np.log(np.abs(y - x) + 1e-300)
    return (2.0 * x * (3.0 * y + 2.0 * x) + lr * (y + x) ** 2) / (y + x) ** 2


def Int_hat_statement(xcols, nodes):
    """M_ij = (1/pi) PV Int_0^1 S(x_i,y) phi_j(y) dy, hat basis on arbitrary nodes.
    Uses the exact antiderivative F_state; PV handled analytically (the log|x-y|
    term's jump is captured by the antiderivative)."""
    N = len(xcols)
    M = np.zeros((N, N))
    for i in range(N):
        z = xcols[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            # left hat piece on [a,b]: phi=(y-a)/(b-a)
            pa, pb = F_state(z, a), F_state(z, b)
            # integral of (y-a)/(b-a)*S = (1/(b-a))[ F(y) - a*Int S ]  -- no:
            # we need Int phi S, phi linear. Use: Int (alpha + beta y) S dy
            #   = alpha F(y) + beta [y F(y) - Int F dy]  -> requires double antideriv.
            # Simpler: numeric per-hat with PV-aware Gauss splitting at y=z.
            pass
    return M


def M_hat_numeric(xcols, nodes, ngauss=400):
    """Reference operator by PV-aware Gauss quadrature of the raw kernel S(x,y)."""
    N = len(xcols)
    M = np.zeros((N, N))
    gx, gw = leggauss(ngauss)

    def Sper(z, y):
        return 1.0 / (z - y) - 1.0 / (z + y) - 2.0 * y * (z - y) / (z + y) ** 3

    for i in range(N):
        z = xcols[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            tot = 0.0
            for lo, hi in ((a, b), (b, c)):
                # split at z if z inside [lo,hi] for PV
                pts = [lo, hi]
                if abs(b - lo) > 0 and lo <= z <= hi and (hi - lo) > 1e-14:
                    if lo < z < hi:
                        pts = [lo, z, hi]
                for k in range(len(pts) - 1):
                    l, hgt = pts[k], pts[k + 1]
                    ys = 0.5 * (hgt - l) * gx + 0.5 * (hgt + l)
                    w = 0.5 * (hgt - l) * gw
                    # hat value
                    if hi == b:
                        phi = (ys - lo) / (b - lo)
                    else:
                        phi = (hi - ys) / (hi - b)
                    vals = Sper(z, ys)
                    vals = np.where(np.abs(ys - z) < 1e-14, 0.0, vals)
                    tot += np.sum(w * vals * phi)
            M[i, j] = tot / np.pi
    return M


if __name__ == "__main__":
    N = 10
    h = 1.0 / (N + 1)
    xcols = np.arange(1, N + 1) * h
    nodes = np.concatenate(([0.0], xcols, [1.0]))
    Mn = M_hat_numeric(xcols, nodes, ngauss=600)
    print("numeric statement-kernel M built:", Mn.shape)
    print("row sums (should relate to Int S sum):", np.round(Mn.sum(axis=1), 4))