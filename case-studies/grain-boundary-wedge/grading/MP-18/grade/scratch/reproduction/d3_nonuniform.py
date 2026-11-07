import numpy as np
from fd_nonuniform import build_D1D2


def build_D3_nonuniform(nodes, ghost_d=1.0):
    """Nonuniform D3; ghost nodes placed symmetrically about x=1 (distance = last
    interior spacing) so that s'(1)=s''(1)=0 -> s_{M+1}=s_{M-1}, s_{M+2}=s_{M-2}."""
    x = np.asarray(nodes, float)
    M = len(x) - 2
    hM = x[M + 1] - x[M]        # spacing just before x=1
    # extended coordinate array: indices 0..M+2 (M+1 is x=1, M+2 ghost)
    xe = np.concatenate((x, [x[M + 1] + hM]))     # length M+3
    D3 = np.zeros((M, M))

    def coeffs(xs, k):
        import math
        n = len(xs)
        A = np.zeros((n, n))
        for p in range(n):
            for q in range(n):
                A[p, q] = (xs[q] - xs[k]) ** p
        b = np.zeros(n)
        b[3] = math.factorial(3)
        return np.linalg.solve(A, b)

    xi = np.arange(M) + 1
    for i in range(M):
        fi = xi[i]
        # 5-node stencil
        nodeidx = [fi - 2, fi - 1, fi, fi + 1, fi + 2]
        if min(nodeidx) < 0:
            nodeidx = list(range(0, 5))
        if max(nodeidx) > M + 2:
            nodeidx = list(range(M - 3, M + 2))
        nodeidx = [k for k in nodeidx if 0 <= k <= M + 2]
        if fi not in nodeidx:
            nodeidx = sorted(set(nodeidx) | {fi})
        xs = np.array([xe[k] for k in nodeidx])
        kk = nodeidx.index(fi)
        c = coeffs(xs, kk)
        for cc, node in zip(c, nodeidx):
            if node == 0:
                continue
            if node == M + 2:
                D3[i, M - 2] += cc     # ghost -> s_{M-2} (index M-2 means node M+... wait)
            elif node == M + 1:
                D3[i, M - 1] += cc     # node M+1 is x=1; treat s(1) via BC? 
            else:
                D3[i, node - 1] += cc
    return D3


if __name__ == "__main__":
    N = 20
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D3 = build_D3_nonuniform(nodes)
    f = x ** 3
    print("D3 x^3 (should be 6 interior):", np.round(D3 @ f, 4))