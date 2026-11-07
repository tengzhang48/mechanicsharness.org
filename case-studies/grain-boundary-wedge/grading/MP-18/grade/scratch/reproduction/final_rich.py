import numpy as np
from graded_D3 import rates_D3


def mode_near(N, ref):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    r = rates_D3(nodes)
    return r[np.argmin(np.abs(r - ref))]


if __name__ == "__main__":
    Ns = np.array([40, 80, 160, 320, 480])
    m60 = np.array([mode_near(int(N), 60) for N in Ns])
    m1070 = np.array([mode_near(int(N), 1070) for N in Ns])
    m2370 = np.array([mode_near(int(N), 2375) for N in Ns])
    print("N        ~60        ~1070      ~2375")
    for N, a, b, c in zip(Ns, m60, m1070, m2370):
        print("%5d  %.6f   %.5f   %.5f" % (N, a, b, c))

    # Richardson lam ~ L + c/N + d/N^2 on finest points
    def rich(Ns, vals, kmax=2):
        # solve for L minimizing; use p=1,2 basis
        A = np.vstack([np.ones(len(Ns)), 1.0 / Ns, 1.0 / Ns ** 2]).T
        sol, *_ = np.linalg.lstsq(A, vals, rcond=None)
        return sol

    print("\n~60  L,c,d =", np.round(rich(Ns[-4:], m60[-4:]), 6))
    print("~1070 L,c,d =", np.round(rich(Ns[-4:], m1070[-4:]), 6))
    print("~2375 L,c,d =", np.round(rich(Ns[-4:], m2370[-4:]), 6))