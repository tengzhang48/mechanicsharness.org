import numpy as np
from converge_all import rates_chain, rates_D3


def get(seq, ref):
    return seq[np.argmin(np.abs(seq - ref))]


if __name__ == "__main__":
    Ns = np.array([40, 80, 160, 320, 640, 960], float)
    lam1_chain = np.array([rates_chain(int(N))[1] for N in Ns])   # 2nd slowest = genuine leading
    lam1_d3 = np.array([get(rates_D3(int(N)), 60) for N in Ns])
    print("N      chain-mode2    D3-mode(~60)")
    for N, a, b in zip(Ns, lam1_chain, lam1_d3):
        print("%4d   %.6f      %.6f" % (N, a, b))

    def rich(Ns, vals, p):
        A = np.vstack([np.ones(len(Ns)), Ns ** (-p)]).T
        sol, *_ = np.linalg.lstsq(A, vals, rcond=None)
        return sol[0]
    for p in [1.0, 1.5, 2.0]:
        print("p=%.1f  chain extrap=%.6f   D3 extrap=%.6f" %
              (p, rich(Ns[-3:], lam1_chain[-3:], p), rich(Ns[-3:], lam1_d3[-3:], p)))