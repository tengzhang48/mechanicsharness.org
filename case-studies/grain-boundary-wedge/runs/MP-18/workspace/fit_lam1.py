import numpy as np
from conv_phys import leading


def lam1_uniform(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    lam, ev, M, D1, D2 = leading(nodes)
    return lam


if __name__ == "__main__":
    Ns = [20, 40, 80, 160, 320, 640]
    lams = []
    for N in Ns:
        lam = lam1_uniform(N)
        lams.append(lam)
    lams = np.array(lams)
    print("N        lambda1        lam*sqrt(N)   lam*N/ln N   lam*N^0.75")
    for N, lam in zip(Ns, lams):
        print("%4d  %.8f  %.5f   %.4f   %.4f" %
              (N, lam, lam * np.sqrt(N), lam * N / np.log(N), lam * N ** 0.75))
    # fit log(lam) vs log(N)
    A = np.polyfit(np.log(Ns), np.log(lams), 1)
    print("fit: lambda1 ~ %.5f * N^(%.4f)" % (np.exp(A[1]), A[0]))