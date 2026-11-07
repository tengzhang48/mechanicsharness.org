import numpy as np
from scipy.optimize import curve_fit
from conv_phys import leading


def lam1_uniform(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    lam, ev, M, D1, D2 = leading(nodes)
    return lam


if __name__ == "__main__":
    Ns = np.array([30, 60, 120, 240, 400], float)
    lams = np.array([lam1_uniform(int(N)) for N in Ns])
    print("data:", list(zip(Ns.astype(int), np.round(lams, 6))))

    # model A: lam = c * N^(-alpha)   (lambda*=0)
    pA, _ = curve_fit(lambda N, c, al: c * N ** (-al), Ns, lams, p0=[3, 0.5])
    print("model c*N^-a   : c=%.5f alpha=%.5f" % tuple(pA))
    residA = lams - pA[0] * Ns ** (-pA[1])
    print("   residual:", np.round(residA, 6))

    # model B: lam = lam* + c*N^(-alpha)
    pB, _ = curve_fit(lambda N, lstar, c, al: lstar + c * N ** (-al), Ns, lams,
                      p0=[0.0, 3, 0.5], maxfev=20000)
    print("model l*+c*N^-a: lam*=%.5f c=%.5f alpha=%.5f" % tuple(pB))
    residB = lams - (pB[0] + pB[1] * Ns ** (-pB[2]))
    print("   residual:", np.round(residB, 6))
    print("   rms A=%.3e  rms B=%.3e" % (np.sqrt(np.mean(residA**2)),
                                         np.sqrt(np.mean(residB**2))))