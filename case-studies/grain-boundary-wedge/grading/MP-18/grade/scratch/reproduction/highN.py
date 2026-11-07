import numpy as np
from converge_all import rates_chain, rates_D3


def get(seq, ref):
    return seq[np.argmin(np.abs(seq - ref))]


if __name__ == "__main__":
    Ns = [640, 960, 1280, 1600]
    chain1 = []
    for N in Ns:
        rc = rates_chain(N)
        chain1.append(rc[1])          # genuine 2nd-slowest
    chain1 = np.array(chain1)
    print("high-N chain mode2:", np.round(chain1, 6))
    # differences
    print("diffs:", np.round(np.diff(chain1), 6))
    # Richardson with p=1 on last 2 and last 3
    def rich(Ns, vals, p):
        A = np.vstack([np.ones(len(Ns)), np.array(Ns, float) ** (-p)]).T
        s, *_ = np.linalg.lstsq(A, vals, rcond=None)
        return s[0]
    for p in [1.0, 1.5, 2.0, 2.5]:
        print("p=%.1f  extrap(last3)=%.6f  extrap(last2)=%.6f" %
              (p, rich(Ns[-3:], chain1[-3:], p), rich(Ns[-2:], chain1[-2:], p)))