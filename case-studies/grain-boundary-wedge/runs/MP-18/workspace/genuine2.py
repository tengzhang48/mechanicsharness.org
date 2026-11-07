import numpy as np
from conv_D3 import leading_D3


def all_rates(N):
    x, ev, A, M, D3 = leading_D3(N)
    return np.sort(-ev.real)


if __name__ == "__main__":
    # Collect mode sequences by matching across N (nearest-neighbour tracking)
    Ns = [30, 60, 120, 240]
    seqs = []
    prev = None
    for N in Ns:
        r = all_rates(N)
        r = r[r > 1e-6]     # drop the zero
        print("N=%3d number of positive rates: %d" % (N, len(r)))
        print("   smallest 8:", np.round(np.sort(r)[:8], 4))
        seqs.append(r)