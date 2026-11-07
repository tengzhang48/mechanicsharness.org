import numpy as np
from conv_D3 import leading_D3


def track_genuine_D3():
    """For -M D3, separate genuine modes: track eigenvalues that are stable
    (grid-independent) across N."""
    print("=== -M D3 (physical kernel): eigenvalues sorted by |value| ===")
    for N in [20, 40, 80, 160]:
        x, ev, A, M, D3 = leading_D3(N)
        vals = np.sort(-ev.real)      # decay rates
        small = vals[np.argsort(np.abs(vals))][:6]
        print("N=%3d  smallest |decay rates|: %s" % (N, np.round(small, 5)))


if __name__ == "__main__":
    track_genuine_D3()