import numpy as np
from chain import lam1_chain
import grids


def test_grading():
    print("=== geometric grading near x=0, physical kernel, chain operator ===")
    for p in [1.0, 1.5, 2.0]:
        print("-- x=eta^%s --" % p)
        for N in [40, 80, 160]:
            nodes = grids.graded_nodes(N, p)
            lam, ev = lam1_chain(nodes)
            # count spurious huge modes
            top = np.sort(-ev.real)[::-1]
            print("  N=%3d lam1=%.6f  biggest=%.3e  (2nd slow=%.4f)" %
                  (N, lam, top[0], top[1]))


if __name__ == "__main__":
    test_grading()