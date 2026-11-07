import numpy as np
from solve_D3n import solve_D3


if __name__ == "__main__":
    print("N      slow(~60)   ~1070      ~2378     #modes")
    for N in [40, 80, 160, 320]:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        r = solve_D3(nodes)
        # genuine: pick rates nearest to reference values
        def nearest(v):
            return r[np.argmin(np.abs(r - v))]
        print("%4d   %.5f    %.5f    %.5f   %d" %
              (N, nearest(60), nearest(1075), nearest(2378), len(r)))