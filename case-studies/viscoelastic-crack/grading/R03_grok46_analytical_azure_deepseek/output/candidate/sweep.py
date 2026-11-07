import time
import numpy as np
import fastroots

NP = (200, 120, 200)
NY = 120


def run(V, e0, a0, L):
    a, eps, K = fastroots.root(V, e0, NP, NY, a0=a0, L=L, steps=18, tol=1e-4, step_tol=1e-5)
    return a, eps, K


if __name__ == "__main__":
    fh = open("sweep.txt", "a")
    # chain through speeds
    pts = [(0.001, 0.576, 0.4496), (0.01, 0.60, 0.46), (0.05, 0.68, 0.49),
           (0.1, 0.7683, 0.5157), (0.5, 1.25, 0.63), (1.0, 1.596, 0.6998),
           (3.0, 3.0, 0.85), (10.0, 5.8, 0.955), (30.0, 10.0, 1.05),
           (100.0, 20.0, 1.15), (300.0, 40.0, 1.2)]
    for (V, e0, a0) in pts:
        L = max(14.0, e0 + 10.0)
        t0 = time.time()
        try:
            a, eps, K = run(V, e0, a0, L)
        except Exception as ex:
            a, eps, K = float('nan'), float('nan'), float('nan')
            print('ERR', V, repr(ex), flush=True)
        R = a * a / 0.2
        line = "V=%-9g L=%-5.1f a=%.6f eps=%.6f R=%.5f K=%.2e t=%.0fs" % (
            V, L, a, eps, R, K, time.time() - t0)
        print(line, flush=True)
        fh.write(line + "\n")
        fh.flush()
    fh.close()