import time
import fastroots
import solve

if __name__ == "__main__":
    fh = open("v10L.txt", "a")
    for L in [40.0, 28.0]:
        for np_, ny in [((400, 250, 400), 250)]:
            t0 = time.time()
            a, eps, K = fastroots.root(10.0, 5.717, np_, ny, a0=0.9612, L=L,
                                       steps=20, tol=2e-4, step_tol=1e-5)
            r = solve.bulk_ratio(10.0, a, eps, np_, ny, L)
            line = "L=%g res=%d a=%.6f eps=%.6f K=%.2e R=%.5f ratio=%.5f t=%.0fs" % (
                L, np_[0], a, eps, K, a * a / 0.2, r, time.time() - t0)
            print(line, flush=True)
            fh.write(line + "\n")
            fh.flush()
    fh.close()