import time
import fastroots

JOBS = [
    (30.0, 10.5, 1.03, (320, 200, 320), 200, 28.0),
    (100.0, 21.0, 1.14, (320, 200, 320), 200, 40.0),
]

if __name__ == "__main__":
    fh = open("highV.txt", "a")
    for (V, e0, a0, np_, ny, L) in JOBS:
        t0 = time.time()
        a, eps, K = fastroots.root(V, e0, np_, ny, a0=a0, L=L, steps=22, tol=2e-4, step_tol=1e-5)
        R = a * a / 0.2
        line = "V=%-8g a=%.6f eps=%.6f R=%.5f K=%.3e t=%.0fs" % (V, a, eps, R, K, time.time() - t0)
        print(line, flush=True)
        fh.write(line + "\n")
        fh.flush()
    fh.close()