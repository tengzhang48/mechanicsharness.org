import time
import fastroots

RES = [
    ((240, 150, 240), 150),
    ((320, 200, 320), 200),
    ((400, 250, 400), 250),
]

if __name__ == "__main__":
    fh = open("v10.txt", "a")
    e0, a0 = 6.0, 0.955
    for (np_, ny) in RES:
        t0 = time.time()
        a, eps, K = fastroots.root(10.0, e0, np_, ny, a0=a0, steps=18, tol=1e-4, step_tol=1e-5)
        e0, a0 = eps, a
        line = "V=10 nx=%-4d ny=%-4d a=%.7f eps=%.7f K=%.2e t=%.0fs" % (
            np_[0], ny, a, eps, K, time.time() - t0)
        print(line, flush=True)
        fh.write(line + "\n")
        fh.flush()
    fh.close()