import time
import fastroots

RES = [
    ((240, 150, 240), 150),
    ((320, 200, 320), 200),
    ((400, 250, 400), 250),
]
GUESS = {
    1e-6: (0.5705, 0.4479),
    0.001: (0.5705, 0.4479),
    0.1: (0.7683, 0.5158),
    10.0: (1.05, 1.04),
}

if __name__ == "__main__":
    fh = open("final.txt", "a")
    for V in [1e-6, 0.001, 0.1, 10.0]:
        e0, a0 = GUESS[V]
        for (np_, ny) in RES:
            t0 = time.time()
            try:
                a, eps, K = fastroots.root(V, e0, np_, ny, a0=a0, steps=22, tol=1e-7)
            except Exception as e:
                a, eps, K = float('nan'), float('nan'), float('nan')
                print('ERR', V, np_, repr(e), flush=True)
            e0, a0 = eps, a
            line = "V=%-8g nx=%-4d ny=%-4d a=%.8f eps=%.8f K=%.2e t=%.0fs" % (
                V, np_[0], ny, a, eps, K, time.time() - t0)
            print(line, flush=True)
            fh.write(line + "\n")
            fh.flush()
    fh.close()