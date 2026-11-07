import time
import numpy as np
import roots_nested as rn

RES = [
    ((240, 150, 240), 150),
    ((320, 200, 320), 200),
    ((400, 250, 400), 250),
    ((480, 290, 480), 290),
]
GUESS = {
    1e-6: (0.4471, 0.5687),
    0.001: (0.4471, 0.5687),
    0.1: (0.5155, 0.7675),
    10.0: (1.04, 1.05),
}

if __name__ == "__main__":
    fh = open("nested.txt", "a")
    for V in [1e-6, 0.001, 0.1, 10.0]:
        g = GUESS[V]
        for (np_, ny) in RES:
            t0 = time.time()
            try:
                a, eps, K = rn.root(V, g[1], np_, ny, a0=g[0], steps=14, tol=1e-7)
            except Exception as e:
                a, eps, K = float('nan'), float('nan'), float('nan')
                print('ERR', V, np_, e)
            g = (a, eps)
            line = "V=%-8g nx=%-4d ny=%-4d a=%.7f eps=%.7f K=%.2e t=%.0fs" % (
                V, np_[0], ny, a, eps, K, time.time() - t0)
            print(line, flush=True)
            fh.write(line + "\n")
            fh.flush()
    fh.close()