import time
import numpy as np
import solve

# resolutions (npanel tuple, ny)
RES = [
    ((240, 150, 240), 150),
    ((320, 200, 320), 200),
    ((400, 250, 400), 250),
    ((480, 290, 480), 290),
]
GUESS = {
    1e-6: (0.4471, 0.5687),
    0.001: (0.4471, 0.5687),
    0.1: (0.5155, 0.7466),
    10.0: (1.04, 1.05),
}


def res2(V, a, eps, np_, ny, L=14.0):
    f, K, T0, W, S2, X, Y, K0 = solve.residuals(V, a, eps, npanel=np_, ny=ny, L=L)
    return f - 0.1, K


def newton_safe(V, a, eps, np_, ny, L=14.0, itmax=15, tol=1e-8):
    for it in range(itmax):
        r0 = np.array(res2(V, a, eps, np_, ny, L))
        if np.linalg.norm(r0) < tol:
            return a, eps, r0, True
        da = 1e-3 * max(1.0, abs(a))
        de = 1e-3 * max(1.0, abs(eps))
        r1 = np.array(res2(V, a + da, eps, np_, ny, L))
        r2 = np.array(res2(V, a, eps + de, np_, ny, L))
        J = np.column_stack([(r1 - r0) / da, (r2 - r0) / de])
        try:
            step = np.linalg.solve(J, -r0)
        except np.linalg.LinAlgError:
            break
        nrm = np.linalg.norm(step)
        if nrm > 0.4:
            step *= 0.4 / nrm
        a += step[0]
        eps += step[1]
    r0 = np.array(res2(V, a, eps, np_, ny, L))
    return a, eps, r0, np.linalg.norm(r0) < 1e-6


if __name__ == "__main__":
    fh = open("conv.txt", "a")
    for V in [1e-6, 0.001, 0.1, 10.0]:
        g = GUESS[V]
        for (np_, ny) in RES:
            t0 = time.time()
            a, eps, rr, ok = newton_safe(V, g[0], g[1], np_, ny)
            g = (a, eps)
            line = ("V=%-8g nx=%-4d ny=%-4d a=%.7f eps=%.7f res=%.2e,%.2e ok=%d t=%.0fs"
                    % (V, np_[0], ny, a, eps, rr[0], rr[1], ok, time.time() - t0))
            print(line, flush=True)
            fh.write(line + "\n")
            fh.flush()
    fh.close()