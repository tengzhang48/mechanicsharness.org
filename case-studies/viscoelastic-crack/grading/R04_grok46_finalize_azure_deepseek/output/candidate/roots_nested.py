import numpy as np
import solve


def fK(V, a, eps, np_, ny, L=14.0):
    f, K, T0, W, S2, X, Y, K0 = solve.residuals(V, a, eps, npanel=np_, ny=ny, L=L)
    return f - 0.1, K


def a_for_eps(V, eps, np_, ny, a0, L=14.0, it=10):
    a = a0
    for _ in range(it):
        f, _ = fK(V, a, eps, np_, ny, L)
        if abs(f) < 1e-9:
            break
        da = 1e-3 * max(1.0, abs(a))
        f1, _ = fK(V, a + da, eps, np_, ny, L)
        slope = (f1 - f) / da
        if abs(slope) < 1e-12:
            break
        a -= f / slope
    return a


def K_of_eps(V, eps, np_, ny, a0, L=14.0):
    a = a_for_eps(V, eps, np_, ny, a0, L)
    _, K = fK(V, a, eps, np_, ny, L)
    return K, a


def root(V, eps0, np_, ny, L=14.0, a0=0.45, steps=30, tol=1e-6):
    # secant on K(eps) with noise averaging (use central difference with sizable step)
    eps = eps0
    a0 = a_for_eps(V, eps, np_, ny, a0, L)
    K, a = K_of_eps(V, eps, np_, ny, a0)
    for _ in range(steps):
        if abs(K) < tol:
            break
        de = 2e-3 * max(1.0, abs(eps))
        K2, _ = K_of_eps(V, eps + de, np_, ny, a, L)
        slope = (K2 - K) / de
        if abs(slope) < 1e-12:
            break
        step = -K / slope
        step = max(-0.15, min(0.15, step))
        eps += step
        a = a_for_eps(V, eps, np_, ny, a, L)
        K, a = K_of_eps(V, eps, np_, ny, a)
    return a, eps, K


if __name__ == "__main__":
    for V, e0 in [(0.1, 0.7466), (1e-6, 0.5687)]:
        for np_, ny in [((240, 150, 240), 150), ((320, 200, 320), 200)]:
            a, eps, K = root(V, e0, np_, ny)
            print('V=%-8g res=%s a=%.7f eps=%.7f K=%.2e' % (V, np_, a, eps, K), flush=True)