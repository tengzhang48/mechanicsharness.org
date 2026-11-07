import numpy as np
import solve

_cache = {}


def get_state(V, eps, npanel, ny, L=14.0):
    key = (V, round(eps, 12), npanel, ny, L)
    st = _cache.get(key)
    if st is None:
        st = solve.State(V, eps, L=L, npanel=npanel, ny=ny)
        _cache[key] = st
        # keep cache tiny: memory-heavy factorizations
        while len(_cache) > 3:
            _cache.pop(next(iter(_cache)))
    return st


def a_for_eps(V, eps, npanel, ny, a0, L=14.0, it=8, tol=1e-9):
    st = get_state(V, eps, npanel, ny, L)
    a = a0
    f = st.residual(a)[0]
    for _ in range(it):
        if abs(f - 0.1) < tol:
            break
        da = 1e-3 * max(1.0, abs(a))
        f1 = st.residual(a + da)[0]
        slope = (f1 - f) / da
        if abs(slope) < 1e-12:
            break
        a = a - (f - 0.1) / slope
        f = st.residual(a)[0]
    return a


def K_of_eps(V, eps, npanel, ny, a0, L=14.0):
    a = a_for_eps(V, eps, npanel, ny, a0, L)
    return get_state(V, eps, npanel, ny, L).residual(a)[1], a


def root(V, eps0, npanel, ny, a0=0.45, L=14.0, steps=12, tol=2e-4, step_tol=1e-5, verbose=False):
    eps = eps0
    a = a_for_eps(V, eps, npanel, ny, a0, L)
    K, a = K_of_eps(V, eps, npanel, ny, a, L)
    for it in range(steps):
        if abs(K) < tol:
            break
        de = 1e-3 * max(1.0, abs(eps))
        K2, _ = K_of_eps(V, eps + de, npanel, ny, a, L)
        slope = (K2 - K) / de
        if abs(slope) < 1e-12:
            break
        step = -K / slope
        step = max(-0.2, min(0.2, step))
        if abs(step) < step_tol:
            break
        eps = eps + step
        a = a_for_eps(V, eps, npanel, ny, a, L)
        K, a = K_of_eps(V, eps, npanel, ny, a, L)
        if verbose:
            print('   it%d eps=%.8f a=%.8f K=%.3e step=%.2e' % (it, eps, a, K, step))
    return a, eps, K