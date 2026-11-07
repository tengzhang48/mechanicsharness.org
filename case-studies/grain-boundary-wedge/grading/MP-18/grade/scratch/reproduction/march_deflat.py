import numpy as np
import scipy.linalg as sla
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast


def march_deflated(N, T=0.5, nsteps=2000):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    D3 = build_D3_nonuniform(nodes)
    A = -(M @ D3)
    n = len(x)
    # spurious null vector
    w, V = np.linalg.eig(A)
    i0 = np.argmax(w.real)
    v0 = V[:, i0].real
    v0 = v0 / np.linalg.norm(v0)
    P = np.eye(n) - np.outer(v0, v0)     # projector removing v0
    Ad = P @ A @ P
    dt = T / nsteps
    Lm = np.eye(n) - 0.5 * dt * Ad
    Rp = np.eye(n) + 0.5 * dt * Ad
    lu = sla.lu_factor(Lm)
    s = np.ones(n)
    s = s - v0 * (v0 @ s)
    s = s / np.linalg.norm(s)
    logs = []
    for k in range(nsteps):
        s = sla.lu_solve(lu, Rp @ s)
        s = s - v0 * (v0 @ s)
        v = np.linalg.norm(s)
        s = s / v
        logs.append(np.log(v))
    rate = -logs[-1] / dt
    return rate, s


if __name__ == "__main__":
    for N in [80, 160, 320]:
        rate, s = march_deflated(N)
        x = np.linspace(0, 1, N + 2)[1:-1]
        vv = s / s[-1]
        print("N=%3d  deflated march rate = %.5f   f(0.5)=%.5f" %
              (N, rate, np.interp(0.5, x, vv)))