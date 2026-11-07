#!/usr/bin/env python3
"""Independent check: time-march  d_tau s = M D3 s  (statement kernel) from s=1 and
read the late-time decay rate; compare with the leading eigenvalue of M D3."""
import numpy as np
from scipy.linalg import expm
import stmt_spec
from d3_uniform import build_D3_uniform


def march(N=200, tmax=0.5, nstep=4000):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = stmt_spec.M_stmt(x, nodes)
    D3 = build_D3_uniform(N)
    A = M @ D3
    s = np.ones(N)
    # project out the null space of D3 (the marginal discrete mode)
    U, sv, Vt = np.linalg.svd(D3)
    null = Vt[-1]                      # right null vector of D3
    s = s - (s @ null) * null
    dt = tmax / nstep
    # implicit Euler (matrix exponential per step is fine for moderate N)
    E = expm(A * dt)
    amp = []
    ts = []
    for k in range(nstep):
        s = E @ s
        ts.append((k + 1) * dt)
        amp.append(np.max(np.abs(s)))
    ts = np.array(ts)
    amp = np.array(amp)
    ok = amp > 1e-12
    ts, amp = ts[ok], amp[ok]
    # fit late-time slope of log|amp|
    n0 = int(0.6 * len(ts))
    p = np.polyfit(ts[n0:], np.log(amp[n0:]), 1)
    return p[0], ts, amp


if __name__ == "__main__":
    for N in [100, 200]:
        rate, ts, amp = march(N)
        print(f"N={N}  measured late-time decay rate (should be -lambda1) = {rate:.4f}")