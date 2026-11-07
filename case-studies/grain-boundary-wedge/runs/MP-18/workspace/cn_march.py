import numpy as np
import scipy.linalg as sla
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast


def march(N, T=1.0, nsteps=4000):
    """Explicit-ish stable march of  s_tau = -M D3 s  via Crank-Nicolson on the
    deviation.  We march  ds/dtau = A s  with A = -(M D3) using CN:
      (I - dt/2 A) s_{k+1} = (I + dt/2 A) s_k .
    Eigen-decay of ||s|| gives the dominant (slowest) mode rate."""
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    D3 = build_D3_nonuniform(nodes)
    A = -(M @ D3)
    n = len(x)
    dt = T / nsteps
    Lm = np.eye(n) - 0.5 * dt * A
    Rp = np.eye(n) + 0.5 * dt * A
    lu = sla.lu_factor(Lm)
    s = np.ones(n)
    s = s / np.linalg.norm(s)
    ts = np.zeros(nsteps + 1)
    ns = np.zeros(nsteps + 1)
    ns[0] = 1.0
    for k in range(nsteps):
        s = sla.lu_solve(lu, Rp @ s)
        v = np.linalg.norm(s)
        s = s / v
        ts[k + 1] = (k + 1) * dt
        ns[k + 1] = np.log(v)     # log growth per step
    return ts, ns


if __name__ == "__main__":
    print("CN march: constant log-decay/step should -> -dt*lambda_slowest")
    for N in [80, 160]:
        ts, ns = march(N, T=0.5, nsteps=2000)
        rate_est = -ns[-1] / (ts[1])   # per unit time
        print("N=%3d  late log-decay/step=%.6e  -> lambda ~ %.5f" %
              (N, ns[-1], rate_est))