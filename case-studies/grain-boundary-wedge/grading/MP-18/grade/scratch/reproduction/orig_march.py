#!/usr/bin/env python3
"""Independent check of lambda1 from the ORIGINAL (dimensional) equations.

Everything set to 1: h=1, sigma0=1, E*=1, Omega=1, delta_gb D_gb=1, kT=1.
Then t0 = 4 pi kT h^3/(E* delta_gb D_gb Omega) = 4 pi.

   sigma_gb(z)=1 - (E*/2pi) Int_0^1 S(z,z') u'(z') dz'      (S = statement kernel)
   j(z) = (delta D/(kT)) sigma_gb'(z)                        (from mu=-sigma Omega)
   u_t = -(Omega/2) j'(z)
   BCs: sigma_gb(0)=0 ; sigma_gb'(1)=sigma_gb''(1)=0 ; u(0)=u'(1)=0 ; sigma_gb(.,0)=1.

March u implicitly; s rebuilt exactly each step from the kernel integral.
Decay rate in physical time = lambda1/t0, so lambda1 = decay_rate * t0 = decay_rate*4pi.
"""
import numpy as np
import stmt_spec


def march(N=200, T=None, nstep=20000):
    h = 1.0 / (N + 1)
    x = np.arange(1, N + 1) * h
    nodes = np.concatenate(([0.0], x, [1.0]))
    # M = (1/pi) Int S phi_j  ; the physical operator is (E*/2pi) Int S = 1 * pi * M
    Mpi = stmt_spec.M_stmt(x, nodes)          # = (1/pi) Int S phi
    Mfull = np.pi * Mpi                        # = Int S phi  (E*/(2pi)*... with E*=1: (1/2pi)Int S)

    # discrete second-derivative D2 on interior nodes with s(0)=0, s'(1)=s''(1)=0
    def D2_op():
        D = np.zeros((N, N))
        for i in range(N):
            D[i, i] = -2.0 / h ** 2
            if i > 0:
                D[i, i - 1] = 1.0 / h ** 2
            if i < N - 1:
                D[i, i + 1] = 1.0 / h ** 2
        # BCs: s_0 = 0 (drop), and s'(1)=s''(1)=0 -> at i=N-1 reflect
        D[N - 1, N - 2] += 1.0 / h ** 2        # s_{N+1}=s_{N-1}
        return D

    D2 = D2_op()

    def sigma(u):
        # sigma(z) = 1 - (1/2pi) Int S u'  ;  u' via hat-difference on u values
        # We solve for w=u' directly (w in hat basis, coefficients w_j)
        return 1.0 - 0.5 * (Mfull @ u)          # u here are the hat-coefficients of w=u'

    # unknowns: w_j (values of u' at interior nodes).  u = integral of w with u(0)=0.
    # s = 1 - (1/2pi) Mfull w .  j = s' , u_t = -j'/2 -> w_t = -j''/2 = -s'''/2 ... 
    # March w with  w_tau = ... ; easier: march u via u_t = -(Omega/2) j'(z), j=(dD/kT) s'
    # s' via D1, then u_t = -(1/2) D1 s' = -(1/2) D2 s (interior), with BCs.
    u = np.zeros(N)

    def D1_op():
        D = np.zeros((N, N))
        for i in range(N):
            if 0 < i < N - 1:
                D[i, i - 1] = -0.5 / h
                D[i, i + 1] = 0.5 / h
        D[0, 1] = 0.5 / h                       # one-sided at i=0 (s_0=0)
        D[N - 1, N - 2] = -0.5 / h              # s'(1)=0 -> reflect
        return D

    D1 = D1_op()

    def rhs(u):
        s = sigma(u)
        j = D1 @ s                              # j = s'
        return -0.5 * (D1 @ j)                  # u_t = -(Omega/2) j'

    if T is None:
        T = 4 * np.pi * 0.5                     # t0/2
    dt = T / nstep
    # Implicit Euler: u <- u + dt*rhs(u)  is  u <- (I - dt*G)^{-1} u  with rhs linear:
    # i.e. rhs(u) = G u , G = -0.5 D1 D1 (1 - 0.5 Mfull) ... build G explicitly.
    G = -0.5 * (D1 @ (D1 @ (1.0 * np.eye(N)) - 0.5 * (D1 @ Mfull)))
    # sanity: rhs(u) = -0.5 D1 D1 (1 - 0.5 Mfull u)? recompute cleanly:
    #   s = 1 - 0.5 Mfull u  ; j = D1 s = D1*1 - 0.5 D1 Mfull u ; u_t = -0.5 D1 j
    #   -> u_t = -0.5 D1(D1*1) + 0.25 D1 D1 Mfull u  -> linear part G = +0.25 D1 D1 Mfull
    G = 0.25 * (D1 @ (D1 @ Mfull))
    I = np.eye(N)
    A = np.linalg.inv(I - dt * G)
    amp = []
    ts = []
    for k in range(nstep):
        u = A @ u
        amp.append(np.max(np.abs(u)))
        ts.append((k + 1) * dt)
    ts = np.array(ts)
    amp = np.array(amp)
    ok = amp > 1e-14
    ts, amp = ts[ok], amp[ok]
    n0 = int(0.5 * len(ts))
    p = np.polyfit(ts[n0:], np.log(amp[n0:]), 1)
    return p[0], ts, amp


if __name__ == "__main__":
    rate, ts, amp = march(N=100)
    print("physical decay rate = %.5f  ; t0 = 4pi = %.5f ; lambda1 = rate*t0 = %.4f"
          % (rate, 4 * np.pi, rate * 4 * np.pi))