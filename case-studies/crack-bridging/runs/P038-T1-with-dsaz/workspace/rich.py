"""Nyström collocation solver for the finite bridged center crack.

Crack |x|<c; closure sigma0 on a<|x|<c (both tips), traction-free on |x|<a.
The unknown is the (climb) dislocation density b(x) on the crack; the
traction on the crack line is

    sigma_yy(x) = -(E'/2pi) f.p. int_{-c}^{c} b(s)/(x-s) ds + remote

(the Cauchy principal value). For a symmetric crack, using b(x) odd:

    sigma_yy(x) = -(E'/(2pi)) f.p. int_0^c b(s)[1/(x-s) - 1/(x+s)] ds + remote

with b(s) the (even) density on the positive half. Closure: for x in the
zone, sigma_yy = sigma0; for x in the slit, sigma_yy = 0. Remote = Sigma
everywhere -> move to rhs.

Unknown b(s), s in (0,c). Tip behavior: b ~ sqrt(c-s) at the outer tip
(s=c) and b ~ (s-a)^{-1/2}-ish at the junction s=a?  Actually the junction
has the log-type behavior; we represent b as

    b(s) = A w_tip(s) + (regular part interpolated)

and use a Gauss-Chebyshev Nyström rule in the variable s = c*cos-ish...

This is more robust if we use the semi-infinite zone kernel (zone.py) for
the interior and match.  For time, we instead use a high-order approach:
substitute s = c * (1 - (1 - a/c) * sin^2(theta)^...) -- too fiddly.

Instead: use the Krenk-Gauss quadrature with the sqrt weight and REMOVE
the log/jump singularity at s=a by adding the analytic 'closure only on a
half-plane strip' base.  Too long for budget.

Pragmatic decision: keep the spectral solver but accelerate convergence by
a Richardson fit over N in {300,600,1000,1500,2200}.
"""
import numpy as np
from final import state, ell, an, K0n, Ep, sigma0, delta_c

if __name__ == '__main__':
    Ns = [300, 600, 1000, 1500, 2200]
    Ls, Ss = [], []
    for N in Ns:
        lo, hi = 0.35, 0.75
        for _ in range(20):
            mid = 0.5 * (lo + hi)
            _, dt = state(mid, N)
            if dt < 1.0:
                lo = mid
            else:
                hi = mid
        Ln = 0.5 * (lo + hi)
        Sig, dt = state(Ln, N)
        Ls.append(Ln * ell * 1e3)
        Ss.append(Sig * 160.0)
        print(f"N={N:5d}: L={Ls[-1]:.6f} mm  Sigma={Ss[-1]:.5f} MPa")
    Ls = np.array(Ls)
    Ss = np.array(Ss)
    # try Richardson: assume L(N)=Linf + A/N^p ; fit p via three N
    Ns = np.array(Ns, dtype=float)
    print("finite-N values:", list(zip(Ns, Ls, Ss)))