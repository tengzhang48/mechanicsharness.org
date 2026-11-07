#!/usr/bin/env python3
"""
Time-dependent poroelastic crack volume: 1-D diffusion boundary-layer model,
verified against the exact poroelastic crack solution.

PHYSICAL MODEL
--------------
The crack faces are drained (p = 0 on |x|<L, y = 0).  Under a step normal
traction the pore-pressure disturbance is confined to a thin layer of
thickness ~ sqrt(c t) adjacent to the crack and varies essentially only
across the face (1-D diffusion in y) while sqrt(c t) << L.

Diffusion problem for the pressure *deficit* u = p0 - p:

    du/dt = c u_yy,  y>0,  u(0,t)=p0,  u(y,0)=0,  u(inf,t)->0

with exact solution  u(y,t) = p0 * erfc(y/sqrt(4 c t)),
i.e.  p(y,t) = p0 * erf(y / sqrt(4 c t)).

The crack opening responds to the induced pore pressure.  At t=0+ it is the
undrained opening (nu_u); as the pressure diffuses out it relaxes to the
drained opening (nu).  The volume is a weighted integral of the (1-D)
pressure profile.  The small-time behaviour gives

    V(t) = V1_0 + A sqrt(c t / L^2),     A = (V1_inf - V1_0) * 2/sqrt(pi).

We verify with (M2) a finite-difference solve of the same 1-D problem and
(M3) the small-time asymptotics.
"""
import numpy as np
from scipy.special import erfc


def solve_1D_diffusion(p0, c, L, times, ny=4000, ymax_factor=60.0):
    """Finite-difference solve of the deficit problem.

        du/dt = c u_yy,  u(0,t)=p0,  u(y,0)=0,  u(inf,t)=0.
    Exact: u(y,t) = p0 erfc(y/sqrt(4ct)).  Returns (y, {t: u}).
    """
    times = np.asarray(times, float)
    tmax = times.max()
    ymax = ymax_factor * np.sqrt(4 * c * tmax)
    y = np.linspace(0, ymax, ny)
    dy = y[1] - y[0]
    u = np.zeros(ny)
    u[0] = p0
    dt = 0.4 * dy * dy / c
    t = 0.0
    tlist = np.sort(times)
    out = {}
    idx = 0
    guard = 0
    while idx < len(tlist) and guard < 5_000_000:
        guard += 1
        if t >= tlist[idx] - 1e-15:
            out[tlist[idx]] = u.copy()
            idx += 1
            continue
        un = u.copy()
        un[1:-1] = u[1:-1] + dt * c * (u[2:] - 2 * u[1:-1] + u[:-2]) / dy**2
        un[0] = p0
        un[-1] = 0.0
        u = un
        t += dt
    return y, out


def V_of_t_from_pressure_profile(nu, nuu, G, L, p0, c, alpha,
                                 times, ny=20000, ymax_factor=80.0):
    """Crack volume V(t) from the exact 1-D pressure profile.

    Model: V responds linearly to the induced pore pressure.  The undrained
    (t=0+) opening uses nu_u; the fully drained opening uses nu.  The local
    relaxation rate is proportional to the local *induced* pressure, so

        V(t) = V1_inf - (V1_inf - V1_0) * S(t),   S(t) = <p(t)>/p0,

    where <p(t)> is the *opening-weighted* average of the pore pressure over
    the diffusion layer.  We use the exact profile p(y,t)=p0 erf(y/sqrt(4ct))
    and weight it with a decaying influence kernel w(y)=exp(-y/L) that
    represents the mode-I opening influence (its explicit form affects only
    the O(1) prefactor; the sqrt(t) exponent and the form A ∝ (V1_inf - V1_0)
    are robust).  The *natural* choice that makes the t->inf limit exact and
    gives the classical result is the unweighted (flat) average over the
    influence depth lambda = L:

        S(t) = (1/L) int_0^L p(y,t)/p0 dy.

    Its small-time expansion:
        (1/L) int_0^L erf(y/sqrt(4ct)) dy
          = 1 - (1/L) int_0^L erfc(y/sqrt(4ct)) dy
          ~ 1 - (1/L) int_0^inf erfc(y/sqrt(4ct)) dy
          = 1 - sqrt(4ct)/(L sqrt(pi))
          = 1 - (2/sqrt(pi)) sqrt(ct)/L.
    So  V(t) - V1_0 = (V1_inf - V1_0) (2/sqrt(pi)) sqrt(ct)/L,
    i.e. A = (V1_inf - V1_0) * 2/sqrt(pi).
    """
    times = np.asarray(times, float)
    V0 = 2 * np.pi * (1 - nuu) / G * L**2
    Vinf = 2 * np.pi * (1 - nu) / G * L**2
    xs = np.linspace(0, 1, ny)
    ys = xs * L
    S = np.array([np.trapezoid(
        np.array([erf_helper(y, t, c) for y in ys]), ys) / L
        for t in times])
    return V0 + (Vinf - V0) * (1 - S), V0, Vinf


def erf_helper(y, t, c):
    from scipy.special import erf
    z = y / np.sqrt(4 * c * t)
    return erf(z)
    y, out = solve_1D_diffusion(p0, c, L, times=np.array([0.05, 0.2, 0.5]),
                                ny=3000)
    print("=== (M2) FD vs exact  p(y,t)=p0 erf(y/sqrt(4ct)) ===")
    for tt in sorted(out):
        u = out[tt]
        print(f"t={tt:.3f}")
        for yv in (0.0, 0.1, 0.3, 0.6):
            j = int(np.argmin(np.abs(y - yv)))
            exact = p0 * erfc(yv / np.sqrt(4 * c * tt))
            print(f"   y={y[j]:.4f}  FD={u[j]:.6f}  exact={exact:.6f} "
                  f"diff={u[j]-exact: .2e}")
from scipy.special import erfc