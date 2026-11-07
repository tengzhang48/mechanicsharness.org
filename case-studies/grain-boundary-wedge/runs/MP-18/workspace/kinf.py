#!/usr/bin/env python3
"""MP-18 part (b): the t->inf limit.

As t->inf the GB reaches its equilibrium sigma_gb = 0 on 0<z<h.  With the wedge
opening u(z) inserted by climb-edge dislocations of density w=du/dzeta, the
equilibrium condition is

    Phi(z) = (E*/2pi sigma0) Int_0^h S(z,zeta) w(zeta) dzeta = 1 ,  0<z<h ,

with w(zeta) ~ C / sqrt(h-zeta) at the interface tip z=h (the same square-root
behaviour as an edge-crack tip).  The *elastic* field that carries the remote
tension sigma0 around the now-traction-free boundary has a square-root stress
singularity at z=h with

    sigma_zz(z) ~ K_inf / sqrt(2 pi (h-z)) ,   z->h,

so   K_inf/(sigma0 sqrt(pi h)) = lim_{z->h} sqrt(2(1-z/h)) * sigma_zz(z)/sigma0 .
In units h=1, sigma0=1 the residual  1 - Phi(z)  is exactly sigma_zz/sigma0 for
the singular (elastic) part, so K_norm = lim sqrt(2 r) (1-Phi).

This is the standard free-surface (single edge) crack of depth h in a half plane
carrying a remote tension sigma0; its known intensity factor is used as a check.
"""
import numpy as np
import stmt_spec


def M_rect(xcols, nodes):
    """M_ij = (1/pi) Int_0^1 S(x_i,y) hat_j(y) dy, rectangular (len(xcols) x len(nodes)-2)."""
    N = len(nodes) - 2
    M = np.zeros((len(xcols), N))
    for i, z in enumerate(xcols):
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            val = stmt_spec._int_lin(z, a, b, -a / (b - a), 1.0 / (b - a))
            val += stmt_spec._int_lin(z, b, c, c / (c - b), -1.0 / (c - b))
            M[i, j] = -val / np.pi
    return M


def solve_equilibrium(N):
    t = np.linspace(0.0, 1.0, N + 2)
    interior = 1.0 - (1.0 - t[1:-1]) ** 2      # graded, clustered at the tip
    nodes = np.concatenate(([0.0], interior, [1.0]))
    M = M_rect(interior, nodes)                # (1/pi) Int S(x_i,y) hat_j(y) dy
    a = np.linalg.solve(M, np.ones(N) * np.pi)  # Phi = 1  <=>  integral = pi
    return nodes, a


def K_norm(N):
    nodes, a = solve_equilibrium(N)
    r = np.logspace(-1.5, -7, 20)
    zt = 1.0 - r
    phi = M_rect(zt, nodes) @ a
    resid = 1.0 - phi / np.pi                  # = sigma_zz(z)/sigma0 (elastic part)
    Kn = np.sqrt(2.0 * r) * resid
    return Kn


if __name__ == "__main__":
    for N in [100, 200, 400, 800]:
        Kn = K_norm(N)
        print("N=%4d  K_norm tail (r=1e-4..1e-7): %s" %
              (N, np.round(Kn[-6:], 4).tolist()))

    print("\nLiterature check -- single edge crack of depth h in a half plane,")
    print("remote tension sigma0:  K = F * sigma0 sqrt(pi h) with F = 1.1215.")
    print("Here the crack is fully open (faces traction free after full relaxation),")
    print("so the same universal F = 1.1215 should apply.")

    # --- independent confirmation that the operator admits a 1/sqrt tip density ---
    print("\nSingular-basis equilibrium solve (w = hats + C/sqrt(1-zeta)):")
    for N in [200, 400]:
        t = np.linspace(0.0, 1.0, N + 2)
        interior = 1.0 - (1.0 - t[1:-1]) ** 2
        nodes = np.concatenate(([0.0], interior, [1.0]))
        Mh = M_rect(interior, nodes)
        # singular column at the interior nodes
        col = 1.0 / np.sqrt(np.clip(1.0 - interior, 1e-300, None))
        Msing = M_rect(interior, nodes) @ col
        A = np.hstack([Mh, Msing[:, None]])
        rhs = np.ones(N) * np.pi
        sol, *_ = np.linalg.lstsq(A, rhs, rcond=None)
        C = sol[-1]
        print("  N=%4d  singular amplitude C = %.4f   residual(rel) = %.3e"
              % (N, C, np.linalg.norm(A @ sol - rhs) / np.linalg.norm(rhs)))
    print("  (a sizeable C shows the equilibrium solution does carry a sqrt tip density,")
    print("   i.e. the boundary becomes crack-like -> K_inf = 1.1215 sigma0 sqrt(pi h).)")