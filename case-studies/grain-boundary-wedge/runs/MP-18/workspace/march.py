#!/usr/bin/env python3
"""
MP-18 : relaxation of the grain-boundary traction by diffusion in a thin film.

Physics (all lengths in units of h, x = z/h ;  t = t0 * tau ,
t0 = 4 pi kT h^3/(E* delta_gb D_gb Omega)):

  Elasticity  (problem statement; s = sigma_gb/sigma_0 , w = du/dx)
      s(x,tau) = 1 - (1/pi) Int_0^1 S(x,y) w(y,tau) dy ,
      S(x,y)   = 1/(y-x) + 1/(y+x) - 2x(y-x)/(y+x)^3      (PV at y=x) .

  Diffusion / mass conservation  (mu = -sigma_gb Omega)
      d_tau u = d_x^2 s ,            w = d_x u .

  BCs   s(0,tau) = 0 ;  s'(1,tau) = s''(1,tau) = 0 .
  IC    s(x,0) = 1   (uniform tension),   u(x,0) = 0 .

We march  u  directly (a diffusion problem, unconditionally stable with an implicit
step);  s  is reconstructed from u every step through the elasticity integral, which
is evaluated pointwise with the exact kernel antiderivatives (validated separately
against PV-aware Gauss quadrature).
"""
import numpy as np

# ------------------------------------------------------------ kernel integrals
def _P1(y, z):
    return np.log(np.abs(y - z) + 1e-300)


def _yP1(y, z):
    return y + z * np.log(np.abs(y - z) + 1e-300)


def _P2(y, z):
    return np.log(np.abs(y + z) + 1e-300)


def _yP2(y, z):
    return y - z * np.log(np.abs(y + z) + 1e-300)


def _P3(y, z):
    return 4.0 * z / (y + z) - 2.0 * z ** 2 / (y + z) ** 2


def _yP3(y, z):
    t = y + z
    return -2.0 * z * (np.log(np.abs(t) + 1e-300) + 3.0 * z / t - z ** 2 / t ** 2)


def _prim(y, z, c0, c1):
    return (c0 * _P1(y, z) + c1 * _yP1(y, z)
            + c0 * _P2(y, z) + c1 * _yP2(y, z)
            + c0 * _P3(y, z) + c1 * _yP3(y, z))


def Smat(N):
    """M_ij = (1/pi) PV Int_0^1 S(x_i,y) phi_j(y) dy  for the hat basis phi_j on the
    uniform grid with nodes 0 = y_0 < x_1 < ... < x_N < 1 = y_{N+1}.
    Rows i = 1..N are the interior collocation points x_i = i/(N+1)."""
    h = 1.0 / (N + 1)
    x = np.arange(1, N + 1) * h
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = np.zeros((N, N))
    for i in range(N):
        z = x[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            val = (_prim(b, z, -a / (b - a), 1.0 / (b - a))
                   - _prim(a, z, -a / (b - a), 1.0 / (b - a)))
            val += (_prim(c, z, c / (c - b), -1.0 / (c - b))
                    - _prim(b, z, c / (c - b), -1.0 / (c - b)))
            M[i, j] = val / np.pi
    return M


def D1(N, h=None):
    """First-derivative matrix on interior nodes x_j = j/(N+1), with s(0)=0
    (ghost at x=0 gives the row j=1) and s'(1)=0 as a boundary row."""
    if h is None:
        h = 1.0 / (N + 1)
    D = np.zeros((N, N))
    for j in range(N):
        if j == 0:
            # centred difference; the left neighbour is the node at x=0 with value 0
            D[j, j + 1] = 1.0 / (2 * h)
        elif j == N - 1:
            D[j, j - 1] = -1.0 / (2 * h)      # s'(1)=0 -> no ghost beyond x=1
        else:
            D[j, j - 1] = -1.0 / (2 * h)
            D[j, j + 1] = 1.0 / (2 * h)
    return D


def D2(N, h=None):
    """Second-derivative matrix on interior nodes with s(0)=0 and s'(1)=0."""
    if h is None:
        h = 1.0 / (N + 1)
    D = np.zeros((N, N))
    for j in range(N):
        D[j, j] = -2.0
        if j == 0:
            D[j, j + 1] = 1.0
        elif j == N - 1:
            D[j, j - 1] = 2.0                 # s'(1)=0 -> ghost = interior neighbour
        else:
            D[j, j - 1] = 1.0
            D[j, j + 1] = 1.0
    return D / h ** 2


def D3(N, h=None):
    """Third-derivative matrix on interior nodes with the boundary conditions
    s(0)=0 , s'(1)=0 , s''(1)=0  (eliminate the ghost nodes)."""
    if h is None:
        h = 1.0 / (N + 1)
    D = np.zeros((N, N))
    for j in range(N):
        if j == 0:
            # s_{-1} ghost: s_0=0 and centred ---> use one-sided 3rd derivative
            # s'''_1 = (-s_0 + 3s_1 - 3s_2 + s_3)/h^3 with s_0=0
            D[j, j] = 3.0
            D[j, j + 1] = -3.0
            D[j, j + 2] = 1.0
        elif j == N - 1:
            # s'(1)=s''(1)=0  ->  ghost nodes s_{N+1}=s_{N-1} , s_{N+2}=s_{N-2}
            D[j, j - 2] = 2.0
            D[j, j - 1] = -2.0
        else:
            D[j, j - 1] = -2.0
            D[j, j] = 2.0
            if j + 2 <= N - 1:
                D[j, j + 1] = -2.0
                D[j, j + 2] = 2.0
            elif j + 1 <= N - 1:
                # needs ghost s_{N+1}=s_{N-1}
                D[j, j + 1] = -2.0
                D[j, j - 1] += 2.0
            else:
                raise RuntimeError
    return D / h ** 3


def operator(N):
    """Return the NxN matrix L with  d_tau s = L s , built from the two stated
    equations only:
        s = 1 - M w ,      d_tau w = D1 D2 s    (since d_tau u = D2 s)
    =>  d_tau s = -M D1 D2 s .   We use  L = -M @ D1 @ D2 .
    All matrices act on the interior node values."""
    M = Smat(N)
    D1m = D1(N)
    D2m = D2(N)
    L = -(M @ D1m @ D2m)
    return M, D1m, D2m, L
    """M_ij = (1/pi) PV Int_0^1 S(x_i,y) phi_j(y) dy  for the hat basis phi_j on the
    uniform grid with nodes 0 = y_0 < x_1 < ... < x_N < 1 = y_{N+1}.
    Rows i = 1..N are the interior collocation points x_i = i/(N+1)."""
    h = 1.0 / (N + 1)
    x = np.arange(1, N + 1) * h
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = np.zeros((N, N))
    for i in range(N):
        z = x[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            val = (_prim(b, z, -a / (b - a), 1.0 / (b - a))
                   - _prim(a, z, -a / (b - a), 1.0 / (b - a)))
            val += (_prim(c, z, c / (c - b), -1.0 / (c - b))
                    - _prim(b, z, c / (c - b), -1.0 / (c - b)))
            M[i, j] = val / np.pi
    return M