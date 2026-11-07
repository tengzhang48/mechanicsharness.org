"""Solver 2 (independent): f-formulation Chebyshev-Lobatto collocation.

Unknown: f at Chebyshev-Lobatto points x[0]=1 ... x[N]=0 on [0,1].
Equation rows (interior i):  -lam f(x_i) = K[f'''](x_i)
  K[g](x) = PV-int_0^1 S(x,xi) g(xi) dxi,  computed by
  sum_j W_j g_j k_img(x_i,x_j)                       (image kernels, regular)
+ sum_{j!=i} W_j (g_j - g_i)/(x_i - x_j) - W_i g'(x_i)  (PV subtraction)
+ g_i ln(x_i/(1-x_i))
with W = Clenshaw-Curtis weights (exact for degree <= N), g = f''' from
stable Chebyshev-coefficient differentiation.
BC rows: f(0)=0, f'(1)=0, f''(1)=0.  Pencil A v = mu B v, mu = -lam.
"""
import numpy as np


def lobatto_grid(N):
    return (1.0 - np.cos(np.arange(N + 1) * np.pi / N)) / 2.0


def cheb_V(x, N):
    y = 2.0 * np.asarray(x) - 1.0
    V = np.zeros((y.size, N + 1))
    V[:, 0] = 1.0
    if N >= 1:
        V[:, 1] = y
    for n in range(1, N):
        V[:, n + 1] = 2 * y * V[:, n] - V[:, n - 1]
    return V


def dcoef_matrix(N):
    D = np.zeros((N + 1, N + 1))
    for n in range(N + 1):
        c = np.zeros(N + 1); c[n] = 1.0
        d = np.zeros(N + 3)
        for k in range(N, 0, -1):
            d[k - 1] = d[k + 1] + 2.0 * k * c[k]
        d[0] *= 0.5
        D[:, n] = 2.0 * d[:N + 1]
    return D


def assemble(N):
    x = lobatto_grid(N)
    V = cheb_V(x, N)
    moments = np.zeros(N + 1)
    moments[0] = 1.0
    for n in range(2, N + 1, 2):
        moments[n] = 1.0 / (1.0 - n * n)
    W = np.linalg.solve(V.T, moments)          # CC weights on [0,1]
    D = dcoef_matrix(N)
    Vinv = np.linalg.inv(V)
    D1v = V @ D @ Vinv
    D2v = V @ (D @ D) @ Vinv
    D3v = V @ (D @ D @ D) @ Vinv

    A = np.zeros((N + 1, N + 1))               # rows: K[f'''] (interior), BC rows
    B = np.zeros((N + 1, N + 1))
    for i in range(1, N):
        xi_ = x[i]
        kimg = -1.0 / (xi_ + x) - 2.0 * x * (xi_ - x) / (xi_ + x) ** 3
        row = W * kimg
        j = np.arange(N + 1)
        mask = j != i
        sub = np.where(mask, W / (xi_ - x), 0.0)
        row = row + sub
        row = row - sub.sum() * np.eye(N + 1)[i]   # -g_i * sum_j W_j/(x_i-x_j)
        row = row - W[i] * D1v[i, :]               # -g'(x_i) weight W_i
        row = row + np.log(xi_ / (1.0 - xi_)) * np.eye(N + 1)[i]
        A[i, :] = (row @ V) @ (D @ D @ D)
        B[i, :] = V[i, :]
    A[0, :] = D1v[0, :]                         # f'(1) = 0
    A[N, :] = np.eye(N + 1)[N]                  # f(0) = 0
    A[N - 1, :] = D2v[0, :]                     # f''(1)=0 (replaces eq at i=N-1)
    B[N - 1, N - 1] = 0.0                       # BC row must not carry eigen-relation
    return x, V, Vinv, W, D1v, D2v, D3v, A, B
