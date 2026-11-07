#!/usr/bin/env python3
"""Part C: loading-1 transient. p diffuses (drained crack faces), eigenstress
trace via 2D FFT -> crack opening; check V1(0+), short-time sqrt coefficient A,
and K1(t) constancy. Plus numerical Eshelby-disk check of the trace kernel."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import time as time_mod

nu, nu_u, alpha = 0.2, 0.4, 0.89
G, L = 1.0, 1.0
kappa = 3 - 4 * nu
gamma = alpha * (1 - 2 * nu) / (1 - nu)
V1_0 = np.pi * (1 - nu_u)
V1_inf = np.pi * (1 - nu)
A_pred = 8 * (nu_u - nu) / (np.sqrt(np.pi) * (1 - nu))
C0 = 2 * (nu_u - nu) / (alpha * (1 - 2 * nu))     # Skempton face value

# ---------- quick numerical Eshelby-disk check of kernel ----------
print("Eshelby disk check of smooth-field kernel (smoothed disk r<1, value 1):")
A_b, Nf = 8.0, 1024
xs = (np.arange(Nf) - Nf // 2) * (2 * A_b / Nf)
Xg, Yg = np.meshgrid(xs, xs, indexing="xy")
R = np.sqrt(Xg**2 + Yg**2)
P = 0.5 * (1 - np.tanh((R - 1.0) / 0.05))
kx = 2 * np.pi * np.fft.fftfreq(Nf, d=2 * A_b / Nf)
KX, KY = np.meshgrid(kx, kx, indexing="xy")
K2 = KX**2 + KY**2
ker = np.zeros_like(K2); ker[K2 > 0] = KX[K2 > 0] ** 2 / K2[K2 > 0]
Syy = np.fft.ifft2(-gamma * ker * np.fft.fft2(P)).real
i0 = Nf // 2
print("  sigma_yy^sm at disk center = %.5f   (exact inclusion value -(gamma/2) = %.5f)"
      % (Syy[i0, i0], -gamma / 2))
print("  sigma_xx^sm (via s=-gamma p)       = %.5f" % (-gamma * P[i0, i0] - Syy[i0, i0]))

# ---------- crack solver (as in verify.py) ----------
def crack_VK(q_nodes, Nc=161, Nm=120):
    j = np.arange(Nc)
    theta = (j + 0.5) * np.pi / Nc
    x = L * np.cos(theta)
    n = np.arange(1, Nm + 1)
    M = np.sin(np.outer(theta, n))
    b = -L * (kappa + 1) * q_nodes * np.sin(theta) / (2 * G)
    a, *_ = np.linalg.lstsq(M, b, rcond=None)
    odd = (n % 2 == 1)
    V = -np.pi * L / 2 * np.sum(a[odd] / n[odd])
    K = np.sum(q_nodes * np.sqrt((L + x) / (L - x)) * np.pi / Nc * np.sin(theta)) / np.sqrt(np.pi * L)
    return V, K

# ---------- FD grid for diffusion ----------
Xmax, Ymax = 4.0, 4.0
Nx, Ny = 800, 500
x = np.linspace(0.0, Xmax, Nx)
y = Ymax * (np.linspace(0.0, 1.0, Ny) ** 2)      # clustered near y=0
dx = x[1] - x[0]
dyv = np.diff(y)
Ntot = Nx * Ny
idx = lambda i, j: j * Nx + i      # i: x-index, j: y-index

rows, cols, vals = [], [], []
def add(r, c, v):
    rows.append(r); cols.append(c); vals.append(v)

for j in range(Ny):
    for i in range(Nx):
        r = idx(i, j)
        if j == 0:
            # y=0 boundary: Dirichlet p=0 for x<L ; Neumann p[0]=p[1] for x>L
            if x[i] < L:
                add(r, r, 1.0)
            else:
                add(r, r, 1.0); add(r, idx(i, 1), -1.0)
            continue
        if i == Nx - 1 or j == Ny - 1:
            add(r, r, 1.0); continue     # far field p=0
        if i == 0:
            add(r, r, 1.0); add(r, idx(1, j), -1.0); continue   # symmetry dp/dx=0
        # interior: d/dy( dy p ) + d2/dx2
        dyl, dyu = dyv[j - 1], dyv[j]
        coef0 = 2 / (dyl * dyu)
        add(r, idx(i, j), -coef0 - 2 / dx**2)
        add(r, idx(i, j - 1), 2 / (dyl * (dyl + dyu)))
        add(r, idx(i, j + 1), 2 / (dyu * (dyl + dyu)))
        add(r, idx(i - 1, j), 1 / dx**2)
        add(r, idx(i + 1, j), 1 / dx**2)
Lap = sp.coo_matrix((vals, (rows, cols)), shape=(Ntot, Ntot)).tocsc()

# initial condition: Skempton field of undrained elastic solution
z = x[None, :] + 1j * y[:, None]
z = x[None, :] + 1j * (y[:, None] + 1e-6)   # tiny offset: tips stay finite (1/sqrt singularity integrable)
sq = np.sqrt(z * z - L * L + 0j)
p_init = C0 * (1 - np.real(z / sq))

dt = 5e-6
Aimp = (sp.identity(Ntot) - dt * Lap).tocsc()
# rows that are BC equations must get their own RHS (0), not p^n
bc_rows = []
for j in range(Ny):
    for i in range(Nx):
        r = idx(i, j)
        if j == 0 or i == Nx - 1 or j == Ny - 1 or i == 0:
            bc_rows.append(r)
bc_rows = np.array(bc_rows)
t0 = time_mod.time()
lu = spla.splu(Aimp)
print("\nFD diffusion: Nx=%d Ny=%d, splu done in %.1f s" % (Nx, Ny, time_mod.time() - t0))

# FFT projection grid (uniform, full plane by mirroring)
Nf2 = 2048
Af = 4.0
xf = (np.arange(Nf2) - Nf2 // 2) * (2 * Af / Nf2)
kxf = 2 * np.pi * np.fft.fftfreq(Nf2, d=2 * Af / Nf2)
KXf, KYf = np.meshgrid(kxf, kxf, indexing="xy")
K2f = KXf**2 + KYf**2
kerf = np.zeros_like(K2f); kerf[K2f > 0] = KXf[K2f > 0] ** 2 / K2f[K2f > 0]
xh = xf[Nf2 // 2:]                      # 0..Af
kcol = kxf
row0 = Nf2 // 2

from scipy.interpolate import RegularGridInterpolator
Yq, Xq = np.meshgrid(xh, xh, indexing="ij")

def process(pfield):
    interp = RegularGridInterpolator((y, x), pfield, bounds_error=False, fill_value=0.0)
    Ph = interp(np.stack([Yq, Xq], axis=-1))
    Pfull = np.zeros((Nf2, Nf2))
    Pfull[Nf2 // 2:, Nf2 // 2:] = Ph
    Pfull[Nf2 // 2:, :Nf2 // 2] = Ph[:, ::-1]
    Pfull[:Nf2 // 2, :] = Pfull[Nf2 // 2:, :][::-1, :]
    Syy = np.fft.ifft2(-gamma * kerf * np.fft.fft2(Pfull)).real
    qrow = Syy[row0, :]         # crack-solver input = sigma_yy^sm on the faces
    return xf, qrow

def p_field_from_vec(pv):
    return pv.reshape(Ny, Nx)

# t=0+ check
xf, q0 = process(p_init)
q_nodes0 = np.interp(L * np.cos((np.arange(161) + 0.5) * np.pi / 161), xf, q0)
Vc0, Kc0 = crack_VK(q_nodes0)
print("pipeline at t=0+: V = pi*(1-nu)+Vcrack = %.5f + (%.5f) = %.5f  [undrained exact %.5f]"
      % (np.pi * (1 - nu), Vc0, np.pi * (1 - nu) + Vc0, V1_0))
print("pipeline at t=0+: K/sqrt(pi L) = %.5f  [exact 1.0]" % Kc0)

pv = p_init.flatten().copy()
pv[idx(0, 0)] = 0.0
out_t, out_V, out_K = [0.0], [np.pi * (1 - nu) + Vc0], [Kc0]
nsteps = 800
for st in range(1, nsteps + 1):
    b = pv.copy(); b[bc_rows] = 0.0
    pv = lu.solve(b)
    if st % 100 == 0 or st == nsteps:
        pf = p_field_from_vec(pv)
        xf, qr = process(pf)
        qn = np.interp(L * np.cos((np.arange(161) + 0.5) * np.pi / 161), xf, qr)
        Vc, Kc = crack_VK(qn)
        out_t.append(st * dt); out_V.append(np.pi * (1 - nu) + Vc); out_K.append(Kc)
        print("t=%9.6f  V=%.5f  K=%.5f   (V1_0+A*sqrt(t) = %.5f)"
              % (st * dt, np.pi * (1 - nu) + Vc, Kc, V1_0 + A_pred * np.sqrt(st * dt)))

out_t = np.array(out_t); out_V = np.array(out_V); out_K = np.array(out_K)
sel = (out_t > 4e-4) & (out_t < 4e-3)
cf = np.polyfit(np.sqrt(out_t[sel]), out_V[sel], 1)
print("\nfit over ct in [4e-4,4e-3]: A_fit = %.5f  (analytic %.5f, diff %.2f%%)"
      % (cf[0], A_pred, 100 * abs(cf[0] - A_pred) / A_pred))
print("fit intercept = %.5f (V1_0 = %.5f)" % (cf[1], V1_0))
print("K1(t)/sqrt(pi L): min %.5f max %.5f over outputs" % (out_K[1:].min(), out_K[1:].max()))
np.savez("transient_out.npz", t=out_t, V=out_V, K=out_K)
