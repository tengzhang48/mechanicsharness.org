"""MP-46 independent poroelastic FEM check (plane strain, Biot).

Backward-Euler, linear triangles, double-noded crack slit.  Unknowns (ux,uy,p).
Run: python3 fem_poro.py
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

nu, nu_u, alpha = 0.2, 0.4, 0.89
G, kap, L = 1.0, 1.0, 1.0
lam = 2.0 * G * nu / (1.0 - 2.0 * nu)
W, h = 3.0, 0.05
N = int(round(2 * W / h))
xs = np.linspace(-W, W, N + 1)
jc = N // 2
iL, iR = int(round((-L + W) / h)), int(round((L + W) / h))
cracked = lambda i: iL <= i <= iR

# ---- nodes
xy = []
up, dn, sid = {}, {}, {}
for i in range(N + 1):
    if cracked(i):
        up[i] = len(xy); xy.append((xs[i], 0.0))
        dn[i] = len(xy); xy.append((xs[i], 0.0))
for j in range(N + 1):
    for i in range(N + 1):
        if j == jc and cracked(i):
            continue
        sid[(i, j)] = len(xy); xy.append((xs[i], xs[j]))
xy = np.array(xy)


def nid(i, jj, row):
    if jj != jc or not cracked(i):
        return sid[(i, jj)]
    return up[i] if row >= jc else dn[i]


tris = []
for j in range(N):
    for i in range(N):
        a = nid(i, j, j); b = nid(i + 1, j, j)
        d = nid(i + 1, j + 1, j); cc = nid(i, j + 1, j)
        tris.append((a, b, d)); tris.append((a, d, cc))
tris = np.array(tris)
nn = len(xy)
print("nodes", nn, "elems", len(tris))

# ---- element data
GRAD = np.zeros((len(tris), 3, 2)); AREA = np.zeros(len(tris))
for e, (n1, n2, n3) in enumerate(tris):
    P = xy[[n1, n2, n3]]
    J = np.array([[P[1, 0] - P[0, 0], P[2, 0] - P[0, 0]],
                  [P[1, 1] - P[0, 1], P[2, 1] - P[0, 1]]])
    AREA[e] = abs(np.linalg.det(J)) / 2.0
    GRAD[e] = np.array([[-1, -1], [1, 0], [0, 1]]) @ np.linalg.inv(J).T

D = np.array([[2 * G + lam, lam, 0], [lam, 2 * G + lam, 0], [0, 0, G]])


def Bmat(e):
    B = np.zeros((3, 6))
    for a in range(3):
        B[0, 2 * a] = GRAD[e, a, 0]
        B[1, 2 * a + 1] = GRAD[e, a, 1]
        B[2, 2 * a] = GRAD[e, a, 1]
        B[2, 2 * a + 1] = GRAD[e, a, 0]
    return B


su = alpha * (1 + nu)
a1 = alpha * (1 - 2 * nu) / (2 * G)
bcoef = -3 * alpha * a1 + alpha**2 * (1 - 2 * nu)**2 / (2 * G * (nu_u - nu))

# ---- assemble constant matrices
Kuu_r, Kuu_c, Kuu_v = [], [], []
Kup_r, Kup_c, Kup_v = [], [], []
Cp_r, Cp_c, Cp_v = [], [], []
Mp_r, Mp_c, Mp_v = [], [], []
Kp_r, Kp_c, Kp_v = [], [], []
for e in range(len(tris)):
    n = tris[e]; A = AREA[e]; B = Bmat(e); g = GRAD[e]
    Ke = A * (B.T @ D @ B)
    for a in range(3):
        for b in range(3):
            for ia in range(2):
                for ib in range(2):
                    Kuu_r.append(2 * n[a] + ia); Kuu_c.append(2 * n[b] + ib)
                    Kuu_v.append(Ke[2 * a + ia, 2 * b + ib])
            for ib in range(2):
                Kup_r.append(2 * n[a] + ib); Kup_c.append(n[b])
                Kup_v.append(-alpha * (A / 3.0) * g[a, ib])
                Cp_r.append(n[a]); Cp_c.append(2 * n[b] + ib)
                Cp_v.append((A / 3.0) * g[b, ib])
    me = (A / 12.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    for a in range(3):
        for b in range(3):
            Mp_r.append(n[a]); Mp_c.append(n[b]); Mp_v.append(me[a, b])
            Kp_r.append(n[a]); Kp_c.append(n[b])
            Kp_v.append(kap * A * (g[a, 0] * g[b, 0] + g[a, 1] * g[b, 1]))

Kuu = sp.csr_matrix((Kuu_v, (Kuu_r, Kuu_c)), shape=(2 * nn, 2 * nn))
Kup = sp.csr_matrix((Kup_v, (Kup_r, Kup_c)), shape=(2 * nn, nn))
Cp = sp.csr_matrix((Cp_v, (Cp_r, Cp_c)), shape=(nn, 2 * nn))
Mp = sp.csr_matrix((Mp_v, (Mp_r, Mp_c)), shape=(nn, nn))
Kp = sp.csr_matrix((Kp_v, (Kp_r, Kp_c)), shape=(nn, nn))
rowsum_p = np.array(Mp.sum(axis=1)).ravel()


def edges(ids):
    return [(ids[i], ids[i + 1]) for i in range(iL, iR)]


def run(loading, dt=2.0e-3, nsteps=400, elastic=False):
    ufix = np.zeros(2 * nn, bool); pfix = np.zeros(nn, bool)
    uval = np.zeros(2 * nn); pval = np.zeros(nn)
    for k, (x, y) in enumerate(xy):
        if abs(x) >= W - 1e-9 or abs(y) >= W - 1e-9:
            ufix[2 * k] = ufix[2 * k + 1] = True
            pfix[k] = True
    cr = [up[i] for i in range(iL, iR + 1)] + [dn[i] for i in range(iL, iR + 1)]
    pval[cr] = 1.0 if loading == 2 else 0.0
    pfix[cr] = True
    if elastic:            # drained-elastic check: p=0 everywhere, single solve
        pfix[:] = True; pval[:] = 0.0; nsteps = 1; dt = 1.0
    F = np.zeros(2 * nn)
    if loading == 1:
        for ids, ty in ((up, 1.0), (dn, -1.0)):
            for (a, b) in edges(ids):
                F[2 * a + 1] += ty * h / 2.0
                F[2 * b + 1] += ty * h / 2.0
    fixed = np.concatenate([ufix, pfix]); vals = np.concatenate([uval, pval])
    un = np.zeros(2 * nn); pn = np.zeros(nn)
    log = []
    A = sp.bmat([[Kuu, Kup], [(su / dt) * Cp, (bcoef / dt) * Mp + kap * Kp]]).tocsr()
    A = A.tolil()
    for d in np.where(fixed)[0]:
        A.rows[d] = [d]; A.data[d] = [1.0]
    A = A.tocsr()
    lu = spla.splu(A.tocsc())
    for step in range(nsteps):
        t = (step + 1) * dt
        du = np.zeros(nn)
        ux = un[0::2]; uy = un[1::2]
        d = np.zeros(len(tris))
        for a in range(3):
            na = tris[:, a]
            d += GRAD[:, a, 0] * ux[na] + GRAD[:, a, 1] * uy[na]
        wq = np.repeat(d * AREA / 3.0, 3)
        du = np.bincount(tris.ravel(), weights=wq, minlength=nn) / rowsum_p
        zeta_n = su * du + bcoef * pn
        rhs = np.concatenate([F, (1.0 / dt) * (Mp @ zeta_n)])
        rhs[fixed] = vals[fixed]
        sol = lu.solve(rhs)
        un = sol[:2 * nn]; pn = sol[2 * nn:]
        if step % 10 == 0 or step == nsteps - 1 or elastic:
            w = np.array([un[2 * up[i] + 1] - un[2 * dn[i] + 1]
                          for i in range(iL, iR + 1)])
            xc = xs[iL:iR + 1]
            V = np.trapezoid(w, xc) if hasattr(np, "trapezoid") else np.trapz(w, xc)
            # K from COD near the right tip:  w(r) = (4(1-nu)/G) K sqrt(r/2pi)
            Kd = Kal = None
            r = L - xc[-2]          # node next to the tip
            wn = w[-2]
            K_d = G / (4 * (1 - nu)) * np.sqrt(2 * np.pi / r) * wn
            K_u = G / (4 * (1 - nu_u)) * np.sqrt(2 * np.pi / r) * wn
            log.append((t, V, K_d, K_u, w[len(w) // 2]))
    return log, un, pn


if __name__ == "__main__":
    print("--- elastic calibration (drained, nu=0.2): expect V=%.5f K=1.0 ---"
          % (np.pi * (1 - nu)))
    for (t, V, Kd, Ku, wm) in run(1, elastic=True):
        print("  V=%9.5f  K(nu)=%9.5f  K(nu_u)=%9.5f" % (V, Kd, Ku))
    for ld in (1, 2):
        print("=== loading", ld, "===")
        for (t, V, Kd, Ku, wm) in run(ld, dt=5.0e-3, nsteps=300):
            print("  t=%7.4f  V=%9.5f  K(nu)=%9.5f  K(nu_u)=%9.5f  w(0)=%9.5f" %
                  (t, V, Kd, Ku, wm))