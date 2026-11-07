exec(open('<round-runs>/runs/MP-40/explore26.py').read().split('prev = {}')[0])

import numpy as np
import mpmath as mp


def build(M):
    basis = make_basis(M)
    n = len(basis); m = n
    nodes = [mp.mpf(0.5)*(1-mp.cos(mp.pi*(i+1)/(m+1))) for i in range(m)]
    On = mp.matrix(m, n); Vn = mp.matrix(m, n)
    for i in range(m):
        xi = nodes[i]
        for k, p in enumerate(basis):
            acc = mp.mpf(0)
            for d, c in enumerate(p):
                if abs(c) > 1e-14:
                    acc += mp.mpf(c)*Nval(d, xi)
            On[i, k] = acc
            Vn[i, k] = mp.polyval(list(p[::-1]), xi)
    return basis, nodes, On, Vn


for M in [20, 24, 28]:
    basis, nodes, On, Vn = build(M)
    n = len(basis)
    L = mp.inverse(Vn)*On
    ev, ER = mp.eig(L)
    # determine convention
    k0 = 0
    # pick eigenvalue with smallest positive -Re
    cands = sorted([(-complex(ev[k]).real, k) for k in range(n)
                    if abs(complex(ev[k]).imag) < 1e-9*max(1, abs(complex(ev[k]))) and -complex(ev[k]).real > 1e-6])
    lam1, k1 = cands[0]
    # test column
    colres = 0.0
    for i in range(n):
        val = sum(L[i, j]*ER[j, k1] for j in range(n))
        colres = max(colres, abs(val - ev[k1]*ER[j if False else i, k1]) if False else abs(val - ev[k1]*ER[i, k1]))
    # test row
    rowres = 0.0
    for i in range(n):
        val = sum(L[i, j]*ER[k1, j] for j in range(n))
        rowres = max(rowres, abs(val - ev[k1]*ER[k1, i]))
    print("M=%2d lam1=%.5f  colres=%.3e rowres=%.3e" % (M, lam1, float(colres), float(rowres)), flush=True)
    # time stepping to confirm
    a = mp.inverse(Vn)*mp.matrix([mp.mpf(1)]*n)
    dt = mp.mpf('0.01')
    Ainv = mp.inverse(mp.eye(n) - dt*L)
    norms = []
    tau = mp.mpf(0)
    for step in range(200):
        a = Ainv*a
        tau += dt
        nrm = mp.sqrt(sum(x*x for x in a))
        norms.append((float(tau), float(nrm)))
    norms = np.array(norms)
    mask = (norms[:, 0] > 0.3) & (norms[:, 0] < 1.5)
    p = np.polyfit(norms[mask, 0], np.log(norms[mask, 1]), 1)
    r = -p[0]
    lam_be = (np.exp(r*0.01)-1)/0.01
    print("   time-step rate=%.5f -> lam1=%.5f" % (r, lam_be), flush=True)
    if M == 24:
        # mode shape using the correct convention
        conv = 'col' if colres < rowres else 'row'
        vec = [ER[j, k1] if conv == 'col' else ER[k1, j] for j in range(n)]
        mx = max(abs(x) for x in vec)
        vec = [x/mx for x in vec]

        def sig(xv):
            return sum(vec[k]*mp.polyval(list(basis[k][::-1]), xv) for k in range(n))
        f1 = sig(mp.mpf(1))
        print("   convention=%s  f(0.5)/f(1)=%.6f f(0.25)=%.6f f(0.75)=%.6f" %
              (conv, float(sig(mp.mpf('0.5'))/f1), float(sig(mp.mpf('0.25'))/f1), float(sig(mp.mpf('0.75'))/f1)))