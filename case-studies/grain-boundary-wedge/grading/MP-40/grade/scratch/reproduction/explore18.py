import numpy as np

exec(open('<round-runs>/runs/MP-40/explore13.py').read().split('print("Correct')[0])

K = 22
G, H, Q, NQ = build(K)
ev, V = np.linalg.eig(np.linalg.solve(H, G))
cands = [( -ev[k].real, k) for k in range(len(ev)) if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k])) and -ev[k].real > 0]
cands.sort()
order = [k for _, k in cands]
xs = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
idx = np.argsort(XI)
for n, k in enumerate(order[:3]):
    a = np.real(V[:, k])
    mode = (Q @ a)
    mode = mode/np.max(np.abs(mode))
    # fix sign so mode(1)>0
    if np.interp(1.0, XI[idx], mode[idx]) < 0:
        mode = -mode
    interp = np.interp(xs, XI[idx], mode[idx])
    f1 = np.interp(1.0, XI[idx], mode[idx])
    print("lam=%.5f  f/f(1) at x=%s : %s" % (cands[n][0], list(xs), np.round(interp/f1, 6)))
# near-boundary behavior of mode 1
a = np.real(V[:, order[0]]); mode = (Q @ a)
if np.interp(1.0, XI[idx], mode[idx]) < 0:
    mode = -mode
f1 = np.interp(1.0, XI[idx], mode[idx])
print("\nnear x=0 (should vanish):")
for x in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]:
    print("  x=%.0e  f/f(1)=%.6e" % (x, np.interp(x, XI[idx], mode[idx])/f1))
xx = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
print("  power fit near 0:", np.polyfit(np.log(xx), np.log(np.abs(np.interp(xx, XI[idx], mode[idx])/f1)), 1))
print("near x=1:")
for x in [0.999, 0.995, 0.99, 0.98, 0.95]:
    print("  1-x=%.0e  1-f/f1=%.6e" % (1-x, 1-np.interp(x, XI[idx], mode[idx])/f1))