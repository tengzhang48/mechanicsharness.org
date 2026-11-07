"""Extrapolate s->inf for the tip quantities a(gamma) (K factor) and chat(gamma)
(tip volume coefficient: V_W - V_u ~ (1-nu_u) chat sqrt(pi)/2 ... in Laplace: (sV_W - V_u) sqrt(s)/(1-nu_u) -> chat).
Compare with first-order (weak-coupling) theory: a = 1 - 3 gamma/4, chat = 2 pi gamma."""
import numpy as np, json, glob
rows = []
for f in sorted(glob.glob('a_of_gamma_*.json')):
    rows += json.load(open(f))
rows.sort(key=lambda r: r['gamma'])
print(' gamma    a(gamma)    1-3g/4     sqrt(1-3g/2)  chat(gamma)   2*pi*g    chat/(2 pi g)')
out = []
for r in rows:
    g = r['gamma']; d = r['data']
    s = np.array([x['s'] for x in d]); W = np.array([x['W'] for x in d]); cW = np.array([x['cW'] for x in d])
    # a: fit a + b s^-1/2 + c s^-1 (3 points)
    A = np.vstack([np.ones(3), s**-0.5, s**-1.0]).T; a = np.linalg.solve(A, W)[0]
    # chat: fit c - d s^-1/4 + e s^-1/2
    B = np.vstack([np.ones(3), s**-0.25, s**-0.5]).T; ch = np.linalg.solve(B, cW)[0]
    out.append(dict(gamma=g, a=a, chat=ch))
    print(f"{g:6.3f}  {a:.7f}  {1-0.75*g:.7f}  {np.sqrt(1-1.5*g):.7f}  {ch:.5f}  {2*np.pi*g:.5f}  {ch/(2*np.pi*g):.5f}")
json.dump(out, open('a_chat_vs_gamma.json', 'w'), indent=1)
