"""Final s->infinity extrapolations for nu=0.2, nu_u=0.4, alpha=0.89 (gamma = 0.25).
Laplace-domain data (converged in N and quadrature) -> t->0+ limits via the initial-value theorem:
  K1(0+) = lim s K1~ ;  A (V1 ~ V_u + A sqrt(t)) = (2/sqrt(pi)) lim (s V1~ - V_u) sqrt(s)."""
import numpy as np, json, glob
from crack_solver2 import Params
P = Params(0.2, 0.4, 0.89)
data = {}
# converged values (highest N) from the runs
for f in glob.glob('large_s_*.json'):
    for r in json.load(open(f)):
        if r['s'] not in data or r['N'] > data[r['s']]['N']:
            data[r['s']] = r
extra = {100.0: (0.8650232778471896, 2.028587164286242, -0.0734811080445404, -0.11336879472227795),
         1e3: (0.83963065, 1.93942880, -0.04170787, -0.04167034),
         1e4: (0.82546865, 1.90398655, -0.02352233, -0.01433188),
         1e5: (0.81759171, 1.89131882, -0.01323970, -0.00475444)}
for s, v in extra.items():
    if s not in data:
        data[s] = dict(s=s, sK1=v[0], sV1=v[1], sK2=v[2], sV2=v[3])
ss = np.array(sorted(data)); D = [data[s] for s in ss]
sK1 = np.array([d['sK1'] for d in D]); sV1 = np.array([d['sV1'] for d in D])
sK2 = np.array([d['sK2'] for d in D]); sV2 = np.array([d['sV2'] for d in D])
b2 = 2*P.beta; Vu = np.pi*(1 - P.nuu)
W = sK1 + b2*sK2; VW = sV1 + b2*sV2
print(' s        sK1        W=s(K1+2bK2)   (sV1-Vu)rt(s)  (sV_W-Vu)rt(s)  sV2 rt(s)   s^1/4 sK2')
for i, s in enumerate(ss):
    print(f"{s:8.0e} {sK1[i]:.8f} {W[i]:.8f} {(sV1[i]-Vu)*np.sqrt(s):.6f} {(VW[i]-Vu)*np.sqrt(s):.6f} {sV2[i]*np.sqrt(s):.6f} {sK2[i]*s**0.25:.6f}")
def fit(x, y, powers):
    A = np.vstack([x**(-p) for p in powers]).T
    return np.linalg.lstsq(A, y, rcond=None)[0]
sel = ss >= 1e3
a = fit(ss[sel], W[sel], [0, 0.5, 1.0])[0]
pw = [0, 0.25, 0.5, 0.75] if sel.sum() >= 5 else [0, 0.25, 0.5]   # V-type data: s^-1/4 corrections
cW = fit(ss[sel], (VW[sel]-Vu)*np.sqrt(ss[sel]), pw)[0]
c2 = fit(ss[sel], sV2[sel]*np.sqrt(ss[sel]), pw)[0]
cK2 = fit(ss[sel], sK2[sel]*ss[sel]**0.25, [0, 0.5, 1.0])[0]
c1 = fit(ss[sel], (sV1[sel]-Vu)*np.sqrt(ss[sel]), pw)[0]
print('powers used for V fits:', pw, ' points:', int(sel.sum()))
from scipy.special import gamma as Gf
res = dict(a=a, K1_0=a, chat=cW/(1 - P.nuu), A_tip=2/np.sqrt(np.pi)*cW, A_face_numeric=-b2*2/np.sqrt(np.pi)*c2,
           A_face_formula=16*P.gam*(1 - P.nuu)/np.sqrt(np.pi), A_total_direct=2/np.sqrt(np.pi)*c1,
           A_total_sum=2/np.sqrt(np.pi)*cW - b2*2/np.sqrt(np.pi)*c2,
           A_exactface_plus_tip=16*P.gam*(1 - P.nuu)/np.sqrt(np.pi) + 2/np.sqrt(np.pi)*cW,
           K2_t14_coeff=cK2/Gf(1.25), K1_t14_coeff=-b2*cK2/Gf(1.25))
for k, v in res.items():
    print(f"{k:16s} {v:.6f}")
json.dump({k: float(v) for k, v in res.items()}, open('final_extrapolation.json', 'w'), indent=1)
