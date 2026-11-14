import numpy as np, sys, time, json
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

res=[]
def run(tag, ridge=3e-3, use_prior=True, **kw):
    t0=time.time()
    p = WavyFrontProblem(**kw)
    prior = p.solve_prior(ridge=1e-4) if use_prior else None
    p.solve(ridge=ridge, prior=prior)
    f = p.fit_partA()
    if 'C1' in f:
        msg = (f"{tag}: relres {p.fitinfo['relres']:.1e} K0err {f['K0_1']-1:+.4f} "
               f"C1 {f['C1']:+.3f} S1 {f['S1']:+.3f} C2 {f['C2']:+.3f} S2 {f['S2']:+.3f} "
               f"C3 {f['C3']:+.3f} S3 {f['S3']:+.3f}")
    else:
        msg = (f"{tag}: relres {p.fitinfo['relres']:.1e} K0err {f['K0_1']-1:+.4f} "
               f"K2err {f['K0_2']-0.6:+.4f} K3err {f['K0_3']-0.8:+.4f}")
    print(msg + f" [{time.time()-t0:.0f}s]", flush=True)
    res.append(dict(tag=tag, relres=p.fitinfo['relres'], **f, ridge=ridge, **kw))
    return p, f

base = dict(harm=1, Lcr=2.0, Lx=3.5, taper=(0.30,0.95), Nxf=4096, Nz=64, M=18, Nzmode=4)
print("--- straight check with new defaults ---")
p0, f0 = run("A=0      ", A=0.0, **base)
print("--- A scan with ridge+prior ---")
for A in [0.01, 0.02, 0.04]:
    run(f"A={A}    ", A=A, **base)
print("--- ridge sensitivity at A=0.02 ---")
for rg in [1e-3, 3e-3, 1e-2]:
    run(f"ridge{rg}", A=0.02, ridge=rg, **base)
json.dump(res, open('/work/sweep3.json','w'), default=float)
