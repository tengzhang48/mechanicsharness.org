import numpy as np, sys, time, json
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

res=[]
def run(tag, **kw):
    opts = dict(cut=kw.pop('cut',1e-9), use_weights=kw.pop('use_weights',False),
                margin_dx=kw.pop('margin_dx',3.0))
    t0=time.time()
    p = WavyFrontProblem(**kw)
    p.solve(**opts)
    f = p.fit_partA()
    print(f"{tag}: relres {p.fitinfo['relres']:.1e} K0err {f['K0_1']-1:+.4f} "
          f"C1 {f['C1']:+.3f} S1 {f['S1']:+.3f} C2 {f['C2']:+.3f} S2 {f['S2']:+.3f} "
          f"C3 {f['C3']:+.3f} S3 {f['S3']:+.3f} [{time.time()-t0:.0f}s]", flush=True)
    res.append(dict(tag=tag, relres=p.fitinfo['relres'], **f, **opts, **kw))

base = dict(harm=1, Lcr=2.0, Lx=3.5, taper=(0.30,0.95), Nxf=4096, Nz=64, M=18, Nzmode=4)
print("--- weights sweep (A=0.02) ---")
run("A.02 noW m3 ", A=0.02, **base)
run("A.02 W m3   ", A=0.02, use_weights=True, **base)
run("A.02 noW m8 ", A=0.02, margin_dx=8.0, **base)
run("A.02 W m8 c-7", A=0.02, use_weights=True, margin_dx=8.0, cut=1e-7, **base)
print("--- A scan with weights m8 ---")
for A in [0.01, 0.02, 0.04]:
    run(f"A={A} W m8  ", A=A, use_weights=True, margin_dx=8.0, **base)
json.dump(res, open('/work/sweep2.json','w'), default=float)
