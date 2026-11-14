import numpy as np, sys, time, json
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

res=[]
def run(tag, **kw):
    o = dict(ridge=kw.pop('ridge',3e-3), margin_dx=kw.pop('margin_dx',15.0),
             cut=kw.pop('cut',1e-9))
    t0=time.time()
    p = WavyFrontProblem(**kw)
    prior = p.solve_prior()
    p.solve(prior=prior, **o)
    f = p.fit_partA()
    if 'C1' in f:
        msg=(f"{tag}: relres {p.fitinfo['relres']:.1e} K0err {f['K0_1']-1:+.4f} "
             f"C1 {f['C1']:+.3f} S1 {f['S1']:+.3f} C2 {f['C2']:+.3f} S2 {f['S2']:+.3f} "
             f"C3 {f['C3']:+.3f} S3 {f['S3']:+.3f}")
    else:
        msg=(f"{tag}: relres {p.fitinfo['relres']:.1e} K0err {f['K0_1']-1:+.4f} "
             f"K2err {f['K0_2']-0.6:+.4f} K3err {f['K0_3']-0.8:+.4f}")
    print(msg+f" [{time.time()-t0:.0f}s]", flush=True)
    res.append(dict(tag=tag, relres=p.fitinfo['relres'], **f, **o, **kw))

base = dict(harm=1, Lcr=2.0, Lx=3.5, taper=(0.30,0.95), A=0.02)
print("--- margin / ridge / resolution on wavy A=0.02 ---")
run("m15 r3e-3 M18 Nzm4 Nx4096", **base, Nxf=4096, Nz=64, M=18, Nzmode=4)
run("m25 r3e-3 M18 Nzm4 Nx4096", margin_dx=25.0, **base, Nxf=4096, Nz=64, M=18, Nzmode=4)
run("m25 r1e-3 M18 Nzm4 Nx4096", margin_dx=25.0, ridge=1e-3, **base, Nxf=4096, Nz=64, M=18, Nzmode=4)
run("m25 r1e-4 M18 Nzm4 Nx4096", margin_dx=25.0, ridge=1e-4, **base, Nxf=4096, Nz=64, M=18, Nzmode=4)
run("m25 r3e-3 M22 Nzm6 Nx8192", margin_dx=25.0, **base, Nxf=8192, Nz=96, M=22, Nzmode=6)
run("m25 r1e-3 M22 Nzm6 Nx8192", margin_dx=25.0, ridge=1e-3, **base, Nxf=8192, Nz=96, M=22, Nzmode=6)
json.dump(res, open('/work/sweep4.json','w'), default=float)
