import numpy as np, time, sys, json
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

res = []
def run(tag, **kw):
    t0=time.time()
    p = WavyFrontProblem(**kw)
    p.solve()
    f = p.fit_partA()
    row = dict(tag=tag, relres=p.fitinfo['relres'], nsv=p.fitinfo['nsv'],
               K0err=f['K0_1']-1.0, K20err=f['K0_2']-0.6, K30err=f['K0_3']-0.8,
               C1=f['C1'], S1=f['S1'], C2=f['C2'], S2=f['S2'],
               C3=f['C3'], S3=f['S3'], t=time.time()-t0, **kw)
    res.append(row)
    print(f"{tag}: relres {row['relres']:.1e} K0err {row['K0err']:+.4f} "
          f"C1 {f['C1']:+.3f} S1 {f['S1']:+.3f} C2 {f['C2']:+.3f} S2 {f['S2']:+.3f} "
          f"C3 {f['C3']:+.3f} S3 {f['S3']:+.3f} [{row['t']:.0f}s]", flush=True)

base = dict(harm=1, Lcr=2.0, Lx=3.5, taper=(0.30,0.95))
for A in [0.02, 0.01, 0.04]:
    run(f"A={A} Nzm4 M18 Nx4096", A=A, Nxf=4096, Nz=64, M=18, Nzmode=4, **base)
run("A=0.02 Nzm8 M18 Nx4096", A=0.02, Nxf=4096, Nz=64, M=18, Nzmode=8, **base)
run("A=0.02 Nzm8 M22 Nx8192 Nz128", A=0.02, Nxf=8192, Nz=128, M=22, Nzmode=8, **base)
run("A=0.01 Nzm8 M22 Nx8192 Nz128", A=0.01, Nxf=8192, Nz=128, M=22, Nzmode=8, **base)
run("A=0.005 Nzm8 M22 Nx8192 Nz128", A=0.005, Nxf=8192, Nz=128, M=22, Nzmode=8, **base)
json.dump(res, open('/work/sweep1.json','w'), default=float)
