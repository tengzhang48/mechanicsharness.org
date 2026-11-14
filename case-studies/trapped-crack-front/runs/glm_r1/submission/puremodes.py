import numpy as np, sys, time, json
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

res=[]
def run(tag, K0, **kw):
    t0=time.time()
    p = WavyFrontProblem(K0=K0, **kw)
    p.solve()
    f = p.fit_partA()
    print(f"{tag}: relres {p.fitinfo['relres']:.1e} "
          f"K01 {f['K0_1']:+.4f} K02 {f['K0_2']:+.4f} K03 {f['K0_3']:+.4f}\n"
          f"    C1 {f['C1']:+.3f} S1 {f['S1']:+.3f} C2 {f['C2']:+.3f} S2 {f['S2']:+.3f} "
          f"C3 {f['C3']:+.3f} S3 {f['S3']:+.3f} [{time.time()-t0:.0f}s]", flush=True)
    res.append(dict(tag=tag, K0=list(K0), **f, **kw))
    return f

base = dict(A=0.01, harm=1, Lcr=2.0, Lx=3.5, taper=(0.30,0.95))
print("--- pure mode I: expect C2=C3=S2=S3=0, S1=0 ---")
run("pure I ", (1.0, 0.0, 0.0), Nxf=4096, Nz=64, M=18, Nzmode=4, **base)
print("--- pure mode II: expect S1=0? C3=S3=0? (no K3^0, no K1^0) ---")
run("pure II", (0.0, 1.0, 0.0), Nxf=4096, Nz=64, M=18, Nzmode=4, **base)
print("--- pure mode III ---")
run("pure III",(0.0, 0.0, 1.0), Nxf=4096, Nz=64, M=18, Nzmode=4, **base)
print("--- mixed (0,0.6,0.8) for reference ---")
run("mixed  ", (1.0, 0.6, 0.8), Nxf=4096, Nz=64, M=18, Nzmode=4, **base)
json.dump(res, open('/work/puremodes.json','w'), default=float)
