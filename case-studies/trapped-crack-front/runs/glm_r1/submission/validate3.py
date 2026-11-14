import numpy as np, time, sys
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

def run(tag, **kw):
    t0=time.time()
    p = WavyFrontProblem(**kw)
    p.solve()
    th, K1, K2, K3 = p.K_along_front()
    print(f"{tag}: nsv {p.fitinfo['nsv']} relres {p.fitinfo['relres']:.2e} "
          f"K1err {K1.mean()-1.0:+.4f} K2err {K2.mean()-0.6:+.4f} "
          f"K3err {K3.mean()-0.8:+.4f} dev1 {np.abs(K1-K1.mean()).max():.1e} "
          f"[{time.time()-t0:.0f}s]")

print("--- straight crack: taper/M/cut sweep ---")
for M in [14, 18, 22]:
    run(f"M={M} taper(0.45,0.95)", A=0.0, harm=1, Nxf=4096, Nz=64, M=M, Nzmode=2, Lcr=2.0, Lx=3.5)
print("--- wider taper (0.3,0.95) ---")
for M in [14, 18]:
    run(f"M={M} taper(0.30,0.95)", A=0.0, harm=1, Nxf=4096, Nz=64, M=M, Nzmode=2, Lcr=2.0, Lx=3.5, taper=(0.30,0.95))
print("--- shift invariance ---")
run("shift=0.05 A=0  ", A=0.0, harm=1, shift=0.05, Nxf=4096, Nz=64, M=18, Nzmode=2, Lcr=2.0, Lx=3.5)
run("shift=-0.05 A=0 ", A=0.0, harm=1, shift=-0.05, Nxf=4096, Nz=64, M=18, Nzmode=2, Lcr=2.0, Lx=3.5)
