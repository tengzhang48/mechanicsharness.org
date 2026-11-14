import numpy as np, time, sys
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem

def run(tag, expect_slope=False, **kw):
    t0=time.time()
    p = WavyFrontProblem(**kw)
    p.solve()
    th, K1, K2, K3 = p.K_along_front()
    out = f"{tag}: relres {p.fitinfo['relres']:.2e} " \
          f"K1err {K1.mean()-1.0:+.4f} K2err {K2.mean()-0.6:+.4f} K3err {K3.mean()-0.8:+.4f}"
    if expect_slope:
        f = p.fit_partA()
        out += (f"\n    C1 {f['C1']:+.3f} S1 {f['S1']:+.3f}  C2 {f['C2']:+.3f} "
                f"S2 {f['S2']:+.3f}  C3 {f['C3']:+.3f} S3 {f['S3']:+.3f}")
    print(out + f" [{time.time()-t0:.0f}s]")
    return p

print("--- straight & shift invariance ---")
run("A=0 shift 0   ", A=0.0, harm=1, Nxf=4096, Nz=64, M=18, Nzmode=2, Lcr=2.0, Lx=3.5, taper=(0.30,0.95))
run("A=0 shift+.05 ", A=0.0, harm=1, shift=0.05, Nxf=4096, Nz=64, M=18, Nzmode=2, Lcr=2.0, Lx=3.5, taper=(0.30,0.95))
run("A=0 shift-.05 ", A=0.0, harm=1, shift=-0.05, Nxf=4096, Nz=64, M=18, Nzmode=2, Lcr=2.0, Lx=3.5, taper=(0.30,0.95))
print("--- first wavy test A=0.02 (predictions: C1~-pi, S1~0, C2~-0.6pi, S2~+1.6pi, C3~-0.8pi, S3~-1.2pi) ---")
p = run("A=0.02       ", A=0.02, harm=1, Nxf=4096, Nz=64, M=18, Nzmode=4, Lcr=2.0, Lx=3.5,
        taper=(0.30,0.95), expect_slope=True)
