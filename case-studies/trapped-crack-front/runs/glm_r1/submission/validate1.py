import numpy as np, time, sys
sys.path.insert(0, '/work')
from crack3d import WavyFrontProblem, material

t0 = time.time()
mat = material()
print("=== V1: straight crack a=0 (harm=1, A=0) ===")
p = WavyFrontProblem(A=0.0, harm=1, Nxf=4096, Nz=64, M=16, Nzmode=4,
                     Lcr=2.0, Lx=3.5)
info = p.solve(verbose=True)
th, K1, K2, K3 = p.K_along_front()
print(f"  K1 = {K1.mean():.6f} (exact {p.K10}), max dev {np.abs(K1-K1.mean()).max():.2e}")
print(f"  K2 = {K2.mean():.6f} (exact {p.K20}), max dev {np.abs(K2-K2.mean()).max():.2e}")
print(f"  K3 = {K3.mean():.6f} (exact {p.K30}), max dev {np.abs(K3-K3.mean()).max():.2e}")
print(f"  elapsed {time.time()-t0:.1f}s")

t0 = time.time()
print("=== V2: constant shift a=0.05 (harm=0) ===")
p = WavyFrontProblem(A=0.05, harm=0, Nxf=4096, Nz=64, M=16, Nzmode=2,
                     Lcr=2.0, Lx=3.5)
info = p.solve(verbose=True)
th, K1, K2, K3 = p.K_along_front()
print(f"  K1 = {K1.mean():.6f} (exact {p.K10}), dev {np.abs(K1-K1.mean()).max():.2e}")
print(f"  K2 = {K2.mean():.6f} (exact {p.K20}), dev {np.abs(K2-K2.mean()).max():.2e}")
print(f"  K3 = {K3.mean():.6f} (exact {p.K30}), dev {np.abs(K3-K3.mean()).max():.2e}")
print(f"  elapsed {time.time()-t0:.1f}s")
