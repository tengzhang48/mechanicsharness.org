import sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigvals
from scipy.special import roots_legendre
sys.path.insert(0, str(Path(__file__).parent / 'reproduction'))
from solver5 import RhoOps, cheb_T
from solver7 import RhoOps32
from solver10 import Res10
from core import kernel_matrix
from crack import solve_crack, stress_ahead

xs = np.array([0.05,0.2,0.4,0.6,0.85,0.99])
for n in [16,24,32]:
    op = RhoOps(n)
    h = op.h_pieces(xs,0)[0]
    print('sqrt_h_vs_core', n, np.max(abs(h-kernel_matrix(xs,n))))
for n in [16,24,32,40,48]:
    a,xc,A = solve_crack(n)
    print('crack_N',n,'K',-np.pi*np.sqrt(2)*sum(a))
    if n==40:
        print('unused_stress_ahead_helper',stress_ahead(np.array([1.001]),a,n)[0]*np.sqrt(0.002))

for nq,kq in [(280,120),(560,240),(1120,480)]:
    r=Res10(32,Mquad=nq,Kleg=kq)
    ls=eigvals(-r.H,r.F)
    ls=sorted(v.real for v in ls if abs(v.imag)<1e-7 and v.real>0)
    print('quadrature',nq,kq,'first_six',ls[:6])

try:
    op = RhoOps(8)
    op.h_pieces(np.array([0.5,1.0]),2)
    print('mixed_endpoint_h_pieces_ok')
except Exception as exc:
    print('mixed_endpoint_h_pieces_error',type(exc).__name__,str(exc))

from solver8 import Resolvent
r=Resolvent(8)
try:
    print('discarded_mode_from_residue',r.mode_from_residue(-11.74,0.5))
except Exception as exc:
    print('discarded_mode_from_residue_error',type(exc).__name__,str(exc))

from solver2 import assemble
x,V,Vinv,W,D1,D2,D3,A,B=assemble(12)
print('solver2_endpoint_coordinates',x[0],x[-1])
print('solver2_replaced_bc_B_row_norm',np.linalg.norm(B[-2]))
