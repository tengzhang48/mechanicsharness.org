import numpy as np
from solve import solve_finite, a_half, ell, sigma0, delta_c, Gamma0
# normalized geometry: alpha=a/ell, c=eta*alpha, L=(eta-1)*alpha
an=a_half/ell; K0n=np.sqrt(Gamma0/(sigma0*delta_c))
eta=2.702337141; Ln=(eta-1)*an; c=eta*an
print(f'alpha={an:.6f} eta={eta:.6f} Ln={Ln:.6f} c={c:.6f} K0n={K0n:.6f}')
N=2000
# closure-only and remote-only unit solves
s0=solve_finite(c,Ln,sigma0_=1.0,remote=0.0,N=N,E_=1.0)
s1=solve_finite(c,Ln,sigma0_=0.0,remote=1.0,N=N,E_=1.0)
Sig=(K0n-s0['K_tip'])/s1['K_tip']
sol=solve_finite(c,Ln,sigma0_=1.0,remote=Sig,N=N,E_=1.0)
th_a=np.arccos((c-Ln)/c)
dt=sol['delta_at'](np.array([th_a]))[0]
print(f'spectral N={N}: Sigma={Sig*160:.5f} MPa  delta_t/dc={dt:.5f}')
# closed forms for reference
Ep=1.0
r=(Sig)-(2/np.pi)*np.arccos(1/eta)
dt_cf=4*np.sqrt(eta**2-1)*r+(8/np.pi)*np.log(eta)   # in units sigma0*a/E' == alpha*... 
# careful: formula gives delta_t*E'/(sigma0*a); a_norm=alpha => delta_t = alpha*[...]
dt_cf=an*(4*np.sqrt(eta**2-1)*r+(8/np.pi)*np.log(eta))
print(f'closed form (45): Sigma={Sig*160:.5f} MPa  delta_t/dc={dt_cf:.5f}')
