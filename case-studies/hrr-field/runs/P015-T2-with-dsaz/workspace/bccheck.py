import numpy as np
from scipy.integrate import solve_ivp
from notch import make_rhs, quadratic_coeff

n=4.0; alpha=2*np.pi/3
s=1.843824739714; pf=-1.263078996343

# integrate dense to alpha, get f,f',f'',f''' at alpha and at 0
sol=solve_ivp(lambda t,y:[y[1],y[2],y[3],make_rhs(n,s)(y[0],y[1],y[2],y[3])],
              [0.0,alpha],[1.0,0.0,pf,0.0],rtol=1e-13,atol=1e-15,max_step=1e-3,dense_output=True)
fa=sol.y[0,-1]; fpa=sol.y[1,-1]
# BC at notch faces
print('f(alpha)         = %.3e   (target 0)'%fa)
print("f'(alpha)        = %.3e   (target 0)"%fpa)
print('sigma_tt(alpha)  = %.3e   (= K r^(s-2) s(s-1) f(alpha))'%(s*(s-1)*fa))
print('sigma_rt(alpha)  = %.3e   (= K r^(s-2) (1-s) f(alpha))'%((1-s)*fpa))

# effective stress at theta=0 should be nonzero (symmetric opening, amplitude>0)
M=quadratic_coeff(n,s)
z0=np.array([1.0,0.0,pf])
Q0=z0@M@z0
print('se_tilde(0)      = %.6f'%np.sqrt(Q0))

# stresses at theta=0 via Airy
sth0=s*(s-1)*1.0
srr0=s*1.0+pf
print('sigma_tt(0)      = %.10f'%sth0)
print('sigma_rr(0)      = %.10f'%srr0)
print('R                = %.12f'%(srr0/sth0))
# also the cartesian ratio at theta=0: sig_yy/sig_xx == sig_tt/sig_rr at theta=0
print('1/R              = %.12f'%(sth0/srr0))
