import numpy as np
from scipy.integrate import solve_ivp
n=3.5; k=1.0/n
def shoot(a_disp):
    a=a_disp; coef=((a-1.0)*k+1.0)*a
    def rhs(theta,y):
        f,fp=y; h2=a*a*f*f+fp*fp
        num=-coef*f-(k-1.0)*a*a*f*fp*fp/h2
        den=1.0+(k-1.0)*fp*fp/h2
        return (fp,num/den)
    s=solve_ivp(rhs,(0,np.pi),(0.0,1.0),rtol=1e-12,atol=1e-14)
    return s.y[1,-1]     # f'(pi); zero <=> traction-free face
for da in (0.0, 1e-3, 1e-2):
    print(f"a = 1/(n+1) + {da:+.0e}:  f'(pi) = {shoot(1.0/(n+1.0)+da):+.6e}")
# diagnostics at requested angles
