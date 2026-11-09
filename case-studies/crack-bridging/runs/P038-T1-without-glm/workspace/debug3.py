import numpy as np
from scipy.integrate import quad
import solve as S

Ep=S.Ep; a=S.a; t=S.t
c = S.c_s
SD = (2*t/np.pi)*S.theta(c)
xa = a/c

def Ffun(xi):
    return -(4.0*c/Ep)*(SD - t*(np.abs(xi)>xa))

# TRUE Chebyshev-U coefficients via quadrature (weight sqrt(1-xi^2))
def btrue(n):
    val = quad(lambda x: Ffun(x)*np.sin(n*np.arccos(x))*np.sqrt(1-x*x)/(np.sin(np.arccos(x))), -1,1, points=[-xa,xa], limit=400)
    # U_{n-1}(x) = sin(n theta)/sin(theta); use direct eval instead:
    return None

def Un(x, n):  # U_{n-1}
    th = np.arccos(x)
    return np.sin(n*th)/np.sin(th)

for n in [1,2,3,4,5,10,20,50,100,200,400]:
    val = quad(lambda x: Ffun(x)*Un(x,n)*np.sqrt(1-x*x), -1, 1, points=[-xa,xa], limit=400)[0]
    print(n, (2/np.pi)*val)
