"""Full check at n=1: for phi=a cos(3t/2)+b cos(t/2), verify all three
constitutive relations (I),(II),(III) hold, i.e. does an Ur,Ut exist?
(I): p Ur = q(A-B)/2 -> defines Ur.
(III): Ur' + (p-1)Ut = 2qC -> defines Ut.
(II): Ur + Ut' = q(B-A)/2 -> must hold automatically.
We test (II) with exact symbolic-style numerics via sympy for reliability."""
import sympy as sp
t = sp.symbols('t')
p = sp.Rational(1,2); k = sp.Rational(3,2)
a,b = sp.symbols('a b')
phi = a*sp.cos(sp.Rational(3,2)*t)+b*sp.cos(sp.Rational(1,2)*t)
A = k*phi+sp.diff(phi,t,2)
B = k*(k-1)*phi
C = -(k-1)*sp.diff(phi,t)
q = sp.Integer(1)  # n=1
Ur = q*(A-B)/(2*p)
Ut = (2*q*C - sp.diff(Ur,t))/(p-1)
resid = sp.simplify(Ur + sp.diff(Ut,t) - q*(B-A)/2)
print("(II) residual for a cos3t/2+b cos t/2:", sp.expand(resid))

# Now test the actual ODE solution space: use the general symmetric Airy
# biharmonic solution at lambda: phi = cos((1+s)t)+... Let's instead directly
# test the ODE with phi=cos(t/2) and phi=cos(3t/2) individually using FINITE
# differences on the RHS via hrr_clean, at n=1.
import numpy as np
import hrr_clean as H
h=1e-4
def fd(f,t,m):
    if m==0: return f(t)
    if m==1: return (f(t+h)-f(t-h))/(2*h)
    if m==2: return (f(t+h)-2*f(t)+f(t-h))/h**2
    if m==3: return (f(t+2*h)-2*f(t+h)+2*f(t-h)-f(t-2*h))/(2*h**3)
    if m==4: return (f(t+2*h)-4*f(t+h)+6*f(t)-4*f(t-h)+f(t-2*h))/h**4
for f in [lambda t: np.cos(t/2), lambda t: np.cos(1.5*t)]:
    print("--- phi=", "cos(t/2)" if f(0.0)>0.99 else "cos(1.5t)")
    for t in [0.3,1.0,2.0]:
        Y=np.array([fd(f,t,0),fd(f,t,1),fd(f,t,2),fd(f,t,3)])
        r=H.rhs(t,Y,n=1.0,p=0.5,k=1.5)
        print(f"   t={t}: d4_fd={fd(f,t,4):.5f} rhs={r[3]:.5f} diff={fd(f,t,4)-r[3]:.2e}")