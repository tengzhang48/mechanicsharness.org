import sympy as sp, numpy as np
import hrr_clean as H
t=sp.symbols('t')
def test(fexpr,name):
    f0=sp.lambdify(t,fexpr,'numpy')
    d=[sp.lambdify(t,sp.diff(fexpr,t,m),'numpy') for m in range(1,5)]
    print("---",name)
    for tt in [0.3,1.0,2.0]:
        Y=np.array([f0(tt),d[0](tt),d[1](tt),d[2](tt)])
        r=H.rhs(tt,Y,n=1.0,p=0.5,k=1.5)
        print(f"   t={tt}: d4={d[3](tt):.5f} rhs={r[3]:.5f} diff={d[3](tt)-r[3]:.2e}")
test(sp.cos(t/2),"cos(t/2)")
test(sp.cos(sp.Rational(3,2)*t),"cos(3t/2)")
test(sp.cos(sp.Rational(3,2)*t)+sp.Rational(1,3)*sp.cos(t/2),"cos3/2+1/3cos1/2")
test(sp.cos(sp.Rational(3,2)*t)+3*sp.cos(t/2),"cos3/2+3cos1/2")