"""Verify that p = r^(-1/2) cos(theta/2) erf(r/(2 sqrt(c t))) solves dp/dt = c lap p,
vanishes on the faces theta = +-pi, tends to r^(-1/2) cos(theta/2) as t->0+, and is regular
(~ r^(1/2)) at the tip; and evaluate the tip-region volume integral used in Sec. 4.2."""
import sympy as sp
r, t, c, th = sp.symbols('r t c theta', positive=True)
p = r**sp.Rational(-1, 2)*sp.cos(th/2)*sp.erf(r/(2*sp.sqrt(c*t)))
lap = sp.diff(p, r, 2) + sp.diff(p, r)/r + sp.diff(p, th, 2)/r**2
print('PDE residual:', sp.simplify(sp.diff(p, t) - c*lap))
print('face value:', sp.simplify(p.subs(th, sp.pi)))
print('near tip:', sp.series(p.subs(th, 0), r, 0, 2))
# int_0^inf erfc(r/(2 sqrt(ct))) dr and int cos^2(th/2) dth
print(sp.integrate(sp.erfc(r/(2*sp.sqrt(c*t))), (r, 0, sp.oo)), sp.integrate(sp.cos(th/2)**2, (th, -sp.pi, sp.pi)))
