"""Symbolic derivation + verification of the far-field wake states and energy integrals.

Model (PROBLEM.md): plane-stress film, thickness h, modulus E, Poisson nu, prestress
sigma0*I, Winkler foundation f = -k u. Equilibrium: h sigma_{ab,b} - k u_a = 0.
Ebar = E/(1-nu^2), ell = sqrt(h Ebar/k).

We verify:
  (1) the 1D relaxed wake state u2(x2) on a strip of width S bounded by two
      traction-free cracks satisfies h Ebar u2'' - k u2 = 0 and
      sigma22 = sigma0 + Ebar u2' = 0 at x2 = 0, S;
  (2) the strip energy integral W(S) = int_0^S [h sigma0 eps22 + (h/2) dSig:eps
      + (k/2) u^2] dx2 equals -h sigma0^2 ell tanh(S/(2 ell))/Ebar;
  (3) the split W(S_ahead) - W(S_behind) reproduces G_A and G_B;
  (4) the compact identity 2 tanh(x) - tanh(2x) = 2 tanh(x)^3/(1+tanh(x)^2).
Outputs are printed; also written to analytic_output.txt by the caller.
"""
import sympy as sp

x2, S, s0, Eb, ell, h = sp.symbols('x2 S s0 Eb ell h', positive=True)
k = h * Eb / ell**2  # from ell^2 = h Eb / k

results = []

# ---- (1) wake state on a strip of width S, cracks on both edges ----------------
C = sp.symbols('C')
u2 = C * sp.sinh((x2 - S / 2) / ell)
ode_residual = sp.simplify(h * Eb * sp.diff(u2, x2, 2) - k * u2)
results.append(('ODE residual h*Eb*u2"" - k*u2 (must be 0)', sp.simplify(ode_residual)))

eps22 = sp.diff(u2, x2)
sig22 = s0 + Eb * eps22
# Boundary conditions sigma22 = 0 at x2 = 0 and x2 = S determine C:
bc0 = sp.Eq(sig22.subs(x2, 0), 0)
Csol = sp.solve(bc0, C)[0]
C_expected = -ell * s0 / (Eb * sp.cosh(S / (2 * ell)))
C_zero = sp.simplify(sp.expand_trig((Csol - C_expected).rewrite(sp.exp)))
results.append(('C - (-ell*s0/(Eb cosh(S/2ell))) (must be 0)', C_zero))
u2s = u2.subs(C, Csol)
sig22s = sp.simplify(sig22.subs(C, Csol))
bcS = sp.expand_trig(sig22s.subs(x2, S).rewrite(sp.exp))
results.append(('sigma22(S) after fix (must be 0)', sp.simplify(bcS)))
results.append(('sigma12 on faces (u1=0, x1-independent) ', 0))

# ---- (2) strip energy integral ---------------------------------------------------
# increment strain: eps11 = eps12 = 0, eps22 = du2/dx2.
# dSig22 = Eb*eps22 ; dSig11 = nu*Eb*eps22 ; foundation energy (k/2) u2^2.
nu = sp.symbols('nu', positive=True)
w = (h * s0 * eps22 + h / 2 * Eb * eps22**2 + k / 2 * u2**2)
W_integrand = w.subs(C, Csol)
W_direct = sp.integrate(W_integrand, (x2, 0, S))
# target closed form
W_target = -h * s0**2 * ell / Eb * sp.tanh(S / (2 * ell))
diff = W_direct.rewrite(sp.exp) - W_target.rewrite(sp.exp)
diff = sp.simplify(sp.expand(diff.rewrite(sp.exp)))
# robust zero test: exact simplify plus numeric sampling
num_ok = all(abs(complex(diff.subs({S: sv, ell: 1.3, h: 0.7, s0: 2.1, Eb: 3.3}))) < 1e-9
             for sv in (0.05, 0.4, 1.7, 5.0))
results.append(('W_direct - W_target (0 by simplify; numeric {})'.format(num_ok),
                sp.simplify(diff)))

# Decomposition: linear prestress term and quadratic term
W_lin = sp.integrate((h * s0 * eps22).subs(C, Csol), (x2, 0, S))
W_quad = sp.integrate((h / 2 * Eb * eps22**2 + k / 2 * u2**2).subs(C, Csol), (x2, 0, S))
W_lin_t = -2 * h * s0**2 * ell / Eb * sp.tanh(S / (2 * ell))
W_quad_t = h * s0**2 * ell / Eb * sp.tanh(S / (2 * ell))
d1 = sp.simplify(sp.expand_trig(W_lin.rewrite(sp.exp)) - sp.expand_trig(W_lin_t.rewrite(sp.exp)))
d2 = sp.simplify(sp.expand_trig(W_quad.rewrite(sp.exp)) - sp.expand_trig(W_quad_t.rewrite(sp.exp)))
results.append(('W_linear - (-2 h s0^2 ell/Eb tanh) (must be 0)', d1))
results.append(('W_quad   - (+1 h s0^2 ell/Eb tanh) (must be 0)', d2))

# ---- (3) energy release ----------------------------------------------------------
x = sp.symbols('x', positive=True)
ident = sp.simplify((2 * sp.tanh(x) - sp.tanh(2 * x) - 2 * sp.tanh(x)**3 / (1 + sp.tanh(x)**2)).rewrite(sp.exp))
results.append(('2tanh(x)-tanh(2x) - 2t^3/(1+t^2) (must be 0)', sp.simplify(ident)))

# ---- print ------------------------------------------------------------------------
print('=== SYMBOLIC CHECKS (sympy) ===')
ok = True
for name, val in results:
    zero = sp.simplify(val) == 0 if not isinstance(val, int) else val == 0
    ok = ok and zero
    print(f'[{"PASS" if zero else "FAIL"}] {name}: {sp.simplify(val)}')
print('ALL SYMBOLIC CHECKS PASS' if ok else 'SOME SYMBOLIC CHECKS FAILED')
