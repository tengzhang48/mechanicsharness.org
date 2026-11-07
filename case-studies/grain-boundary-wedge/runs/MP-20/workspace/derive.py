"""Symbolic derivation / verification of closed-form principal-value integrals
for the grain-boundary kernel.

We need, for a collocation point a in (0,1) and a polynomial weight P(x):

  I1(a) = p.v. int_0^1 P(x)/(a-x)^3 dx
  I2(a) =       int_0^1 P(x)/(a+x)^3 dx
  I3(a) =       int_0^1 P(x)/(a+x)^4 dx
  I4(a) =       int_0^1 P(x)/(a+x)^5 dx
  I5(a) =       int_0^1 P(x)/(a+x)^6 dx

because  S_ss(a,x) = 2/(a-x)^3 - 2/(a+x)^3 + 12x(a-x)/(a+x)^5 - 12x/(a+x)^4
                   = 2/(a-x)^3 - 2/(a+x)^3 - 12x/(a+x)^4 + 12x/(a+x)^5 - 24x^2/(a+x)^6

(the last step expands 12x(a-x)/(a+x)^5 = 12x/(a+x)^4 - 12x^2/(a+x)^5 - ... let
sympy do the algebra so no sign slips through).

Closed forms for a monomial P(x) = x^n:
  int x^n/(a+x)^k dx                : elementary
  p.v. int x^n/(a-x)^3 dx           : logarithmic in (a-1)/(a-0) plus rational
"""

import sympy as sp


def pv_int_xn_over_amx3(n):
    """p.v. int_0^1 x^n / (a-x)^3 dx  as a function of the symbol a."""
    a, x = sp.symbols("a x", positive=True)
    expr = x ** n / (a - x) ** 3
    # use the antiderivative and take the finite part (drop the 1/(a-x)^2 pole)
    F = sp.integrate(expr, x)
    F = sp.simplify(F)
    return sp.simplify(F.subs(x, 1) - F.subs(x, 0))


def int_xn_over_apx_k(n, k):
    """int_0^1 x^n / (a+x)^k dx."""
    a, x = sp.symbols("a x", positive=True)
    F = sp.integrate(x ** n / (a + x) ** k, x)
    return sp.simplify(F.subs(x, 1) - F.subs(x, 0))


if __name__ == "__main__":
    a, x = sp.symbols("a x", positive=True)
    print("pv int x^n/(a-x)^3 for n=0..4:")
    for n in range(5):
        print(" n=%d :" % n, sp.simplify(pv_int_xn_over_amx3(n)))
    print()
    print("int x^n/(a+x)^k:")
    for k in (3, 4, 5, 6):
        for n in range(k + 2):
            print("  k=%d n=%d :" % (k, n), int_xn_over_apx_k(n, k))
    print()
    print("check S_ss expansion:")
    S = (2 / (a - x) ** 3 - 2 / (a + x) ** 3
         + 12 * x * (a - x) / (a + x) ** 5 - 12 * x / (a + x) ** 4)
    print(sp.simplify(sp.together(sp.expand(S))))
    target = (2 / (a - x) ** 3 - 2 / (a + x) ** 3 - 12 * x / (a + x) ** 4
              + 12 * x / (a + x) ** 5 - 24 * x ** 2 / (a + x) ** 6)
    print("difference:", sp.simplify(S - target))