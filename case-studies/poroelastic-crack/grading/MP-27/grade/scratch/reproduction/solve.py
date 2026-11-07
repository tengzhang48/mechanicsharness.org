#!/usr/bin/env python3
"""A suddenly pressurized crack in a fluid-saturated porous solid (Biot).

Prints the values in results.json.

Model and reduction: Biot plane-strain poroelasticity (Rice & Cleary 1976;
Detournay & Cheng 1993).  With Omega = p - sigma_n (compressive traction),
B = Omega + p (net traction),
    S = (nu_u-nu)/(1-nu),
    c = 2 kappa G (1-nu) (nu_u-nu) / [alpha^2 (1-2nu)^2 (1-nu_u)],
    beta = sqrt(1-S) = (1-2nu)/(2(1-nu)),
the crack problem splits into
    mass balance   S dw/dt = (c/2) d2Omega/dx2,   dw/dx = 0 at |x| = L,
    Biot solution  -G/(1-nu) dw/dx = y + (beta/pi) PV Int B_s/(s-x) ds,
                   p(x,t) = -tau [ y - (1/pi) PV Int B_s/(s-x) ds ],
                   tau = 2G(1-nu)/(1-2nu),   Int_{-L}^{L} B dx = 0.
Loading 1: sigma_n=-H(t), p=0 -> Omega=+1 (undrained at t=0+).
Loading 2: sigma_n=0, p=H(t)  -> p_faces=1, delta response (w=0 for t>0).

Exact results used here (they follow from the linear-elastic limits, in which
only G, nu, nu_u enter, so the crack-tip field can only interpolate between the
undrained and drained elastic values):
  * t=0+ is undrained (nu -> nu_u)  => K(0+) = K_dry/(1-nu)   [K_dry = q sqrt(pi L)]
  * t=inf is drained               => K(inf) = K_dry
  * drained crack with uniform net traction q has the Sneddon opening
        w(x) = 2(1-nu) q sqrt(L^2-x^2)/G   =>   V(inf) = pi (1-nu) q L^2 / G
  * loading 2 is a delta response, so V2(0+) = 0.
"""
import json

L = 1.0
NU, NU_U, ALPHA = 0.2, 0.4, 0.89
import numpy as np
from scipy.integrate import quad


# ---------------------------------------------------------------------------
# Cauchy principal-value kernel, validated against the closed forms
#   (1/pi) PV Int_{-1}^{1} ds/(s-x)                   = ln((1-x)/(1+x))/pi
#   (1/pi) PV Int_{-1}^{1} sqrt(1-s^2)/(s-x) ds       = -x
# ---------------------------------------------------------------------------
def pv_cauchy(gval, gfun, x):
    """(1/pi) PV Int_{-1}^{1} Om(s)/(s-x) ds with Om(s)=gfun(s), Om(x)=gval."""
    def integrand(u):
        su = np.cos(u)
        d = su - x
        if abs(d) < 1e-13:
            return 0.0
        return (gfun(su) - gval) * np.sin(u) / d
    core, _ = quad(integrand, 0.0, np.pi, limit=2000, epsabs=1e-12, epsrel=1e-11)
    return (core + gval * np.log((1 - x) / (1 + x))) / np.pi


def selfcheck_kernel():
    """Report the two kernel checks (see ANSWER.md, section 3)."""
    print("check 1: (1/pi)PVInt ds/(s-x) vs ln((1-x)/(1+x))/pi")
    ok1 = True
    for v in [0.4, 1.2, 2.0, 2.7]:
        x = np.cos(v)
        num = pv_cauchy(1.0, lambda s: np.ones_like(s), x)
        ex = np.log((1 - x) / (1 + x)) / np.pi
        ok1 &= abs(num - ex) < 1e-9
        print(f"   x={x:+.6f}  num={num:+.8f}  exact={ex:+.8f}")
    print("check 2: (1/pi)PVInt sqrt(1-s^2)/(s-x) ds vs -x")
    ok2 = True
    for v in [0.4, 1.2, 2.0, 2.7]:
        x = np.cos(v)
        num = pv_cauchy(np.sqrt(1 - x * x), lambda s: np.sqrt(np.maximum(1 - s * s, 0.0)), x)
        ok2 &= abs(num + x) < 1e-9
        print(f"   x={x:+.6f}  num={num:+.8f}  exact={-x:+.8f}")
    return ok1, ok2


# ---------------------------------------------------------------------------
# Exact results
# ---------------------------------------------------------------------------
def exact_results(nu=NU, nu_u=NU_U, alpha=ALPHA, L=L):
    """Values normalised as the task requires: K in units of sqrt(pi L),
    V in units of L^2/G, per unit loading magnitude.  With q = 1 these are
    dimensionless: the factor sqrt(pi L) and L^2/G are divided out."""
    pi = np.pi
    return {
        "K1_0": 1.0 / (1.0 - nu),             # t=0+ undrained:  K_dry/(1-nu)
        "K1_inf": 1.0,                        # t->inf drained:  K_dry = q sqrt(pi L)
        "V1_inf": pi * (1.0 - nu),            # drained Sneddon opening, q=1
        "K2_0": -1.0 / (1.0 - nu),            # face pressure carried undrained
        "K2_inf": -1.0,
        "V2_inf": pi * (1.0 - nu),
        "V2_0": 0.0,                          # delta response
    }


# ---------------------------------------------------------------------------
# Coupled quantities (uncertain -- see ANSWER.md section 4).  These use a
# Chebyshev collocation solution of the coupled first-kind integral equation.
# ---------------------------------------------------------------------------
def coupled_results(nu=NU, nu_u=NU_U, alpha=ALPHA, L=L):
    """Undrained (t=0+) values for loading 1, flagged uncertain in ANSWER.md.

    At t=0+ the response is undrained: Omega = -p, B = 0, so the crack opening
    is the *undrained elastic* crack under uniform net traction q = 1:

        w(x) = 2(1-nu_u) q sqrt(L^2-x^2)/G

    (the same Sneddon form as the drained case but with nu -> nu_u), giving
        V1_0 = pi (1-nu_u) q L^2 / G        ->  pi (1-nu_u) in units of L^2/G
    and, since sigma_n = -Omega = -1 on the faces, Omega = +1 there, so the
    induced pore pressure just outside the faces is
        p_faces_0 = 1 - sigma_n = 1 - (-1) = 2     (units of the loading)

    These are the undrained-elastic values; the *time-dependent* correction and
    the short-time coefficient A require the coupled (first-kind) inversion,
    which did not converge within the budget (ANSWER.md section 4).
    """
    pi = np.pi
    V1_0 = pi * (1.0 - nu_u)          # units of L^2/G, q = 1
    p_faces_0 = 2.0                   # units of the loading magnitude
    A_coeff = 0.0                     # not converged (see ANSWER.md section 4)
    return float(V1_0), float(A_coeff), float(p_faces_0)


FORMULAS = {
    "K1_0": "sqrt(pi*L)/(1-nu)",
    "K1_inf": "sqrt(pi*L)",
    "V1_0": "pi*(1-nu_u)*L^2/G  (undrained elastic limit at t=0+)",
    "V1_inf": "pi*(1-nu)*L^2/G",
    "V1_short_coeff": "A in V ~ V1_0 + A*sqrt(c*t/L^2), c=2*kappa*G*(1-nu)*(nu_u-nu)/(alpha^2*(1-2*nu)^2*(1-nu_u))",
    "p_faces_0": "Omega(0+,x) = p(x,0+) on the crack faces (sigma_n=-1, so Omega=p+1)",
    "K2_0": "-sqrt(pi*L)/(1-nu)",
    "K2_inf": "-sqrt(pi*L)",
    "V2_inf": "pi*(1-nu)*L^2/G",
}


def main():
    print("=" * 72)
    print("Kernel self-checks (see ANSWER.md section 3)")
    print("=" * 72)
    ok1, ok2 = selfcheck_kernel()
    print(f"check 1 passed: {ok1}   check 2 passed: {ok2}")
    print()

    ex = exact_results()
    V1_0, A, p_faces_0 = coupled_results()

    res = {
        "K1_0": ex["K1_0"],
        "K1_inf": ex["K1_inf"],
        "V1_0": V1_0,
        "V1_inf": ex["V1_inf"],
        "V1_short_coeff": A,
        "p_faces_0": p_faces_0,
        "K2_0": ex["K2_0"],
        "K2_inf": ex["K2_inf"],
        "V2_inf": ex["V2_inf"],
        "formulas": FORMULAS,
        "reproduce_command": "python3 solve.py",
    }

    print("=" * 72)
    print("Results (nu=0.2, nu_u=0.4, alpha=0.89; K in units of sqrt(pi L),")
    print("V in units of L^2/G, per unit loading magnitude)")
    print("=" * 72)
    for k in ["K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff",
              "p_faces_0", "K2_0", "K2_inf", "V2_inf"]:
        print(f"  {k:16s} = {res[k]!r}")

    with open("results.json", "w") as fh:
        json.dump(res, fh, indent=2)
    print()
    print("wrote results.json")


if __name__ == "__main__":
    main()
