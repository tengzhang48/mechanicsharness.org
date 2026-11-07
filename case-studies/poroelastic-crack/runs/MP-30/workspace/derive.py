#!/usr/bin/env python3
"""
Derivation checks for the poroelastic suddenly-loaded crack.

Constitutive relations (given):
  2G e_ij = sigma_ij - nu sigma_kk d_ij + alpha (1-2nu) d_ij p          (1)
  2G zeta = alpha (1-2nu) sigma_kk + alpha^2 (1-2nu)^2 p /(nu_u - nu)   (2)
  sigma_ij,j = 0                                                        (3)
  q_i = -kappa p,i                                                      (4)
  dzeta/dt + q_i,i = 0                                                  (5)

Invert (1):
  sigma_ij = 2G e_ij + { 2 G nu/(1-2nu) e_kk - alpha (1-2nu)... } p ...

We use the standard Biot identifications.  Let us derive several useful
combinations symbolically (by hand) and verify them numerically against the
given relations.

KEY IDENTITIES to verify numerically:
  (I1) undrained limit: setting zeta=0 in (2) and substituting sigma_kk from
       trace of (1) must reproduce Hooke's law with nu_u and same G.
  (I2) drained limit (p=0): Hooke's law with nu.
  (I3) Biot modulus / M and Skempton B from the given constants.
"""
import sympy as sp


def main():
    nu, nuu, alpha, G = sp.symbols('nu nu_u alpha G', positive=True)

    # From trace of (1):  2G e_kk = (1-2nu) sigma_kk + 2 alpha (1-2nu) p
    # => sigma_kk = (2G e_kk + ... ) -- instead invert (1) directly.
    # Solve (1) for sigma_ij:
    #   2G e_ij = sigma_ij - nu sigma_kk d_ij + alpha(1-2nu) p d_ij
    # Let s = sigma_kk. Trace: 2G e = (1-2nu) s + 2 alpha(1-2nu) p   [e=e_kk]
    #  => s = (2G e - 2 alpha (1-2nu) p)/(1-2nu)
    e = sp.symbols('e')
    pp = sp.symbols('pp')
    s = sp.simplify((2*G*e - 2*alpha*(1-2*nu)*pp)/(1-2*nu))
    print("sigma_kk =", s)

    # Undrained: zeta = 0 -> solve (2) for p:
    #   0 = alpha(1-2nu) s + alpha^2 (1-2nu)^2 p/(nu_u-nu)
    #   => p = -(nu_u-nu)/(alpha(1-2nu)) * s
    p_und = sp.simplify(-(nuu-nu)/(alpha*(1-2*nu))*s)
    print("p (undrained) =", sp.simplify(p_und))

    # Now plug p_und back into (1) to get effective undrained stress-strain.
    # 2G e_ij = sigma_ij - nu s d_ij + alpha(1-2nu) p_und d_ij
    # Compare with Hooke: 2G e_ij = sigma_ij - nu_u sigma_kk d_ij
    #   -> require: -(nu_u) s = -nu s + alpha(1-2nu) p_und
    rhs = -nu*s + alpha*(1-2*nu)*p_und
    print("coeff of s in (1) with p_und =", sp.simplify(rhs/s), " (should equal -nu_u =", -nuu, ")")
    return


def surface_pore_pressure():
    """Loading 1, t=0+ : induced pore pressure just next to the crack faces.

    Just outside the crack face the solid sees a traction-free boundary with
    applied normal traction sigma_yy = -1 (loading 1 drives sigma_n=-H) and,
    for a smooth crack face, sigma_xy = 0.  Plane strain (e_zz = 0).
    Undrained at t=0+ (zeta = 0).

    Unknowns: sigma_xx, p   (sigma_yy = -1, sigma_xy = 0, e_zz = 0, zeta = 0)
    Equations:
      2G e_zz = sigma_zz - nu sigma_kk + alpha(1-2nu) p       with e_zz=0
      2G zeta = alpha(1-2nu) sigma_kk + alpha^2(1-2nu)^2 p/(nu_u-nu)  with zeta=0
    where sigma_kk = sigma_xx + sigma_yy.
    """
    import sympy as sp
    G, nu, nuu, alpha, sxx, p = sp.symbols('G nu nu_u alpha sxx p', positive=True)
    syy = sp.Integer(-1)
    _ = sp  # sigma_zz = nu sigma_kk... in plane strain sigma_zz != 0 but e_zz=0
    # Need sigma_zz: from the third constitutive eq with e_zz=0 we get
    # 0 = sigma_zz - nu*(sxx+syy+sigma_zz) + alpha(1-2nu)p
    # => sigma_zz*(1-nu) = nu*(sxx+syy) - alpha(1-2nu)p
    szz = (nu * (sxx + syy) - alpha * (1 - 2 * nu) * p) / (1 - nu)
    skk = sxx + syy + szz
    eq1 = sp.Eq(0, alpha * (1 - 2 * nu) * skk + alpha**2 * (1 - 2 * nu)**2 * p / (nuu - nu))
    # The remaining equation: the imposed field is the mode-I plane-strain field
    # with sigma_yy = -1 at the surface.  For the *uniform* (zero-frequency)
    # mode-I field the in-plane stress at the surface of a pressurized crack
    # satisfies sigma_xx = -sigma_yy?  We use the plane-strain constraint that
    # the far field is uniaxial strain in x?  Instead, we note that just outside
    # the crack face along the extension the state is "uniaxial stress" normal
    # to the face for the *edge* of a half-space: for a half-space with a
    # uniformly loaded strip (crack face) the horizontal stress is not imposed.
    # The physically correct BVP: this is a *uniform* loading over the whole
    # line? No -- over the crack only.  Near the surface, however, by the
    # boundary-layer argument (uniform pressure over the crack, the relevant
    # local state at a surface point is that of a half-space under uniform
    # normal traction -1): for a half-space with traction -1 applied on the
    # free surface, sigma_xx = -1 as well (hydrostatic in-plane), and e_zz=0
    # (plane strain) => a well defined state.  Take sxx = -1.
    sol = sp.solve(eq1.subs({sxx: sp.Integer(-1)}), p)
    return sp.simplify(sol[0])


def surface_pore_pressure_correct():
    """Loading 1, t=0+ : pore pressure next to the crack face.

    Exact mode-I solution gives, on the face (see mode1_stress.py):
        sigma_yy = sigma_xx = -P   (P = applied face traction = 1)
        sigma_xy = 0
    Plane strain: e_zz = 0.  Undrained at t=0+: zeta = 0.

    Constitutive (given):
      2G e_ij = sigma_ij - nu sigma_kk d_ij + alpha(1-2nu) d_ij p
    Third equation (e_zz = 0), delta_zz=1:
      0 = sigma_zz - nu sigma_kk + alpha(1-2nu) p        (*)
    Trace: 2G e_kk = (1-2nu) sigma_kk + 2 alpha (1-2nu) p
    zeta eq:
      0 = alpha(1-2nu) sigma_kk + alpha^2 (1-2nu)^2 p/(nu_u-nu)   (**)

    Unknowns sigma_zz, p ;  sigma_xx = sigma_yy = -1.
    """
    import sympy as sp
    nu, nuu, alpha, szz, p = sp.symbols('nu nu_u alpha szz p', positive=True)
    skk = sp.Integer(-2) + szz            # sigma_xx + sigma_yy + sigma_zz
    eq_plane = sp.Eq(0, szz - nu * skk + alpha * (1 - 2 * nu) * p)
    eq_zeta = sp.Eq(0, alpha * (1 - 2 * nu) * skk
                    + alpha**2 * (1 - 2 * nu)**2 * p / (nuu - nu))
    sol = sp.solve([eq_plane, eq_zeta], [szz, p], dict=True)[0]
    return sp.simplify(sol[szz]), sp.simplify(sol[p])


if __name__ == "__main__":
    import sympy as sp
    pv = surface_pore_pressure()
    print("p_faces_0 (old, wrong sigzz) =", pv)
    szz, pc = surface_pore_pressure_correct()
    print("sigma_zz (undrained, plane strain) =", sp.simplify(szz))
    print("p_faces_0 (correct) =", sp.simplify(pc))
    nu, nuu, alpha = sp.symbols('nu nu_u alpha')
    val = pc.subs({nu: sp.Rational(1, 5), nuu: sp.Rational(2, 5), alpha: sp.Rational(89, 100)})
    print("p_faces_0 at (.2,.4,.89) =", float(val))