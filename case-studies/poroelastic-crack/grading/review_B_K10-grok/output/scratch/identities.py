#!/usr/bin/env python3
"""Independent identities for the tip layer and small-time volume.

1. Energy identity (configurational / J-integral with Darcy term)
2. Face-only volume coefficient (paper)
3. Weak-coupling expansions
4. Paper mesh vs diffusion length
"""
import math
import numpy as np

nu, nu_u, alpha = 0.2, 0.4, 0.89
d = nu_u - nu
a = alpha * (1 - 2 * nu)
eta_mod = d / (1 - nu_u)          # (1-nu)/(1-nu_u) - 1
eta_paper = a / (2 * (1 - nu))
r = (1 - nu_u) / (1 - nu)         # paper K1_0 / sqrt(pi L)

print("=== parameters ===")
print(f"eta_mod={eta_mod}  eta_paper={eta_paper}  r={r}  1+eta_mod={1+eta_mod}")
print(f"eta_mod vs eta_paper relative {(eta_mod-eta_paper)/eta_paper:.6e}")
print("NOTE: K1_0 depends only on nu, nu_u through eta_mod, not on alpha.")
print("Paper eta is the long-time loading-2 backstress coefficient.")

# Elastic energy release: G_elastic = (1-nu_*) K^2 / (2G)  in plane strain
# (standard Irwin: G = K^2 (1-nu^2)/E = K^2 (1-nu)/(2G) )
print("\n=== energy candidates ===")
# Outer undrained flux ~ (1-nu_u) K_outer^2 / (2G)
# Inner drained flux ~ (1-nu) K_tip^2 / (2G)
# If they were equal: K_tip/K_outer = sqrt( (1-nu_u)/(1-nu) ) = sqrt(r) ~ 0.866
# Paper sets K_tip/K_outer = r = 0.75 by substituting undrained COD into drained Irwin
# GPT energy: n11^2 + n21^2 = 1+eta, K_tip/K_outer = n11/(1+eta)
# so (1+eta)^2 (K_ratio)^2 + n21^2 = 1+eta
# (1-nu) K_tip^2 / (2G)  vs (1-nu_u) K_outer^2 / (2G)
# K_outer in units sqrt(pi L) is 1, K_tip = n11/(1+eta)
# inner/outer elastic G ratio = (1-nu)/(1-nu_u) * K_ratio^2 = (1+eta) * n11^2 / (1+eta)^2
# = n11^2 / (1+eta)
# remaining hydraulic share = 1 - n11^2/(1+eta) = n21^2/(1+eta)
# so hydraulic dissipation fraction = n21^2 / (1+eta)

n11 = 1.0767906703323102
n21 = 0.4169594526473146
print("n11^2+n21^2", n11**2+n21**2, "1+eta", 1+eta_mod, "err", n11**2+n21**2-(1+eta_mod))
Kgpt = n11/(1+eta_mod)
print("K_gpt", Kgpt)
print("sqrt(r) energy-equal", math.sqrt(r))
print("paper r", r)
print("G_inner/G_outer", (1-nu)/(1-nu_u) * Kgpt**2)
print("hydraulic fraction of outer G", 1 - (1-nu)/(1-nu_u)*Kgpt**2)

# Weak coupling: n11 = 1 + eta/4 + O(eta^2), n21=sqrt(eta/2)+O(eta^{3/2})
# K = n11/(1+eta) = (1+eta/4)/(1+eta) +... = 1 - 3eta/4 + O(eta^2)
print("\n=== weak coupling ===")
print("K weak", 1 - 0.75*eta_mod, "vs gpt", Kgpt, "vs paper", r, "vs 1-eta", 1-eta_mod)
print("paper r = 1/(1+eta) = 1 - eta + eta^2 - ...")
print("weak K = 1 - 0.75 eta + ...  so paper over-reduces at small eta")
print("energy-equal sqrt(1/(1+eta)) = 1 - 0.5 eta + ...")

# Volume face-only
print("\n=== volume coefficients ===")
# Paper eq 18: 16 L^2 (1-nu_u)(nu_u-nu) / [G sqrt(pi) (1-nu)] * sqrt(ct/L^2)
A_paper = 16 * (1-nu_u) * d / ((1-nu) * math.sqrt(math.pi))
# MP-37 face: 16 d / ((1+eta) sqrt(pi))
# 1+eta = (1-nu)/(1-nu_u) so 16 d / ((1+eta)sqrt(pi)) = 16 d (1-nu_u)/((1-nu)sqrt(pi))
A_face = 16 * d / ((1+eta_mod)*math.sqrt(math.pi))
print("A_paper", A_paper)
print("A_face MP37", A_face, "match", abs(A_paper-A_face)<1e-15)

# Reciprocity identity of paper:
# mode2 small-time volume for unit p = - (something)
# Paper: V2(t) ~ - 8 L^2 (1-nu_u)(1-2nu) alpha / [G something wait]
# From KEY: A = 16 (1-nu_u)(nu_u-nu)/[(1-nu)sqrt(pi)]
# Induced p = 2d/a, times mode-2 short-time volume coefficient.

# Tip contribution MP-37:
# A_tips = 4 d I / ((1+eta) sqrt(pi)), I = pi sqrt(2/eta) n21
A_tips = 4 * d * math.pi * math.sqrt(2/eta_mod) * n21 / ((1+eta_mod)*math.sqrt(math.pi))
A_tot = A_face + A_tips
print("A_tips", A_tips)
print("A_tot", A_tot)
print("A_tips/A_face", A_tips/A_face)
print("formula boxed", 4*d/((1+eta_mod)*math.sqrt(math.pi))*(4+math.pi*math.sqrt(2/eta_mod)*n21))

# Why tips contribute at same order:
# Face layer: zeta ~ O(1) over area ~ L * sqrt(ct)  => int zeta ~ L sqrt(ct)
# Tip layer: two tips, s^a ~ sqrt(L/r), zeta ~ O(1) over area ~ ct, but weighted
# int zeta s^a ~ sqrt(L) * int zeta /sqrt(r) dA over region size sqrt(ct)
# dA ~ ct, 1/sqrt(r) ~ (ct)^{-1/4}? wait size delta=sqrt(ct), int_0^delta (1/sqrt(r)) r dr ~ int r^{1/2} dr ~ delta^{3/2} = (ct)^{3/4} too small?
# 2D: dA = r dr dtheta, f~ r^{-1/2}, zeta ~ Z(R) with R=r/delta, dA_phys = delta^2 dA_nd
# int zeta s^a dA = const * sqrt(L) * delta * int Z R^{-1/2} dA_nd
# because s^a ~ sqrt(L/r) = sqrt(L/delta) R^{-1/2}, dA = delta^2 dA_nd
# product sqrt(L/delta)*delta^2 = sqrt(L) delta^{3/2}  Hmm that's (ct)^{3/4} not (ct)^{1/2}!
#
# Wait: Laplace domain. delta_z = sqrt(c/z), Laplace of sqrt(t) is ~ z^{-3/2}
# The identity is V-V0 = (d/a) int zeta s^a dA
# zeta Laplace ~ 1/z * Z(X) because step loading, length scale delta=sqrt(c/z)
# int zeta s^a ~ (1/z) * sqrt(L) * delta^{3/2}? still wrong order.
#
# s^a ~ sqrt(2L/r) which is O(sqrt(L/delta)) = O( (L^2 z/c)^{1/4} ) in the tip cell.
# zeta_hat ~ O(1/z) * (fluid content of order the undrained p * storativity)
# area ~ delta^2 = c/z
# int ~ (1/z)*(c/z)*sqrt(L)*sqrt(z/c) = (c/z^2) * sqrt(L z/c) wait:
# sqrt(2L/r)=sqrt(2L/delta)*R^{-1/2}
# int dA = delta^2
# => (1/z) * sqrt(L/delta) * delta^2 * I_nd = (1/z) sqrt(L) delta^{3/2} I_nd
# delta=sqrt(c/z), delta^{3/2}=(c/z)^{3/4}, times 1/z = z^{-1} z^{-3/4}=z^{-7/4}  WRONG
#
# Unless zeta itself is O(1/sqrt(r)) near the tip so extra 1/sqrt(delta):
# p ~ B_p sqrt(r) so bounded*sqrt, zeta related to Z=P+H
# H is harmonic, like the stress trace, H ~ 1/sqrt(r) (K-field)!
# Z = P + H, P is Helmholtz, P ~ B_p sqrt(r) near origin? They said P ~ B_p sqrt(R) cos(theta/2)
# H ~ far undrained ~ const / sqrt(r) * delta scaling.
# If Z ~ O(R^{-1/2}) then extra 1/sqrt(delta):
# int zeta s^a ~ (1/z)* (delta^{-1/2}) * sqrt(L/delta) * delta^2 = (1/z) sqrt(L) delta
# = (1/z) sqrt(L) sqrt(c/z) = sqrt(c L) z^{-3/2}
# inverse Laplace ~ sqrt(c t) * L^{1/2}  wait we need L * sqrt(ct) for dimensionless A * L^2/G * sqrt(ct/L^2) = (L/G) sqrt(ct)*L
# V ~ L^2/G * A * sqrt(ct/L^2) = (L/G) A sqrt(c t)
# int zeta s^a has units of zeta * stress * area. This matches sqrt(t) if Z ~ r^{-1/2}.
print("\nTip volume is O(sqrt(t)) because the inner fluid-content/stress-trace field Z is K-like, r^{-1/2},")
print("so int zeta s^a ~ delta * sqrt(L) / z  -> z^{-3/2} -> sqrt(t). Face layer is also z^{-3/2}.")

# Paper numerics
print("\n=== paper mesh vs tip layer ===")
# 6 quadratic elements over 2L, equal size => element length = 2L/6 = L/3
# smallest t* in Fig 9: the text mentions t* down to 0.0001 in Fig 7 axis
# Fig 6: t*<0.001 is where they say small-time volume is excellent
# Fig 9: "limiting values shown as dashed lines"
elem = 1.0/3  # in units of L
for tstar in (1e-4, 1e-3, 1e-2, 0.1):
    delta = math.sqrt(tstar)  # sqrt(ct)/L
    print(f"t*={tstar:g}  delta/L={delta:.5f}  delta/elem={delta/elem:.4f}  elem/delta={elem/delta:.1f}")

print("""
Stehfest with 6 terms: the inversion samples Laplace parameters
s_n = n ln(2)/t, n=1..6. Smallest physical time t*=1e-4 gives
diffusion length 0.01 L, while the smallest element is 0.33 L,
and the shape function phi is only constrained as a quadratic
over that element. The drained-tip layer is therefore not resolved
whenever t* << 0.1. Fig 9's apparent approach to 0.75 sqrt(pi L)
is the method returning the outer COD-based intensity
K = G sqrt(pi) phi(L) / (2(1-nu)), i.e. exactly the paper's
heuristic (undrained opening inserted into drained Irwin),
not a resolved inner K.
""")

# Rice-Simons high-speed limit
print("=== Rice-Simons correspondence ===")
print("""
A steadily moving semi-infinite crack at speed V has a diffusion length c/V.
As V -> inf the process zone shrinks. Rice & Simons 1976 / Simons 1977:
the tip is drained (r^{-1/2} singularity with drained moduli) while the
outer field is undrained. The ratio K_tip/K_outer is NOT the modulus
ratio (1-nu_u)/(1-nu); it is a function of a speed parameter that
approaches a nontrivial high-speed limit.
The sudden stationary crack at t=0+ is the same inner problem after
Laplace transform (diffusion length sqrt(c/z), z->inf), or equivalently
the high-speed moving-crack inner problem with the permeable-face
boundary condition of loading 1 (p=0 on the crack, no flux on the
ligament). The paper invoked Rice & Simons only to keep the r^{-1/2}
singularity and then replaced the whole COD by the undrained ellipse,
which is the wrong matching. That is why their Eq. 29 is a heuristic
and underestimates K(0+).
""")

# Face-only volume misses tips of same order, ratio A_tips/A_face ~ 0.8
print("A_tot/A_paper", A_tot/A_paper)
