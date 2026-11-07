"""Independent derivation of the plane-strain poroelastic surface symbol.

Laplace in t with parameter z>0 (or s), Fourier in x with wavenumber xi.
k = |xi|, m = sqrt(k^2 + z/c). Decaying as y -> +infty.

Airy-like representation used by MP-37 is re-derived from first principles:

  Phi = A0 e^{-k y} - (1+eta) H /(2k) y e^{-k y} - b P /(m^2-k^2) e^{-m y}
  p   = P e^{-m y} - beta H e^{-k y}

with H the surface amplitude of h = s + b p, and P = p0 + beta H.

We compute T = sigma_yy(y=0), Q = p_y(y=0), and opening from strain.
"""
import numpy as np

# Constitutive in plane strain (tension positive):
# 2G e_ij = sigma_ij - nu s delta + a p delta
# a = alpha(1-2nu)
# 2G zeta = a s + a^2 (1-2nu)/((nu_u-nu)) wait: PROBLEM: 2G zeta = a s + a^2 p / (nu_u-nu)
# a = alpha(1-2nu)
# c = 2 kappa G (1-nu)(nu_u-nu) / [a^2 (1-nu_u)]
# eta_paper = a / (2(1-nu))
# MP-37: b = a/(1-nu), eta = d/(1-nu_u), d=nu_u-nu, D_f = a^2 (1-nu_u)/(d(1-nu)), beta=a/D_f, b*beta=eta

nu, nu_u, alpha = 0.2, 0.4, 0.89
G = 1.0
a = alpha * (1 - 2 * nu)
d = nu_u - nu
b = a / (1 - nu)
eta_37 = d / (1 - nu_u)
Df = a**2 * (1 - nu_u) / (d * (1 - nu))
beta = a / Df
print("a,b,eta_37,Df,beta,b*beta", a, b, eta_37, Df, beta, b * beta)
print("paper eta", a / (2 * (1 - nu)))
print("1+eta_37", 1 + eta_37, "(1-nu)/(1-nu_u)", (1 - nu) / (1 - nu_u))

# Check MP-37 eta vs paper eta:
# paper eta = a/(2(1-nu)) = 0.33375
# MP-37 eta = d/(1-nu_u) = 0.2/0.6 = 1/3
# They differ! MP-37 used eta = d/(1-nu_u) which is NOT the paper's eta.
# For these numbers they are close: 1/3 vs 0.33375

print("ratio paper/37", (a / (2 * (1 - nu))) / eta_37)

# Coupling parameter that appears in drained vs undrained:
# G/(1-nu) vs G/(1-nu_u)
# (1-nu)/(1-nu_u) = 0.8/0.6 = 4/3 = 1 + eta_37
# paper's 1+eta_paper = 1.33375, not exactly 4/3.

# The mechanical stiffness at high k (drained, m~k):
# T = - G k / (2(1-nu)) * w   if p=0? wait opening w=2 u_y
# Elastic: sigma = - G/(2(1-nu)) |d/dx| * opening for Griffith
# Fourier: T_hat = - (G k /(2(1-nu))) w_hat for drained p=0? 

# Undrained: T_hat = - (G k /(2(1-nu_u))) w_hat
# ratio undrained/drained stiffness = (1-nu)/(1-nu_u) = 1+eta_37

# So MP-37's eta is the modulus ratio parameter, not paper's eta.
# Paper eta = alpha(1-2nu)/(2(1-nu))
# Relation: paper_eta = (nu_u-nu)/(1-nu_u) * ? let's see
# Skempton-like: actually
# (1-nu)/(1-nu_u) = 1 + (nu_u-nu)/(1-nu_u)
# The coupling in Biot that produces this is related to alpha.
# From PROBLEM, undrained from zeta=0:
# 2G zeta = a s + a^2 p / d = 0 => p = - d s / a
# 2G e_ij = sigma_ij - nu s delta + a p delta
#         = sigma_ij - nu s + a (-d/a) s = sigma_ij - (nu+d) s = sigma_ij - nu_u s
# Yes independent of alpha! The modulus ratio is independent of alpha.
# Paper's eta involves alpha. For K2_inf paper uses eta_paper.
# For K1_0 the modulus ratio is (1-nu_u)/(1-nu) independent of alpha.

print("K paper heuristic", (1-nu_u)/(1-nu))
