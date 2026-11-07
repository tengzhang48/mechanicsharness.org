#!/usr/bin/env python3
"""Check that MP-37's D equals the mechanical drained-stiffness factor
after the k_+ k_- scaling, and that off-diagonal product is 2 eta t^2/(m k).
"""
import numpy as np

eta = 1.0 / 3.0
xi = np.array([0.01, 0.1, 1.0, 3.0, 10.0, 100.0])
k = np.abs(xi)
m = np.hypot(1.0, xi)
t = k / (m + k)
D = 1 + eta - 2 * eta * k * t / m
dmech = 1 + eta * (m - k) / (m + k)
# After scaling S_s = 2(1-nu) S / k_-, W_s = k_+ W, N11 should be Dmech * k / (k_+ k_-)
# but k_+ k_- = |xi| = k, so N11_naive = dmech. MP-37 uses D instead.
# Alternative standard scaling uses k / m factors from mixed hydraulic/mechanical.
# det N = D - |N12|^2
m_plus = np.sqrt(1 - 1j * xi)
k_plus = np.sqrt(-1j * xi)
N12 = np.sqrt(2 * eta) * t / (m_plus * k_plus)
det = D - np.abs(N12) ** 2
print("D", D)
print("dmech", dmech)
print("|N12|^2", np.abs(N12) ** 2)
print("det N", det)
print("D(0) should be 1+eta", 1 + eta)
# det N at 0 is (1+eta)*1
# At inf D->1, N12->0, det->1
# Product of off-diagonals:
print("|N12|^2 vs 2 eta t^2 /(m k)", 2 * eta * t**2 / (m * k))

# Is D = dmech * m*(m+k)/(k*(something))?
# Direct: D = 1+eta - 2 eta k^2 /(m(m+k))
# dmech = 1 + eta (m-k)/(m+k) = [ (m+k) + eta(m-k) ]/(m+k)
# Try D / dmech
print("D/dmech", D / dmech)

# For energy identity, what matters is N(0)=diag(1+eta,1) and N(inf)=I
# and Hermitian plus/minus factorization. det N on the real line:
# If det has a known scalar factorization, n11^2+n21^2 = det N(0) / something.

# Scalar check: at xi=0, N=diag(1+eta,1), det=1+eta.
# Canonical factors -> N_-(0)^* N_-(0) = N(0) because N_+(0)=N_-(0)^* for
# Hermitian kernels with N(inf)=I (plus/minus at 0 related by conjugation
# through the real axis / infinity identification).
# Then n11^2 + n21^2 = N(0)_11 = 1+eta. This is THE energy identity.

print("energy identity is ||first column of N_-(0)||^2 = N11(0) = 1+eta")
print("paper heuristic would require n11=1, K=1/(1+eta)=0.75")
print("energy-equal no Darcy would require n21=0, n11=sqrt(1+eta), K=1/sqrt(1+eta)")
