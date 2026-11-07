"""Compute the equilibrium wedge w_inf and the stress intensity factor at z=h.

At t -> infinity the traction vanishes identically,
    sigma_gb(z, inf) = 0,
and the wedge w_inf solves   A1 w_inf = 1   (sigma0 normalised to 1),
i.e. the wedge whose elastic field cancels the imposed tension everywhere.

The stress intensity factor at the blocked interface z = h is the coefficient
of the square-root singular part of the opening *derivative* there.  For a
wedge whose flanks are stress free and whose tip is wedged against a rigid
blocked interface the opening vanishes at the tip with a square-root profile

    u(h - eps) ~ (K/E*) sqrt(eps) / sqrt(2 pi) * const ,

so K can be read from a fit of u against sqrt(h-z) near z=h.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
import gb3

xq, wq = leggauss(2000); xq = 0.5 * (xq + 1); wq = 0.5 * wq

# Solve A1 w = 1 on interior Gauss nodes, gauge-fixed to zero mean.
n = 30
sg, wg = leggauss(n)
si = 0.5 * (sg + 1.0)
E = gb3.elasticity_matrix(si, xq, wq)
Aug = np.vstack([E, np.ones(n)[None, :]])
rhs = np.concatenate([np.ones(n), [0.0]])
w, *_ = np.linalg.lstsq(Aug, rhs, rcond=None)

# The physical opening u is proportional to w; fix the constant of
# proportionality by the source: the applied tension sigma0 == 1 is exactly
# cancelled, so no further normalisation is needed for the SHAPE.
# Slope near the tip (si -> 1): fit w(si) = w1 + c*sqrt(1-si) + ...
sel = si > 0.85
X = np.sqrt(1.0 - si[sel])
A = np.vstack([np.ones(len(X)), X]).T
coef, *_ = np.linalg.lstsq(A, w[sel], rcond=None)
print('tip behaviour  w ~ %.5f + %.5f*sqrt(1-s)' % (coef[0], coef[1]))
# stress intensity factor: K/(sigma0 sqrt(pi h)) = -(E*/sigma0) * d/dz[u(h-z)] ...
# For the standard wedge/crack relation the dimensionless intensity is
#   K_norm = -c * sqrt(2)   (from u ~ c sqrt(1-s) = c sqrt((h-z)/h) = (c/sqrt h) sqrt(h-z))
c = coef[1]
print('c (opening sqrt-coefficient) = %.5f' % c)
for conv, lab in ((1.0, 'K = c'), (np.sqrt(2), 'K = c*sqrt2'), (np.sqrt(2 * np.pi), 'K = c*sqrt(2pi)')):
    print('  %-16s -> K_norm = %+.5f' % (lab, conv * c))
print('\nequilibrium opening (E*-scaled) at a few depths:')
for k in range(0, n, 5):
    print('  s=%.4f  w=%+.5f' % (si[k], w[k]))