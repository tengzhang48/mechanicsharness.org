import numpy as np
import sympy as sp

# =================== Coordinate-transform formulation ===================
# x = eta^p,  eta in [0,1].  Unknown v(eta) = u(x(eta)).
# We keep the operator in the form  u_tau = D2_x (1 - M s)  with s = 1 - M D1_x u.
# In eta coordinates:
#   d/dx = (1/x'(eta)) d/deta
#   d^2/dx^2 = (1/x'^2) d^2/deta^2 - (x''/x'^3) d/deta
# and M is applied in x space (its columns correspond to x-nodes).
# We therefore discretize the *differential* parts in eta (uniform grid) and map
# node values back and forth: this resolves the x=0 layer via clustering.

p, eta = sp.symbols('p eta', positive=True)
xof = eta ** p
xp = sp.diff(xof, eta)
xpp = sp.diff(xof, eta, 2)
print("x' =", xp)
print("x'' =", xpp)

# operator on u in eta:
# u_tau = D2_x s ; s = 1 - M D1_x u
# We build:  d/dx on eta-grid, then M in x-space weighted appropriately.
# Since M acts on values of w=d_x u at x-nodes, we can compute w at x-nodes via
# d_x = (1/x') d_eta applied to u-at-eta-nodes. This is a nodal operation in eta.
print("transform ready")