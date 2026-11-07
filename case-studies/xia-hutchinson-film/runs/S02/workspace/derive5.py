import sympy as sp

# ======================================================================
# Energy release, Case A, with ell = sqrt(h Ebar/k), Ebar=E/(1-nu^2).
# ======================================================================
# Far ahead (x1 -> +inf): intact film, uniform equibiaxial. Pick gauge
# u=0. sigma = s0 I everywhere (no increment). DP=0.
#
# Far behind (x1 -> -inf): wake fully split. Traction-free faces.
# The film relaxes by u1 = U e^{x1/ell}. Energies per unit x2-length,
# per half-period (one crack line), region of length dx1 moving with tip.
#
# We use the local (near-tip) representation: just behind the tip the
# release is a 1D uniaxial relaxation over the strip of half-width... 
# Better: use the mode-I local near tip. Let's derive G via energy of
# a "unit advance" of all tips by da, comparing far fields.
#
# Far behind state: relaxed. Compute its energy relative to reference in
# the region swept by the advance (a strip of width H in x2, length da).
# During advance by da, material points in [tip, tip+da] transition from
# intact/far-ahead state to relaxed wake. Energy release per tip:
#   dPi/tip = -[stored energy difference per period H per unit advance]
# The swept strip per period H is a rectangle da x H. But the transition
# is NOT local: the wake extends over length ~ell behind, so we must use
# the far-field energy densities and the "J-integral"/energy release =
# (far-ahead energy density - far-behind energy density) * H * da.
#
# Actually G*da per tip = change in potential energy per period H when
# the tip advances da, evaluated via far-field fluxes (J-integral).
# For mode I with the given bookkeeping, the correct statement:
#   G = (1/h) * d/dA (release). Per tip advance da, area advanced = h da
#   -> G = release per unit da divided by h.
# Release per unit da per period H = (W_ahead - W_behind)*H where W are
# stored energy per unit area.  So G = (W_ahead - W_behind)*H/h.
#
# W_ahead = 0 (increment field zero).
# W_behind = energy density of the relaxed wake relative to reference:
#   from DP = h s0 eps_aa + h/2 Dsig_ab eps_ab + k/2 u_a u_a
# relaxed wake: u1 = U e^{x1/ell}, sigma11 -> 0 means ...
#   sigma11 = s0 + Ebar u1' = 0  => u1' = -s0/Ebar.
#   u1 = -(s0/Ebar) x1 + const  (linear ramp, because traction-free
#   uniaxial requires the stress increment to cancel prestress).
# Wait: traction-free means sigma11=0 -> Ebar eps11 = -s0.
#   eps11 = u1' = -s0/Ebar. So u1 grows linearly; the foundation term
#   then requires decay to zero far away -> combine: u1'' = u1/ell^2.
#   General: u1 = -s0/Ebar + A e^{x1/ell}; as x1->-inf, u1 must stay
#   bounded? No: foundation limits it. Solve h Ebar u1'' - k u1=0 gives
#   u1 = C e^{x1/ell}. That gives eps11 = C/ell e^{x1/ell} -> 0, so
#   sigma11 -> s0, NOT zero. Contradiction with traction-free.
print("contradiction -> need care")