import sympy as sp

# ======================================================================
# Wiener-Hopf formulation, Case A (matched asymptotic / strip).
# ======================================================================
# Work with mode I dominant. Strip x2 in [0,H], cracks at x2=0,x2=H for
# x1<0. Use the mode-I scalar reduction? For a periodic crack array under
# remote tension normal to crack plane, the exact reduction is to a single
# function. Let us instead use the "cosine" eigenfunctions in x2 that
# vanish-traction at crack faces (x2=0,H) for x1<0 and satisfy symmetry
# for x1>0.
#
# Take u1=0 and w=u2 (mode I). Then:
#   sig22 = Ebar w_,2 (wait eps22=w,2 only if u1 const)
# Correction: u1 not zero. But for mode I symmetric about each crack
# PLANE? The cracks are at x2=0,H; the mid-plane is x2=H/2, a symmetry
# plane (w even about H/2, u1 odd). So consider x2 in [0,H/2].
# For mode I (opening), w(x2) even around H/2 => w' =0 at H/2.
#
# Set u1=0 and w=w(x1,x2): eps22=w,2 ; eps11 = 0? No, even about H/2
# means u1=0 on that plane but not everywhere.
# The problem DOES couple. Use full vector Wiener-Hopf.
print("Do full vector symbolic WH in a separate script")