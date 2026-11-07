import numpy as np
import scipy.sparse as spa
import scipy.sparse.linalg as spla

# ======================================================================
# CASE A: finite-difference solution of the film-foundation strip problem.
# ======================================================================
# Domain: x1 in [-Lf, La] (Lf large behind, La ahead), x2 in [-H/2, H/2].
# Period H in x2 (periodic BC). Cracks on x2 = ±H/2 for x1 < 0.
# Crack faces traction free: sigma_22 = 0, sigma_12 = 0 on those faces.
# (Note crack at x2=+H/2 and -H/2 are the SAME line by periodicity:
#  identify x2=H/2 with -H/2. So we have period H with a crack at the
#  seam x2=0 if we shift.)  Simplify: x2 in [0,H], crack at x2=0 and H,
#  periodic. Symmetric about x2=H/2 and about x2=0? Crack at 0 disconnects.
# Use x2 in [0,H], crack at x2=0 (identified with x2=H).
# ----------------------------------------------------------------------
# Governing (plane stress): h*div(sigma) - k u = 0,
#   sigma = s0 I + Ebar[(1-nu)eps + nu tr(eps) I].
# Far ahead (x1->+inf): u->0.
# Far behind (x1->-inf): relaxed.
#
# Weak form for FE is easier. Use linear triangles on a structured mesh.
# We'll build with vector P1 elements, handle crack faces as free surfaces.
# ----------------------------------------------------------------------
# To keep it robust, use a virtual-work FEM with the periodic strip and
# a very long x1 domain; extract G via the J-integral or via K from the
# near-tip opening. Let's implement.
print("FEM")