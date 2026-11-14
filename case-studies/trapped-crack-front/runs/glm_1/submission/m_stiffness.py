"""
Fourier-domain crack-plane stiffness matrix M(p,q) for an isotropic elastic solid.

A displacement-discontinuity (DD) distribution B=(Bx,By,Bz) on plane y=0 with Fourier
factor exp(i(px+qz)) generates traction t = M(p,q) B on that plane.

Half-space solution (y>0):  u = (a + y*b) exp(-r y) exp(i(px+qz)),  r^2=p^2+q^2,
with b = c (i p, -r, i q)  (only form that solves Navier, incl. the Jordan branch).
Navier fixes  c = (lam+mu) D0 /(r (lam+3 mu)),  D0 = i p a_x + i q a_z - r a_y.
Sector BCs (from gluing the two half-spaces by the y-/xz-reflections):
  opening sector (B=(0,By,0)):  t_x = t_z = 0,  a_y = By/2
  shear sector   (B=(Bx,0,Bz)): t_y = 0,  a_x = Bx/2, a_z = Bz/2
t = sigma . e_y at y=0 from the analytic fields.

This module: solves those 3x3 complex systems exactly (numpy) and compares with the
closed forms
  M_yy = -mu r /(2(1-nu))
  M_ij = -mu/(2(1-nu)) [ (1-nu) r delta_ij + nu k_i k_j / r ]   (i,j = x,z; k=(p,q))
which the hand derivation gives.
"""
import numpy as np

NU = 0.3
MU = 1.0

def fields(a, c, p, q, r, x, y, z):
    """Return u and needed derivative products of u=(a+y b)e^{-ry}e^{i(px+qz)}, b=c(ip,-r,iq)."""
    b = c*np.array([1j*p, -r, 1j*q])
    E = np.exp(-r*y)*np.exp(1j*(p*x+q*z))
    u = (a + y*b)*E
    return u, b, E

def traction(a, c, p, q, r, mu=MU, nu=NU):
    """t = sigma.e_y on y=0 from the half-space solution. lam = lame."""
    lam = 2*mu*nu/(1-2*nu)
    b = c*np.array([1j*p, -r, 1j*q])
    D0 = 1j*p*a[0] + 1j*q*a[2] - r*a[1]
    Da = D0 + b[1]
    tx = mu*(b[0] - r*a[0] + 1j*p*a[1])
    ty = 2*mu*(b[1] - r*a[1]) + lam*Da
    tz = mu*(b[2] - r*a[2] + 1j*q*a[1])
    return np.array([tx, ty, tz])

def c_of_a(a, p, q, r, mu=MU, nu=NU):
    lam = 2*mu*nu/(1-2*nu)
    D0 = 1j*p*a[0] + 1j*q*a[2] - r*a[1]
    return (lam+mu)*D0/(r*(lam+3*mu))

def navier_residual(a, c, p, q, r, mu=MU, nu=NU, y=0.7, n=5):
    """Numerically verify Navier for u=(a+yb)e^{-ry}e^{i(px+qz)} via spectral check:
    plug analytic Laplacian and gradient-divergence (all analytic here)."""
    lam = 2*mu*nu/(1-2*nu)
    b = c*np.array([1j*p, -r, 1j*q])
    # Delta u = -2 r b e  (using r^2=p^2+q^2);  div u = (Da + y*Db) e, Db = i p b_x + i q b_z - r b_y = 0
    # grad(div u) = e * ( i p, d/dy, i q ) (Da + y Db) ;  d/dy(div u) = -r Da e (Db=0)
    Da = 1j*p*a[0] + 1j*q*a[2] - r*a[1] + b[1]
    lap = -2*r*b
    gdiv = np.array([1j*p*Da, -r*Da, 1j*q*Da])
    return mu*lap + (lam+mu)*gdiv

def M_numeric(p, q, mu=MU, nu=NU):
    """Exact numeric M(p,q) (3x3, complex) from solving the sector mixed problems."""
    r = np.hypot(p, q)
    M = np.zeros((3,3), complex)
    # opening sector: unknowns a; equations: t_x=0, t_z=0, a_y=1/2
    A = np.zeros((3,3), complex); rhs = np.zeros(3, complex)
    a = np.array([1,0,0], complex)
    # t_x is linear in (a, c(a)); write t = T(a) with c linear in a.
    # Build matrices: for basis vector e_j, c_j = c_of_a(e_j), traction column = traction(e_j, c_j)
    cols = []
    for j in range(3):
        ej = np.zeros(3, complex); ej[j] = 1.0
        cj = c_of_a(ej, p, q, r, mu, nu)
        cols.append(traction(ej, cj, p, q, r, mu, nu))
    Tm = np.array(cols).T   # Tm[:,j] = traction from a=e_j
    # opening: t_x=0 (row0), t_z=0 (row2), a_y=1/2 (row for a_y)
    A[0,:] = Tm[0,:]; rhs[0] = 0
    A[1,1] = 1.0;      rhs[1] = 0.5
    A[2,:] = Tm[2,:]; rhs[2] = 0
    a = np.linalg.solve(A, rhs)
    M[:,1] = Tm @ a     # traction; B = (0,1,0)
    # shear sector: t_y=0 (row1), a_x=1/2, a_z=1/2 (for each column)
    for j in (0,2):
        A2 = np.zeros((3,3), complex); rhs2 = np.zeros(3, complex)
        A2[1,:] = Tm[1,:]; rhs2[1] = 0
        A2[0,0] = 1.0; rhs2[0] = 0.0
        A2[2,2] = 1.0; rhs2[2] = 0.5 if j==2 else 0.0
        if j==0:
            A2[0,0] = 1.0; rhs2[0] = 0.5
            A2[2,2] = 1.0; rhs2[2] = 0.0
        a2 = np.linalg.solve(A2, rhs2)
        M[:,j] = Tm @ a2
    return M

def M_closed(p, q, mu=MU, nu=NU):
    r = np.hypot(p, q)
    M = np.zeros((3,3), complex)
    M[1,1] = -mu*r/(2*(1-nu))
    k = np.array([p, q])
    Min = -mu/(2*(1-nu))*((1-nu)*r*np.eye(2) + (nu/r)*np.outer(k,k))
    M[np.ix_([0,2],[0,2])] = Min
    return M

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    print("closed-form verification (should be ~1e-15):")
    worst = 0.0
    for _ in range(6):
        p = rng.uniform(-3,3); q = rng.uniform(0.2,3)
        Mn = M_numeric(p,q); Mc = M_closed(p,q)
        err = np.max(np.abs(Mn-Mc))/max(1.0,np.max(np.abs(Mc)))
        worst = max(worst, err)
        print(f"  p={p:+.3f} q={q:.3f}  rel.err={err:.2e}")
        # structure checks
        assert np.max(np.abs(Mn-Mn.T)) < 1e-12*max(1,np.max(np.abs(Mn))), "symmetry"
        assert abs(Mn[0,1])+abs(Mn[2,1]) < 1e-12, "sector decoupling"
    print("worst rel.err:", worst)
    # 2D limits
    print("\n2D limits q->0 (r=|p|), p=1.7:")
    p0=1.7
    M0 = M_closed(p0, 1e-9)
    print("  M_yy/(mu|p|) =", M0[1,1].real/-(MU*p0), " expect 1/(2(1-nu)) =", 1/(2*(1-NU)))
    print("  M_xx/(mu|p|) =", M0[0,0].real/-(MU*p0), " expect 1/(2(1-nu)) =", 1/(2*(1-NU)))
    print("  M_zz/(mu|p|) =", M0[2,2].real/-(MU*p0), " expect 1/2")
    # Navier residual of the ansatz (with c fixed by Navier itself)
    a = np.array([.3,-.2,.5])
    c = c_of_a(a, 1.1, 0.7, np.hypot(1.1,0.7))
    print("\nNavier residual check:", np.max(np.abs(navier_residual(a, c, 1.1, 0.7, np.hypot(1.1,0.7)))))
    # isotropy check: M(Rk) = R3 M(k) R3^T  with R3 a rotation about y
    th = 0.6; ct, st = np.cos(th), np.sin(th)
    R3 = np.array([[ct,0,st],[0,1,0],[-st,0,ct]])
    pk = ct*1.3 + st*0.8; qk = -st*1.3 + ct*0.8
    err = np.max(np.abs(M_closed(pk,qk) - R3@M_closed(1.3,0.8)@R3.T))
    print("isotropy residual:", err)
