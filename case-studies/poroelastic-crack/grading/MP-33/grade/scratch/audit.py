import math
import numpy as np

nu, nu_u, alpha, G = .2, .4, .89, 1.
p0 = 2*(nu_u-nu)/(alpha*(1-2*nu))
A = 8*(nu_u-nu)/((1-nu)*math.sqrt(math.pi))
conversion = 4*(1-2*nu)*alpha/(2*G*(1-nu))
print("Finite-difference integral audit: original right sum versus trapezoidal sum")
rows = []
for Ny in (1000, 2000, 4000, 8000):
    dy = 1/Ny
    dt = .4*dy**2
    nsteps = math.ceil(5e-4/dt)
    p = np.full(Ny+1, p0)
    p[0] = 0.
    for step in range(nsteps):
        p[1:-1] += .4*(p[2:]-2*p[1:-1]+p[:-2])
    t = nsteps*dt
    right = np.sum(p0-p[1:])*dy
    trap = right + .5*p0*dy
    original = conversion*right/math.sqrt(t)
    corrected = conversion*trap/math.sqrt(t)
    rows.append((Ny, original-A, corrected-A))
    print(f"Ny={Ny}: original_error={original-A:.12g}, trapezoid_error={corrected-A:.12g}, missing_endpoint_term={conversion*.5*p0*dy/math.sqrt(t):.12g}")
for label, col in (("original",1),("trapezoid",2)):
    print(label, "error ratios:", [rows[i+1][col]/rows[i][col] for i in range(len(rows)-1)])

eta = alpha*(1-2*nu)/(2*(1-nu))
print("\nExact plane-strain expanding circular pressure patch: p=1 inside r=R, p=0 outside")
eps = alpha*(1-2*nu)/(4*G*(1-nu))
lam = 2*G*nu/(1-2*nu)
sigma_inner = 2*(G+lam)*eps-alpha
sigma_outer_at_R = -2*G*eps
print("inner radial/hoop strain:",eps)
print("inner radial/hoop total stress:", sigma_inner)
print("outer radial traction at R:", sigma_outer_at_R)
print("expected inner total stress -eta:",-eta)
print("free swelling strain proposed by solver:",alpha*(1-2*nu)/(2*G))
print("traction-free crack correction K2_inf:", -eta)
print("traction-free crack correction V2_inf:", -eta*math.pi*(1-nu)/G)
assert abs(sigma_inner-sigma_outer_at_R)<1e-14
assert abs(sigma_inner+eta)<1e-14
