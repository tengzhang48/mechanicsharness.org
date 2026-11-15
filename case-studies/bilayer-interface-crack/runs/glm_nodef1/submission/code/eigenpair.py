"""Exact G-K relation for a bimaterial interface crack via the J-integral of
the crack-tip eigenpair (all quantities computed numerically, analytically
exact inputs).

Construction (material 1 above the interface, cut along bonded side x>0,
branch arg z in (0, 2pi)), lam = -1/2 + i*eps_m with
    e^{2 pi i lam} = -(mu2 k1 + mu1)/(mu1 k2 + mu2) =: -num/den,
    c* = conj(e^{2 pi i lam}) = -num/den (real for isotropic pairs),

    F1 =  a z^lam (UHP)      F2 = -a z^lam (LHP)          [phi' = F]
    Om1'= -conj(a) c* z^{lam~} (UHP)   Om2' = +conj(a) c* z^{lam~} (LHP)
    psi' = Om' - phi' - z phi''.

Interface checks built in: traction continuity ahead, displacement
continuity ahead, traction-free crack faces, K and delta(K) amplitudes.
"""
import numpy as np


def zpow(z, lam):
    """z^lam with branch arg(z) in (0, 2pi) (cut on the positive real axis)."""
    return np.exp(lam*(np.log(np.abs(z)) + 1j*np.mod(np.angle(z), 2*np.pi)))


class Eig:
    def __init__(self, E1, nu1, E2, nu2):
        self.mu1, self.mu2 = E1/(2*(1+nu1)), E2/(2*(1+nu2))
        self.k1, self.k2 = 3-4*nu1, 3-4*nu2
        num = self.mu2*self.k1 + self.mu1
        den = self.mu1*self.k2 + self.mu2
        self.lam = 0.5 + np.log(num/den)/(2j*np.pi) - 1.0     # singular
        self.cstar = np.conj(np.exp(2j*np.pi*self.lam))
        self.E, self.nu = np.array([E1, E2]), np.array([nu1, nu2])
        self.eps_m = self.lam.imag

    def fields(self, z, a):
        up = (np.imag(z) >= 0)
        lam, lc = self.lam, np.conj(self.lam)
        s = np.where(up, 1.0, -1.0)
        F = s*a*zpow(z, lam)
        Om = -s*np.conj(a)*self.cstar*zpow(z, lc)
        phi = s*a*zpow(z, lam+1)/(lam+1)
        phi_pp = lam*s*a*zpow(z, lam-1)
        psi_p = Om - F - z*phi_pp
        # (only psi' is needed: stresses and displacement derivatives)
        mu = np.where(up, self.mu1, self.mu2)
        kap = np.where(up, self.k1, self.k2)
        zc = np.conj(z)
        tr = 4*F.real
        dev = 2*(zc*phi_pp + psi_p)
        sxx = 0.5*(tr - dev.real)
        syy = 0.5*(tr + dev.real)
        sxy = 0.5*dev.imag
        # du/dx = (kap phi' - conj(phi') - z conj(phi'') - conj(psi')) / (2 mu)
        dudx = (kap*F - (np.conj(F) + z*np.conj(phi_pp) + np.conj(psi_p)))/(2*mu)
        return dict(sxx=sxx, syy=syy, sxy=sxy, ux=dudx.real, vx=dudx.imag, up=up)


def interface_checks(eig, a):
    print("--- interface checks")
    x = np.array([0.001, 0.01, 0.1])
    fU = eig.fields(x+1e-14j, a)
    fL = eig.fields(x-1e-14j, a)
    print("traction ahead (up):", fU["syy"]+1j*fU["sxy"])
    print("traction ahead (lo):", fL["syy"]+1j*fL["sxy"])
    xm = -x
    fUm = eig.fields(xm+1e-14j, a)
    fLm = eig.fields(xm-1e-14j, a)
    print("crack-face tractions (up):", fUm["syy"]+1j*fUm["sxy"],
          " (lo):", fLm["syy"]+1j*fLm["sxy"])
    # displacement derivative continuity ahead
    print("du/dx ahead up:", fU["ux"]+1j*fU["vx"])
    print("du/dx ahead lo:", fL["ux"]+1j*fL["vx"])
    # traction amplitude -> K
    tU = fU["syy"]+1j*fU["sxy"]
    K = tU*np.sqrt(2*np.pi*x)*np.exp(+1j*eig.eps_m*np.log(x))  # t ~ x^{-1/2 - i eps_m}
    print("K from traction:", K)


def jint(eig, a, rho):
    th = np.linspace(0, 2*np.pi, 200001)[:-1] + 1e-7
    z = rho*np.exp(1j*th)
    f = eig.fields(z, a)
    mu = np.where(f["up"], eig.mu1, eig.mu2)
    iE = np.where(f["up"], 0, 1)
    E, nu = eig.E[iE], eig.nu[iE]
    sxx, syy, sxy = f["sxx"], f["syy"], f["sxy"]
    ee11 = (1+nu)/E*((1-nu)*sxx - nu*syy)
    ee22 = (1+nu)/E*((1-nu)*syy - nu*sxx)
    ee12 = (1+nu)/E*sxy
    W = 0.5*(sxx*ee11 + syy*ee22 + 2*sxy*ee12)
    nx, ny = np.cos(th), np.sin(th)
    tx, ty = sxx*nx + sxy*ny, sxy*nx + syy*ny
    integrand = (W*nx - (tx*f["ux"] + ty*f["vx"]))*rho
    return np.trapz(integrand, th)


if __name__ == "__main__":
    E1, nu1, E2, nu2 = 9.0, 1/3, 1.0, 1/3
    eig = Eig(E1, nu1, E2, nu2)
    print(f"lam = {eig.lam:.8f}  eps_m = {eig.eps_m:.8f}  c* = {eig.cstar:.6f}")
    a = 0.7+0.5j
    interface_checks(eig, a)
    for rho in [0.01, 0.03, 0.1, 0.3]:
        J = jint(eig, a, rho)
        print(f"J(rho={rho}) = {J:.8f}")
    K = np.sqrt(2*np.pi)*np.conj(a)*(1-eig.cstar)
    print("K formula value:", K, " |K|^2 =", abs(K)**2)
    print("c = J/|K|^2:")
    for rho in [0.01, 0.1]:
        print("   ", jint(eig, a, rho)/abs(K)**2)
    # homogeneous limit
    eigh = Eig(1.0, 1/3, 1.0, 1/3)
    print("homogeneous lam:", eigh.lam)
    Jh = jint(eigh, a, 0.05)
    Kh = np.sqrt(2*np.pi)*np.conj(a)*(1-eigh.cstar)
    print(f"homogeneous c = {Jh/abs(Kh)**2:.8f}   1/Ebar = {1/(1/(1-nu1**2)):.8f}"
          f"  S/2 = {(1/(E1/(1-nu1**2)) + 1/(E2/(1-nu1**2)))/2:.8f}")
