"""Validate the L-operator quadrature against closed forms / reference quadrature."""
import numpy as np
from scipy.integrate import quad
from kernel import LOperator, S0_analytic, s_reg


def s_full(Z, Zp):
    return 1.0/(Z-Zp) - 1.0/(Z+Zp) - 2.0*Zp*(Z-Zp)/(Z+Zp)**3


def L_reference(Z, g, gp):
    """Trustworthy reference: PV by subtraction + adaptive quad."""
    gZ = g(Z)
    f1 = lambda zp: (g(zp)-gZ)/(Z-zp) if abs(zp-Z) > 1e-13 else -gp(Z)
    cauchy, _ = quad(f1, 0.0, 1.0, points=[Z], epsabs=1e-12, epsrel=1e-11, limit=500)
    cauchy += gZ*np.log(Z/(1.0-Z))
    f2 = lambda zp: s_reg(Z, zp)*g(zp)
    reg, _ = quad(f2, 0.0, 1.0, epsabs=1e-12, epsrel=1e-11, limit=500)
    return cauchy + reg


def test_g1():
    Lop = LOperator(M=256)
    Gfun = lambda y: 2.0*y            # g = 1
    Gpfun = lambda y: 2.0*np.ones_like(y)
    Zs = np.array([1e-3, 0.05, 0.2, 0.5, 0.8, 0.95, 0.999])
    got = Lop.apply(Zs, Gfun, Gpfun)
    want = S0_analytic(Zs)
    print("test g=1 vs S0 analytic:")
    for z, a, b in zip(Zs, got, want):
        print(f"  Z={z:8.4f}  quad={a: .10f}  analytic={b: .10f}  diff={a-b: .2e}")
    return


def test_singular():
    # g(Z') = (1-Z')^{-1/2}:  G(y) = 2y*(1/y) = 2
    Lop = LOperator(M=256)
    Gfun = lambda y: 2.0*np.ones_like(y)
    Gpfun = lambda y: np.zeros_like(y)
    Zs = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0])
    print("test g=(1-Z')^{-1/2} vs reference quad:")
    for z in Zs:
        a = Lop.apply(np.array([z]), Gfun, Gpfun)[0]
        if z < 1.0:
            b = L_reference(z, lambda zp: 1.0/np.sqrt(1.0-zp),
                            lambda zp: 0.5/(1.0-zp)**1.5)
            yz = np.sqrt(1.0-z)
            cana = np.log((1-yz)/(1+yz))/yz
            print(f"  Z={z:6.3f}  quad={a: .10f}  ref={b: .10f}  diff={a-b: .2e}   (analytic cauchy part={cana: .8f})")
        else:
            reg, _ = quad(lambda zp: s_reg(1.0, zp)/np.sqrt(1.0-zp), 0, 1,
                          epsabs=1e-12, epsrel=1e-11, limit=500)
            print(f"  Z={z:6.3f}  quad={a: .10f}  ref(cauchy -2 + reg {reg: .8f})={-2.0+reg: .10f}  diff={a-(-2+reg): .2e}")


def test_poly():
    # g(Z') = Z'^2  -> G(y) = 2y (1-y^2)^2 ; compare with reference quad
    Lop = LOperator(M=256)
    Gfun = lambda y: 2.0*y*(1.0-y**2)**2
    Gpfun = lambda y: 2.0*(1.0-y**2)**2 - 8.0*y**2*(1.0-y**2)
    Zs = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    print("test g=Z'^2 vs reference quad:")
    for z in Zs:
        a = Lop.apply(np.array([z]), Gfun, Gpfun)[0]
        b = L_reference(z, lambda zp: zp**2, lambda zp: 2.0*zp)
        print(f"  Z={z:6.3f}  quad={a: .10f}  ref={b: .10f}  diff={a-b: .2e}")


if __name__ == "__main__":
    test_g1()
    test_singular()
    test_poly()
