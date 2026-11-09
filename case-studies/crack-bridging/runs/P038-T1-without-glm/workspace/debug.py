import numpy as np
import solve as S  # reuse constants

# debug spectral solver on anchor cases
def solve_SIE_dbg(Sr, c, tt, M=2048):
    a = S.a; Ep = S.Ep
    xa = a/c
    k  = np.arange(1, M+1)
    xi = np.cos(k*np.pi/(M+1))
    s2 = np.sin(k*np.pi/(M+1))**2
    g  = np.where(np.abs(xi) > xa, Sr - tt, Sr)
    F  = -(4.0*c/Ep)*g
    N  = M
    b  = np.zeros(N+1)
    U0 = np.ones(M); U1 = 2*xi
    b[1] = (2.0/(M+1))*np.dot(s2*F, U0)
    b[2] = (2.0/(M+1))*np.dot(s2*F, U1)
    for n in range(3, N+1):
        U2 = 2*xi*U1 - U0
        b[n] = (2.0/(M+1))*np.dot(s2*F, U2)
        U0, U1 = U1, U2
    K = -(Ep/4)*np.sqrt(np.pi/c)*b[1:].sum()
    def delt(xi0):
        s = 0.0; V0 = 1.0; V1 = 2*xi0
        s = b[1]*V0
        for n in range(2, N+1):
            V2 = 2*xi0*V1 - V0
            s += b[n]*V2/n
            V0, V1 = V1, V2
        return -np.sqrt(1-xi0*xi0)*s
    return K, delt, b

c = S.c_s
# case A: pure remote (tt=0) -> K = Sr*sqrt(pi c), delta(x)=4 Sr sqrt(c^2-x^2)/Ep
KA, dA, bA = solve_SIE_dbg(100e6, c, 0.0, M=1024)
print("A: K spec/exact:", KA/1e6, 100e6*np.sqrt(np.pi*c)/1e6)
print("A: delta(0) spec/exact:", dA(0.0)*1e6, 4*100e6*c/S.Ep*1e6)
print("A: delta(a) spec/exact:", dA(S.a/c)*1e6, 4*100e6*np.sqrt(c*c-S.a*S.a)/S.Ep*1e6)

# case B: Dugdale (Sigma_D, t) -> K=0, delta(a)= 6.8914 um
SD = (2*S.t/np.pi)*S.theta(c)
KB, dB, bB = solve_SIE_dbg(SD, c, S.t, M=1024)
print("B: K spec (should be ~0):", KB/1e6)
print("B: delta(a) spec vs lit:", dB(S.a/c)*1e6, (8*S.t*S.a/(np.pi*S.Ep))*np.log(1/np.cos(S.theta(c)))*1e6)
print("B: b_n tail size:", np.abs(bB[-5:]).max(), " b1..b3:", bB[1:4])
