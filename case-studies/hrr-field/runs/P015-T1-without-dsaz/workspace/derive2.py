"""Check constitutive relations (I),(II),(III) directly for Williams field at n=1."""
import sympy as sp
t = sp.symbols('t')
p = sp.Rational(1,2)
k = sp.Rational(3,2)

def make(phiW):
    A = k*phiW + sp.diff(phiW,t,2)
    B = k*(k-1)*phiW
    C = -(k-1)*sp.diff(phiW,t)
    se2 = sp.Rational(3,4)*(A-B)**2 + 3*C**2
    se = sp.sqrt(sp.simplify(se2))
    return A,B,C,se

# generic Williams combination
a,b = sp.symbols('a b')
phiW = a*sp.cos(sp.Rational(3,2)*t)+b*sp.cos(sp.Rational(1,2)*t)
A,B,C,se = make(phiW)
print("se^2 =", sp.simplify(sp.expand(se**2)))

# the 3 constitutive relations at n=1 (q=1):
#   p Ur = (A-B)/2 ; Ur+Ut'=(B-A)/2 ; Ur'+(p-1)Ut = 2C
# Eliminate Ur,Ut to get the compatibility condition that phi must satisfy:
Ur = (A-B)/(2*p)
# from III: Ut = (2C - Ur')/(p-1)
Ut = (2*C - sp.diff(Ur,t))/(p-1)
resid = sp.simplify(Ur + sp.diff(Ut,t) - (B-A)/2)
print("compat residual (II):", sp.expand(sp.simplify(resid)))

# Also check the Williams Mode I "true" eigenfunction (stress intensity form)
# For the standard K-field: sigtilde_rr = cos(t/2)(1-sin(t/2)sin(3t/2)) etc.
# Build phi from B=k(k-1)phi = sigtt:
srr = sp.cos(t/2)*(1-sp.sin(t/2)*sp.sin(3*t/2))
stt = sp.cos(t/2)*(1+sp.sin(t/2)*sp.sin(3*t/2))
srt = sp.cos(t/2)*sp.sin(t/2)*sp.cos(3*t/2)
# B = k(k-1) phi, k(k-1)=3/2*1/2=3/4 => phi = stt/(3/4)
phi_K = stt/sp.Rational(3,4)
AK,BK,CK,seK = make(phi_K)
print("check A matches srr:", sp.simplify(AK-srr))
print("check C matches srt:", sp.simplify(CK-srt))
UrK = (AK-BK)/(2*p)
UtK = (2*CK - sp.diff(UrK,t))/(p-1)
print("compat residual for K-field:", sp.simplify(sp.expand(UrK + sp.diff(UtK,t) - (BK-AK)/2)))