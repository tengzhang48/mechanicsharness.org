"""Test corrected ODE at n=1: should reproduce a biharmonic eigenfunction.
At n=1 the material is linear (se^{n-1}=1, q=1), the ODE linearizes.
Mode I Williams: sigma ~ r^{-1/2}, i.e. p=1/(n+1)=1/2. Check that
phi = cos((1)theta/2)? Actually with p=1/2, k=2-1/2=3/2. The symmetric
biharmonic eigenfunction for opening mode: Phi = r^{3/2} phi with phi(0)=1,
phi'(0)=0 and A(pi)=C(pi)=0. Solution: phi = c1 cos(3th/2)+c2 cos(th/2)
(with the right combination)."""
import numpy as np
# For n=1 the ODE reduces to biharmonic: (d^2/dth^2 + ... ) ; verify the
# candidate phi = cos(3 th/2) - 3 cos(th/2)? Let's check which satisfies BC.
# General symmetric Airy eigenfunc r^k phi, k=3/2: A=k phi+phi'', C=-(k-1)phi'
# Try phi = a*cos(3th/2)+b*cos(th/2). phi(0)=a+b=1. phi'(0)=0 auto.
# A = k phi + phi'' : phi'' = a*(-9/4)cos(3th/2)+b*(-1/4)cos(th/2)
# At pi: cos(3pi/2)=0, cos(pi/2)=0 -> A(pi)=0 automatically?! 
# C(pi)=-(k-1)phi'(pi); phi'= a*(-3/2)sin(3th/2)+b*(-1/2)sin(th/2)
# at pi: sin(3pi/2)=-1, sin(pi/2)=1 -> phi'(pi)=(3/2)a-(1/2)b
# set =0 => 3a=b. With a+b=1 => 4a=1 => a=1/4,b=3/4.
a=0.25;b=0.75
print("n=1 candidate phi=0.25 cos(3t/2)+0.75 cos(t/2); phi(0)=%.3f phi'(0)=%.3f"%(a+b,0))
# This is the classic HRR n=1 (=Williams) angular function (up to normalization).
print("phi''(0) =", a*(-9/4)+b*(-1/4))
print("phi'''(0) =", a*0+b*0)  # third deriv: a*(27/8)sin(3t/2)+b*(1/8)sin(t/2)|0 =0