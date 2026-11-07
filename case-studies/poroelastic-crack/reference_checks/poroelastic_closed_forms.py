# Independent arithmetic of D&C 1991 closed forms and the quoted correction factors.
import math
nu,nuu,al=0.2,0.4,0.89
print('K1_0 heuristic', (1-nuu)/(1-nu))
print('V1_0', math.pi*(1-nuu), 'V1_inf', math.pi*(1-nu))
A_face=16*(1-nuu)*(nuu-nu)/(math.sqrt(math.pi)*(1-nu)); print('A faces', A_face, 'A total', A_face+1.086165)
print('p_faces_0', 2*(nuu-nu)/(al*(1-2*nu)))
eta=al*(1-2*nu)/(2*(1-nu)); print('eta DC', eta, 'V2_inf', -al*math.pi*(1-2*nu)/2)
eta_t=(nuu-nu)/(1-nuu); print('eta_t', eta_t)
print('0.75*1.07679067 =', 0.75*1.07679067)
print('exponent B', math.log(0.807593)/math.log(0.75), ' energy-matching (1+eta_t)^-1/2 =', (1+eta_t)**-0.5)
print('Wilson eta_t=1 exponent', -math.log(0.601850)/math.log(2))
# confined pressurized region (plane strain, thermal analogy): free 3D strain = al(1-2nu)p/E, interior stress = -E eps/(2(1-nu))
E=1.0; eps=al*(1-2*nu)/E; print('interior stress/p', -E*eps/(2*(1-nu)), '(should equal -eta)')
