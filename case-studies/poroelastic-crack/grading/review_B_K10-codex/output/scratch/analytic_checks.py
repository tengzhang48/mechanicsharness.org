"""Algebra, weak-coupling tip solution, and six-term inversion checks."""
import json
import math
from fractions import Fraction
from pathlib import Path
from datetime import datetime,timezone
import sympy as sp
import numpy as np
from scipy.integrate import quad

k,m,y,eta,b,H,P,nu,G=sp.symbols('k m y eta b H P nu G',positive=True)
q2=m*m-k*k
B=-(1+eta)*H/(2*k)
C=-b*P/q2
A=(B-m*C)/k
Phi=(A+B*y)*sp.exp(-k*y)+C*sp.exp(-m*y)
p=P*sp.exp(-m*y)-eta/b*H*sp.exp(-k*y)
trace=sp.diff(Phi,y,2)-k*k*Phi
h=sp.simplify(trace+b*p)
checks={
 'harmonic_stress_trace':sp.simplify(h-H*sp.exp(-k*y)),
 'zero_surface_shear':sp.simplify(sp.diff(Phi,y).subs(y,0)),
 'pressure_diffusion':sp.simplify(sp.diff(p,y,2)-k*k*p-q2*(p+eta/b*h)),
}
# Integrate strain modes algebraically from 0 to infinity.
intPhi=A/k+B/k**2+C/m
intTrace=(1+eta)*H/k-b*P/m
intP=P/m-eta/b*H/k
uy0=(k*k*intPhi+nu*intTrace-b*(1-nu)*intP)/(2*G)
checks['surface_displacement']=sp.simplify(uy0+(1-nu)*H/(2*G*k))
checks['normal_stress']=sp.simplify(-k*k*Phi.subs(y,0)-((1+eta)*H/2-b*k/(m+k)*P))
checks['normal_pressure_gradient']=sp.simplify(sp.diff(p,y).subs(y,0)-(-m*P+eta*k/b*H))
print(json.dumps({'utc':datetime.now(timezone.utc).strftime('%H:%M:%S'),'symbolic_residuals':{name:str(value) for name,value in checks.items()}}))
assert all(value==0 for value in checks.values())
# Weak coupling: H=f=r^(-1/2)cos(theta/2), Z=e^(-r)f,
# P=(e^(-r)-1)f, so B_pressure=-1 and I=pi.
radial=quad(lambda R:math.exp(-R),0,np.inf,epsabs=1e-12)[0]
angular=quad(lambda th:math.cos(th/2)**2,-math.pi,math.pi,epsabs=1e-12)[0]
print(json.dumps({'weak_coupling_tip_integral':radial*angular,'exact_pi':math.pi,
                  'A_over_d_exact_weak_limit':(16+4*math.pi)/math.sqrt(math.pi),
                  'A_over_d_paper_weak_limit':16/math.sqrt(math.pi),
                  'K_weak_slope_exact':-.75,'K_weak_slope_paper':-1}))
# The same weak-coupling cell directly in physical time (c=1).
r,t=sp.symbols('r t',positive=True)
pr=-r**(-sp.Rational(1,2))*sp.erf(r/(2*sp.sqrt(t)))
heat=sp.simplify(sp.diff(pr,t)-sp.diff(pr,r,2)-sp.diff(pr,r)/r+pr/(4*r*r))
print(json.dumps({'weak_cell_physical_time_heat_residual':str(heat),
                  'weighted_tip_integral_over_sqrt_ct':2*math.sqrt(math.pi)}))
assert heat==0
# Independent six-term Stehfest weights, as in Appendix B.
N=6; weights=[]
for n in range(1,N+1):
    total=Fraction(0)
    for j in range((n+1)//2,min(n,N//2)+1):
        total+=Fraction(j**(N//2)*math.factorial(2*j),
                        math.factorial(N//2-j)*math.factorial(j)*math.factorial(j-1)*
                        math.factorial(n-j)*math.factorial(2*j-n))
    weights.append((-1)**(n+N//2)*total)
const=sum(w/Fraction(i) for i,w in enumerate(weights,1))
sqrt_coef=sum(float(w)/(i*math.log(2))**1.5 for i,w in enumerate(weights,1))*math.log(2)
print(json.dumps({'Stehfest_N6_weights':[str(w) for w in weights],
                  'step_factor':str(const),'sqrt_t_coefficient':sqrt_coef,
                  'exact_sqrt_t_coefficient':2/math.sqrt(math.pi),
                  'relative_error_sqrt_t':sqrt_coef/(2/math.sqrt(math.pi))-1}))
tau=1e-4
deltas=[math.sqrt(tau/(j*math.log(2))) for j in range(1,7)]
print(json.dumps({'tau':tau,'element_length_over_L':1/3,'nearest_stress_node_from_tip_over_L':1/12,
                  'nearest_nonendpoint_density_node_from_tip_over_L':1/6,
                  'Laplace_tip_layer_lengths_over_L':deltas,
                  'stress_node_to_layer_ratios':[(1/12)/delta for delta in deltas]}))
nu=.2;nuu=.4;d=nuu-nu;eta=d/(1-nuu)
K=.807593003;n21=.4169594526473146
face=16*d/((1+eta)*math.sqrt(math.pi))
tip=4*d/((1+eta)*math.sqrt(math.pi))*math.pi*math.sqrt(2/eta)*n21
print(json.dumps({'volume_face_term':face,'volume_tip_term':tip,'volume_total':face+tip,
                  'energy_hydraulic_fraction':n21*n21/(1+eta),
                  'paper_K_error_percent':100*(K/.75-1)}))
# MP-39's own finite-crack volumes: adjacent differences cancel the unknown
# periodic-box initial volume. The claimed resolvent slope is only 1.2.
ell=np.array([.1,.05,.025,.0125])
vol=np.array([2.0429251423911685,1.9794399541751084,1.9422132554242526,1.921049745437891])
print(json.dumps({'MP39_adjacent_volume_slopes':((vol[:-1]-vol[1:])/(ell[:-1]-ell[1:])).tolist(),
                  'paper_resolvent_slope':face*math.sqrt(math.pi)/2,
                  'correct_resolvent_slope':(face+tip)*math.sqrt(math.pi)/2}))
