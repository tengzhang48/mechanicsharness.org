"""Recheck the half-space reduction directly from the task's Biot equations.
python3 airy_mode_check.py; symbolic rational material values, general m>k>0.
No archived solver or factorization is imported.  Prints exact zero residuals.
"""
import json
import sympy as s
k,m,G,w,p0=s.symbols('k m G w p0',positive=True)
nu=s.Rational(1,5); nuu=s.Rational(2,5); alpha=s.Rational(89,100)
d=nuu-nu; a=alpha*(1-2*nu); b=a/(1-nu)
eta=d/(1-nuu); beta=eta/b; C=G/(2*(1-nu))
q2=m*m-k*k
kappa=a*a*(1-nuu)/(2*G*(1-nu)*d)  # c=1
H=-2*C*k*w; P=p0+beta*H
B=-(1+eta)*H/(2*k); D=-b*P/q2; A=(B-m*D)/k
# Phi=(A+B*y)*exp(-k*y)+D*exp(-m*y), p=P*exp(-m*y)-beta*H*exp(-k*y).
trace_k=-2*k*B; trace_m=q2*D
integral_phi=A/k+B/k**2+D/m
integral_trace=trace_k/k+trace_m/m
integral_p=P/m-beta*H/k
opening=-(-k*k*integral_phi-nu*integral_trace+a*integral_p)/G
normal=-k*k*(A+D)
gradient=-m*P+beta*k*H
normal_ref=-C*k*(1+eta*(m-k)/(m+k))*w-b*k/(m+k)*p0
gradient_ref=-m*p0+2*beta*C*k*(m-k)*w
residuals={
 'compatibility_diffusive':q2*q2*D+b*q2*P,
 'zero_shear':B-k*A-m*D,
 'harmonic_h':trace_k-b*beta*H-H,
 'mass_harmonic':a/(2*G)*(trace_k-a/d*beta*H),
 'mass_diffusive':q2*a/(2*G)*(trace_m+a/d*P)-kappa*q2*P,
 'opening':opening-w,
 'normal_surface_operator':normal-normal_ref,
 'pressure_gradient_surface_operator':gradient-gradient_ref,
}
out={name:str(s.factor(value)) for name,value in residuals.items()}
assert all(value=='0' for value in out.values()),out
print(json.dumps(out,indent=2))
