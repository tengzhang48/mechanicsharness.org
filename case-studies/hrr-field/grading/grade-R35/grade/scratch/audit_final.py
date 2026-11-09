"""Grader audit: evaluate the submitted functions with the reported tip curvature.
Uses original function bodies; excludes top-level scans and demonstrations.
"""
import ast
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

source = ast.parse(Path('scratch/final.py').read_text())
functions = ast.Module(body=[node for node in source.body if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))], type_ignores=[])
ns = {'__name__': 'submission_functions'}
exec(compile(functions, 'scratch/final.py', 'exec'), ns)
n = 6.5
s = (2*n+1)/(n+1)
ns.update(n=n, s=s, alpha=0.02, sig0=400e6, E=200e9, eps0=0.002,
          p=1/(n+1), aexp=n*(s-2)+1)
rhs = ns['make_rhs'](n, s)
sol = solve_ivp(lambda t,y: [y[1],y[2],y[3],rhs(*y)], [0,np.pi],
                [1,0,-0.689256761930,0], rtol=1e-12, atol=1e-15,
                dense_output=True, max_step=0.002)
M = ns['Mcoeff'](s)
y = sol.sol(np.linspace(0,np.pi,20001))
qmax = np.sqrt(np.einsum('ij,jk,ik->i', y[:3].T,M,y[:3].T)).max()
ns.update(sol=sol, M=M, CNORM=1/qmax)
print('Reported-curvature field: original strain-check ratio vs correct tensor check')
errors=[]
for x,y in [(1e-3,1e-4),(5e-4,3e-4),(2e-3,1e-3),(1e-3,5e-4)]:
    x,y=np.array([x]),np.array([y])
    exx,eyy,exy=ns['strain_cart'](x,y,1.0)
    em=(exx+eyy)/3
    corrected=np.sqrt((2/3)*((exx-em)**2+(eyy-em)**2+em**2+2*exy**2))
    expected=ns['eps_e_from_stress'](x,y,1.0)
    bad=ns['eps_e_from_u'](x,y,1.0)
    sxx,syy,sxy,*_=ns['sigma_cart'](x,y,1.0)
    se=ns['se_cart'](x,y,1.0)
    rho=0.02*0.002*(se/400e6)**(n-1)/400e6
    predicted=np.array([0.75*rho*(sxx-syy), -0.75*rho*(sxx-syy), 1.5*rho*sxy])
    actual=np.array([exx,eyy,exy])
    rel=np.max(np.abs(actual-predicted))/np.max(np.abs(predicted))
    errors.append(float(rel))
    print(f'x={x[0]:.1e},y={y[0]:.1e}: original={bad[0]/expected[0]:.9f}, correct={corrected[0]/expected[0]:.9f}, tensor_relative_error={rel:.3e}')
print('Maximum tensor constitutive error:',max(errors))
assert max(errors)<1e-6
th=np.linspace(0,np.pi,40001)
sr,st,srt,se,ur,urp,ut,utp=ns['angular'](th)
a=ns['aexp']
I=2*np.trapezoid((n/(n+1))*se**(n+1)*np.cos(th)
                 -np.sin(th)*(sr*(ut-urp)-srt*(ur+utp))
                 -a*np.cos(th)*(sr*ur+srt*ut),th)
Jn=ns['J_numeric'](1.0,1e-3,nth=40001)
Jcorrect=ns['alpha']*ns['eps0']*ns['sig0']*I
print(f'Original Cartesian J={Jn:.10e}, correct amplitude formula J={Jcorrect:.10e}, ratio={Jn/Jcorrect:.10f}')
assert abs(Jn/Jcorrect-1)<1e-6
