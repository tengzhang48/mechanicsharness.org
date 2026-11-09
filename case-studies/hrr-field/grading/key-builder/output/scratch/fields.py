import numpy as np
from scipy.integrate import solve_bvp, solve_ivp, simpson
from scipy.optimize import root


def fourth_rhs(theta, y, n, s):
    f, v, w, z = y
    h = s - 1
    c = s * (2 - s)
    k = n * (s - 2)
    d = w + c * f
    dp = z + c * v
    a, b = 0.75, 3 * h * h
    vv = a * d * d + b * v * v
    vv = np.maximum(vv, 1e-26)
    vp = 2 * (a * d * dp + b * v * w)
    vpp_rest = 2 * (a * (dp * dp + d * c * w) + b * (w * w + v * z))
    m = (n - 1) / 2
    qp = m * vp / vv
    qpp_rest = m * vpp_rest / vv + m * (m - 1) * (vp / vv) ** 2
    coefficient = 1 + (n - 1) * a * d * d / vv
    other = c * w + 2 * qp * dp + d * qpp_rest
    f4 = (k * (k + 2) * d - 4 * h * (k + 1) * (w + qp * v) - other) / coefficient
    return np.array([v, w, z, f4])


def initial(theta, beta):
    t = theta / beta
    # Smooth symmetric positive trial function with a nonzero normal stress on faces.
    return np.array([1 - 2*t*t + t**4,
                     (-4*t+4*t**3)/beta,
                     (-4+12*t*t)/beta**2,
                     24*t/beta**3])


def solve_field(n, beta=np.pi, tol=1e-8, guess=None):
    theta = np.linspace(0, beta, 241)
    if guess is None:
        y = initial(theta, beta)
        s0 = 2 - 1/(n+1)
    else:
        y = guess.sol(theta)
        s0 = guess.p[0]
    def bc(a, b, par):
        return np.array([a[0]-1, a[1], a[3], b[0], b[1]])
    sol = solve_bvp(lambda th, y, p: fourth_rhs(th,y,n,p[0]), bc,
                    theta, y, p=[s0], tol=tol, max_nodes=30000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def shoot_field(n, beta, start):
    def ivp(par, tight=False):
        s, curvature = par
        return solve_ivp(lambda th,y: fourth_rhs(th,y,n,s), [0,beta],
                         [1,0,curvature,0], method='DOP853',
                         rtol=1e-11 if tight else 5e-10,
                         atol=1e-12 if tight else 5e-11, dense_output=True,
                         max_step=beta/160)
    def residual(par):
        out=ivp(par)
        if not out.success:
            return np.array([1e4,1e4])
        return out.y[:2,-1]
    opt = root(residual, start, tol=1e-9)
    if not opt.success and np.max(np.abs(residual(opt.x)))>1e-8:
        raise RuntimeError(opt.message)
    return opt.x, ivp(opt.x, True)


def properties(n, s, eval_y, beta=np.pi, samples=20001):
    th=np.linspace(0,beta,samples)
    f,v,w,z=eval_y(th)
    h=s-1
    d=w+s*(2-s)*f
    dp=z+s*(2-s)*v
    rr=s*f+w
    tt=s*h*f
    rt=-h*v
    q2=0.75*d*d+3*rt*rt
    q=np.sqrt(q2)
    multiplier=q**(n-1)
    er=0.75*multiplier*d
    et=-er
    ert=1.5*multiplier*rt
    # u = r^a [g(theta), h(theta)] with engineering shear 2*eps_rt.
    exponent=n*(s-2)+1
    ur=er/exponent
    qprime=(0.75*d*dp+3*h*h*v*w)/q
    erprime=0.75*((n-1)*q**(n-2)*qprime*d+multiplier*dp)
    urprime=erprime/exponent
    ut=(2*ert-urprime)/(exponent-1)
    utprime=et-ur
    # At r=1, Cartesian x derivative expressed in the rotating polar basis.
    ux_r=exponent*ur*np.cos(th)-(urprime-ut)*np.sin(th)
    ux_t=exponent*ut*np.cos(th)-(utprime+ur)*np.sin(th)
    energy=n/(n+1)*q**(n+1)
    j_density=energy*np.cos(th)-rr*ux_r-rt*ux_t
    j=2*simpson(j_density,x=th)
    energy_angle=2*simpson(energy,x=th)
    normalization=q.max()
    radial_power=(n+1)*(s-2)+2
    return dict(s=s, p=2-s, curvature=w[0], ratio=rr[0]/tt[0],
                qmax=normalization, qmin=q.min(),
                I=j/normalization**(n+1),
                H=tt[0]/j**(1/(n+1)),
                energy_ratio=energy_angle/(radial_power*j),
                j_unnormalized=j, radial_energy_power=radial_power,
                traction_residual=max(abs(tt[-1]),abs(rt[-1])))


if __name__=='__main__':
    import json
    for n, beta in [(3,np.pi),(5,np.pi),(6.5,np.pi),(4,2*np.pi/3),(3,3*np.pi/4)]:
        try:
            sol=solve_field(n,beta)
            result=properties(n, sol.p[0], sol.sol, beta)
            print(n,beta,json.dumps(result,indent=2))
            par, check=shoot_field(n,beta,[sol.p[0],sol.y[2,0]])
            print('shoot',par,properties(n,par[0],check.sol,beta))
        except Exception as e:
            print(n,beta,'ERROR',e)
