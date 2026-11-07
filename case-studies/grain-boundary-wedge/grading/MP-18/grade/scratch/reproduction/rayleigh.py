import numpy as np
from phys_oper import M_phys
from fd_nonuniform import build_D1D2


def local_rate(N, shapes):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D1, D2 = build_D1D2(nodes)
    M = M_phys(x, nodes)
    A = -(M @ D1 @ D2)
    print("N=%d" % N)
    for name, f in shapes:
        Lf = A @ f
        # Rayleigh quotient: lambda = -<f,Lf>/<f,f>
        rq = -np.dot(f, Lf) / np.dot(f, f)
        print("   %-16s Rayleigh -<f,Lf>/<f,f> = %+.6e" % (name, rq))


if __name__ == "__main__":
    N = 200
    x = np.linspace(0, 1, N + 2)[1:-1]
    shapes = [
        ("x", x),
        ("1-(1-x)^2", 1 - (1 - x) ** 2),
        ("sin(pi x/2)", np.sin(np.pi * x / 2)),
        ("x(2-x)", x * (2 - x)),
    ]
    local_rate(N, shapes)