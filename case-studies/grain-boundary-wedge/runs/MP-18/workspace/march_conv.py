import numpy as np, march


def march_slowest(N):
    M, D1m, D2m, L = march.operator(N)
    ev = np.linalg.eigvals(L).real
    return np.sort(ev)[::-1]        # descending; slowest decay = max (closest to 0)


if __name__ == "__main__":
    print("=== march.operator(N): top real eigenvalues (slowest first) ===")
    for N in [30, 60, 120, 240, 480]:
        r = march_slowest(N)
        print("N=%3d  top: %s" % (N, np.round(r[:4], 6)))