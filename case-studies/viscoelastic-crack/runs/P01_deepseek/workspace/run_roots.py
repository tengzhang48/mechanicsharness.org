import numpy as np
import newton

RES = [(160, 100, 160, 160), (240, 150, 240, 190), (320, 200, 320, 230)]
GUESS = {0.001: (0.44, 0.52), 0.1: (0.516, 0.77), 10.0: (0.95, 1.0)}

if __name__ == "__main__":
    out = {}
    for V, g in GUESS.items():
        for (n1, n2, n3, ny) in RES:
            a, eps = newton.newton(V, g[0], g[1], (n1, n2, n3), ny, verbose=False)
            print('V=%-8g nx=%-4d ny=%-4d a=%.7f eps=%.7f' % (V, n1, ny, a, eps), flush=True)
            g = (a, eps)