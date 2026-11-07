import fastroots

NP = (200, 120, 200)
NY = 120
L = 40.0

if __name__ == "__main__":
    eps_list = [6.0, 10.0, 14.0, 18.0, 24.0, 30.0, 36.0]
    fh = open("scan.txt", "a")
    for V in [12.0, 15.0, 18.0, 20.0, 22.0, 25.0, 30.0, 40.0]:
        vals = []
        for eps in eps_list:
            a = fastroots.a_for_eps(V, eps, NP, NY, 1.0, L=L)
            K = fastroots.get_state(V, eps, NP, NY, L).residual(a)[1]
            vals.append((eps, a, K))
        s = "V=%-5g " % V + " | ".join("%.0f:%.4f/%.1e" % (e, a, K) for e, a, K in vals)
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()
    fh.close()