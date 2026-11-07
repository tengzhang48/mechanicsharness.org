import fastroots

NP = (300, 180, 300)
NY = 180
L = 70.0
EPS = [8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0]

if __name__ == "__main__":
    fh = open("crit.txt", "a")
    for V in [12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0]:
        vals = []
        for eps in EPS:
            a = fastroots.a_for_eps(V, eps, NP, NY, 1.0, L=L)
            K = fastroots.get_state(V, eps, NP, NY, L).residual(a)[1]
            vals.append((eps, a, K))
        s = "V=%-5g " % V + " ".join("%.0f:%.3f/%.1e" % (e, a, K) for e, a, K in vals)
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()
    fh.close()