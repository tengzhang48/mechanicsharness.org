import fastroots

NP = (400, 250, 400)
NY = 250
L = 100.0

if __name__ == "__main__":
    fh = open("crit2.txt", "a")
    for V in [13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8]:
        row = []
        for eps in [40.0, 60.0, 80.0]:
            a = fastroots.a_for_eps(V, eps, NP, NY, 1.0, L=L)
            K = fastroots.get_state(V, eps, NP, NY, L).residual(a)[1]
            row.append("%.0f:%.5f/%.2e" % (eps, a, K))
        s = "V=%-5g " % V + " ".join(row)
        print(s, flush=True)
        fh.write(s + "\n")
        fh.flush()
    fh.close()