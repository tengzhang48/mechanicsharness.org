"""Collect results/*.json into convergence tables and the final table (printed as Markdown)."""
import json
import glob
import numpy as np

SETS = ["V1", "S1", "S2", "S3", "S4"]
WIN = "0.0003-0.03"


def load(s, m, L):
    try:
        with open(f"results/{s}_m{m:g}_L{L:g}.json") as f:
            return json.load(f)
    except FileNotFoundError:
        return None


def mint(c):
    ps = np.array([d["psi"] for d in c["rings"]])
    ks = np.array([d["absK"] for d in c["rings"]])
    return ps.mean(), ps.max() - ps.min(), ks.mean()


def main():
    final = {}
    print("## Convergence tables (psi in degrees; M-int = mean over 4 rings, spread = max-min over rings)\n")
    for s in SETS:
        print(f"### {s}\n")
        print("| m | L/(h+H) | dofs | psi_A M-int | ring spread | psi_A traction fit | psi_B M-int | ring spread | "
              "psi_B traction fit | abs(K_A)/abs(K)_G - 1 | abs(K_B)/abs(K)_G - 1 | J_A/G-1 | J_B/G-1 |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for (m, L) in [(1, 4), (2, 4), (3, 4), (4, 4), (2, 2), (2, 8), (3, 8)]:
            r = load(s, m, L)
            if r is None:
                continue
            A, B = r["cases"]["A"], r["cases"]["B"]
            pa, sa, ka = mint(A)
            pb, sb, kb = mint(B)
            ja = np.mean([d["J_over_G"] for d in A["rings"]]) - 1
            jb = np.mean([d["J_over_G"] for d in B["rings"]]) - 1
            print(f"| {m} | {L} | {r['ndofs']} | {pa:.4f} | {sa:.4f} | {A['traction_fit'][WIN]['psi']:.4f} | "
                  f"{pb:.4f} | {sb:.4f} | {B['traction_fit'][WIN]['psi']:.4f} | {ka / A['absK_from_G'] - 1:+.1e} | "
                  f"{kb / B['absK_from_G'] - 1:+.1e} | {ja:+.1e} | {jb:+.1e} |")
        print()
        # final values
        rec = {}
        for cn in ("A", "B"):
            vals = {}
            for key in [(3, 4), (4, 4), (3, 8), (2, 4), (2, 8), (2, 2)]:
                r = load(s, *key)
                if r is not None:
                    vals[key] = r["cases"][cn]
            best = vals[(4, 4)] if (4, 4) in vals else vals[(3, 4)]
            p_best, spread, k_best = mint(best)
            d_mesh = abs(mint(vals[(4, 4)])[0] - mint(vals[(3, 4)])[0]) if (4, 4) in vals else abs(
                mint(vals[(3, 4)])[0] - mint(vals[(2, 4)])[0])
            d_trunc = abs(mint(vals[(3, 8)])[0] - mint(vals[(3, 4)])[0])
            # traction fit, Richardson extrapolation in m^-2 using the two finest meshes at L=4
            if (4, 4) in vals:
                t3, t4 = vals[(3, 4)]["traction_fit"][WIN]["psi"], vals[(4, 4)]["traction_fit"][WIN]["psi"]
                t_rich = t4 - (t3 - t4) * (1 / 16) / (1 / 9 - 1 / 16)
                k3, k4 = vals[(3, 4)]["traction_fit"][WIN]["absK"], vals[(4, 4)]["traction_fit"][WIN]["absK"]
                k_rich = k4 - (k3 - k4) * (1 / 16) / (1 / 9 - 1 / 16)
            else:
                t2, t3 = vals[(2, 4)]["traction_fit"][WIN]["psi"], vals[(3, 4)]["traction_fit"][WIN]["psi"]
                t_rich = t3 - (t2 - t3) * (1 / 9) / (1 / 4 - 1 / 9)
                k2, k3 = vals[(2, 4)]["traction_fit"][WIN]["absK"], vals[(3, 4)]["traction_fit"][WIN]["absK"]
                k_rich = k3 - (k2 - k3) * (1 / 9) / (1 / 4 - 1 / 9)
            rec[cn] = dict(psi=p_best, ring_spread=spread, d_mesh=d_mesh, d_trunc=d_trunc, psi_traction_rich=t_rich,
                           absK=k_best, absK_G=best["absK_from_G"], absK_traction_rich=k_rich,
                           K=np.mean([d["K"] for d in best["rings"]], axis=0).tolist(), eps=best["eps"],
                           G=best["G_beam"])
        final[s] = rec
    print("## Final values (M-integral on the finest mesh, L = 4(h+H); traction fit Richardson-extrapolated in m^-2)\n")
    print("| set | eps | case | psi M-int | ring spread | mesh change | truncation change (L 4->8) | psi traction (extrap.) | "
          "K1 | K2 | abs(K) M-int | abs(K) from G | rel. diff | abs(K) traction (extrap.) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in SETS:
        for cn in ("A", "B"):
            d = final[s][cn]
            print(f"| {s} | {d['eps']:.6f} | {cn} | {d['psi']:.4f} | {d['ring_spread']:.4f} | {d['d_mesh']:.4f} | "
                  f"{d['d_trunc']:.4f} | {d['psi_traction_rich']:.4f} | {d['K'][0]:.6f} | {d['K'][1]:.6f} | {d['absK']:.6f} | "
                  f"{d['absK_G']:.6f} | {d['absK'] / d['absK_G'] - 1:+.1e} | {d['absK_traction_rich']:.6f} |")
    with open("results/final.json", "w") as f:
        json.dump(final, f, indent=1)


if __name__ == "__main__":
    main()
