"""Which checks catch each planted error (log P240). Run: python3 mutant_detection.py > mutant_detection.txt"""
import numpy as np
import test_checks as T
from crack_tip_checks import check_field

def main():
    th, a, b, c = T.hrr_field(3, True)
    ths, as_, bs, cs = T.hrr_field(3, False)
    wn = T.williams_notch(2 * np.pi / 3)
    wc = T.williams_crack()
    keep = np.abs(wc[0]) <= 0.9 * np.pi + 1e-12
    cases = {
        "wrong p (HRR n=3, p=0.26)": check_field(th, a, b, c, p=0.26, n=3),
        "notch, p off by 0.02 (R23 type)": check_field(*wn[:4], p=wn[4] + 0.02, law="linear", alpha=2 * np.pi / 3),
        "S_rt scaled by 1.05": check_field(th, a, b, 1.05 * c, p=0.25, n=3),
        "declared n=5 for an n=3 field": check_field(th, a, b, c, p=0.25, n=5),
        "plane stress field declared plane strain": check_field(ths, as_, bs, cs, p=0.25, n=3, plane="strain"),
        "reported I = 5.6476 (dropped basis derivatives)": check_field(ths, as_, bs, cs, p=0.25, n=3, plane="stress",
                                                                       reported_J=5.6476),
        "faces at 0.9 pi not traction-free": check_field(wc[0][keep], wc[1][keep], wc[2][keep], wc[3][keep], p=0.5,
                                                         law="linear", alpha=wc[0][keep][-1]),
        "broken symmetry": check_field(wc[0], wc[1] + 0.01 * np.max(np.abs(wc[1])) * np.sin(wc[0]), wc[2], wc[3],
                                       p=0.5, law="linear"),
        "infinite-energy exponent p=0.6, n=3": check_field(th, a, b, c, p=0.6, n=3),
    }
    print("Planted error -> checks that FAIL")
    for k, r in cases.items():
        print(f"- {k}: " + "; ".join(sorted(x["check"] for x in r if x["passed"] is False)))


if __name__ == "__main__":
    main()
