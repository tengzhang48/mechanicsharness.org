"""Robustness: change the mesh topology (polar-zone radius R0 = 0.25 instead of 0.4, and tip-ring size rmin)
at m = 3, L = 4(h+H), and compare psi with the production runs."""
import json
import sys
import bilayer

out = {}
for s in sys.argv[1:]:
    r = bilayer.run(s, m=3.0, Lfac=4.0, R0=0.25, rmin=1e-5, verbose=False)
    out[s] = {c: [d["psi"] for d in r["cases"][c]["rings"]] for c in ("A", "B")}
    out[s]["traction"] = {c: r["cases"][c]["traction_fit"]["0.0003-0.03"]["psi"] for c in ("A", "B")}
    print(s, json.dumps(out[s]))
with open("results/robust_R0.json", "w") as f:
    json.dump(out, f, indent=1)
