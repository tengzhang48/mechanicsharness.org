"""Size of the complete law's sigma_y term relative to the leading field on theta = 0 (audit 2026-11-07)."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recompute_V2 import full_law
res = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "recompute_V2_results.json")))["values"]
T, K = res["T"], res["K"]
for r in (0.01, 0.0613252, 0.2):
    s = K / np.sqrt(2 * np.pi * r)
    ezz_full = full_law(s * np.diag([1.0, 1.0, 0.0]), T)[2, 2]
    ezz_lead = -0.005 / np.sqrt(r)
    print(f"r = {r:.6g} mm: sigma_e = {s:.4g} MPa, sigma_y/sigma_e = {10/s:.4f}, eps_zz lead {ezz_lead:.5g}, "
          f"complete law at leading stress {ezz_full:.5g}, relative difference {ezz_full/ezz_lead-1:+.4f}")
