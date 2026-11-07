import runpy
from pathlib import Path
import sys

path = Path(__file__).parent / "submission"
sys.path.insert(0, str(path))
state = runpy.run_path(str(path / "t_time.py"))
sol = state["sol"]
print("INTEGRATOR STATUS:", sol.success, sol.message,
      "last integrated time:", sol.t[-1], "requested end:", 5.0)
