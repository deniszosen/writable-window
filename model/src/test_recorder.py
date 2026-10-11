"""Unit tests for the model's components.

`verify_numbers.py` checks the manuscript against a stored results file and
`check_keys.py` checks that the file can be produced. Neither tests whether the
model is correct, which is what a reader auditing the code wants to know. These
assertions test properties that must hold by construction, so each one fails if
the corresponding piece of the model breaks.

    python test_recorder.py        # exits non-zero on the first failure
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate
import analysis as A

BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000, T=150)
ok, bad = 0, []


def check(name, cond, note=""):
    global ok
    if cond:
        ok += 1
        print(f"  PASS  {name}{'  ' + note if note else ''}")
    else:
        bad.append(name)
        print(f"  FAIL  {name}{'  ' + note if note else ''}")


s = simulate(Params(**BASE))

# --- the array is written in order and never overfilled ---------------------
check("occupancy never exceeds array length",
      int(s["slot"].max()) <= Params(**BASE).L,
      f"max {int(s['slot'].max())} of {Params(**BASE).L}")
check("occupancy is non-negative", int(s["slot"].min()) >= 0)

# --- writing is bounded by divisions ----------------------------------------
check("no cell holds more symbols than it had divisions",
      bool(np.all(s["slot"] <= s["ndiv"] + 1e-9)),
      f"worst surplus {float(np.max(s['slot'] - s['ndiv'])):.3g}")

# --- occupancy matches the binomial expectation ------------------------------
p = Params(**BASE)
exp = float(np.mean(np.minimum(s["ndiv"] * p.q_div, p.L)))
got = float(s["slot"].mean())
check("mean occupancy matches n_div x q_div", abs(got - exp) / exp < 0.05,
      f"expected {exp:.3f}, got {got:.3f}")

# --- a full array cannot record anything further -----------------------------
full = simulate(Params(**{**BASE, "L": 1, "q_div": 1.0}))
check("a one-monomer array at q = 1 fills and stops",
      float(full["slot"].max()) == 1.0)

# --- no writing without writing ----------------------------------------------
silent = simulate(Params(**{**BASE, "q_div": 0.0, "q_nondiv": 0.0}))
check("q = 0 writes nothing", float(silent["slot"].max()) == 0.0)

# --- the gradient is 1 when the tissue never stops cycling -------------------
flat = simulate(Params(**{**BASE, "pp_start": 0.3, "pp_end": 0.3,
                          "nn_start": 0.0, "nn_end": 0.0,
                          "tc_start": 30.0, "tc_end": 30.0}))
g = A.epoch_gradient(flat)
# With no cycle exit the population grows exponentially and every terminal
# lineage passes through the same number of divisions in any window, so the
# two trends that produce front-loading both vanish and the gradient collapses
# to 1. This is the null the figure 2 axis draws and no mechanism arm reaches.
check("a population that never exits the cycle has no gradient",
      0.75 < g < 1.3, f"gradient {g:.2f} against 5.28 with cycle exit")

# --- determinism --------------------------------------------------------------
a = simulate(Params(**BASE))
check("the same seed gives the same occupancy",
      float(a["slot"].mean()) == float(s["slot"].mean()))
b = simulate(Params(**{**BASE, "seed": BASE["seed"] + 1}))
check("a different seed gives a different realisation",
      float(b["slot"].mean()) != float(s["slot"].mean()))

# --- cycling fraction is a fraction ------------------------------------------
cf = np.asarray(s["cyc_frac"], float)
check("cycling fraction stays in [0, 1]",
      bool(np.all((cf >= -1e-12) & (cf <= 1 + 1e-12))))
check("cycling fraction starts at 1", abs(float(cf[0]) - 1.0) < 1e-9)
check("cycling fraction decreases over the culture", float(cf[-1]) < float(cf[0]),
      f"{float(cf[0]):.2f} -> {float(cf[-1]):.3f}")

print(f"\n{ok} passed, {len(bad)} failed")
sys.exit(1 if bad else 0)
