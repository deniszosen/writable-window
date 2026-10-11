"""Does the sensitivity grid show interaction, and how much could it have seen?

`run_sensitivity.py` reports the gradient over a grid of writing rate against
fate-switch midpoint. "No interaction" means the surface is additive,
g(q, mid) = mu + a(q) + b(mid); the interaction is whatever an additive fit
leaves behind. That residual is never zero, because each cell carries Monte
Carlo noise, so it is compared here against the residual that the cells' own
standard errors produce on a surface built to be exactly additive.

Testing separability as a rank-1
approximation. That was wrong twice over: an additive surface is rank 2, not
rank 1, and the perturbation used to measure power was itself an outer
product and so contained no interaction at all.

    python interaction_test.py            report, and check the stored block
    python interaction_test.py --write    merge into results.json
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
R = json.loads((HERE / "results.json").read_text())
S = R["sensitivity_2d"]
g = np.array(S["gradient"], float)
sem = np.array([[c["gradient_sem"] for c in S["cells"] if c["fate_midpoint"] == m]
                for m in S["fate_midpoint"]], float)


def interaction_ss(m):
    """Residual sum of squares after the additive model mu + a_i + b_j."""
    resid = (m - m.mean(axis=1, keepdims=True) - m.mean(axis=0, keepdims=True)
             + m.mean())
    return float(np.sum(resid ** 2))


additive = (g.mean(axis=1, keepdims=True) + g.mean(axis=0, keepdims=True)
            - g.mean())                                  # the no-interaction fit
observed = interaction_ss(g)

rng = np.random.default_rng(0)
draws = np.array([interaction_ss(additive + rng.normal(0.0, sem))
                  for _ in range(4000)])
crit = float(np.percentile(draws, 95))
p = float(np.mean(draws >= observed))

# An interaction pattern orthogonal to every row and column effect, so it is
# pure interaction and invisible to the additive model.
pat = np.outer(np.linspace(-1, 1, g.shape[0]), np.linspace(-1, 1, g.shape[1]))
pat = pat - pat.mean(axis=1, keepdims=True) - pat.mean(axis=0, keepdims=True) + pat.mean()
pat = pat / np.sqrt(np.mean(pat ** 2))                   # unit RMS


def power_at(delta, n=600, seed=1):
    r = np.random.default_rng(seed)
    return float(np.mean([interaction_ss(additive + delta * pat + r.normal(0.0, sem))
                          >= crit for _ in range(n)]))


lo, hi = 0.0, 3.0
for _ in range(16):
    mid = 0.5 * (lo + hi)
    if power_at(mid) < 0.80:
        lo = mid
    else:
        hi = mid
mde = 0.5 * (lo + hi)

out = dict(observed_interaction_ss=observed,
           noise_only_median=float(np.median(draws)),
           noise_only_p95=crit,
           p_value=p,
           detectable_interaction_rms=mde,
           grid_range=float(g.max() - g.min()),
           n_bootstrap=int(draws.size), rng_seed=0,
           verdict=("no detectable interaction" if p > 0.05
                    else "interaction beyond Monte Carlo noise"))

print(f"interaction sum of squares, observed : {observed:7.3f}")
print(f"  from Monte Carlo noise alone, median: {out['noise_only_median']:7.3f}")
print(f"  95th percentile                     : {crit:7.3f}")
print(f"  p                                   : {p:7.3f}")
print(f"smallest detectable interaction       : {mde:7.2f} gradient units RMS")
print(f"  grid spans                          : {out['grid_range']:7.2f}")
print(f"verdict                               : {out['verdict']}")

if "--write" in sys.argv:
    R["sensitivity_2d"]["interaction_test"] = out
    R["sensitivity_2d"].pop("noise_comparison", None)
    (HERE / "results.json").write_text(json.dumps(R, indent=1, default=float))
    print("\nwrote sensitivity_2d.interaction_test")
else:
    st = S.get("interaction_test")
    if st:
        same = abs(st["p_value"] - p) < 1e-12 and abs(st["detectable_interaction_rms"] - mde) < 1e-12
        print("\nstored block matches this run" if same else "\nSTORED BLOCK DIFFERS")
        sys.exit(0 if same else 1)
