"""Daily cell number and cycling fraction, extracted from the simulator.

results.json stores only the endpoint expansion (~1,650-fold), not a daily
population curve, so figure 1 cannot be drawn from it. This runs the released
simulator and records, for every day, the number of tracked cells (scaled back
up through the subsampling factor) and how many of them are still cycling.
No new model and no new parameters: Params() defaults are the central design.

Writes growth_trajectory.json beside results.json.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
# The simulator sits at model/src in the deposit and at release/model/src in
# the authoring tree. Pointing outside either would make the released pipeline
# unrunnable for anyone who downloads it, which is how it was broken before.
for _cand in (HERE.parent / "model" / "src",
              HERE.parent / "release" / "model" / "src"):
    if (_cand / "recorder.py").exists():
        SRC = _cand
        break
else:
    sys.exit("cannot find recorder.py in model/src beside this directory")
sys.path.insert(0, str(SRC))
import recorder as RC                                               # noqa: E402

SEEDS = (0, 1, 2)


def run(seed):
    """Call the released simulator and read its own per-day series.

    simulate() records cyc_frac (cycling fraction of the tracked set) and pool
    (scale * number of cycling cells) for every day, so the total tracked
    population is pool / cyc_frac. Reimplementing the division loop by hand
    got the expansion wrong, because asymmetric divisions also add a daughter.
    """
    sim = RC.simulate(RC.Params(seed=seed))
    cf = np.asarray(sim["cyc_frac"], float)
    pool = np.asarray(sim["pool"], float)
    tot = np.where(cf > 0, pool / np.where(cf > 0, cf, 1.0), np.nan)
    return tot, pool


def main():
    T = RC.Params().T
    tots, cycs = [], []
    for s in SEEDS:
        a, b = run(s)
        tots.append(a); cycs.append(b)
    tot = np.mean(tots, axis=0); cycg = np.mean(cycs, axis=0)
    out = dict(
        day=list(range(T)),
        n_cells=[float(x) for x in tot],
        n_cycling=[float(x) for x in cycg],
        cycling_fraction=[float(c / t) if t else 0.0 for c, t in zip(cycg, tot)],
        n_founder=RC.Params().n_founder,
        expansion=float(tot[-1] / tot[0]),
        seeds=list(SEEDS),
        note=("Daily totals from the released simulator at Params() defaults, "
              "averaged over three seeds. Counts are scaled back up through the "
              "subsampling factor, so they are the tracked lineage population, "
              "not a separate model."),
    )
    (HERE / "growth_trajectory.json").write_text(json.dumps(out, indent=1))
    print(f"day 0: {tot[0]:,.0f} cells, {out['cycling_fraction'][0]:.0%} cycling")
    print(f"day 150: {tot[-1]:,.0f} cells, {out['cycling_fraction'][-1]:.2%} cycling")
    print(f"expansion {out['expansion']:.0f}-fold "
          f"(results.json endpoint: ~1,650-fold)")


if __name__ == "__main__":
    main()
