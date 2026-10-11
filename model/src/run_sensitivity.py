"""Two analyses answering the interaction and benchmarking objections.

J. A two-dimensional sensitivity surface over the two parameters that rank
   highest in section 2.7, the writing rate and the fate-switch midpoint. Every
   other sweep in this paper moves one parameter at a time, which cannot show
   interaction. Ten populations per cell: a surface needs structure, not the
   precision the headline arms carry.

K. Writing moved off division entirely. Setting the per-division rate to zero
   leaves only the per-day rate, which this model applies to post-mitotic cells,
   so the writer records after cycle exit instead of at division. This is NOT
   the polymerase ticker tape of the Fisher-information analyses, which writes
   in every cell irrespective of cycle state and which this model cannot
   express without modification. What it does isolate is the coupling itself:
   the same population, the same schedule, the same readout, with only the
   event that triggers a write changed.
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate
import analysis as A

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
R = json.load(open(os.path.join(OUT, "results.json")))
BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000)
GRID_SEEDS = tuple(range(int(os.environ.get("N_GRID", "10"))))
SEEDS = tuple(range(int(os.environ.get("N_POPULATIONS", "50"))))
t0 = time.time()


def tick(m):
    print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)


def sem(v):
    v = np.asarray(v, float)
    return float(np.nanstd(v, ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0


def stats(seeds, **kw):
    g, occ = [], []
    for sd in seeds:
        s = simulate(Params(**{**BASE, **kw, "seed": BASE["seed"] + sd}))
        g.append(A.epoch_gradient(s))
        occ.append(float(s["slot"].mean()))
    return (float(np.mean(g)), sem(g), float(np.mean(occ)))


# ------------------------------------- J. two-dimensional sensitivity -------
tick(f"J. gradient over writing rate x fate-switch midpoint, "
     f"{len(GRID_SEEDS)} populations per cell")
qs = [0.06, 0.12, 0.18, 0.30, 0.50]
mids = [30.0, 40.0, 45.0, 55.0, 65.0]
cells = []
for mid in mids:
    row = []
    for q in qs:
        g, s_, _ = stats(GRID_SEEDS, q_div=q, fate_midpoint=mid)
        row.append(g)
        cells.append(dict(q_div=q, fate_midpoint=mid, gradient=g, gradient_sem=s_))
    tick("   midpoint %4.0f d: " % mid + "  ".join(f"{v:5.2f}" for v in row))
grid = np.array([[c["gradient"] for c in cells if c["fate_midpoint"] == m]
                 for m in mids], float)
# Is the surface separable? If the gradient were a product of a q effect and a
# midpoint effect there would be no interaction, and the rank-1 approximation
# would be exact.
u, sv, vt = np.linalg.svd(grid - grid.mean())
interaction = float(1.0 - sv[0] ** 2 / float(np.sum(sv ** 2)))
R["sensitivity_2d"] = dict(
    description="early-to-late gradient over a grid of writing rate and "
                "fate-switch midpoint, the two parameters section 2.7 ranks "
                "highest. residual_after_rank1 is the share of variance about "
                "the mean that a separable (non-interacting) surface leaves "
                "unexplained.",
    q_div=qs, fate_midpoint=mids, n_seeds=len(GRID_SEEDS),
    gradient=[[float(x) for x in r] for r in grid],
    gradient_min=float(grid.min()), gradient_max=float(grid.max()),
    residual_after_rank1=interaction, cells=cells)
tick(f"   range {grid.min():.2f} to {grid.max():.2f}; variance left by a "
     f"separable surface {interaction*100:.1f}%")

# ------------------------------------------- K. a time-coupled writer -------
tick(f"K. writing moved off division (q_div = 0), {len(SEEDS)} populations")
rows = []
for qn in (0.005, 0.01, 0.02, 0.04):
    g, s_, occ = stats(SEEDS, q_div=0.0, q_nondiv=qn)
    rows.append(dict(q_nondiv=qn, q_div=0.0, gradient=g, gradient_sem=s_,
                     occupancy=occ))
    tick(f"   {qn:.3f} per day: gradient {g:5.2f} +-{s_:.2f}, occupancy {occ:.2f}")
R["time_coupled"] = dict(
    description="the same model with writing moved off division entirely. The "
                "per-day rate applies to post-mitotic cells, so this writer "
                "records after cycle exit rather than at division; it is not "
                "the polymerase ticker tape of the Fisher-information "
                "analyses, which writes in every cell regardless of cycle "
                "state and which this model cannot express. It isolates the "
                "coupling: population, schedule and readout are unchanged and "
                "only the event that triggers a write differs.",
    division_coupled_gradient=R["multiplicity"]["reference_gradient"], sweep=rows)

# ------------------------------- M. the ticker-tape benchmark --------------
tick(f"M. a true time-coupled writer (q_time, every cell), {len(SEEDS)} populations")
rows_tt = []
for qt in (0.005, 0.01, 0.02, 0.04):
    g, s_, occ = stats(SEEDS, q_div=0.0, q_time=qt)
    rows_tt.append(dict(q_time=qt, gradient=g, gradient_sem=s_, occupancy=occ))
    tick(f"   {qt:.3f} per day in every cell: gradient {g:5.2f} +-{s_:.2f}, "
         f"occupancy {occ:.2f}")
R["ticker_tape"] = dict(
    description="writing at a constant rate per unit time in every cell, "
                "cycling or not, which is the regime the Fisher-information "
                "analyses of a polymerase ticker tape assume. Each terminal "
                "lineage occupies exactly one cell at every instant, so the "
                "expected number of ancestral write opportunities inside a "
                "window of fixed width is the same at every epoch and the "
                "gradient is 1 by construction. This is the benchmark the "
                "division-coupled 5.28 should be read against.",
    division_coupled_gradient=R["multiplicity"]["reference_gradient"],
    predicted_gradient=1.0, sweep=rows_tt)

# ------------------------------------------------- L. the null control -----
tick(f"L. a population that never exits the cycle, {len(SEEDS)} populations")
g, s_, occ = stats(SEEDS, pp_start=0.30, pp_end=0.30, nn_start=0.0, nn_end=0.0,
                   tc_start=30.0, tc_end=30.0)
R["no_exit_null"] = dict(
    description="every division keeps both daughters cycling and the cycle "
                "never lengthens, so nothing leaves the cycle. The population "
                "grows exponentially and every terminal lineage passes through "
                "the same number of divisions in any window, so both trends "
                "that produce front-loading vanish. This is the null the "
                "figure 2 axis draws at 1 and which no mechanism arm reaches.",
    gradient=g, gradient_sem=s_, occupancy=occ,
    reference_gradient=R["multiplicity"]["reference_gradient"])
tick(f"   gradient {g:.2f} +-{s_:.2f} against "
     f"{R['multiplicity']['reference_gradient']:.2f} with cycle exit")

with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=float)
tick("merged sensitivity_2d and time_coupled")
