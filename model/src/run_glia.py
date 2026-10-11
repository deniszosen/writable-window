"""Three analyses added in response to review.

G. A dividing glial compartment, in two forms. Section 2.2 had tested a
   persistent compartment that does not expand and transit-amplifying
   progenitors that amplify a division already taken; neither is gliogenesis.
   In a 150-day culture the neurogenic-to-gliogenic switch occurs and
   astrocyte precursors and OPCs divide, which is the remaining route by which
   late recording capacity could be restored. The arms reported in the paper
   switch the compartment on at day 80, 100 or 120. A variant that cycles from
   t = 0 is kept as a bound on the proxy: it blunts the window MORE, because by
   day 150 it has had the whole culture to expand in, so it is the optimistic
   bound on gliogenesis rather than the pessimistic one.

H. Heterogeneity in fate-switch timing between lineages. This is a spread
   over lineages, not over spatial domains: a region is a contiguous patch in
   which many lineages share a schedule, which a well-mixed model cannot
   represent. The gradient was never reported across this sweep, only
   clone-size statistics.

I. The confound, decomposed. The recorded per-epoch signal is the product of a
   true biological susceptibility and the probability that an exposure at that
   epoch is recorded at all. Section 2.1 gives the second factor; this block
   puts a true susceptibility through it, which is what makes the confound a
   demonstration rather than an inference. Closed form from the released
   per-epoch recording probabilities, no new simulation.
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate
import analysis as A

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
R = json.load(open(os.path.join(OUT, "results.json")))
BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000)
SEEDS = tuple(range(int(os.environ.get("N_POPULATIONS", "50"))))
t0 = time.time()


def tick(m):
    print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)


def sem(v):
    v = np.asarray(v, float)
    return float(np.nanstd(v, ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0


def half_day(sim):
    return int(np.argmax(np.cumsum(A.write_day_density(sim)) >= 0.5))


def pop_stats(**kw):
    g, h, e, late_cells, late_writes = [], [], [], [], []
    for sd in SEEDS:
        s = simulate(Params(**{**BASE, **kw, "seed": BASE["seed"] + sd}))
        g.append(A.epoch_gradient(s))
        h.append(half_day(s))
        e.append(A.expansion(s))
        lf = s.get("late")
        if lf is not None:
            late_cells.append(float(np.mean(lf)))
            w = s["slot"]
            tot = float(np.sum(w))
            late_writes.append(float(np.sum(w[lf]) / tot) if tot > 0 else 0.0)
    out = dict(gradient=float(np.mean(g)), gradient_sem=sem(g),
               half_day=float(np.mean(h)), half_day_sem=sem(h),
               expansion=float(np.mean(e)), n_seeds=len(SEEDS))
    if late_cells:
        out["late_share_endpoint"] = float(np.mean(late_cells))
        out["late_share_writes"] = float(np.mean(late_writes))
    return out


# ---------------------------------------------- G. glial compartment --------
tick("G. proliferative glial compartment")
ref = R["multiplicity"]["reference_gradient"]
target = float(np.mean(R["model"]["expansion"]))
rows = []
# Proliferation is pushed until the compartment either moves the gradient or
# takes total expansion outside the calibrated range, which is the same bound
# section 2.2 applies to post-mitotic death and to a late-expanding compartment.
for frac, pp, nn, tcf in ((0.20, 0.06, 0.02, 2.0),
                          (0.20, 0.12, 0.02, 2.0),
                          (0.20, 0.20, 0.02, 2.0),
                          (0.20, 0.30, 0.02, 2.0),
                          (0.40, 0.20, 0.02, 2.0)):
    st = pop_stats(late_fraction=frac, late_pp=pp, late_nn=nn,
                   late_tc_factor=tcf, late_mid_shift=0.0)
    st.update(late_fraction=frac, late_pp=pp, late_nn=nn, late_tc_factor=tcf,
              expansion_ratio=st["expansion"] / target)
    rows.append(st)
    tick(f"   frac {frac:.2f} pp {pp:.2f}: gradient {st['gradient']:5.2f} "
         f"+-{st['gradient_sem']:.2f}, expansion {st['expansion']:7.0f} "
         f"({st['expansion_ratio']:.2f}x reference), glia are "
         f"{st['late_share_endpoint']*100:4.1f}% of endpoint cells and "
         f"{st['late_share_writes']*100:4.1f}% of writes")
# The same compartment, but switched on at the gliogenic transition instead of
# cycling from the start. This is the arm the second round of review asked for:
# a pool that writes early props up the numerator of the very ratio being
# tested, so an always-on compartment bounds the mechanism from the wrong side.
rows_sw = []
for start in (80.0, 100.0, 120.0):
    st = pop_stats(late_fraction=0.20, late_pp=0.30, late_nn=0.02,
                   late_tc_factor=2.0, late_start_day=start)
    st.update(late_fraction=0.20, late_pp=0.30, late_start_day=start,
              expansion_ratio=st["expansion"] / target)
    rows_sw.append(st)
    tick(f"   switch on day {start:.0f}: gradient {st['gradient']:5.2f} "
         f"+-{st['gradient_sem']:.2f}, expansion {st['expansion']:7.0f} "
         f"({st['expansion_ratio']:.2f}x), glia are "
         f"{st['late_share_endpoint']*100:4.1f}% of endpoint cells and "
         f"{st['late_share_writes']*100:4.1f}% of writes")
R["gliogenesis_switched"] = dict(
    description="the glial compartment switched on at a gliogenic transition "
                "rather than cycling from t = 0. Before the switch day the "
                "flagged cells are indistinguishable from the main population, "
                "so the compartment contributes no early writing and the arm "
                "bounds the mechanism from the correct side.",
    reference_gradient=ref, reference_expansion=target, sweep=rows_sw)

R["gliogenesis"] = dict(
    description="a slowly cycling compartment with net-positive proliferation, "
                "standing in for astrocyte precursors and OPCs. The compartment "
                "cycles throughout rather than switching on at the gliogenic "
                "transition, which the model cannot express; it therefore writes "
                "early as well as late and understates how far a purely late "
                "glial wave would blunt the gradient.",
    reference_gradient=ref, reference_expansion=target, sweep=rows)

# ------------------------------- H. fate-timing (regional) heterogeneity ----
tick("H. heterogeneity in fate-switch timing between lineages")
rows = []
for sd_days in (0.0, 5.0, 10.0, 15.0, 20.0):
    st = pop_stats(fate_mid_sd=sd_days)
    st["fate_mid_sd"] = sd_days
    rows.append(st)
    tick(f"   sd {sd_days:4.0f} d: gradient {st['gradient']:5.2f} "
         f"+-{st['gradient_sem']:.2f}, half-capacity day {st['half_day']:.0f}, "
         f"expansion {st['expansion']:7.0f}")
R["fate_heterogeneity"] = dict(
    description="per-lineage spread on the fate-switch midpoint, standing in "
                "for regional heterogeneity between domains on different "
                "schedules.", reference_gradient=ref, sweep=rows)

# ------------------------------------------- I. the confound, decomposed ----
tick("I. decomposition of the recorded signal")
p_rec = np.asarray(R["writable"]["per_epoch_record"], float)
epochs = list(R["writable"]["epochs"])
rel = p_rec / p_rec[0]
scen = {}
# a true susceptibility profile, relative to the first epoch
profiles = {
    "no biological window": np.ones(4),
    "true early window, 3-fold": np.array([3.0, 1.0, 1.0, 1.0]),
    "true late window, 3-fold": np.array([1.0, 1.0, 1.0, 3.0]),
    "true late window, 5.3-fold": np.array([1.0, 1.0, 1.0, float(ref)]),
}
for name, s in profiles.items():
    obs = s * p_rec                       # what the recorder returns
    obs_rel = obs / obs[0]
    corrected = obs / p_rec               # divide out the writable population
    scen[name] = dict(
        true_relative=[float(x) for x in s / s[0]],
        observed=[float(x) for x in obs],
        observed_relative=[float(x) for x in obs_rel],
        observed_early_to_late=float(obs[0] / obs[3]),
        corrected_relative=[float(x) for x in corrected / corrected[0]],
        recovers_truth=bool(np.allclose(corrected / corrected[0], s / s[0], rtol=1e-9)))
    tick(f"   {name:28s} -> recorder reports early:late = "
         f"{scen[name]['observed_early_to_late']:5.2f}, "
         f"correction recovers truth: {scen[name]['recovers_truth']}")
R["confound_decomposition"] = dict(
    description="recorded signal = true susceptibility x per-epoch recording "
                "probability. The second factor is the writable window.",
    epochs=epochs, p_record=[float(x) for x in p_rec],
    p_record_relative=[float(x) for x in rel], scenarios=scen)

with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=float)
tick("merged gliogenesis, fate_heterogeneity and confound_decomposition")
