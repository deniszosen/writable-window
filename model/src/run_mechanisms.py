"""The mechanism arms of section 2.2, and the corrected power calculation.

  A  The power calculation on the right footing. The between-culture variance
     is estimated from the same handful of cultures, so the contrast is a
     two-sample t rather than a z, and at one culture per arm there is no test
     at all. Both corrections move the recommendation.
  B  Readers richer than "does this array carry an exposure symbol", which
     measures how much the binary observable costs and whether a richer reader
     lifts the miscalling ceiling.
  C  Prime editor silencing over 150 days.
  D  Cell death, both uniform and biased towards early-born neurons.
  E  Transit-amplifying intermediate progenitors, matched on total expansion.
     This is a different object from the persistent compartment: that one
     self-renews, an intermediate progenitor amplifies and then exits.
  F  Position-dependent editing efficiency along the array.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate, readout, EXPOSURE_A
import analysis as A
import power as P

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000)
# Populations averaged over. The released results.json, and every
# 50-population quantity the manuscript reports from this script, were
# produced at 50. Override with N_POPULATIONS for a faster pass, which
# will NOT reproduce the published numbers.
SEEDS = tuple(range(int(os.environ.get("N_POPULATIONS", "50"))))
t0 = time.time()


def tick(m):
    print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)


def _enc(o):
    if isinstance(o, (np.floating, np.integer)):
        o = o.item()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return [_enc(x) for x in o.tolist()]
    return float(o) if isinstance(o, (int, float)) else str(o)


def sem(v):
    v = np.asarray(v, float)
    return float(np.nanstd(v, ddof=1) / np.sqrt(np.sum(~np.isnan(v)))) if v.size > 1 else 0.0


R = json.load(open(os.path.join(OUT, "results.json")))


def flush():
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=_enc)


def half_day(sim):
    """Day by which half of a lineage's writes have happened."""
    return int(np.argmax(np.cumsum(A.write_day_density(sim)) >= 0.5))


def pop_stats(**kw):
    """Gradient, half-capacity day, expansion and occupancy over seeds."""
    g, h, e, occ = [], [], [], []
    for sd in SEEDS:
        s = simulate(Params(**{**BASE, **kw, "seed": BASE["seed"] + sd}))
        g.append(A.epoch_gradient(s))
        h.append(half_day(s))
        e.append(A.expansion(s))
        occ.append(float(s["slot"].mean()))
    return dict(gradient=float(np.mean(g)), gradient_sem=sem(g),
                half_day=float(np.mean(h)), expansion=float(np.mean(e)),
                occupancy=float(np.mean(occ)))


# ------------------------------------------ A. the power calculation, fixed --
tick("A. power calculation with estimated between-culture variance")
pi = R["variance"]["between_culture"]["pi_mean"]
net = R["variance"]["between_culture"]["net_mean"]
m_cult = R["variance"]["between_culture"]["n_cultures"]
curves = []
for nf, sd_b in ((150, R["founder_scaling"]["sweep"][0]["pi_sd"]),
                 (600, R["founder_scaling"]["sweep"][2]["pi_sd"])):
    for n_arrays in (500, 1000, 5000):
        c = P.curve(sd_b, pi, net, n_arrays, m_cultures=m_cult)
        c["n_founder"] = nf
        curves.append(c)
p3 = [c for c in curves
      if c["n_founder"] == 150 and c["arrays_per_culture"] == 1000][0]
need = {f"{f:.2f}": P.cultures_needed(p3["sd_between"], pi, net, 1000, f,
                                      m_cultures=m_cult)
        for f in (1.2, 1.33, 1.5, 2.0)}
R["power_curves_t"] = dict(
    pi=pi, net=net, m_cultures=m_cult, curves=curves, cultures_needed=need,
    note="two-sided alpha 0.05, 80 per cent power, two-sample t on culture-level "
         "mark rates with 2(k-1) degrees of freedom, solved exactly against the "
         "noncentral t. One culture per arm has no test. The _sd_low and "
         "_sd_high fields carry the chi-squared interval on sigma_between.")
r3 = p3["points"][2]
tick(f"   150 founders, 1000 arrays: 3 cultures/arm detect "
     f"{r3['min_detectable_fold']:.2f}-fold exact against "
     f"{r3['min_detectable_fold_z']:.2f}-fold under the z form; "
     f"2 cultures/arm {p3['points'][1]['min_detectable_fold']:.2f}-fold; "
     f"1 culture/arm has no test")
tick(f"   sigma_between 90 per cent interval {p3['sd_between_ci'][0]:.4f} to "
     f"{p3['sd_between_ci'][1]:.4f} from {m_cult} cultures")
tick(f"   cultures/arm for 1.33-fold: {need['1.33']['point']} at the point "
     f"estimate, {need['1.33']['sigma_high']} at the pessimistic end")
flush()

# --------------------------------------- B. what a better reader buys -------
tick("B. binary, count and joint readers across dynamic range")
rows = []
for D in (100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0):
    ri = A.reader_information(Params(**{**BASE, "dynamic_range": D}), seeds=SEEDS)
    rows.append(dict(dynamic_range=D, **{k: ri[k] for k in ri}))
    tick(f"   D={D:7.0f}  binary {ri['binary']['bits']:.4f}  "
         f"count {ri['count']['bits']:.4f}  joint {ri['joint']['bits']:.4f}  "
         f"(permutation floor {ri['joint']['permutation_floor']:.5f})")
best = {k: max(r[k]["bits"] for r in rows) for k in ("binary", "count", "joint")}
R["readers"] = dict(
    sweep=rows, ceiling_bits=best,
    gain_joint_over_binary=float(best["joint"] / best["binary"]),
    note="Mutual information in bits per recovered cell about a binary exposure, "
         "for three observables of the same simulated arrays. Miller-Madow "
         "corrected; the permutation floor is the residual estimation bias after "
         "shuffling condition labels and is the scale a gain must be judged "
         "against.")
tick(f"   saturating information: binary {best['binary']:.4f}, "
     f"count {best['count']:.4f}, joint {best['joint']:.4f} bits "
     f"({best['joint']/best['binary']:.2f}-fold over binary)")
flush()

# ------------------------------------------------- C. editor silencing ------
tick("C. editor silencing")
rows = []
for hl in (np.inf, 100.0, 75.0, 50.0, 30.0):
    st = pop_stats(q_halflife=hl)
    rows.append(dict(q_halflife=(None if not np.isfinite(hl) else hl), **st))
    tick(f"   half-life {str(int(hl)) if np.isfinite(hl) else 'none':>4} d: "
         f"gradient {st['gradient']:5.2f} +- {st['gradient_sem']:.2f}, "
         f"half-capacity day {st['half_day']:.0f}, "
         f"occupancy {st['occupancy']:.2f}")
R["silencing"] = dict(sweep=rows,
                      note="writing rate multiplied by 2**(-t/half-life). The "
                           "direction is obvious; this is the size.")
flush()

# ----------------------------------------------------------- D. death ------
tick("D. death of post-mitotic cells")
rows = []
for dpm in (0.0, 0.002, 0.005, 0.010):
    st = pop_stats(death_pm=dpm)
    rows.append(dict(mode="uniform", death_pm=dpm, death_pm_early=0.0, **st))
    tick(f"   uniform {dpm:.3f}/d: gradient {st['gradient']:5.2f} "
         f"+- {st['gradient_sem']:.2f}, half-capacity day {st['half_day']:.0f}, "
         f"expansion {st['expansion']:.0f}")
for de in (0.005, 0.010, 0.020):
    st = pop_stats(death_pm_early=de, death_split=40)
    rows.append(dict(mode="early_biased", death_pm=0.0, death_pm_early=de,
                     death_split=40, **st))
    tick(f"   early-biased {de:.3f}/d before day 40: gradient {st['gradient']:5.2f} "
         f"+- {st['gradient_sem']:.2f}, half-capacity day {st['half_day']:.0f}, "
         f"expansion {st['expansion']:.0f}")
R["death"] = dict(sweep=rows,
                  note="per-day hazard on post-mitotic cells. 'early_biased' "
                       "applies only to cells that left the cycle before "
                       "death_split, the case whose direction is not obvious "
                       "in advance.")
flush()

# --------------------------------- E. transit-amplifying progenitors --------
tick("E. transit-amplifying intermediate progenitors, matched on expansion")
ref = pop_stats()
target = ref["expansion"]
tick(f"   reference expansion {target:.0f}-fold, gradient {ref['gradient']:.2f}")
rows = [dict(ipc_fraction=0.0, ipc_rounds=0, pp_start=Params().pp_start,
             matched=True, **ref)]
for frac, rounds in ((0.3, 2), (0.5, 2), (0.5, 3), (0.8, 2)):
    raw = pop_stats(ipc_fraction=frac, ipc_rounds=rounds)
    # an intermediate progenitor amplifies a division that already happened, so
    # it adds cells. Comparing its gradient against the central design at a
    # different total expansion compares two cultures that cannot both be the
    # calibrated one, so the proliferative fraction is reduced to restore it.
    pp_m, exp_m = A.match_expansion(
        Params(**{**BASE, "ipc_fraction": frac, "ipc_rounds": rounds}), target)
    st = pop_stats(ipc_fraction=frac, ipc_rounds=rounds, pp_start=pp_m)
    rows.append(dict(ipc_fraction=frac, ipc_rounds=rounds, pp_start=float(pp_m),
                     matched=True, expansion_unmatched=raw["expansion"],
                     gradient_unmatched=raw["gradient"], **st))
    tick(f"   fraction {frac:.1f}, {rounds} rounds: unmatched expansion "
         f"{raw['expansion']:.0f}-fold gradient {raw['gradient']:.2f}; "
         f"matched at p_PP {pp_m:.3f} expansion {st['expansion']:.0f}-fold "
         f"gradient {st['gradient']:5.2f} +- {st['gradient_sem']:.2f}, "
         f"half-capacity day {st['half_day']:.0f}")
R["ipc"] = dict(sweep=rows, target_expansion=float(target),
                note="a fraction of the daughters that would have left the cycle "
                     "instead divide ipc_rounds further times on a 0.6-fold "
                     "cycle. pp_start is bisected to restore the reference "
                     "expansion, because the compartment adds cells.")
flush()

# ------------------------- F. position-dependent editing efficiency --------
tick("F. position-dependent editing efficiency")
rows = []
for hl in (np.inf, 12.0, 6.0, 3.0):
    st = pop_stats(q_pos_halflife=hl)
    rows.append(dict(q_pos_halflife=(None if not np.isfinite(hl) else hl), **st))
    tick(f"   monomer half-life {str(int(hl)) if np.isfinite(hl) else 'none':>4}: "
         f"occupancy {st['occupancy']:5.2f}, gradient {st['gradient']:5.2f} "
         f"+- {st['gradient_sem']:.2f}")
# the interaction that matters is with array length, because a gradient in
# per-position efficiency is what would make extra length costly
length = []
for hl in (np.inf, 6.0):
    key = None if not np.isfinite(hl) else hl
    for L in (5, 10, 12, 20, 50):
        bits = []
        for sd in SEEDS:
            q = Params(**{**BASE, "L": L, "q_pos_halflife": hl,
                          "seed": BASE["seed"] + sd})
            bits.append(A.contrast(q, seeds=(0,))["bits"])
        length.append(dict(q_pos_halflife=key, L=L, bits=float(np.mean(bits)),
                           bits_sem=sem(bits)))
    sub = [r for r in length if r["q_pos_halflife"] == key]
    top = max(r["bits"] for r in sub)
    tick(f"   half-life {str(int(hl)) if np.isfinite(hl) else 'none':>4}: bits by "
         f"length " + ", ".join(f"L{r['L']}={r['bits']:.4f}" for r in sub)
         + f"  (best {top:.4f})")
R["position_efficiency"] = dict(
    sweep=rows, by_length=length,
    note="write probability at monomer s scaled by 2**(-s/half-life). One of "
         "the two mechanisms that penalise array length.")
flush()

tick("merged into out/results.json")
