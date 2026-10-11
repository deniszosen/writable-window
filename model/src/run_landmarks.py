"""Landmark sections and the array-length loss, recomputed over several seeds.

Two corrections. The conservative reader now bounds an exposure with the last
day of the following pulse window rather than its first day, which is the valid
bound once pulses last more than a day. And the late-epoch loss caused by array
exhaustion is averaged over seeds, because a single realisation cannot resolve
it: records are clonally correlated, so the effective sample size for a
per-epoch rate is the number of independent writing events.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate, readout
import analysis as A

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000)
# The released results.json was produced at 5 for this script.
SEEDS = (0, 1, 2, 3, 4)
t0 = time.time()
def tick(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)
def _enc(o):
    if isinstance(o, (np.floating, np.integer)): o = o.item()
    if isinstance(o, float) and not np.isfinite(o): return None
    if isinstance(o, (np.bool_, bool)): return bool(o)
    if isinstance(o, np.ndarray): return [_enc(x) for x in o.tolist()]
    return float(o) if isinstance(o, (int, float)) else str(o)
def sem(v):
    v = np.asarray(v, float)
    return float(np.nanstd(v, ddof=1) / np.sqrt(np.sum(~np.isnan(v)))) if v.size > 1 else 0.0

R = json.load(open(os.path.join(OUT, "results.json")))

def bracket_over_seeds(over, reader, seeds=SEEDS, max_rows=5000):
    keys = ("acc_all", "acc_two_sided", "median_width_all", "median_width_two_sided",
            "frac_two_sided", "right_censored")
    acc = {k: [] for k in keys}
    for sd in seeds:
        p = Params(**{**BASE, **over, "seed": BASE["seed"] + sd})
        s = simulate(p)
        sym, got, wd, eid = readout(s, p)
        b = A.brackets(s, p, sym, got, wd, reader=reader, max_rows=max_rows)
        for k in keys:
            acc[k].append(b[k])
    out = {k: float(np.nanmean(acc[k])) for k in keys}
    out.update({k + "_sem": sem(acc[k]) for k in keys})
    return out

tick("landmark alphabet, both readers, 5 seeds")
rows = []
for k in (1, 2, 3, 4, 8):
    for reader in ("naive", "conservative"):
        rows.append(dict(alphabet=k, n_orthogonal_inducers=k, reader=reader,
                         **bracket_over_seeds(dict(fid_alphabet=k), reader)))
    tick(f"   k={k}: naive acc2 {rows[-2]['acc_two_sided']:.3f} w2 "
         f"{rows[-2]['median_width_two_sided']:.0f}d | conservative acc2 "
         f"{rows[-1]['acc_two_sided']:.3f}+-{rows[-1]['acc_two_sided_sem']:.3f} w2 "
         f"{rows[-1]['median_width_two_sided']:.0f}d")
R["landmark_alphabet"] = rows

tick("pulse duration; brackets only where pulses do not overlap the exposures")
rows = []
for dur in (1, 2, 4, 8, 16, 24):
    caps = []
    for sd in SEEDS:
        p = Params(**{**BASE, "fid_duration": dur, "seed": BASE["seed"] + sd})
        s = simulate(p)
        sym, got, wd, eid = readout(s, p)
        caps.append(A.pulse_capture(s, p, sym, got))
    # a pulse window that covers an exposure day makes the bracket degenerate
    fd = Params(**{**BASE, "fid_duration": dur}).fiducial_days()
    cover = {int(d) + k for d in fd for k in range(dur)}
    clean = not (cover & {15, 35, 55, 75})
    row = dict(duration=dur, pulse_capture=float(np.mean(caps)),
               pulse_capture_sem=sem(caps), pulses_clear_of_exposures=clean)
    if clean:
        row.update(bracket_over_seeds(dict(fid_duration=dur), "conservative"))
    rows.append(row)
    tick(f"   {dur:2d}d: capture {row['pulse_capture']:.3f}+-{row['pulse_capture_sem']:.3f}"
         + (f" two-sided {row['frac_two_sided']:.3f} w2 {row['median_width_two_sided']:.0f}d "
            f"acc2 {row['acc_two_sided']:.3f}" if clean else "  (pulses overlap an exposure window)"))
R["fiducial_duration"] = rows

tick("landmark channel leak")
rows = []
for Dl in (np.inf, 1000.0, 100.0):
    b = bracket_over_seeds(dict(fid_dyn_range=Dl), "conservative")
    rows.append(dict(fid_dyn_range=(None if Dl == np.inf else float(Dl)), **b))
    tick(f"   landmark D={Dl}: acc2 {b['acc_two_sided']:.3f}+-{b['acc_two_sided_sem']:.3f}")
R["fiducial_leak"] = rows

tick("array length: information and late-epoch loss over 5 seeds")
rows = []
for L in (5, 10, 12, 20, 50, 200):
    pe, sat, occ = [], [], []
    for sd in SEEDS:
        p = Params(**{**BASE, "L": L, "seed": BASE["seed"] + sd})
        s = simulate(p)
        pe.append(A.per_epoch_recording(s))
        sat.append(float((s["slot"] >= L).mean()))
        occ.append(float(s["slot"].mean()))
    pe = np.array(pe)
    c = A.contrast(Params(**{**BASE, "L": L}), seeds=(0, 1, 2))
    rows.append(dict(L=L, bits=c["bits"], bits_sem=c["bits_sem"], pi=c["pi"], phi=c["phi"],
                     per_epoch=[float(x) for x in pe.mean(0)],
                     per_epoch_sem=[sem(pe[:, i]) for i in range(pe.shape[1])],
                     frac_saturated=float(np.mean(sat)), occupancy=float(np.mean(occ))))
    tick(f"   L={L:3d}: bits {c['bits']:.4f}+-{c['bits_sem']:.4f} occ {rows[-1]['occupancy']:.2f} "
         f"saturated {rows[-1]['frac_saturated']:.3f} late epoch {pe.mean(0)[-1]:.4f}")
ref = rows[-1]
for r in rows:
    r["bits_rel_unbounded"] = float(r["bits"] / ref["bits"])
    r["late_epoch_loss"] = float(1 - r["per_epoch"][-1] / ref["per_epoch"][-1])
    r["late_epoch_loss_sem"] = float(np.hypot(r["per_epoch_sem"][-1], ref["per_epoch_sem"][-1])
                                     / ref["per_epoch"][-1])
    tick(f"   L={r['L']:3d}: information {r['bits_rel_unbounded']:.3f} of unbounded, "
         f"late-epoch loss {r['late_epoch_loss']:+.3f}+-{r['late_epoch_loss_sem']:.3f}")
R["array_length_multiseed"] = rows

with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=_enc)
tick("merged into out/results.json")
