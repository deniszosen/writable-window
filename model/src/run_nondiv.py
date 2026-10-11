"""Section J, recomputed with enough seeds to be quotable.

The early-to-late recording ratio is a ratio of two small probabilities, and
because records are clonally correlated its effective sample size is the number
of independent writing events, not the number of nuclei. One seed per point is
not enough. This recomputes the division-independent writing sweep over five
population realisations and reports a standard error for every quantity, then
merges the result into out/results.json.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate, readout, EXPOSURE_A
import analysis as A

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
    return float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else 0.0


R = json.load(open(os.path.join(OUT, "results.json")))
rows = []
for qn in [0.0, 0.0005, 0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.05]:
    pes, grads, hms, occ, pis, phis, bits, share, fons, reqD = [], [], [], [], [], [], [], [], [], []
    for sd in SEEDS:
        p = Params(**{**BASE, "q_nondiv": qn, "seed": BASE["seed"] + sd})
        s = simulate(p)
        sym, got, wd, eid = readout(s, p)
        pi = float((sym[got] == EXPOSURE_A).any(1).mean())
        phi = A.vehicle_mark_rate(p, seed=p.seed)
        pe = A.per_epoch_recording(s)
        m = A.write_day_density(s)
        f, Dr = A.required_dynamic_range(s, p, 0.1)
        pes.append(pe)
        grads.append(pe[0] / pe[-1] if pe[-1] > 0 else np.nan)
        hms.append(int(np.searchsorted(np.cumsum(m), 0.5)))
        occ.append(float(s["slot"].mean()))
        share.append(float(s["nondiv_writes"].sum() / max((s["wday"] >= 0).sum(), 1)))
        pis.append(pi); phis.append(phi)
        bits.append(float(A.mutual_information(pi, phi)))
        fons.append(f); reqD.append(Dr)
    pes = np.array(pes)
    rows.append(dict(
        q_nondiv=qn,
        per_epoch=[float(x) for x in pes.mean(0)],
        per_epoch_sem=[sem(pes[:, i]) for i in range(pes.shape[1])],
        gradient=float(np.nanmean(grads)), gradient_sem=sem(grads),
        gradient_seeds=[float(g) for g in grads],
        half_mass_day=float(np.mean(hms)), half_mass_day_sem=sem(hms),
        occupancy=float(np.mean(occ)),
        postmitotic_write_share=float(np.mean(share)),
        pi=float(np.mean(pis)), phi=float(np.mean(phis)),
        bits=float(np.mean(bits)), bits_sem=sem(bits),
        f_on=float(np.mean(fons)), required_D=float(np.mean(reqD)),
    ))
    r = rows[-1]
    tick(f"   q_nondiv={qn}: gradient {r['gradient']:.2f}+-{r['gradient_sem']:.2f} "
         f"half-mass d{r['half_mass_day']:.0f} occ {r['occupancy']:.2f} "
         f"postmitotic writes {r['postmitotic_write_share']*100:.0f}% "
         f"bits {r['bits']:.4f}+-{r['bits_sem']:.4f} reqD {r['required_D']:.0f}")

# where does the gradient fall below twofold
g = np.array([r["gradient"] for r in rows])
q = np.array([r["q_nondiv"] for r in rows])
below = np.flatnonzero(g < 2.0)
breakeven = float(q[below[0]]) if below.size else None

R["nondivision"] = dict(
    sweep=rows, seeds=list(SEEDS), gradient_breakeven_twofold=breakeven,
    anchors=dict(
        davis_2024_bulk_cortex_adult=dict(edited_fraction=0.14, days=21,
                                          per_day=float(1 - (1 - 0.14) ** (1 / 21))),
        davis_2024_neonatal_cortex=dict(edited_fraction=0.42, days=21,
                                        per_day=float(1 - (1 - 0.42) ** (1 / 21))),
        davis_2024_sorted_neurons=dict(edited_fraction=0.81, days=21,
                                       per_day=float(1 - (1 - 0.81) ** (1 / 21))),
    ))

# the headline gradient, with the same number of seeds
pe = []
for sd in SEEDS:
    s = simulate(Params(**{**BASE, "seed": BASE["seed"] + sd}))
    pe.append(A.per_epoch_recording(s))
pe = np.array(pe)
gr = pe[:, 0] / pe[:, -1]
R["writable"]["per_epoch_record"] = [float(x) for x in pe.mean(0)]
R["writable"]["per_epoch_sem"] = [sem(pe[:, i]) for i in range(pe.shape[1])]
R["writable"]["per_epoch_sd"] = [float(x) for x in pe.std(0, ddof=1)]
R["writable"]["gradient"] = float(gr.mean())
R["writable"]["gradient_sem"] = sem(gr)
R["writable"]["gradient_seeds"] = [float(x) for x in gr]
R["writable"]["n_seeds"] = len(SEEDS)
tick(f"   headline gradient {gr.mean():.2f} +- {sem(gr):.2f} (seeds "
     f"{np.round(gr,2)}), per-epoch {np.round(pe.mean(0),4)}")

with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=_enc)
tick("merged into out/results.json")
