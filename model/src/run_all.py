"""Produce every number the manuscript quotes, into out/results.json.

The central design uses the full readout model (symbol miscalls can be read as
other symbols; array contraction excises a monomer and shifts the remainder),
the vehicle culture is paired to the exposed one, and every headline quantity
carries a standard error over seeds. Three quantities that a naive analysis
leaves out are computed here: the number of independent writing events behind a
set of marked nuclei, the between-culture variance of the mark rate, and the
effect of division-independent writing.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate, readout, EXPOSURE_A, FID_BASE, CLOCK
import analysis as A

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
os.makedirs(OUT, exist_ok=True)
BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000)
SEEDS = (0, 1, 2)
HEAD_SEEDS = (0, 1, 2, 3, 4, 5)
NUCLEI = 15500          # nuclei assayed per condition
COMPETENT = 0.50        # share descended from recorder-carrying competent progenitors
R = {}
t0 = time.time()
# "The numbers reproduce" means nothing without an environment to be relative
# to, so the interpreter and library versions are recorded in the file.
import platform as _plat
import sys as _sys
R["env"] = dict(python=_sys.version.split()[0], numpy=np.__version__,
                platform=_plat.platform(), machine=_plat.machine())


def tick(msg):
    print(f"[{time.time()-t0:6.1f}s] {msg}", flush=True)


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


def msd(v):
    v = np.asarray(v, float)
    n = v.size
    return float(v.mean()), (float(v.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0), \
        (float(v.std(ddof=1)) if n > 1 else 0.0)


# ---------------------------------------------------------------- A. model --
tick("A. core model, validation, and seed-to-seed spread")
# FIX 7 WITHDRAWN. I briefly changed this to HEAD_SEEDS on the belief that
# writable's gradient block was computed from three seeds while verify.py
# indexed five. That was wrong: run_nondiv.py:108-114 OVERWRITES
# writable["gradient"], ["gradient_sem"], ["gradient_seeds"], ["n_seeds"],
# ["per_epoch_record"], ["per_epoch_sem"] and ["per_epoch_sd"] with five-seed
# values later in the documented run order. verify.py's five-index is therefore
# correct, and an IndexError here comes from running run_all.py in isolation,
# which is not the documented order. Changing this to six seeds would alter
# half_mass_day and the other population quantities that run_nondiv does NOT
# overwrite, for no reason. Left as it is.
sims = [simulate(Params(**{**BASE, "seed": BASE["seed"] + s})) for s in SEEDS]
p0 = Params(**BASE)
per_epoch = np.array([A.per_epoch_recording(s) for s in sims])
R["model"] = dict(
    params={k: (list(v) if isinstance(v, tuple) else (None if v == np.inf else v))
            for k, v in p0.__dict__.items()},
    expansion=[float(s["sym"].shape[0] * s["scale"] / p0.n_founder) for s in sims],
    cycling_fraction={str(d): float(np.mean([s["cyc_frac"][d] for s in sims]))
                      for d in (0, 10, 20, 30, 45, 60, 80, 100, 120, 149)},
    mean_divisions=[float(s["ndiv"].mean()) for s in sims],
    sd_divisions=[float(s["ndiv"].std()) for s in sims],
    occupancy_sim=[float(s["slot"].mean()) for s in sims],
    occupancy_analytic=[float(np.minimum(s["ndiv"] * p0.q_div, p0.L).mean()) for s in sims],
    fiducial_days=[int(d) for d in p0.fiducial_days()],
    alphabet_size=int(p0.alphabet().size),
    confusion_rate=float(A.confusion_rate(p0)),
    # Written here because a consumer reads it. A key that is read but never
    # written makes the released code unable to regenerate the released file,
    # which is the failure check_keys.py exists to catch.
    cycling_fraction_daily=[float(np.mean([s["cyc_frac"][d] for s in sims]))
                            for d in range(p0.T)],
)

wp = np.mean([A.writable_profile(s) for s in sims], axis=0)
mass = wp / wp.sum()
wdd = np.mean([A.write_day_density(s) for s in sims], axis=0)
hm = [int(np.searchsorted(np.cumsum(A.write_day_density(s)), 0.5)) for s in sims]
R["writable"] = dict(
    profile=[float(x) for x in wp], mass=[float(x) for x in mass],
    write_day_density=[float(x) for x in wdd],
    mass_first_30=float(mass[:30].sum()), mass_first_60=float(mass[:60].sum()),
    mass_first_90=float(mass[:90].sum()), mass_after_100=float(mass[100:].sum()),
    half_mass_day=int(np.searchsorted(np.cumsum(mass), 0.5)),
    half_mass_day_writes=int(np.searchsorted(np.cumsum(wdd), 0.5)),
    half_mass_day_seeds=hm,
    per_epoch_record=[float(x) for x in per_epoch.mean(0)],
    per_epoch_sd=[float(x) for x in per_epoch.std(0, ddof=1)],
    per_epoch_cv=[float(a / b) if b else np.nan
                  for a, b in zip(per_epoch.std(0, ddof=1), per_epoch.mean(0))],
    gradient=float(per_epoch.mean(0)[0] / per_epoch.mean(0)[-1]),
    gradient_seeds=[float(r[0] / r[-1]) for r in per_epoch],
    # FIX 5 WITHDRAWN for the same reason as FIX 7: gradient_sem, per_epoch_sem
    # and n_seeds are written into `writable` by run_nondiv.py:108-114, later in
    # the run order. They are not orphaned keys. Only model.cycling_fraction_daily
    # is (FIX 1), and that one is confirmed absent after the FULL pipeline.
    epochs=list(p0.epochs),
    f_on=float(np.mean([A.f_on(s, p0) for s in sims])),
    f_on_division_based=float(np.mean([A.f_on(s, p0, realised=False) for s in sims])),
)
tick(f"   half-mass day {R['writable']['half_mass_day']}, gradient "
     f"{R['writable']['gradient']:.2f} (seeds {np.round(R['writable']['gradient_seeds'],2)})")

# headline contrast with a standard error
c0 = A.contrast(p0, seeds=HEAD_SEEDS)
R["central_contrast"] = c0
tick(f"   central: pi {c0['pi']:.4f} phi {c0['phi']:.4f} net {c0['net']:.4f}+-{c0['net_sem']:.4f} "
     f"bits {c0['bits']:.4f}+-{c0['bits_sem']:.4f}")

# ------------------------------------------------------ B. dynamic range ----
tick("B. dynamic range sweep, both error models")
Ds = [4, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
sweep = []
for D in Ds:
    row = dict(D=D)
    for em in ("confusion", "sentinel"):
        p = Params(**{**BASE, "dynamic_range": float(D), "error_model": em})
        c = A.contrast(p, seeds=SEEDS)
        row[em] = dict(pi=c["pi"], phi=c["phi"], bits=c["bits"], bits_sem=c["bits_sem"],
                       ratio=(c["phi"] / c["pi"] if c["pi"] else np.nan))
    s = simulate(Params(**{**BASE, "dynamic_range": float(D)}))
    ph, fh = A.expected_marks(s, Params(**{**BASE, "dynamic_range": float(D)}))
    row["pi_hat"], row["phi_hat"] = ph, fh
    sweep.append(row)
    tick(f"   D={D}: confusion pi {row['confusion']['pi']:.3f} phi {row['confusion']['phi']:.3f} "
         f"bits {row['confusion']['bits']:.4f} | sentinel phi {row['sentinel']['phi']:.3f} "
         f"bits {row['sentinel']['bits']:.4f}")
f_on_, Dreq = A.required_dynamic_range(sims[0], p0, alpha=0.1)
R["dynamic_range"] = dict(
    sweep=sweep, f_on=f_on_,
    required={str(a): float(A.required_dynamic_range(sims[0], p0, a)[1])
              for a in (0.2, 0.1, 0.05)},
    ceiling=float(A.dynamic_range_ceiling(sims[0], p0)),
    ceiling_by_error={str(e): float(A.dynamic_range_ceiling(
        sims[0], Params(**{**BASE, "symbol_error": e})))
        for e in (0.02, 0.005, 0.002, 0.0005)},
    bits_ceiling_confusion=float(max(r["confusion"]["bits"] for r in sweep)),
    bits_ceiling_sentinel=float(max(r["sentinel"]["bits"] for r in sweep)),
)
tick(f"   required D(0.1) {Dreq:.0f}; miscalling ceiling {R['dynamic_range']['ceiling']:.0f}")

# ----------------------------------------------------- C. writing rate -----
tick("C. writing rate sweep")
rows = []
for q in [0.02, 0.05, 0.10, 0.18, 0.30, 0.45, 0.60, 0.80]:
    p = Params(**{**BASE, "q_div": q})
    c = A.contrast(p, seeds=SEEDS)
    s = simulate(p)
    rows.append(dict(q=q, pi=c["pi"], phi=c["phi"], bits=c["bits"], bits_sem=c["bits_sem"],
                     occupancy=float(s["slot"].mean()),
                     frac_saturated=float((s["slot"] >= p.L).mean()),
                     occ_p95=float(np.percentile(s["slot"], 95))))
    tick(f"   q={q}: bits {c['bits']:.4f} occ {rows[-1]['occupancy']:.1f}")
R["writing_rate"] = rows

tick("C2. array length requirement, and multiple writes per cycle")
rowsL = []
for q in [0.10, 0.18, 0.30, 0.45, 0.60]:
    for wpd in ("bernoulli", "poisson"):
        p = Params(**{**BASE, "q_div": q, "L": 200, "writes_per_div": wpd})
        s = simulate(p)
        need = s["slot"].astype(float)
        rowsL.append(dict(q=q, model=wpd, mean=float(need.mean()),
                          p95=float(np.percentile(need, 95)),
                          p99=float(np.percentile(need, 99)),
                          L_for_99pct=int(np.ceil(np.percentile(need, 99)))))
    tick(f"   q={q}: bernoulli p99 {rowsL[-2]['p99']:.0f}  poisson p99 {rowsL[-1]['p99']:.0f}")
R["array_length"] = rowsL

# ------------------------------------------------------- D. landmarks ------
tick("D. landmark alphabet, reader, and leak")
rows = []
for k in (0, 1, 2, 3, 4, 8):
    p = Params(**{**BASE, "fid_alphabet": k})
    s = simulate(p)
    sym, got, wd, eid = readout(s, p)
    for reader in ("naive", "conservative"):
        b = A.brackets(s, p, sym, got, wd, reader=reader, max_rows=6000)
        rows.append(dict(alphabet=k, n_orthogonal_inducers=(len(p.fiducial_days()) if k == 0 else k),
                         **b))
    tick(f"   alphabet k={k}: naive acc2 {rows[-2]['acc_two_sided']:.3f} w2 "
         f"{rows[-2]['median_width_two_sided']:.0f}d | conservative acc2 "
         f"{rows[-1]['acc_two_sided']:.3f} w2 {rows[-1]['median_width_two_sided']:.0f}d")
R["landmark_alphabet"] = rows

tick("D2. landmark count, both readers")
rows = []
for nf in [2, 4, 6, 8, 12, 16, 24]:
    for k in (0, 1):
        p = Params(**{**BASE, "n_fiducial": nf, "fid_alphabet": k})
        s = simulate(p)
        sym, got, wd, eid = readout(s, p)
        for reader in ("naive", "conservative"):
            rows.append(dict(n_fiducial=nf, alphabet=k, **A.brackets(
                s, p, sym, got, wd, reader=reader, max_rows=6000)))
R["fiducials"] = rows

tick("D3. pulse duration")
rows = []
for dur in [1, 2, 4, 8, 16, 24]:
    p = Params(**{**BASE, "fid_duration": dur})
    s = simulate(p)
    sym, got, wd, eid = readout(s, p)
    cap = A.pulse_capture(s, p, sym, got)
    b = A.brackets(s, p, sym, got, wd, reader="conservative", max_rows=6000)
    rows.append(dict(duration=dur, pulse_capture=cap,
                     occupancy=float(s["slot"].mean()), **b))
    tick(f"   dur={dur}d: capture {cap:.3f} two-sided {b['frac_two_sided']:.3f} "
         f"w2 {b['median_width_two_sided']:.0f}d")
R["fiducial_duration"] = rows

tick("D4. landmark channel leak")
rows = []
for Dl in (np.inf, 10000.0, 3000.0, 1000.0, 300.0, 100.0):
    p = Params(**{**BASE, "fid_dyn_range": Dl})
    s = simulate(p)
    sym, got, wd, eid = readout(s, p)
    b = A.brackets(s, p, sym, got, wd, reader="conservative", max_rows=6000)
    rows.append(dict(fid_dyn_range=(None if Dl == np.inf else float(Dl)), **b))
    tick(f"   landmark D={Dl}: conservative acc_two_sided {b['acc_two_sided']:.3f} "
         f"acc_all {b['acc_all']:.3f}")
R["fiducial_leak"] = rows

s0 = sims[0]
tc_days = float(np.mean([30.0, 110.0])) / 24.0
s_f = p0.fid_share / (p0.clock_share + p0.fid_share + p0.exp_share / p0.dynamic_range)
R["pulse_rule"] = dict(
    s_f=float(s_f), q=float(p0.q_div), tc_days_mean=tc_days,
    divisions_for_90pct=float(np.log(0.1) / np.log(1 - p0.q_div * s_f)),
    days_for_90pct=float(np.log(0.1) / np.log(1 - p0.q_div * s_f) * tc_days))

# ------------------------------------------------------- E. order ----------
tick("E. order accuracy, both contraction models")
rows = []
for dt in [3, 7, 14, 21, 35, 49]:
    r = dict(delta=dt)
    for cm in ("delete", "erase"):
        acc, n = A.order_accuracy(Params(**{**BASE, "contraction_mode": cm}),
                                  day0=20, delta=dt, seeds=(0, 1))
        r[cm] = dict(accuracy=acc, n=n)
    rows.append(r)
    tick(f"   dt={dt}: delete {r['delete']['accuracy']:.3f} (n {r['delete']['n']}) "
         f"erase {r['erase']['accuracy']:.3f}")
R["order"] = rows

# --------------------------------------------- F. clock variability --------
tick("F. clock variability between landmarks")
R["clock"] = A.clock_variability(sims[0], p0)

# ------------------------------- G. independence, variance, sample size ----
tick("G. independence of records, between-culture variance, sample size")
ev = []
for s in sims:
    p = s["params"]
    sym, got, wd, eid = readout(s, p)
    ev.append(A.event_stats(s, p, sym, got, eid))
R["independence"] = dict(
    per_seed=ev,
    pseudoreplication=float(np.mean([e["pseudoreplication"] for e in ev])),
    n_marked=float(np.mean([e["n_marked"] for e in ev])),
    n_events=float(np.mean([e["n_events"] for e in ev])),
    top_decile_share=float(np.mean([e["top_decile_share"] for e in ev])),
    max_clone=float(np.mean([e["max_clone"] for e in ev])),
)
tick(f"   {R['independence']['n_marked']:.0f} marked nuclei carry "
     f"{R['independence']['n_events']:.0f} independent events "
     f"({R['independence']['pseudoreplication']:.1f}x)")

bc = A.between_culture_variance(p0, seeds=range(8))
# The default cap of 40,000 subsamples the tracked population, and
# Methods 4.1 already warns that subsampling "thins subclones, so it
# attenuates estimates of clonal correlation". run_extra.py's founder-scaling
# sweep raises the cap to 400,000 for exactly that reason. This calculation -
# the one that feeds the headline sample-size recommendation - did not, so the
# variance driving "five cultures per arm" was inflated by a simulation
# artefact. Both are now computed and reported side by side.
bc_uncapped = A.between_culture_variance(
    Params(**{**BASE, "cap": 400000}), seeds=range(8))
tick(f"   sigma at cap 40,000: {bc['pi_sd']:.5f};  at cap 400,000: "
     f"{bc_uncapped['pi_sd']:.5f}")
n_rec = NUCLEI * p0.recovery * COMPETENT
R["variance"] = dict(
    between_culture=bc,
    between_culture_cap400k=bc_uncapped,
    cap_note="between_culture uses the default cap of 40,000; "
             "between_culture_cap400k raises the cap to the "
             "value run_extra.py uses for the founder-scaling sweep. The "
             "second is the one to plan with, because subsampling at the "
             "default cap thins the subclones the variance depends on.",
    n_recovered_per_condition=float(n_rec),
    n_recovered_components=dict(
        nuclei=NUCLEI, recovery=float(p0.recovery), competent=COMPETENT,
        note="COMPETENT = 0.50 is the share descended from recorder-carrying "
             "competent progenitors. It is a third multiplier on top of "
             "nuclei and recovery, and is stated here so that a figure like "
             "'15,500 nuclei at 60 per cent recovery give about 4,650' "
             "cannot be read as the whole calculation."),
    design_effect=float(A.design_effect(bc["pi_mean"], bc["pi_var"], n_rec)),
    design_effect_by_n={str(int(n)): float(A.design_effect(bc["pi_mean"], bc["pi_var"], n))
                        for n in (500, 1000, 2000, 5000, 10000)},
)
tick(f"   between-culture SD of pi {bc['pi_sd']:.4f}; design effect at n={n_rec:.0f} "
     f"is {R['variance']['design_effect']:.1f}")

pi, phi = c0["pi"], c0["phi"]
net = pi - phi
folds = (2.0, 1.5, 1.25)
R["power"] = dict(
    pi=pi, phi=phi, net=net,
    n_needed={str(f): dict(
        binomial=float(A.n_per_group(pi, phi + net / f)),
        with_design_effect=float(A.n_per_group(
            pi, phi + net / f, deff=R["variance"]["design_effect"])))
        for f in folds},
    n_needed_presence=float(A.n_per_group(pi, phi)),
    cultures_needed={str(f): float(A.n_cultures_needed(net, bc["net_sd"],
                                                       effect_fraction=1 - 1 / f))
                     for f in folds},
    legacy_usable_cells={str(r): float(A.usable_cells(net, r, NUCLEI, COMPETENT))
                         for r in (0.10, 0.25, 0.40, 0.60, 0.80)},
)
for f in folds:
    tick(f"   detect {f}x change in mark rate: {R['power']['n_needed'][str(f)]['binomial']:.0f} "
         f"nuclei per arm (binomial), "
         f"{R['power']['n_needed'][str(f)]['with_design_effect']:.0f} with the design effect, "
         f"{R['power']['cultures_needed'][str(f)]:.1f} cultures per arm")

# --------------------------------------------- H. feasibility surface ------
tick("H. feasibility over writing rate and recovery")
qs = [r["q"] for r in R["writing_rate"]]
needed = [float(A.n_per_group(r["pi"], r["phi"] + max(r["pi"] - r["phi"], 1e-9) / 1.5))
          for r in R["writing_rate"]]
recs = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
R["feasibility"] = dict(
    q=qs, recovery=recs, n_needed_15fold=needed,
    n_available=[float(NUCLEI * r * COMPETENT) for r in recs],
    nuclei_per_condition=NUCLEI, competent=COMPETENT,
    feasible=[[bool(NUCLEI * r * COMPETENT >= n) for r in recs] for n in needed],
)

# --------------------------------------------- I. robustness ---------------
tick("I. robustness to the population model and to the readout model")
variants = {
    "central": {},
    "slower cycle": dict(tc_start=40.0, tc_end=140.0),
    "faster cycle": dict(tc_start=22.0, tc_end=80.0),
    "later fate switch": dict(fate_midpoint=65.0),
    "earlier fate switch": dict(fate_midpoint=30.0),
    "shallower switch": dict(fate_steep=0.035),
    "more proliferative": dict(pp_start=0.38, nn_start=0.04),
}
rows = []
for name, over in variants.items():
    p = Params(**{**BASE, **over})
    s = simulate(p)
    c = A.contrast(p, seeds=SEEDS)
    m = A.write_day_density(s)
    sym, got, wd, eid = readout(s, p)
    b = A.brackets(s, p, sym, got, wd, reader="conservative", max_rows=5000)
    f, Dr = A.required_dynamic_range(s, p, 0.1)
    rows.append(dict(variant=name, mean_div=float(s["ndiv"].mean()),
                     occupancy=float(s["slot"].mean()),
                     mass_first_60=float(m[:60].sum()),
                     half_mass_day=int(np.searchsorted(np.cumsum(m), .5)),
                     f_on=f, required_D=Dr, pi=c["pi"], phi=c["phi"], bits=c["bits"],
                     frac_two_sided=b["frac_two_sided"],
                     per_epoch=[float(x) for x in A.per_epoch_recording(s)],
                     gradient=float(A.epoch_gradient(s))))
    tick(f"   {name}: half-mass d{rows[-1]['half_mass_day']} reqD {Dr:.0f} "
         f"bits {c['bits']:.4f} gradient {rows[-1]['gradient']:.2f}")
R["robustness"] = rows

tick("I2. robustness to the model's structural switches")
switch_rows = []
for name, over in {
        "central": {},
        "sentinel errors, erasing contraction": dict(
            error_model="sentinel", contraction_mode="erase"),
        "edit after replication (one daughter)": dict(write_timing="post"),
        "multiple writes per cycle": dict(writes_per_div="poisson"),
        "founders 150, no subsampling": dict(n_founder=150, cap=400000),
        "founders 6000": dict(n_founder=6000),
}.items():
    p = Params(**{**BASE, **over})
    s = simulate(p)
    c = A.contrast(p, seeds=SEEDS)
    m = A.write_day_density(s)
    sym, got, wd, eid = readout(s, p)
    es = A.event_stats(s, p, sym, got, eid)
    switch_rows.append(dict(variant=name, pi=c["pi"], phi=c["phi"], bits=c["bits"],
                            mean_div=float(s["ndiv"].mean()),
                            occupancy=float(s["slot"].mean()),
                            half_mass_day=int(np.searchsorted(np.cumsum(m), .5)),
                            gradient=float(A.epoch_gradient(s)),
                            per_epoch=[float(x) for x in A.per_epoch_recording(s)],
                            pseudoreplication=es["pseudoreplication"],
                            n_events=es["n_events"], n_marked=es["n_marked"]))
    tick(f"   {name}: bits {c['bits']:.4f} half-mass d{switch_rows[-1]['half_mass_day']} "
         f"gradient {switch_rows[-1]['gradient']:.2f} pseudorep "
         f"{switch_rows[-1]['pseudoreplication']:.1f}x")
R["switches"] = switch_rows

# ------------------------- J. division-independent writing (main result) ---
tick("J. division-independent writing")
rows = []
for qn in [0.0, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.007, 0.01, 0.02]:
    p = Params(**{**BASE, "q_nondiv": qn})
    s = simulate(p)
    c = A.contrast(p, seeds=SEEDS)
    m = A.write_day_density(s)
    f, Dr = A.required_dynamic_range(s, p, 0.1)
    pe = A.per_epoch_recording(s)
    rows.append(dict(q_nondiv=qn, half_mass_day=int(np.searchsorted(np.cumsum(m), .5)),
                     mass_first_60=float(m[:60].sum()),
                     mass_after_100=float(m[100:].sum()),
                     gradient=float(pe[0] / pe[-1]) if pe[-1] > 0 else None,
                     per_epoch=[float(x) for x in pe],
                     occupancy=float(s["slot"].mean()),
                     postmitotic_write_share=float(s["nondiv_writes"].sum() /
                                                   max((s["wday"] >= 0).sum(), 1)),
                     f_on=f, required_D=Dr, pi=c["pi"], phi=c["phi"], bits=c["bits"]))
    tick(f"   q_nondiv={qn}: half-mass d{rows[-1]['half_mass_day']} "
         f"gradient {rows[-1]['gradient']:.2f} occ {rows[-1]['occupancy']:.1f} "
         f"reqD {Dr:.0f} bits {c['bits']:.4f}")
R["nondivision"] = dict(
    sweep=rows,
    anchors=dict(
        davis_2024_bulk_cortex_adult=dict(edited_fraction=0.14, days=21,
                                          per_day=float(1 - (1 - 0.14) ** (1 / 21))),
        davis_2024_neonatal_cortex=dict(edited_fraction=0.42, days=21,
                                        per_day=float(1 - (1 - 0.42) ** (1 / 21))),
    ))

# ------------------------------- K. exact lineage tree, clone sizes --------
# Six seeds, not one. A single realisation of this block gives 34,909 nuclei,
# 2,908 events, a twelvefold factor and an 800-nucleus maximum subclone, and
# plotting that maximum against a measured 801 would look like validation when
# it is coincidence. Every headline quantity carries a standard error. The
# maximum of a
# heavy-tailed distribution at n = 1 is the least stable statistic available;
# over six seeds it ranges 662 to 1,755. Now looped, with mean and SEM.
tick("K. exact lineage tree without subsampling, over six seeds")
_ex_rows = []
for _sd in HEAD_SEEDS:
    p_ex = Params(**{**BASE, "n_founder": 150, "cap": 10 ** 9,
                     "seed": BASE["seed"] + _sd})
    s_ex = simulate(p_ex)
    sym, got, wd, eid = readout(s_ex, p_ex)
    es = A.event_stats(s_ex, p_ex, sym, got, eid)
    is_exp = sym == EXPOSURE_A
    ids = eid[got & is_exp.any(1)][is_exp[got & is_exp.any(1)]]
    ids = ids[ids > 0]
    uq, ct = np.unique(ids, return_counts=True)
    _ex_rows.append(dict(
        seed=int(p_ex.seed),
        n_founder=int(p_ex.n_founder), n_cells=int(s_ex["sym"].shape[0]),
        expansion=float(s_ex["sym"].shape[0] / p_ex.n_founder),
        n_founders_surviving=int(np.unique(s_ex["founder"]).size),
        **{k: v for k, v in es.items()},
        clone_size_quantiles={str(q): float(np.percentile(ct, q))
                              for q in (50, 75, 90, 99)}))
    tick(f"   seed {p_ex.seed}: {_ex_rows[-1]['n_cells']} cells, "
         f"{es['n_marked']} marked, {es['n_events']} events, "
         f"{es['pseudoreplication']:.2f}x, max clone {es['max_clone']}")

def _agg(key):
    v = np.asarray([r[key] for r in _ex_rows], float)
    m, se, sd = msd(v)
    return dict(mean=m, sem=se, sd=sd, per_seed=[float(x) for x in v])

R["exact_tree"] = dict(
    per_seed=_ex_rows, n_seeds=len(_ex_rows),
    n_founder=150,
    **{k: _agg(k) for k in ("n_cells", "expansion", "n_marked", "n_events",
                            "pseudoreplication", "max_clone", "median_clone",
                            "mean_clone", "top_decile_share")},
    note="Six seeds. The pseudoreplication factor is stable across "
         "them; the maximum subclone is not, and should not be compared "
         "against a single measured clone size.",
)
tick(f"   {R['exact_tree']['n_cells']} cells from {p_ex.n_founder} founders; "
     f"{es['n_marked']} marked nuclei carry {es['n_events']} events "
     f"({es['pseudoreplication']:.1f}x); largest subclone {es['max_clone']}")

with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=_enc)
tick(f"wrote {os.path.join(OUT,'results.json')}")
