"""Compartments, decoding and the readout diagnostics.

Six blocks, each answering a question the headline sweeps cannot:

  A  A persistent, asymmetrically dividing progenitor compartment (outer radial
     glia and similar late-cycling progenitors). Does it rescue late recording
     capacity? It does not, and the reason is structural.
  B  How the miscalling ceiling depends on the size of the symbol alphabet.
  C  Sensitivity, specificity, predictive value and likelihood ratio across
     dynamic range, for readers who do not work in bits.
  D  The distribution of landmark bracket widths and the three-way split of
     records into two-sided, lower-bound-only and unanchored.
  E  The full distribution of nuclei per independent writing event.
  F  Detectable fold change against the number of cultures, which is the form
     of the power calculation an experimenter actually plans with.
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
    return float(np.nanstd(v, ddof=1) / np.sqrt(np.sum(~np.isnan(v)))) if v.size > 1 else 0.0


R = json.load(open(os.path.join(OUT, "results.json")))


def flush():
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=_enc)


# ---------------------------------- A. a persistent progenitor compartment --
tick("A. persistent asymmetrically dividing compartment")
rows = []
variants = [("none", {}),
            ("5 per cent", dict(late_fraction=0.05)),
            ("10 per cent", dict(late_fraction=0.10)),
            ("20 per cent", dict(late_fraction=0.20)),
            ("35 per cent", dict(late_fraction=0.35)),
            ("20 per cent, cycling at the main rate",
             dict(late_fraction=0.20, late_tc_factor=1.0)),
            ("20 per cent, 10 per cent symmetric proliferative",
             dict(late_fraction=0.20, late_pp=0.10))]
for name, over in variants:
    grads, hms, exps, shares, late_writes, divs, m100 = [], [], [], [], [], [], []
    for sd in SEEDS:
        p = Params(**{**BASE, **over, "seed": BASE["seed"] + sd})
        s = simulate(p)
        pe = A.per_epoch_recording(s)
        m = A.write_day_density(s)
        grads.append(float(pe[0] / pe[-1]) if pe[-1] > 0 else np.nan)
        hms.append(int(np.searchsorted(np.cumsum(m), 0.5)))
        exps.append(float(s["sym"].shape[0] * s["scale"] / p.n_founder))
        divs.append(float(s["ndiv"].mean()))
        m100.append(float(m[100:].sum()))
        if s["late"] is not None:
            shares.append(float(s["late"].mean()))
            # share of all writes that happened in the late compartment
            wl = (s["wday"] >= 0) & s["late"][:, None]
            late_writes.append(float(wl.sum() / max((s["wday"] >= 0).sum(), 1)))
        else:
            shares.append(0.0)
            late_writes.append(0.0)
    rows.append(dict(variant=name, **over,
                     gradient=float(np.mean(grads)), gradient_sem=sem(grads),
                     half_mass_day=float(np.mean(hms)), half_mass_day_sem=sem(hms),
                     expansion=float(np.mean(exps)), mean_div=float(np.mean(divs)),
                     mass_after_100=float(np.mean(m100)),
                     terminal_share_late=float(np.mean(shares)),
                     write_share_late=float(np.mean(late_writes))))
    r = rows[-1]
    tick(f"   {name}: gradient {r['gradient']:.2f}+-{r['gradient_sem']:.2f} "
         f"half-mass d{r['half_mass_day']:.0f} expansion {r['expansion']:.0f} "
         f"late compartment is {r['terminal_share_late']*100:.2f}% of endpoint cells "
         f"and {r['write_share_late']*100:.2f}% of writes")
# the alternative construction: a late compartment that keeps the main fate
# probabilities and simply switches later. It expands without bound.
p = Params(**{**BASE, "late_fraction": 0.15, "late_mid_shift": 80.0,
              "late_pp": Params(**BASE).pp_start, "late_nn": Params(**BASE).nn_start,
              "late_tc_factor": 1.0})
s = simulate(p)
R["persistent_compartment"] = dict(
    sweep=rows,
    expanding_alternative=dict(
        description="a late compartment that keeps the main fate probabilities and "
                    "switches 80 days later, rather than dividing asymmetrically",
        expansion=float(s["sym"].shape[0] * s["scale"] / p.n_founder),
        mean_div=float(s["ndiv"].mean()),
        cycling_fraction_day60=float(s["cyc_frac"][60]),
        target_expansion=float(np.mean(R["model"]["expansion"])),
        target_cycling_fraction_day60=float(R["model"]["cycling_fraction"]["60"])))
tick(f"   expanding alternative: expansion {R['persistent_compartment']['expanding_alternative']['expansion']:.0f} "
     f"against a target of {R['persistent_compartment']['expanding_alternative']['target_expansion']:.0f}")
flush()

# ------------------------------------- B. alphabet size and the error floor --
tick("B. how the miscalling ceiling depends on alphabet size")
rows = []
for name, over in (("3 symbols (one shared landmark)", dict(fid_alphabet=1)),
                   ("10 symbols (one per pulse, 8 pulses)", dict(fid_alphabet=0)),
                   ("18 symbols (one per pulse, 16 pulses)",
                    dict(fid_alphabet=0, n_fiducial=16))):
    p = Params(**{**BASE, **over, "dynamic_range": 30000.0})
    a_size = int(p.alphabet().size)
    eps = A.confusion_rate(p)
    s = simulate(p)
    c = A.contrast(p, seeds=(0, 1, 2))
    n_w = float(s["slot"].mean())
    rows.append(dict(variant=name, alphabet_size=a_size,
                     confusion_rate=float(eps),
                     ceiling_analytic=float(A.dynamic_range_ceiling(s, p)),
                     floor_analytic=float(1 - np.exp(-n_w * eps)),
                     floor_simulated=c["phi"], floor_simulated_sem=c["phi_sem"],
                     bits_at_high_D=c["bits"]))
    r = rows[-1]
    tick(f"   {name}: A={a_size} ceiling {r['ceiling_analytic']:.0f} "
         f"false-mark floor analytic {r['floor_analytic']:.4f} simulated "
         f"{r['floor_simulated']:.4f}+-{r['floor_simulated_sem']:.4f}")
R["alphabet_ceiling"] = dict(
    sweep=rows,
    law="ceiling = e (A - 1) / (c u), floor = 1 - exp(-N_w u / (A - 1))",
    note="the ceiling is linear in alphabet size minus one, so a design with a "
         "single shared landmark symbol has a ceiling six times lower than one "
         "with a distinguishable symbol per pulse")
flush()

# ------------------------------------------------ C. familiar diagnostics ---
tick("C. sensitivity, specificity and predictive value across dynamic range")
rows = []
for r in R["dynamic_range"]["sweep"]:
    pi, phi = r["confusion"]["pi"], r["confusion"]["phi"]
    rows.append(dict(D=r["D"], sensitivity=pi, specificity=1 - phi,
                     likelihood_ratio_positive=(pi / phi if phi > 0 else None),
                     ppv_at_even_prior=(pi / (pi + phi) if (pi + phi) > 0 else None),
                     bits=r["confusion"]["bits"]))
    tick(f"   D={r['D']:6d}: sensitivity {pi:.3f} specificity {1-phi:.3f} "
         f"LR+ {rows[-1]['likelihood_ratio_positive']:.1f} "
         f"PPV {rows[-1]['ppv_at_even_prior']:.3f}")
R["diagnostics"] = dict(
    sweep=rows,
    note="a single binary observable gives one operating point, not a curve, so "
         "a receiver operating characteristic would be a single point plus the "
         "two trivial corners; sensitivity, specificity and the positive "
         "likelihood ratio are the informative summaries")
flush()

# ------------------------------------- D. bracket widths and record classes --
tick("D. bracket width distribution and the three-way split of records")
rows = []
for k in (1, 2, 3, 4, 8):
    acc = {}
    for sd in SEEDS:
        p = Params(**{**BASE, "fid_alphabet": k, "seed": BASE["seed"] + sd})
        s = simulate(p)
        sym, got, wd, eid = readout(s, p)
        b = A.brackets(s, p, sym, got, wd, reader="conservative", max_rows=5000)
        for key in ("frac_two_sided", "frac_lower_bound_only", "frac_no_landmark",
                    "acc_two_sided", "median_width_two_sided"):
            acc.setdefault(key, []).append(b[key])
        if b["width_quantiles_two_sided"]:
            for q, v in b["width_quantiles_two_sided"].items():
                acc.setdefault("w" + q, []).append(v)
    row = dict(alphabet=k)
    row.update({key: float(np.nanmean(v)) for key, v in acc.items()})
    row.update({key + "_sem": sem(v) for key, v in acc.items()
                if key in ("frac_two_sided", "acc_two_sided", "median_width_two_sided")})
    rows.append(row)
    tick(f"   k={k}: two-sided {row['frac_two_sided']:.3f} lower-only "
         f"{row['frac_lower_bound_only']:.3f} unanchored {row['frac_no_landmark']:.3f} | "
         f"width 10/50/90 {row.get('w10',0):.0f}/{row.get('w50',0):.0f}/{row.get('w90',0):.0f}d")
R["bracket_distribution"] = rows
flush()

# ------------------------------------------- E. nuclei per writing event ----
tick("E. distribution of nuclei per independent writing event")
p = Params(**{**BASE, "n_founder": 150, "cap": 10 ** 9})
s = simulate(p)
sym, got, wd, eid = readout(s, p)
is_exp = sym == EXPOSURE_A
rowsel = got & is_exp.any(1)
ids = eid[rowsel][is_exp[rowsel]]
ids = ids[ids > 0]
uq, ct = np.unique(ids, return_counts=True)
ct_sorted = np.sort(ct)[::-1]
cum = np.cumsum(ct_sorted) / ct_sorted.sum()
R["event_size_distribution"] = dict(
    n_events=int(uq.size), n_marked=int(rowsel.sum()),
    quantiles={str(q): float(np.percentile(ct, q))
               for q in (10, 25, 50, 75, 90, 95, 99)},
    max=int(ct.max()), mean=float(ct.mean()),
    # complementary cumulative distribution, for plotting
    ccdf_sizes=[int(x) for x in np.unique(ct)],
    ccdf_frac=[float((ct >= x).mean()) for x in np.unique(ct)],
    share_in_top_decile=float(cum[max(int(0.1 * ct.size) - 1, 0)]),
    share_in_top_percent=float(cum[max(int(0.01 * ct.size) - 1, 0)]),
    events_holding_half_the_nuclei=int(np.searchsorted(cum, 0.5) + 1))
e = R["event_size_distribution"]
tick(f"   {e['n_marked']} nuclei from {e['n_events']} events; median {e['quantiles']['50']:.0f}, "
     f"99th {e['quantiles']['99']:.0f}, max {e['max']}; "
     f"{e['events_holding_half_the_nuclei']} events hold half the nuclei")
flush()

# --------------------------------------- F. detectable effect vs cultures ---
tick("F. detectable fold change against the number of cultures")
net = R["variance"]["between_culture"]["net_mean"]
pi = R["variance"]["between_culture"]["pi_mean"]
zsum = A.Z_ALPHA_2 + A.Z_BETA_80
curves = []
for nf, sd_between in ((150, R["founder_scaling"]["sweep"][0]["pi_sd"]),
                       (600, R["founder_scaling"]["sweep"][2]["pi_sd"])):
    for n_arrays in (500, 1000, 5000):
        pts = []
        for k in range(1, 9):
            se_arm = np.sqrt(sd_between ** 2 / k + pi * (1 - pi) / (k * n_arrays))
            delta_min = zsum * np.sqrt(2) * se_arm
            fold = (net / (net - delta_min)) if delta_min < net else None
            pts.append(dict(cultures=k, min_detectable_delta=float(delta_min),
                            min_detectable_fold=(float(fold) if fold else None)))
        curves.append(dict(n_founder=nf, arrays_per_culture=n_arrays,
                           sd_between=float(sd_between), points=pts))
        ok = [p_["min_detectable_fold"] for p_ in pts if p_["min_detectable_fold"]]
        tick(f"   founders {nf}, {n_arrays} arrays/culture: 3 cultures detect "
             f"{pts[2]['min_detectable_fold'] or float('nan'):.2f}-fold, "
             f"8 cultures {pts[7]['min_detectable_fold'] or float('nan'):.2f}-fold")
R["power_curves_normal"] = dict(net=net, pi=pi, curves=curves,
                         note="two-sided alpha 0.05, 80 per cent power, comparing "
                              "net mark rates between two arms")
flush()
tick("merged into out/results.json")
