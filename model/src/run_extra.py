"""Array length, amplicon span, and the published anchors.

It is tempting to call arrays over-provisioned, reasoning from a 50-monomer array
as the reference. Two things in the literature make that the wrong reference.
Current multi-channel designs use a five-unit array (Chen et al. 2024; Nathans
et al. 2026), which is shorter than the mean occupancy this model predicts over
150 days, and recovery of a long tandem array is incomplete: 12.5 of 20 and 8.4
of 12 monomers recovered (Choi et al. 2022), so length is not free. This script
computes the length sweep under both assumptions and writes the result into
out/results.json alongside the main sections.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from recorder import Params, simulate, readout, EXPOSURE_A
import analysis as A

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000)
# The released results.json was produced at 3 for this script.
SEEDS = (0, 1, 2)
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


R = json.load(open(os.path.join(OUT, "results.json")))


def flush():
    with open(os.path.join(OUT, "results.json"), "w") as fh:
        json.dump(R, fh, indent=1, default=_enc)

# ------------------------------------------------- array length, no span cost
tick("L1. array length with whole-array recovery")
rows = []
for q in (0.10, 0.18, 0.30):
    ref = None
    for L in (5, 10, 12, 20, 50, 200):
        p = Params(**{**BASE, "q_div": q, "L": L})
        s = simulate(p)
        c = A.contrast(p, seeds=SEEDS)
        pe = A.per_epoch_recording(s)
        if ref is None and L == 5:
            pass
        rows.append(dict(q=q, L=L, pi=c["pi"], phi=c["phi"], bits=c["bits"],
                         occupancy=float(s["slot"].mean()),
                         frac_saturated=float((s["slot"] >= L).mean()),
                         per_epoch=[float(x) for x in pe],
                         gradient=float(pe[0] / pe[-1]) if pe[-1] > 0 else None))
    base = [r for r in rows if r["q"] == q and r["L"] == 200][0]
    for r in [r for r in rows if r["q"] == q]:
        r["bits_rel_unbounded"] = float(r["bits"] / base["bits"]) if base["bits"] else None
        r["late_epoch_loss"] = float(1 - r["per_epoch"][-1] / base["per_epoch"][-1]) \
            if base["per_epoch"][-1] else None
    tick(f"   q={q}: L=5 keeps {rows[-6]['bits_rel_unbounded']:.2f} of the information, "
         f"saturated {rows[-6]['frac_saturated']:.2f}, late-epoch loss "
         f"{rows[-6]['late_epoch_loss']:.2f}")
R["array_length_sweep"] = rows
flush()

# ------------------------------------------- array length with a recovery cost
tick("L2. array length when amplicon recovery of a long array is incomplete")
# moments matched to Choi et al. 2022: 12.5 +- 4.3 of 20 monomers recovered
span_mean, span_sd = 12.5 / 20.0, 4.3 / 20.0
rows = []
for L in (5, 10, 12, 20, 50, 100):
    p = Params(**{**BASE, "L": L, "span_mean": span_mean, "span_sd": span_sd})
    s = simulate(p)
    c = A.contrast(p, seeds=SEEDS)
    sym, got, wd, eid = readout(s, p)
    rows.append(dict(L=L, pi=c["pi"], phi=c["phi"], bits=c["bits"],
                     occupancy_written=float(s["slot"].mean()),
                     monomers_read=float((sym != 0).sum(1).mean()),
                     frac_saturated=float((s["slot"] >= L).mean())))
    tick(f"   L={L}: bits {c['bits']:.4f} written {rows[-1]['occupancy_written']:.1f} "
         f"read {rows[-1]['monomers_read']:.1f}")
R["array_span"] = dict(span_mean=span_mean, span_sd=span_sd, sweep=rows,
                       optimum_L=int(max(rows, key=lambda r: r["bits"])["L"]))
flush()

# ------------------------------------------------------- published anchors ---
tick("L3. published anchors used to set and to check the parameters")
R["anchors"] = dict(
    editing_rate=dict(
        choi2022_first_site_edited_fraction=0.181,
        choi2022_edits_per_cell_division=2.0,
        choi2022_arrays_per_cell=13,
        implied_per_array_per_division=2.0 / 13,
        model_q_div=0.18,
        note="Choi et al. 2022 report about two edits per cell division across "
             "13 DNA Array arrays in HEK293T, which is 0.15 per array per "
             "division. The model's q = 0.18 is that quantity, not a "
             "population edited fraction."),
    array_length=dict(
        chen2024_engram_array_units=5,
        nathans2026_engram_array_units=5,
        choi2022_common_array_units=5,
        choi2022_extended_array_units=20,
        choi2022_recovered_of_20=12.5, choi2022_recovered_of_20_sd=4.3,
        choi2022_recovered_of_12=8.4, choi2022_recovered_of_12_sd=3.3,
        yu2026_sites_in_vivo=66, yu2026_written_by_e135=55),
    channel_dynamic_range=dict(
        chen2024_teton=11.5, chen2024_nfkb=19.0, chen2024_wnt=22.6,
        note="measured fold-change of stimulated over unstimulated writing for "
             "the three ENGRAM channels reported in Chen et al. 2024"),
    division_independent=dict(
        davis2024_adult_bulk_cortex=0.14, davis2024_neonatal_cortex=0.42,
        davis2024_sorted_cortical_neurons=0.81, davis2024_days=21,
        per_day_adult=float(1 - (1 - 0.14) ** (1 / 21)),
        per_day_neonatal=float(1 - (1 - 0.42) ** (1 / 21)),
        note="dual-AAV prime editing in mouse brain, three weeks after "
             "injection, including postmitotic cortical neurons"),
)

# --------------------------------- L4. clonality, and how it compares with
# --------------------------------- what has been measured in organoids
tick("L4. lineage-size heterogeneity and the independence of records")
rows = []
for sd in (0.0, 8.0, 15.0):
    p = Params(**{**BASE, "n_founder": 150, "cap": 400000, "fate_mid_sd": sd})
    s = simulate(p)
    sym, got, wd, eid = readout(s, p)
    es = A.event_stats(s, p, sym, got, eid)
    f = s["founder"]
    _, ct = np.unique(f, return_counts=True)
    ct = np.sort(ct)[::-1]
    top5 = float(ct[:max(1, int(round(0.05 * ct.size)))].sum() / ct.sum())
    # He et al. 2022 assayed 44,275 nuclei from organoids of order a million
    # cells, so the comparable quantity is clone sizes in a sample of a few per
    # cent of the tissue, not in the whole simulated tree
    rng = np.random.default_rng(5)
    take = max(int(round(0.03 * f.size)), 200)
    pick = rng.choice(f.size, take, replace=False)
    _, sct = np.unique(f[pick], return_counts=True)
    bc = A.between_culture_variance(
        Params(**{**BASE, "fate_mid_sd": sd, "n_founder": 150, "cap": 400000}),
        seeds=range(4))
    n_rec = 15500 * 0.6 * 0.5
    rows.append(dict(fate_mid_sd=sd, n_cells=int(f.size), sampled=int(take),
                     expansion=float(f.size / p.n_founder),
                     n_lineages=int(ct.size), top5_share=top5,
                     clone_orders_of_magnitude=float(np.log10(ct.max() / max(ct.min(), 1))),
                     sampled_clone_mean=float(sct.mean()), sampled_clone_max=int(sct.max()),
                     pseudoreplication=es["pseudoreplication"],
                     n_marked=es["n_marked"], n_events=es["n_events"],
                     max_marked_clone=es["max_clone"],
                     top_decile_share=es["top_decile_share"],
                     between_culture_sd=bc["pi_sd"],
                     design_effect=float(A.design_effect(bc["pi_mean"], bc["pi_var"], n_rec))))
    tick(f"   fate_mid_sd={sd}: top5% share {top5:.2f} (measured 0.80); "
         f"sampled clone mean {rows[-1]['sampled_clone_mean']:.1f} max "
         f"{rows[-1]['sampled_clone_max']} (measured 11, max 801); "
         f"pseudoreplication {es['pseudoreplication']:.1f}x; "
         f"design effect {rows[-1]['design_effect']:.1f}")
# how the culture-to-culture variance depends on the number of founder lineages
tick("L5. between-culture variance against the number of recording founders")
fr = []
for nf in (150, 300, 600):
    p = Params(**{**BASE, "n_founder": nf, "cap": 400000})
    bc = A.between_culture_variance(p, seeds=range(5))
    s = simulate(p)
    sym, got, wd, eid = readout(s, p)
    es = A.event_stats(s, p, sym, got, eid)
    n_star = bc["pi_mean"] * (1 - bc["pi_mean"]) / bc["pi_var"] if bc["pi_var"] > 0 else None
    fr.append(dict(n_founder=nf, pi_mean=bc["pi_mean"], pi_sd=bc["pi_sd"],
                   pi_var=bc["pi_var"], net_sd=bc["net_sd"],
                   arrays_where_noise_equals_floor=(float(n_star) if n_star else None),
                   pseudoreplication=es["pseudoreplication"],
                   n_events=es["n_events"], max_clone=es["max_clone"]))
    tick(f"   founders {nf}: between-culture SD {bc['pi_sd']:.4f}; sequencing more than "
         f"{n_star:.0f} arrays per culture stops helping; pseudoreplication "
         f"{es['pseudoreplication']:.1f}x")
# the variance of a mean over independent lineages scales as 1/n_founder
sc = np.polyfit(np.log([r["n_founder"] for r in fr]),
                np.log([r["pi_var"] for r in fr]), 1)
flush()
R["founder_scaling"] = dict(sweep=fr, log_slope=float(sc[0]),
                            note="slope near -1 means the culture-to-culture variance "
                                 "is set by the number of independently recording "
                                 "founder lineages")
tick(f"   variance scales with founder number to the power {sc[0]:.2f}")

R["clonality"] = dict(
    sweep=rows,
    measured=dict(lindenhofer2024_top5_share_of_progeny=0.80,
                  lindenhofer2024_orders_of_magnitude=5,
                  lindenhofer2024_barcoded_founders=24000,
                  he2022_cells_assayed=44275, he2022_cells_with_readout=22489,
                  he2022_readout_fraction=22489 / 44275,
                  he2022_barcode_family_mean=11, he2022_barcode_family_max=801,
                  he2022_scarred_share_of_barcoded=0.17),
)
flush()

with open(os.path.join(OUT, "results.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=_enc)
tick("merged into out/results.json")
