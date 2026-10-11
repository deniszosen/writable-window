"""Confront the recording model with data it was not calibrated on.

Every number in versions 1 to 10 is internal to the model. The writing rate
q = 0.18 per array per division is anchored to HEK293T (Choi et al. 2022); the
cycling trajectory is an assumption. Reviewers asked, across three versions,
whether the model reproduces an independent published measurement.

One exists. Yu et al. ran DNA Typewriter in a mouse from the zygote to E13.5
and report, for an animal carrying 11 Array integrations of 6 monomers each
(66 writable sites):

    55 +- 5 of 66 sites written per cell (84%)
    median 7 of 11 arrays exhausted (all six monomers used)
    only 0.5% of cells had written all 66 sites; every cell had written >= 33
    per-site editing fell from 96% at monomer 1 to 53% at monomer 6
    "by E13.5 the recording capacity had been >80% used but not saturated"

This is the quantity the recording module produces, in a system the module was
never parameterised for. Nothing here is fitted: q stays at its HEK293T anchor
and the only quantity taken from mouse embryology is the number of divisions
from zygote to E13.5, which is bounded below by log2(cells at E13.5).

Two tests:

  1. OCCUPANCY. Over a grid of (divisions, q), where does mean occupancy fall
     inside the observed 55 +- 5 of 66? Does the HEK293T anchor land in that
     band at a division count mouse embryology allows?

  2. PROFILE SHAPE. A strictly sequential writer makes monomer s reachable
     only once 1..s-1 are filled, so per-site occupancy declines along the
     array from ordering alone, with no efficiency decay. Does ordering alone
     explain the observed 96% -> 53%, or is a per-position term needed?

Test 2 is the one that can fail, and it is the one that bears on section 2.2:
the per-monomer efficiency arm is introduced there as an unmeasured property
of the construct.

Writes a `mouse_validation` block into results.json. Run from code/.

    python mouse_validation.py            recompute and print
    python mouse_validation.py --write    recompute and store
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import special

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"

# ---- what Yu et al. measured, transcribed once ------------------------------
OBS = dict(
    n_arrays=11,
    monomers_per_array=6,
    sites_total=66,
    sites_written_mean=55.0,
    sites_written_sd=5.0,
    per_site_first=0.96,          # monomer 1, 95-100% across the 11 integrations
    per_site_last=0.53,           # monomer 6, 26-77% across the 11 integrations
    median_arrays_exhausted=7,
    frac_cells_complete=0.005,
    min_sites_written=33,         # "every cell had written at least half"
    source=("Yu Q, Kim H, Seidel S, et al. In vivo reconstruction of the cell "
            "lineage history of a developing mouse with DNA Typewriter, from "
            "zygote to late organogenesis. bioRxiv 2026. "
            "doi:10.64898/2026.07.29.741625 (PMC13484211)"),
)

Q_ANCHOR = 0.18       # HEK293T anchor used throughout the manuscript
L = OBS["monomers_per_array"]

# Divisions from zygote to E13.5. A mouse embryo at E13.5 carries of order
# 1e7 to 1e8 cells, so a lineage has undergone at least log2(1e7) = 23 and at
# most about log2(1e8) plus a few asynchronous rounds. The grid spans wider
# than that so the consistency region is visible rather than assumed.
DIVISIONS = np.arange(14, 45)
Q_GRID = np.linspace(0.04, 0.60, 141)


def binom_sf(s, n, q):
    """P(Binomial(n, q) >= s) as the regularised incomplete beta I_q(s, n-s+1).

    scipy.stats.binom.sf needs an integer n; the beta form is the same tail
    and is continuous in n, which the division-count solve below needs.
    """
    s = np.atleast_1d(np.asarray(s, float))
    out = np.where(s <= 0, 1.0,
                   np.where(s > n, 0.0,
                            special.betainc(np.maximum(s, 1e-9),
                                            np.maximum(n - s + 1, 1e-9), q)))
    return out


def occupancy(n_div, q, L=L):
    """Mean monomers filled on one array: E[min(Binomial(n_div, q), L)].

    E[min(K, L)] = sum_{s=1..L} P(K >= s), which is also the per-site
    occupancy profile, so one call gives both tests.
    """
    p_at_least = binom_sf(np.arange(1, L + 1), n_div, q)
    return float(p_at_least.sum()), p_at_least


def divisions_for(target_sites, q, L=L, n_arrays=OBS["n_arrays"]):
    """Fractional division count at which the model reaches `target_sites`/66."""
    target_per_array = target_sites / n_arrays
    lo, hi = 1.0, 400.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if occupancy(mid, q, L)[0] < target_per_array:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def per_cell_distribution(n_div, q, n_cells=200_000, seed=11,
                          n_arrays=OBS["n_arrays"], L=L):
    """Sites written per cell, summing 11 independent arrays in one cell.

    The model treats arrays within a cell as independent given the division
    history, so a cell's 11 arrays are 11 draws conditioned on the same n_div.
    """
    rng = np.random.default_rng(seed)
    k = rng.binomial(n_div, q, size=(n_cells, n_arrays))
    filled = np.minimum(k, L)
    total = filled.sum(1)
    return dict(
        mean=float(total.mean()), sd=float(total.std(ddof=1)),
        median=float(np.median(total)),
        min=int(total.min()), max=int(total.max()),
        frac_complete=float((total == n_arrays * L).mean()),
        frac_at_least_half=float((total >= n_arrays * L / 2).mean()),
        median_arrays_exhausted=float(np.median((k >= L).sum(1))),
        mode=int(np.bincount(total).argmax()),
        hist=[int(c) for c in np.bincount(total, minlength=n_arrays * L + 1)],
    )


def compute():
    lo_sites = OBS["sites_written_mean"] - OBS["sites_written_sd"]
    hi_sites = OBS["sites_written_mean"] + OBS["sites_written_sd"]

    # ---- test 1: where does the model match the observed occupancy? --------
    grid = np.zeros((DIVISIONS.size, Q_GRID.size))
    for i, n in enumerate(DIVISIONS):
        for j, q in enumerate(Q_GRID):
            grid[i, j] = occupancy(n, q)[0] * OBS["n_arrays"]

    n_at_anchor = divisions_for(OBS["sites_written_mean"], Q_ANCHOR)
    n_lo = divisions_for(lo_sites, Q_ANCHOR)
    n_hi = divisions_for(hi_sites, Q_ANCHOR)

    # the reciprocal view: at the division count mouse embryology allows,
    # what writing rate does the observation imply?
    q_implied = {}
    for n in (24, 28, 32):
        lo, hi = 0.001, 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if occupancy(n, mid)[0] * OBS["n_arrays"] < OBS["sites_written_mean"]:
                lo = mid
            else:
                hi = mid
        q_implied[str(n)] = 0.5 * (lo + hi)

    # ---- test 2: does ordering alone explain the per-site profile? ---------
    # Evaluate at the division count that matches the observed mean, so the
    # profile test is not confounded by getting the mean wrong.
    n_match = n_at_anchor
    _, profile = occupancy(n_match, Q_ANCHOR)
    pred_first, pred_last = float(profile[0]), float(profile[-1])
    obs_ratio = OBS["per_site_first"] / OBS["per_site_last"]
    pred_ratio = pred_first / pred_last

    # what per-position half-life would be needed on top of ordering, if the
    # residual decline is attributed to efficiency rather than to ordering?
    # Solve 2**(-(L-1)/h) = pred_ratio / obs_ratio for h.
    resid = obs_ratio / pred_ratio
    h_needed = (float((L - 1) * np.log(2) / np.log(resid))
                if resid > 1.0 else None)

    dist = per_cell_distribution(int(round(n_match)), Q_ANCHOR)

    anchor_ok = 23 <= n_at_anchor <= 34
    profile_ok = pred_ratio >= obs_ratio          # ordering alone suffices

    return dict(
        observed=OBS,
        q_anchor=Q_ANCHOR,
        divisions_grid=[int(d) for d in DIVISIONS],
        q_grid=[float(q) for q in Q_GRID],
        occupancy_grid=[[float(x) for x in row] for row in grid],
        occupancy_band=[lo_sites, hi_sites],
        divisions_required_at_anchor=float(n_at_anchor),
        divisions_required_band=[float(n_lo), float(n_hi)],
        divisions_plausible_range=[23, 34],
        divisions_plausible_note=(
            "A mouse embryo at E13.5 carries of order 1e7 to 1e8 cells, so a "
            "terminal lineage has undergone at least log2(1e7) = 23 divisions; "
            "asynchrony and cell loss raise the mean above that bound."),
        q_implied_at_divisions=q_implied,
        anchor_consistent=bool(anchor_ok),
        profile=dict(
            site=list(range(1, L + 1)),
            predicted_fraction=[float(x) for x in profile],
            observed_first=OBS["per_site_first"],
            observed_last=OBS["per_site_last"],
            predicted_first=pred_first,
            predicted_last=pred_last,
            observed_decline_fold=float(obs_ratio),
            predicted_decline_fold=float(pred_ratio),
            ordering_alone_sufficient=bool(profile_ok),
            residual_fold=float(resid),
            per_position_halflife_needed=h_needed,
        ),
        predicted_cell_distribution=dist,
        verdict=(
            f"At the HEK293T-anchored writing rate the model reaches the "
            f"observed {OBS['sites_written_mean']:.0f} of "
            f"{OBS['sites_total']} sites after {n_at_anchor:.0f} divisions, "
            f"{'inside' if anchor_ok else 'outside'} the 23 to 34 divisions "
            f"mouse embryology allows between the zygote and E13.5."),
        note=("Nothing is fitted. q is held at the value anchored to HEK293T "
              "in Table S1 and the only mouse-specific input is the division "
              "count, which is bounded below by the cell number at E13.5. The "
              "population model of section 4.1 is not used: occupancy depends "
              "on the population only through the division count, so this "
              "tests the recording module alone."),
    )


def main():
    b = compute()
    o = b["observed"]
    print("observed (Yu et al., mouse zygote to E13.5)")
    print(f"  {o['sites_written_mean']:.0f} +- {o['sites_written_sd']:.0f} of "
          f"{o['sites_total']} sites written "
          f"({o['sites_written_mean'] / o['sites_total']:.0%})")
    print(f"  per-site editing {o['per_site_first']:.0%} at monomer 1 -> "
          f"{o['per_site_last']:.0%} at monomer {o['monomers_per_array']}")
    print(f"  median {o['median_arrays_exhausted']} of {o['n_arrays']} arrays "
          f"exhausted; {o['frac_cells_complete']:.1%} of cells complete\n")

    print(f"TEST 1  occupancy, at the HEK293T anchor q = {b['q_anchor']}")
    lo, hi = b["divisions_required_band"]
    print(f"  divisions needed to reach 55 of 66: "
          f"{b['divisions_required_at_anchor']:.1f}  ({lo:.1f} to {hi:.1f} "
          f"across the observed +- 1 sd)")
    print(f"  mouse embryology allows {b['divisions_plausible_range'][0]} to "
          f"{b['divisions_plausible_range'][1]}")
    print(f"  consistent: {b['anchor_consistent']}")
    print("  writing rate the observation implies, by division count: " +
          ", ".join(f"{k} div -> q = {v:.3f}"
                    for k, v in b["q_implied_at_divisions"].items()))

    p = b["profile"]
    print(f"\nTEST 2  per-site profile at "
          f"{b['divisions_required_at_anchor']:.0f} divisions")
    print("  monomer      1     2     3     4     5     6")
    print("  predicted " + " ".join(f"{x:5.2f}" for x in p["predicted_fraction"]))
    print(f"  observed   {p['observed_first']:5.2f}     .     .     .     .  "
          f"{p['observed_last']:5.2f}")
    print(f"  decline across the array: observed "
          f"{p['observed_decline_fold']:.2f}-fold, model "
          f"{p['predicted_decline_fold']:.2f}-fold")
    print(f"  ordering alone explains the observed decline: "
          f"{p['ordering_alone_sufficient']}")
    if p["per_position_halflife_needed"]:
        print(f"  residual needs a per-position half-life of "
              f"{p['per_position_halflife_needed']:.1f} monomers")

    d = b["predicted_cell_distribution"]
    print(f"\n  predicted per-cell sites written: {d['mean']:.1f} +- "
          f"{d['sd']:.1f}, mode {d['mode']}, min {d['min']}")
    print(f"  predicted arrays exhausted (median): "
          f"{d['median_arrays_exhausted']:.0f} of {o['n_arrays']}  "
          f"[observed {o['median_arrays_exhausted']}]")
    print(f"  predicted cells complete: {d['frac_complete']:.2%}  "
          f"[observed {o['frac_cells_complete']:.1%}]")
    print(f"  predicted cells over half written: "
          f"{d['frac_at_least_half']:.2%}  [observed 100%]")
    print(f"\n{b['verdict']}")

    if "--write" in sys.argv:
        R = json.loads(RESULTS.read_text())
        R["mouse_validation"] = b
        RESULTS.write_text(json.dumps(R, indent=1))
        print(f"\nwrote `mouse_validation` into {RESULTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
