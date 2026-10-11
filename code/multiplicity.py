"""Family-wise and false-discovery correction across the Fig. 3 mechanism arms.

The 14 arms are each tested against the same
reference at an uncorrected two-sided alpha = 0.05 and left the correction to
the reader. This computes it.

Each arm is a Welch two-sample t on seed-averaged gradients: the arm's own
standard error against the reference's, with Welch-Satterthwaite degrees of
freedom. Three decision rules are reported side by side, because they do not
agree and the disagreement is informative:

  p_raw     uncorrected two-sided significance
  p_holm    Holm-Bonferroni, controls the family-wise error rate
  p_bh      Benjamini-Hochberg, controls the false discovery rate
  mde_80    the figure's own rule: is the arm outside the difference it had
            80% power to detect. This is STRICTER than p < 0.05, because it
            needs |t| > z(0.975) + z(0.80) ~ 2.80 rather than 1.96, so an arm
            can be significant and still sit inside its bar.

Writes a `multiplicity` block into results.json. Run from code/.

    python multiplicity.py            recompute and check against the stored block
    python multiplicity.py --write    recompute and store
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"
N_SEEDS = 50  # Methods: every Fig. 3 arm is run over fifty populations
ALPHA = 0.05
POWER = 0.80


def welch(g1, se1, n1, g0, se0, n0):
    """Welch t on two means given their standard errors and seed counts."""
    se = float(np.hypot(se1, se0))
    t = (g1 - g0) / se
    df = se**4 / (se1**4 / (n1 - 1) + se0**4 / (n0 - 1))
    p = 2 * stats.t.sf(abs(t), df)
    return float(t), float(df), float(p), se


def mde(se_diff, df):
    """Smallest difference detectable at POWER and two-sided ALPHA.

    Matches figures.mde: the normal-approximation form the figure draws.
    """
    return float(se_diff * (stats.norm.ppf(1 - ALPHA / 2) + stats.norm.ppf(POWER)))


def holm(p):
    """Holm-Bonferroni adjusted p-values, order-preserving (monotone)."""
    p = np.asarray(p, float)
    m = p.size
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(np.argsort(p)):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(running, 1.0)
    return adj


def benjamini_hochberg(p):
    """BH adjusted p-values (q-values), order-preserving."""
    p = np.asarray(p, float)
    m = p.size
    order = np.argsort(p)
    adj = np.empty(m)
    running = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        running = min(running, m * p[i] / (rank + 1))
        adj[i] = min(running, 1.0)
    return adj


def arms(R):
    """The 14 arms figure 2 draws, in the order figures.py collects them."""
    out = []
    nd = {r["q_nondiv"]: r for r in R["nondivision"]["sweep"]}
    for q, lab in ((0.0072, "writing without division, adult-cortex rate"),
                   (0.02, "writing without division, 0.02 per day")):
        k = min(nd, key=lambda z: abs(z - q))
        out.append((lab, nd[k]))
    for r in R["death"]["sweep"]:
        if r.get("death_pm") == 0.010:
            out.append(("death of neurons, 0.010 per day", r))
    pc = R["persistent_compartment"]["sweep"][3]
    out.append((f"persistent progenitor compartment, "
                f"{pc['late_fraction']:.0%} of founders", pc))
    for r in R["ipc"]["sweep"][1:]:
        out.append((f"transit-amplifying IPCs, {r['ipc_fraction']:.0%} of "
                    f"daughters, {r['ipc_rounds']} rounds", r))
    for r in R["position_efficiency"]["sweep"]:
        if r["q_pos_halflife"] == 6:
            out.append(("per-monomer efficiency, 6-monomer half-life", r))
    for r in R["silencing"]["sweep"]:
        if r["q_halflife"] in (50, 30):
            out.append((f"editor silencing, {int(r['q_halflife'])}-day half-life", r))
    for r in R.get("gliogenesis_switched", {}).get("sweep", []):
        out.append((f"gliogenesis from day {r['late_start_day']:.0f}", r))
    return out


def compute(R):
    ref = R["nondivision"]["sweep"][0]
    g0, se0 = ref["gradient"], ref["gradient_sem"]

    rows = []
    for lab, r in arms(R):
        t, df, p, se = welch(r["gradient"], r["gradient_sem"], N_SEEDS,
                             g0, se0, N_SEEDS)
        d = mde(se, df)
        rows.append(dict(
            arm=lab, gradient=r["gradient"], gradient_sem=r["gradient_sem"],
            delta=r["gradient"] - g0, t=t, df=df, p_raw=p, se_diff=se,
            mde_80=d, resolvable_mde=bool(abs(r["gradient"] - g0) > d),
            direction="deepens" if r["gradient"] > g0 else "blunts"))

    p = [r["p_raw"] for r in rows]
    for r, ph, pb in zip(rows, holm(p), benjamini_hochberg(p)):
        r["p_holm"], r["p_bh"] = float(ph), float(pb)
        r["sig_raw"] = bool(r["p_raw"] < ALPHA)
        r["sig_holm"] = bool(ph < ALPHA)
        r["sig_bh"] = bool(pb < ALPHA)

    rows.sort(key=lambda r: r["p_raw"])
    n_mde = sum(r["resolvable_mde"] for r in rows)
    return dict(
        reference_gradient=g0, reference_sem=se0, n_seeds=N_SEEDS,
        n_arms=len(rows), alpha=ALPHA, power=POWER,
        n_sig_raw=sum(r["sig_raw"] for r in rows),
        n_sig_holm=sum(r["sig_holm"] for r in rows),
        n_sig_bh=sum(r["sig_bh"] for r in rows),
        n_resolvable_mde=n_mde,
        mde_and_holm_agree=bool(all(r["resolvable_mde"] == r["sig_holm"]
                                    for r in rows)),
        expected_false_positives_uncorrected=ALPHA * len(rows),
        arms=rows,
        note=("Welch two-sample t on seed-averaged gradients, each arm against "
              "the same reference, fifty populations per arm. Holm-Bonferroni "
              "controls the family-wise error rate over the eleven arms; "
              "Benjamini-Hochberg controls the false discovery rate. mde_80 is "
              "the figure's own power-based rule, which is stricter than "
              "uncorrected significance because it requires |t| > 2.80 rather "
              "than 1.96."),
    )


def main():
    R = json.loads(RESULTS.read_text())
    block = compute(R)

    w = max(len(r["arm"]) for r in block["arms"])
    print(f"reference gradient {block['reference_gradient']:.3f} "
          f"+- {block['reference_sem']:.3f}, {N_SEEDS} populations\n")
    print(f"{'arm':{w}}  {'gradient':>14}  {'p':>9} {'Holm':>9} {'BH':>9}  "
          f"{'|t|':>5}  outside bar")
    for r in block["arms"]:
        print(f"{r['arm']:{w}}  {r['gradient']:6.2f} +- {r['gradient_sem']:4.2f}  "
              f"{r['p_raw']:9.2e} {r['p_holm']:9.2e} {r['p_bh']:9.2e}  "
              f"{abs(r['t']):5.2f}  {'yes' if r['resolvable_mde'] else 'no'}")
    print(f"\nresolvable: {block['n_sig_raw']} uncorrected, "
          f"{block['n_sig_bh']} under BH, {block['n_sig_holm']} under Holm, "
          f"{block['n_resolvable_mde']} outside their own 80%-power bar, "
          f"of {block['n_arms']} arms")
    print("figure rule and Holm agree on every arm: "
          f"{block['mde_and_holm_agree']}")

    if "--write" in sys.argv:
        R["multiplicity"] = block
        RESULTS.write_text(json.dumps(R, indent=1))
        print(f"\nwrote `multiplicity` into {RESULTS}")
        return 0
    if "multiplicity" in R:
        stored = R["multiplicity"]
        same = (len(stored["arms"]) == len(block["arms"]) and
                all(abs(a["p_raw"] - b["p_raw"]) < 1e-12
                    for a, b in zip(block["arms"], stored["arms"])))
        print(f"\nstored block {'matches' if same else 'DIFFERS FROM'} this run")
        return 0 if same else 1
    print("\nno stored block yet; rerun with --write")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
