"""Figures 2 to 5 and the two supplementary figures.

One claim per figure, and nothing that is model description rather than
result. Figure 1 is drawn by fig1_timeline.py, which needs the per-day
trajectory rather than the results file.

  Fig 2  every mechanism on one gradient axis, each null drawn as a bound
  Fig 3  records are clonally correlated, so cultures are the unit of replication
  Fig 4  the exposure schedule sets the channel requirement; readout error caps it
  Fig 5  the recording model against the mouse embryo it was not fitted to
  S1     array length
  S2     temporal anchoring

Run: python figures.py
"""

import json
import os

import matplotlib as mpl

mpl.use("Agg")
import numpy as np                                                  # noqa: E402
from matplotlib import pyplot as plt                                # noqa: E402
from scipy import optimize, stats                                   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
# Read the results file shipped beside this script, so the release is
# self-contained: extracting the deliverable alone must be enough to redraw
# every figure. Reaching into the author's working tree instead would make
# worked in place and failed for anyone who received only the release.
RES = os.path.join(HERE, "results.json")
OUT = os.path.join(os.path.dirname(HERE), "figures")
os.makedirs(OUT, exist_ok=True)
R = json.load(open(RES))

INK, BLUE, ORANGE = "#1a202c", "#2b6cb0", "#c05621"
GREEN, GREY, LIGHT = "#2f855a", "#718096", "#cbd5e0"
PURPLE = "#6b46c1"
# Founder count is an ORDERED quantity, so it gets a single-hue sequential
# ramp rather than three categorical hues. The light step is #5a95cc rather
# than anything paler because a paler blue falls below 3:1 against white.
RAMP = ("#5a95cc", "#2b6cb0", "#17375c")

# Figure 2's classes are two hues on a polarity axis plus marker SHAPE for the
# third, which survives greyscale printing and every form of colour blindness.
# A green/grey/orange encoding was rejected: the green-grey pair separates by
# only dE 12.7 for normal vision and 8.2 under deuteranopia, and that is the
# pair carrying the figure's main distinction.
#   blunts the window   filled blue
#   deepens the window  filled orange
#   not resolvable      open circle, ink outline
BLUNT, DEEPEN = BLUE, ORANGE

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "lines.linewidth": 1.6,
    "figure.dpi": 400, "savefig.dpi": 400, "savefig.bbox": "tight",
    "savefig.facecolor": "white", "svg.fonttype": "none",
    # IOP restricts embedded fonts to Times, Helvetica, Courier and Symbol.
    # The default mathtext set is DejaVu, which an italic $q$ would embed, so
    # map mathtext onto the same Helvetica the rest of the figure uses.
    "mathtext.fontset": "custom",
    "mathtext.rm": "Helvetica", "mathtext.it": "Helvetica:italic",
    "mathtext.bf": "Helvetica:bold", "mathtext.sf": "Helvetica",
    "mathtext.cal": "Helvetica", "mathtext.tt": "Courier",
})


def panel(ax, letter, title):
    ax.set_title(f"{letter}   {title}", loc="left", fontsize=9,
                 fontweight="bold", pad=4)


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}")


def mde(sems, n_seeds=50, alpha=0.05, power=0.80):
    """Smallest difference two seed-averaged means could resolve.

    Reporting 'no resolvable change' without stating what change the seed
    count could have detected says nothing. Figure 2 draws that bound, so a
    null reads as a bound rather than as an absence.
    """
    se = float(np.sqrt(np.sum(np.square(sems))))
    df = 2 * (n_seeds - 1)
    crit = stats.t.ppf(1 - alpha / 2, df)

    def f(d):
        nc = d / se
        upper = stats.nct.sf(crit, df, nc)
        lower = stats.nct.cdf(-crit, df, nc)
        # The opposite-tail term is the probability of rejecting in the wrong
        # direction. scipy returns NaN for it at moderate and large positive
        # noncentrality instead of underflowing to zero, which aborts the
        # solver. A noncentral t is (Z + nc)/sqrt(V/df) with a positive
        # denominator, so T < 0 requires Z < -nc and the term is bounded
        # above by Phi(-nc); every nc at which scipy returns NaN here puts
        # that bound below 1e-17, so zero is exact in double precision.
        if not np.isfinite(lower):
            lower = 0.0
        return upper + lower - power

    return float(optimize.brentq(f, 1e-12, max(20 * se, 1.0)))


# ======================== Fig 2: every mechanism, with each bound drawn ======
def fig2():
    """Every mechanism against the gradient, with each null honestly bounded.

    Each arm is classified against ITS OWN detectable difference, computed from
    its own standard error and the reference's. A single shaded band would be
    wrong, because an arm with a tight SEM can resolve a smaller change than one
    with a loose SEM, and drawing one band made three arms look resolved when
    they are not.
    """
    ref = R["nondivision"]["sweep"][0]
    g0, g0e = ref["gradient"], ref["gradient_sem"]
    rows = []

    nd = {r["q_nondiv"]: r for r in R["nondivision"]["sweep"]}
    # At fifty seeds both the adult-cortex rate and 0.02/day resolvably blunt
    # the gradient, so neither needs the "outside anchored range" caveat the
    # five-seed analysis required.
    for q, lab in ((0.0072, "writing without division,\nadult-cortex rate"),
                   (0.02, "writing without division,\n0.02 per day")):
        k = min(nd, key=lambda z: abs(z - q))
        rows.append((lab, nd[k]["gradient"], nd[k]["gradient_sem"]))
    for r in R["death"]["sweep"]:
        if r.get("death_pm") == 0.010:
            rows.append(("death of neurons,\n0.010 per day", r["gradient"],
                         r["gradient_sem"]))
    pc = R["persistent_compartment"]["sweep"][3]
    rows.append((f"persistent progenitor compartment,\n"
                 f"{pc['late_fraction']:.0%} of founders",
                 pc["gradient"], pc["gradient_sem"]))
    for r in R["ipc"]["sweep"][1:]:
        rows.append((f"transit-amplifying progenitors,\n{r['ipc_fraction']:.0%} of "
                     f"daughters, {r['ipc_rounds']} rounds",
                     r["gradient"], r["gradient_sem"]))
    for r in R["position_efficiency"]["sweep"]:
        if r["q_pos_halflife"] == 6:
            rows.append(("per-monomer efficiency,\n6-monomer half-life",
                         r["gradient"], r["gradient_sem"]))
    for r in R["silencing"]["sweep"]:
        if r["q_halflife"] in (50, 30):
            rows.append((f"editor silencing,\n{int(r['q_halflife'])}-day half-life",
                         r["gradient"], r["gradient_sem"]))
    for r in R.get("gliogenesis_switched", {}).get("sweep", []):
        rows.append((f"gliogenesis,\nswitching on at day {r['late_start_day']:.0f}",
                     r["gradient"], r["gradient_sem"]))

    # classify each arm against its own detectable difference
    out = []
    for lab, g, s_ in rows:
        d = mde([g0e, s_])
        if g - g0 > d:
            fam = "worse"
        elif g0 - g > d:
            fam = "rescue"
        else:
            fam = "null"
        out.append((lab, g, s_, fam, d))
    out.sort(key=lambda z: z[1])

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    col = {"rescue": BLUNT, "null": INK, "worse": DEEPEN}
    ax.axvline(g0, color=INK, lw=1.3, zorder=2)
    ax.axvspan(g0 - g0e, g0 + g0e, color=LIGHT, alpha=0.6, lw=0, zorder=0)

    # Holm-adjusted p-values across the same eleven arms, so the figure can
    # say which verdicts survive a family-wise correction instead of leaving
    # the reader to assume none of them do.
    holm = {a["arm"].replace(", ", ",\n"): a for a in
            R.get("multiplicity", {}).get("arms", [])}

    for i, (lab, g, s_, fam, d) in enumerate(out):
        # the span this arm could not have resolved
        ax.plot([g0 - d, g0 + d], [i, i], color=LIGHT, lw=5, alpha=0.75,
                solid_capstyle="butt", zorder=1)
        # A null is an OPEN marker, not a third hue: shape carries the class
        # in greyscale and under every colour-vision deficiency.
        face = "white" if fam == "null" else col[fam]
        ax.errorbar(g, i, xerr=s_, fmt="o", ms=5.8, color=col[fam],
                    mfc=face, mec=col[fam], mew=1.3,
                    ecolor=col[fam], elinewidth=1.3, capsize=2.6, zorder=3)

    ax.set_yticks(np.arange(len(out)))
    ax.set_yticklabels([r[0] for r in out], fontsize=7.2)
    ax.set_xlabel("early-to-late recording gradient   (day 15 \u00f7 day 75)")
    ax.set_ylim(-0.8, len(out) - 0.3)
    # Keep gradient = 1 on the axis: it is the value at which there is no
    # writable-window effect at all, and truncating it off-scale makes every
    # arm look further from "no effect" than it is.
    ax.set_xlim(0.6, 21.0)
    ax.axvline(1.0, color=INK, lw=0.8, ls=":", alpha=0.65, zorder=1)
    ax.text(1.12, -0.55, "no gradient", fontsize=6.6, color=INK, alpha=0.8,
            va="center")
    ax.text(g0 + 0.25, len(out) - 0.55,
            f"no mechanism added: {g0:.1f} \u00b1 {g0e:.1f}", fontsize=7.5,
            va="center")

    ax.scatter([], [], color=BLUNT, s=34, label="blunts the window")
    ax.scatter([], [], facecolor="white", edgecolor=INK, linewidth=1.3, s=34,
               label="not resolvable at fifty populations")
    ax.scatter([], [], color=DEEPEN, s=34, label="deepens the window")
    ax.plot([], [], color=LIGHT, lw=5, label="range this arm could not resolve")
    ax.legend(loc="lower right", frameon=False, handletextpad=0.4,
              fontsize=7.2)

    mp = R.get("multiplicity")
    if mp:
        worst = max((a["p_holm"] for a in mp["arms"] if a["sig_holm"]),
                    default=None)
        if worst is not None:
            # Mid-right: the only region of the panel with no marks in it.
            # Top-right collides with the 30-day silencing arm at 19.2.
            ax.text(0.985, 0.56,
                    f"all {mp['n_sig_holm']} resolvable arms survive\n"
                    f"Holm correction over {mp['n_arms']} comparisons\n"
                    f"(largest adjusted p = {worst:.0e})",
                    transform=ax.transAxes, ha="right", va="center",
                    fontsize=7.0, color=INK,
                    bbox=dict(boxstyle="round,pad=0.42", facecolor="white",
                              edgecolor=LIGHT, linewidth=0.7))
    fig.tight_layout()
    save(fig, "Fig2_mechanisms")


# ============================================ Fig 3: what sets precision =====
def fig3():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7))

    # The two-bar version of this panel encoded two numbers and a ratio the
    # caption already stated, and could not show the heavy tail the text claims.
    # The complementary cumulative distribution of event size does both: the
    # ratio is annotated, and the tail is the shape of the curve.
    ax = axs[0]
    et = R["exact_tree"]
    pr = et["pseudoreplication"]
    esd = R["event_size_distribution"]
    sizes = np.asarray(esd["ccdf_sizes"], dtype=float)
    frac = np.asarray(esd["ccdf_frac"], dtype=float)
    ax.step(sizes, frac, where="post", color=ORANGE, lw=1.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("nuclei carrying one writing event")
    ax.set_ylabel("fraction of events\nat least this large")

    # the mean is the pseudoreplication factor: one event, this many nuclei
    ax.axvline(pr["mean"], color=INK, lw=0.9, ls=":", alpha=0.75)
    ax.text(pr["mean"] * 1.25, frac.min() * 1.6,
            f"mean {pr['mean']:.0f} nuclei\nper event (SE {pr['sem']:.1f})",
            fontsize=6.4, color=INK, va="bottom")
    # the heavy tail, as a number the bar chart could not carry
    half = esd["events_holding_half_the_nuclei"]
    ax.text(0.97, 0.97, f"{half} of {esd['n_events']:,} events\n"
                        f"hold half of all records",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.4,
            color=INK)
    panel(ax, "a", "Records are clonal")

    ax = axs[1]
    n = np.logspace(2, 4.2, 200)
    for r, c in zip(R["founder_scaling"]["sweep"], RAMP):
        pi, sd = r["pi_mean"], r["pi_sd"]
        ax.plot(n, np.sqrt(sd ** 2 + pi * (1 - pi) / n), color=c,
                label=f"{r['n_founder']} founders")
        ax.plot([r["arrays_where_noise_equals_floor"]], [np.sqrt(2) * sd],
                "o", color=c, ms=5)
    ax.set_xscale("log")
    ax.set_xlabel("recovered arrays sequenced per culture")
    ax.set_ylabel("standard error of the mark rate")
    ax.legend(frameon=False, loc="upper right")
    # Sits in the open wedge between the descending curves and the legend;
    # at the previous lower-left position it overlapped the 600-founder curve.
    ax.text(0.46, 0.46, "markers: beyond here, more\nnuclei stop helping",
            transform=ax.transAxes, fontsize=7, color=GREY)
    panel(ax, "b", "Depth saturates")

    ax = axs[2]
    cur = [c for c in R["power_curves_t"]["curves"]
           if c["arrays_per_culture"] == 1000]
    for c, col in zip(sorted(cur, key=lambda z: z["n_founder"]),
                      (RAMP[0], RAMP[2])):
        pts = [p for p in c["points"] if p.get("min_detectable_fold")]
        ax.plot([p["cultures"] for p in pts],
                [p["min_detectable_fold"] for p in pts], "o-", color=col,
                ms=4, label=f"{c['n_founder']} founders")
    ax.axvspan(0.5, 1.5, color="#f8d7da", alpha=0.75, zorder=0)
    ax.text(1.0, 2.33, "no test\npossible", ha="center", fontsize=7,
            color="#9b2c2c")
    ax.axhline(1.33, color=GREY, lw=0.9, ls=":")
    ax.text(8.4, 1.37, "1.33-fold", fontsize=7, color=GREY, ha="right")
    ax.set_xlabel("cultures per arm")
    ax.set_ylabel("smallest detectable fold change\nin mark rate")
    ax.set_xlim(0.4, 8.6)
    ax.set_ylim(1.0, 2.7)
    ax.legend(frameon=False, loc="upper right")
    panel(ax, "c", "Cultures do")

    fig.tight_layout(w_pad=2.3)
    save(fig, "Fig3_precision")


# ================================= Fig 4: the channel and its ceiling ========
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))
    c_act, e_act = 1.0, 6.0
    ceiling = R["dynamic_range"]["ceiling"]

    ax = axs[0]
    f = np.logspace(-2.3, -0.1, 200)
    for alpha, ls in ((0.2, ":"), (0.1, "-"), (0.05, "--")):
        ax.plot(f * 100, (1 - f) * (c_act + e_act) / (c_act * f * alpha),
                color=BLUE, ls=ls, label=f"α = {alpha}")
    ax.axhline(ceiling, color=ORANGE, lw=1.4)
    # The ceiling label sits at the left edge on its own line and the callout
    # drops below the line: placed any higher the two collide over "2700x".
    ax.text(0.55, ceiling * 1.30, f"readout-error ceiling  {ceiling:.0f}×",
            fontsize=7.2, color=ORANGE, ha="left", va="bottom")
    ax.axhspan(11.5, 22.6, color=GREEN, alpha=0.25, lw=0)
    ax.text(3.2, 16.5, "measured ENGRAM channels, 11.5× to 22.6×",
            fontsize=7.2, color=GREEN, ha="left", va="center")
    f0 = R["dynamic_range"]["f_on"] * 100
    D0 = R["dynamic_range"]["required"]["0.1"]
    ax.plot([f0], [D0], "*", color=INK, ms=13, zorder=5)
    ax.annotate(f"this schedule: four one-day\nwindows in 150 days,\n"
                f"needing {D0:.0f}×",
                xy=(f0, D0), xytext=(11.0, 5200), fontsize=7.2, ha="left",
                arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("share of a cell's writes falling inside\nan exposure window (%)")
    ax.set_ylabel("required channel dynamic range")
    ax.set_ylim(5, 3e4)
    leg = ax.legend(frameon=False, loc="lower left", fontsize=7,
                    title="false : true marks", title_fontsize=7,
                    borderaxespad=0.2)
    leg.get_title().set_color(INK)
    panel(ax, "a", "The schedule sets the requirement")

    ax = axs[1]
    sw = R["dynamic_range"]["sweep"]
    D = np.array([r["D"] for r in sw], float)
    ax.plot(D, [r["sentinel"]["bits"] for r in sw], "o--", color=GREY, ms=4,
            label="if misreads were detectable")
    ax.plot(D, [r["confusion"]["bits"] for r in sw], "o-", color=BLUE, ms=4,
            label="real misreads")
    ax.axvline(ceiling, color=ORANGE, lw=1.4)
    ax.text(ceiling * 0.80, 0.003, f"ceiling\n{ceiling:.0f}×", fontsize=7.5,
            color=ORANGE, ha="right")
    ax.set_xscale("log")
    ax.set_xlabel("channel dynamic range")
    ax.set_ylabel("bits per recovered cell")
    ax.legend(frameon=False, loc="upper left")
    ax.text(0.97, 0.42, "a misread clock symbol reads\nas an exposure symbol, so the\n"
            "curves part above the ceiling",
            transform=ax.transAxes, fontsize=7, color=GREY, ha="right", va="top")
    panel(ax, "b", "Readout error caps what tightening buys")

    fig.tight_layout(w_pad=2.8)
    save(fig, "Fig4_channel")


# ================================================= supplementary =============
def figS1():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7))
    al = sorted(R["array_length_multiseed"], key=lambda r: r["L"])
    L = [r["L"] for r in al]

    ax = axs[0]
    ax.plot(L, [r["frac_saturated"] for r in al], "o-", color=ORANGE,
            label="array full by day 150")
    ax.plot(L, [r["late_epoch_loss"] for r in al], "s--", color=BLUE,
            label="day-75 records lost")
    ax.axvline(12, color=GREY, lw=0.9, ls=":")
    ax.text(12.7, 0.6, "12 monomers", fontsize=7, color=GREY)
    ax.set_xscale("log")
    ax.set_xticks(L)
    ax.set_xticklabels([str(x) for x in L])
    ax.set_xlabel("array length (monomers)")
    ax.set_ylabel("fraction")
    ax.legend(frameon=False)
    panel(ax, "a", "Five-monomer arrays exhaust, and the loss is late")

    ax = axs[1]
    ax.errorbar(L, [r["bits"] for r in al], yerr=[r["bits_sem"] for r in al],
                fmt="o-", color=BLUE, ms=4, capsize=3)
    ax.axvline(12, color=GREY, lw=0.9, ls=":")
    ax.set_xscale("log")
    ax.set_xticks(L)
    ax.set_xticklabels([str(x) for x in L])
    ax.set_xlabel("array length (monomers)")
    ax.set_ylabel("bits per recovered cell")
    ax.text(0.45, 0.12, "no gain beyond twelve", transform=ax.transAxes,
            fontsize=7.5, color=GREY)
    panel(ax, "b", "And no gain from going longer")
    fig.tight_layout(w_pad=2.4)
    save(fig, "FigS1_array_length")


def figS2():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7))
    bd = sorted(R["bracket_distribution"], key=lambda r: r["alphabet"])

    ax = axs[0]
    ax.plot([r["alphabet"] for r in bd],
            [r.get("median_width_two_sided", np.nan) for r in bd], "o-",
            color=BLUE, ms=5)
    ax.axhline(150, color=GREY, lw=0.9, ls=":")
    ax.text(1.3, 142, "length of the whole experiment", fontsize=7, color=GREY)
    ax.set_xlabel("distinguishable landmark symbols\n(= orthogonal inducers needed)")
    ax.set_ylabel("median two-sided bracket (days)")
    ax.set_ylim(0, 168)
    ax.axvspan(3, 4, color=GREEN, alpha=0.18, lw=0)
    ax.text(3.5, 18, "what current\npractice offers", fontsize=7, color=GREEN,
            ha="center")
    panel(ax, "a", "Width is set by inducer count")

    ax = axs[1]
    r = bd[-1]
    parts = [r["frac_two_sided"], r["frac_lower_bound_only"],
             r["frac_no_landmark"]]
    ax.bar(range(3), parts, color=[GREEN, BLUE, GREY], width=0.6)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["landmark on\nboth sides", "lower bound\nonly",
                        "no landmark\nat all"], fontsize=7.4)
    ax.set_ylabel("share of exposure records")
    ax.set_ylim(0, max(parts) * 1.26)
    for i, v in enumerate(parts):
        ax.text(i, v + 0.015, f"{v:.2f}", ha="center", fontsize=7.5)
    panel(ax, "b", "Half of all records are unanchored")
    fig.tight_layout(w_pad=2.4)
    save(fig, "FigS2_anchoring")


# ============ Fig 5: the model against data it was not calibrated on ========
def fig5():
    """The one external test available: mouse DNA Typewriter, zygote to E13.5.

    Panel a is the consistency region, not a fit. Panel b is the shape test
    that could have failed and did not. Panel c is the predicted distribution
    against the four summary statistics Yu et al. report.
    """
    mv = R["mouse_validation"]
    o = mv["observed"]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.75))

    # ---- a: where (divisions, q) reproduces the observed occupancy --------
    ax = axs[0]
    nd = np.array(mv["divisions_grid"], float)
    qg = np.array(mv["q_grid"], float)
    grid = np.array(mv["occupancy_grid"], float)
    lo, hi = mv["occupancy_band"]

    ax.contourf(qg, nd, grid, levels=[lo, hi], colors=[ORANGE], alpha=0.30)
    cs = ax.contour(qg, nd, grid, levels=[o["sites_written_mean"]],
                    colors=[ORANGE], linewidths=1.6)
    ax.clabel(cs, fmt={o["sites_written_mean"]: "55 of 66"}, fontsize=6.6,
              inline=True)

    plo, phi = mv["divisions_plausible_range"]
    ax.axhspan(plo, phi, color=BLUE, alpha=0.12, lw=0)
    ax.axhline(plo, color=BLUE, lw=0.7, ls=":")
    ax.axhline(phi, color=BLUE, lw=0.7, ls=":")
    ax.text(0.585, phi - 1.2, "divisions a mouse embryo\nallows by E13.5",
            fontsize=6.4, color=BLUE, ha="right", va="top")

    ax.axvline(mv["q_anchor"], color=INK, lw=1.0, ls="--")
    ax.plot([mv["q_anchor"]], [mv["divisions_required_at_anchor"]], "*",
            color=INK, ms=13, zorder=5)
    ax.text(mv["q_anchor"] + 0.018, mv["divisions_required_at_anchor"] + 2.4,
            "HEK293T-anchored\nwriting rate", fontsize=6.6, color=INK)

    ax.set_xlabel("writing rate $q$ (per array per division)")
    ax.set_ylabel("divisions since the zygote")
    ax.set_xlim(qg.min(), qg.max())
    ax.set_ylim(nd.min(), nd.max())
    panel(ax, "a", "The anchor lands in the band")

    # ---- b: per-site profile, ordering alone vs the measurement ----------
    ax = axs[1]
    p = mv["profile"]
    s = np.array(p["site"], float)
    ax.plot(s, p["predicted_fraction"], "o-", color=BLUE, ms=5,
            label="model, ordering only")
    ax.plot([1, o["monomers_per_array"]],
            [p["observed_first"], p["observed_last"]], "s", color=ORANGE,
            ms=6.5, label="measured, endpoints only")
    ax.annotate("", xy=(1, p["observed_first"]),
                xytext=(o["monomers_per_array"], p["observed_last"]),
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=1.0,
                                ls=(0, (3, 2))))
    ax.set_xlabel("monomer along the array")
    ax.set_ylabel("fraction of arrays\nwritten at this monomer")
    ax.set_ylim(0, 1.09)
    ax.set_xticks(range(1, o["monomers_per_array"] + 1))
    ax.legend(frameon=False, loc="lower left", fontsize=6.8)
    ax.text(0.97, 0.97,
            f"decline along the array\nmeasured "
            f"{p['observed_decline_fold']:.1f}×, model "
            f"{p['predicted_decline_fold']:.1f}×",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.6,
            color=INK)
    panel(ax, "b", "Ordering explains the decline")

    # ---- c: predicted per-cell distribution vs reported statistics -------
    ax = axs[2]
    d = mv["predicted_cell_distribution"]
    h = np.array(d["hist"], float)
    h = h / h.sum()
    x = np.arange(h.size)
    keep = x >= 25
    ax.bar(x[keep], h[keep], width=0.92, color=BLUE, alpha=0.75, lw=0,
           label="model")
    ax.axvline(o["sites_written_mean"], color=ORANGE, lw=1.5)
    ax.axvspan(o["sites_written_mean"] - o["sites_written_sd"],
               o["sites_written_mean"] + o["sites_written_sd"],
               color=ORANGE, alpha=0.22, lw=0)
    ax.set_xlabel(f"sites written per cell (of {o['sites_total']})")
    ax.set_ylabel("share of cells")
    ax.set_xlim(28, o["sites_total"] + 1)
    ax.text(0.03, 0.97,
            f"measured {o['sites_written_mean']:.0f} ± "
            f"{o['sites_written_sd']:.0f}\nmodel {d['mean']:.0f} ± "
            f"{d['sd']:.0f}, mode {d['mode']}",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.6,
            color=INK)
    panel(ax, "c", "And the whole distribution")

    fig.tight_layout(w_pad=2.4)
    save(fig, "Fig5_external_validation")


if __name__ == "__main__":
    print("figures")
    fig2()
    fig3()
    fig4()
    fig5()
    figS1()
    figS2()
    print("done")
