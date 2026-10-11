"""Figure 1: the population grows, and the writable window closes inside it.

Four panels, every one of them a model output:

  a  cell number in the tracked founder lineages, over time
  b  the cycling fraction against the cumulative share of all writing
  c  the three division outcomes
  d  what the closing window costs an exposure given at each of four days

No geometry. The model is well mixed and has no positions, so it has no
diameter to report: converting its cell count into the diameter of a sphere
would describe the recorder-carrying clone mass packed alone, which is not
the organoid and is not a quantity anyone measures. Where the physical size
of an organoid matters, in the Discussion's argument about delivery into the
interior, it is cited from the literature instead.

Where each quantity comes from:
  cell number         released simulator, three seeds (growth_trajectory.json)
  cycling fraction    released simulator, same run
  cumulative writing  results.json, writable.write_day_density
  per-epoch rates     results.json, writable.per_epoch_record

Run: python growth_trajectory.py && python fig1_timeline.py
"""
import json
import os

import matplotlib as mpl

mpl.use("Agg")
import numpy as np                                                  # noqa: E402
from matplotlib import pyplot as plt                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "figures")
R = json.load(open(os.path.join(HERE, "results.json")))
G = json.load(open(os.path.join(HERE, "growth_trajectory.json")))

INK, BLUE, ORANGE = "#1a202c", "#2b6cb0", "#c05621"
GREY = "#718096"

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
})


def main():
    day = np.array(G["day"])
    n = np.array(G["n_cells"], float)
    cf = np.array(G["cycling_fraction"], float)
    cum = np.cumsum(np.array(R["writable"]["write_day_density"], float))
    half = int(np.searchsorted(cum, 0.5))
    epochs = list(R["writable"]["epochs"])

    fig = plt.figure(figsize=(7.2, 4.4))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.30)

    # ---- (a) the population grows ------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(day, n, color=INK, lw=1.8)
    ax.set_yscale("log")
    ax.set_xlim(0, 149); ax.set_xlabel("day of differentiation")
    ax.set_ylabel("cells in the tracked\nfounder lineages")
    ax.text(0.97, 0.10, f"{G['expansion']:.0f}-fold, from "
            f"{int(n[0]):,} founders\nto {n[-1]/1e6:.1f} million cells",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8,
            color=GREY, linespacing=1.35)
    ax.set_title("a   The population keeps growing", loc="left",
                 fontsize=9, fontweight="bold", pad=5)

    # ---- (b) but it stops being able to write ------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.fill_between(day, 0, cf, color=BLUE, alpha=0.85, lw=0)
    ax.plot(day, cum, color=ORANGE, lw=1.8)
    ax.set_xlim(0, 149); ax.set_ylim(0, 1.0)
    ax.set_xlabel("day of differentiation"); ax.set_ylabel("fraction")
    for d in epochs:
        ax.axvline(d, color=INK, lw=0.7, ls=(0, (2, 2)), alpha=0.5, zorder=4)
    ax.plot([half, half], [0, 0.5], color=ORANGE, lw=0.9, ls=":")
    ax.annotate(f"half of all writing\nby day {half}", xy=(half, 0.5),
                xytext=(half + 20, 0.66), fontsize=6.8, color=ORANGE,
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.7))
    ax.text(96, 0.055, "still cycling", color="white", fontsize=6.8,
            fontweight="bold", ha="center")
    ax.text(0.97, 0.90, "cumulative writing", transform=ax.transAxes,
            ha="right", fontsize=6.8, color=ORANGE, fontweight="bold")
    ax.set_title("b   Writing finishes before the growth does",
                 loc="left", fontsize=9, fontweight="bold", pad=5)

    # ---- (c) why: divisions turn neurogenic --------------------------------
    ax = fig.add_subplot(gs[1, 0])
    p = R["model"]["params"]
    lg = lambda lo, hi: lo + (hi - lo) / (1 + np.exp(-p["fate_steep"] *
                                                     (day - p["fate_midpoint"])))
    pp = lg(p["pp_start"], p["pp_end"]); nn = lg(p["nn_start"], p["nn_end"])
    ax.plot(day, pp, color=BLUE, lw=1.5, label="both daughters cycle")
    ax.plot(day, nn, color=ORANGE, lw=1.5, ls="--", label="both exit")
    ax.plot(day, 1 - pp - nn, color=GREY, lw=1.5, ls=":", label="one of each")
    ax.set_xlim(0, 149); ax.set_ylim(0, 1.0)
    ax.set_xlabel("day of differentiation")
    ax.set_ylabel("probability\nper division")
    ax.legend(frameon=False, fontsize=6.4, loc="upper right",
              bbox_to_anchor=(1.0, 0.80), handlelength=1.6, labelspacing=0.22)
    ax.text(0.03, 0.44, "cell cycle lengthens\n%.0f h to %.0f h"
            % (p["tc_start"], p["tc_end"]),
            transform=ax.transAxes, fontsize=6.4, color=GREY, va="center")
    ax.set_title("c   Divisions turn neurogenic", loc="left", fontsize=9,
                 fontweight="bold", pad=5)

    # ---- (d) what it costs an exposure -------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    pe = np.array(R["writable"]["per_epoch_record"], float)
    se = np.array(R["writable"]["per_epoch_sem"], float)
    x = np.arange(len(pe))
    ax.bar(x, pe, yerr=se, capsize=2.5, width=0.62, color=BLUE,
           error_kw=dict(lw=0.9, ecolor=INK))
    ax.set_xticks(x); ax.set_xticklabels([f"d{d}" for d in epochs], fontsize=7)
    ax.set_ylabel("P(exposure recorded)")
    ax.set_xlabel("day the agent is given")
    ax.set_ylim(0, max(pe) * 1.40)
    g = R["writable"]["gradient"]
    y = max(pe) * 1.16
    ax.annotate("", xy=(0, y), xytext=(3, y),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.0))
    ax.text(1.5, y * 1.05, f"{g:.1f}-fold", ha="center", fontsize=7.4,
            fontweight="bold")
    ax.set_title("d   And an exposure given late is mostly missed", loc="left",
                 fontsize=9, fontweight="bold", pad=5)

    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(OUT, f"Fig1_timeline.{ext}"))
    plt.close(fig)
    print("wrote Fig1_timeline")
    print(f"  day 0: {n[0]:,.0f} cells, {cf[0]*100:.0f}% cycling")
    print(f"  day {len(day)-1}: {n[-1]:,.0f} cells, {cf[-1]*100:.1f}% cycling")
    print(f"  half of all writing by day {half}")


if __name__ == "__main__":
    main()
