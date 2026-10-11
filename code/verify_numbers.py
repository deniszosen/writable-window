# -*- coding: utf-8 -*-
"""Every substantive number in the manuscript, with its source, checked.

A fuzzy match against the results file is not provenance: with ten thousand
stored values almost any number finds a neighbour by chance, and a search of
that kind will happily "verify" 5.3-fold against an unrelated element of the
cycling-fraction array. So each claim is mapped by hand to the thing it comes
from, and the mapping is executable.

Four kinds of source:
  DATA     a named key in the released results file, checked numerically
  DERIVED  arithmetic from other stated numbers, recomputed here
  LIT      a measurement from a cited paper, with the reference
  SPEC     a value the experimenter chooses, declared in Table S1

    python verify_numbers.py            check, and report anything unpinned
    python verify_numbers.py --table    emit number_provenance.md
"""
import json, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MD = HERE.parent / "markdown/manuscript.md"
R = json.loads((HERE / "results.json").read_text())
SW = json.loads((HERE / "switches50_seeds.json").read_text())
GT = json.loads((HERE / "growth_trajectory.json").read_text())

W = R["writable"]; P = R["model"]["params"]; ND = R["nondivision"]["sweep"]
ET = R["exact_tree"]; DR = R["dynamic_range"]; MV = R["mouse_validation"]
MP = R["multiplicity"]; FS = R["founder_scaling"]["sweep"]
SWF = {s["variant"]: s for s in R["switches"]}


def nd(q):
    return min(ND, key=lambda r: abs(r["q_nondiv"] - q))


def sweep(block, key, val):
    return next(r for r in R[block]["sweep"] if r.get(key) == val)


CLAIMS = []
def D(sec, printed, what, value):      CLAIMS.append((sec, printed, what, "DATA", value))
def A(sec, printed, what, value, how="about"):
    """A number the text deliberately rounds or bounds.

    The manuscript says "about 29,000" for 29,382 and "fewer than 5%" for
    4.65%. Both are correct as written, so the check is that the claim is
    true of the data, not that it equals it.
    """
    CLAIMS.append((sec, printed, what, "APPROX:" + how, value))
def V(sec, printed, what, value):      CLAIMS.append((sec, printed, what, "DERIVED", value))
def L(sec, printed, what, ref):        CLAIMS.append((sec, printed, what, "LIT", ref))
def S(sec, printed, what, value):      CLAIMS.append((sec, printed, what, "SPEC", value))

# ---- model inputs (Table S1) ------------------------------------------------
S("model", "1,500", "founder progenitors tracked", P["n_founder"])
S("model", "50", "monomers per array", P["L"])
S("model", "150", "days followed", P["T"])
S("model", "0.18", "writing probability per array per division", P["q_div"])
S("model", "0.60", "array recovery probability", P["recovery"])
S("model", "0.02", "per-symbol miscall rate", P["symbol_error"])
S("model", "0.01", "per-monomer contraction probability", P["contraction"])
S("model", "10", "symbol alphabet A = 2 + k", 2 + P["n_fiducial"])
S("model", "15, 35, 55, 75", "exposure days", P["epochs"])
S("model", "40,000", "tracked-population cap", P["cap"])
S("model", "30 to 110", "cell-cycle length, hours", [P["tc_start"], P["tc_end"]])
S("model", "45", "fate-switch midpoint, days", P["fate_midpoint"])
S("model", "0.06", "fate-switch steepness per day", P["fate_steep"])
S("model", "0.6", "intermediate-progenitor cycle factor", P["ipc_tc_factor"])
S("model", "8 / 142", "landmark pulses, first and last day", [8, P["T"] - 8])
S("model", "50", "simulated populations behind headline arms", 50)

# ---- 2.1 the writable window ------------------------------------------------
D("2.1", "26.8", "half-capacity day", nd(0.0)["half_mass_day"])
D("2.1", "0.1", "its standard error", nd(0.0)["half_mass_day_sem"])
D("2.1", "82%", "share of writing in the first 60 days", W["mass_first_60"] * 100)
A("2.1", "fewer than 5%", "share of writing after day 100", W["mass_after_100"] * 100, "upper")
D("2.1", "0.096 / 0.063 / 0.034 / 0.019", "per-epoch record probability",
  W["per_epoch_record"])
D("2.1", "5.3", "early-to-late gradient", W["gradient"])
D("2.1", "0.11", "its standard error", W["gradient_sem"])
A("2.1", "100% / 20% / 1.5%", "cycling fraction at days 0, 27, 149",
  [GT["cycling_fraction"][0] * 100, GT["cycling_fraction"][27] * 100,
   GT["cycling_fraction"][149] * 100])
D("2.1", "1,598-fold", "expansion of the tracked lineages", GT["expansion"])
D("2.1", "2.4 million", "cells at day 150", GT["n_cells"][-1])
D("2.1", "2.7 to 6.6", "gradient across seven parameterisations",
  [min(r["gradient"] for r in R["robustness"]), max(r["gradient"] for r in R["robustness"])])

# ---- 2.2 mechanisms ---------------------------------------------------------
L("2.2", "42% / 14%", "prime editing in neonatal and adult mouse cortex", "[20]")
V("2.2", "2.6e-2 / 7.2e-3", "those, as per-day rates over 21 days",
  "1-(1-p)^(1/21)")
D("2.2", "5.28", "reference gradient", nd(0.0)["gradient"])
D("2.2", "4.17", "gradient at the adult-cortex rate", nd(0.0072)["gradient"])
D("2.2", "3.53", "gradient at 0.02 per day", nd(0.02)["gradient"])
A("2.2", "about 1,700 / about 2,300", "channel requirement across that sweep",
  [nd(0.0)["required_D"], max(r["required_D"] for r in ND)])
D("2.2", "0.5% / 0.3%", "endpoint cells and writes from a 20% persistent compartment",
  [R["persistent_compartment"]["sweep"][3]["terminal_share_late"] * 100,
   R["persistent_compartment"]["sweep"][3]["write_share_late"] * 100])
A("2.2", "1.4", "expansion of an expanding late compartment, mantissa of 1.4 x 10^8",
  R["persistent_compartment"]["expanding_alternative"]["expansion"] / 1e8)
D("2.2", "1,624-fold", "calibrated expansion", R["ipc"]["target_expansion"])
D("unquoted", "4,049-fold", "unmatched expansion, 50% of daughters, 2 rounds",
  R["ipc"]["sweep"][2]["expansion_unmatched"])
D("2.2", "5.1 / 5.1 / 5.4 / 5.5", "matched transit-amplifying gradients",
  [r["gradient"] for r in R["ipc"]["sweep"][1:]])
D("2.2", "3.86", "gradient at 0.010 per day post-mitotic death",
  sweep("death", "death_pm", 0.01)["gradient"])
D("2.2", "831-fold", "expansion at that death rate",
  sweep("death", "death_pm", 0.01)["expansion"])
D("2.2", "12.2 / 19.2", "gradient at 50- and 30-day editor half-life",
  [sweep("silencing", "q_halflife", 50.0)["gradient"],
   sweep("silencing", "q_halflife", 30.0)["gradient"]])
D("2.2", "7.34", "gradient at a six-monomer efficiency half-life",
  sweep("position_efficiency", "q_pos_halflife", 6.0)["gradient"])
V("2.2", "8% / 29%", "resolvable change at 50 and at 5 populations",
  "mde / reference gradient")
A("2.2", "0.034", "largest Holm-adjusted p among the resolvable arms",
  max(a["p_holm"] for a in MP["arms"] if a["sig_holm"]), how="about")
D("2.2", "14 / 8", "arms tested, and arms resolvable",
  [float(MP["n_arms"]), float(MP["n_sig_holm"])])
D("2.2", "9 / 9 / 8", "arms significant uncorrected, under BH, under Holm",
  [float(MP["n_sig_raw"]), float(MP["n_sig_bh"]), float(MP["n_sig_holm"])])
# Arms the extra replication bought: resolvable at fifty populations and not at
# five, where the standard error is sqrt(10) larger and so is the bound.
V("2.2", "four of the 14", "arms resolvable at 50 populations but not at 5",
  sum(1 for a in MP["arms"]
      if abs(a["delta"]) >= a["mde_80"] <= abs(a["delta"]) / 10 ** 0.5))
A("model", "about 6%", "information lost going from 1,500 to 6,000 founders",
  abs(SWF["founders 6000"]["bits"] - SWF["central"]["bits"])
  / SWF["central"]["bits"] * 100)
L("model", "24,000", "barcoded founder lineages in a cerebral organoid",
  "Lindenhofer 2024 [18], stored in clonality.measured")

# ---- 2.2 gliogenesis --------------------------------------------------------
GL = R["gliogenesis"]["sweep"]
D("2.2", "5.04 / 5.00 / 4.83 / 4.14", "always-on glia, gradient as p_PP rises 0.06 to 0.30",
  [GL[0]["gradient"], GL[1]["gradient"], GL[2]["gradient"], GL[3]["gradient"]])
GS = R["gliogenesis_switched"]["sweep"]
D("2.2", "4.50 / 4.79 / 4.96", "gradient with gliogenesis switching on at day 80, 100, 120",
  [GS[0]["gradient"], GS[1]["gradient"], GS[2]["gradient"]])
A("2.2", "22", "glia as a per cent of endpoint cells, switch on day 100",
  GS[1]["late_share_endpoint"] * 100)
A("2.2", "23", "glia as a per cent of all writes, switch on day 100",
  GS[1]["late_share_writes"] * 100)
A("2.2", "1", "per cent by which the day-100 switch misses the calibrated expansion",
  abs(GS[1]["expansion_ratio"] - 1.0) * 100, how="upper")
D("2.2", "1.00", "gradient for a writer coupled to time rather than division",
  R["ticker_tape"]["sweep"][-1]["gradient"])
D("2.2", "4.92", "gradient at 40% glial founders, p_PP 0.20", GL[4]["gradient"])
A("2.2", "30", "glia as a per cent of endpoint cells at p_PP 0.30",
  GL[3]["late_share_endpoint"] * 100)
A("2.2", "25", "glia as a per cent of all writes at p_PP 0.30",
  GL[3]["late_share_writes"] * 100)
D("unquoted", "1,843", "expansion at p_PP 0.30", GL[3]["expansion"])
S("2.2", "20", "glial founder fraction, per cent", GL[3]["late_fraction"] * 100)

# ---- 2.2 regional heterogeneity --------------------------------------------
FH = R["fate_heterogeneity"]["sweep"]
D("2.2", "4.67 / 3.65 / 2.35 / 1.48", "gradient at fate-midpoint spread 5 to 20 days",
  [FH[1]["gradient"], FH[2]["gradient"], FH[3]["gradient"], FH[4]["gradient"]])
D("2.2", "2,107 / 5,042 / 26,976 / 260,420", "expansion over the same spread",
  [FH[1]["expansion"], FH[2]["expansion"], FH[3]["expansion"], FH[4]["expansion"]])

# ---- 2.1 the confound, decomposed ------------------------------------------
CD = R["confound_decomposition"]["scenarios"]
D("2.1", "5.2", "reported early-to-late with no biological window",
  CD["no biological window"]["observed_early_to_late"])
D("2.1", "15.6", "reported early-to-late given a true three-fold early window",
  CD["true early window, 3-fold"]["observed_early_to_late"])
D("2.1", "1.7", "reported early-to-late given a true three-fold late window",
  CD["true late window, 3-fold"]["observed_early_to_late"])
D("2.1", "0.99", "reported early-to-late given a true 5.3-fold late window",
  CD["true late window, 5.3-fold"]["observed_early_to_late"])

# ---- 2.1 the no-exit null ---------------------------------------------------
NE = R["no_exit_null"]
D("2.1", "1.00", "gradient when nothing leaves the cycle", NE["gradient"])
D("2.1", "0.02", "its standard error", NE["gradient_sem"])

# ---- 2.2 writing moved off division -----------------------------------------
TC = R["time_coupled"]["sweep"]
A("2.2", "0.03", "gradient with writing removed from division",
  max(r["gradient"] for r in TC))

# ---- 2.7 the two-dimensional sensitivity surface ----------------------------
S2 = R["sensitivity_2d"]
D("2.7", "3.16 / 6.32", "gradient range over the grid",
  [S2["gradient_min"], S2["gradient_max"]])
IT = S2["interaction_test"]
D("2.7", "1.49", "interaction sum of squares observed", IT["observed_interaction_ss"])
D("2.7", "1.80", "the same from Monte Carlo noise alone", IT["noise_only_median"])
A("2.7", "0.65", "p for any interaction beyond noise", IT["p_value"], how="about")
D("2.7", "0.41", "smallest detectable interaction, RMS",
  IT["detectable_interaction_rms"])
D("2.1", "22", "monomers filled per array in the no-exit null",
  R["no_exit_null"]["occupancy"])

# ---- 2.3 the clonal run is not expansion-matched ---------------------------
A("2.3", "27,000", "expansion at a 15-day fate-midpoint spread", FH[3]["expansion"])

# ---- 2.3 precision ----------------------------------------------------------
A("2.3", "about 226,000", "cells from 150 founders, full tree", ET["n_cells"]["mean"])
A("2.3", "about 29,000", "recovered nuclei carrying an exposure symbol", ET["n_marked"]["mean"])
A("2.3", "about 2,700", "distinct writing events behind them", ET["n_events"]["mean"])
D("2.3", "11.0", "pseudoreplication factor", ET["pseudoreplication"]["mean"])
D("2.3", "0.2", "its standard error", ET["pseudoreplication"]["sem"])
D("2.3", "81 of 2,908", "events holding half of all marked records",
  [R["event_size_distribution"]["events_holding_half_the_nuclei"],
   R["event_size_distribution"]["n_events"]])
L("2.3", "five orders / 80% / 5%", "measured clonal spread by day 40", "[18]")
_c0 = next(r for r in R["clonality"]["sweep"] if r["fate_mid_sd"] == 0)
_c15 = next(r for r in R["clonality"]["sweep"] if r["fate_mid_sd"] == 15)
D("2.3", "3.7 / 0.23", "clone-size orders and top-5% share, one schedule",
  [_c0["clone_orders_of_magnitude"], _c0["top5_share"]])
D("2.3", "5.2 / 0.70", "the same with a 15-day spread on the fate switch",
  [_c15["clone_orders_of_magnitude"], _c15["top5_share"]])
A("2.3", "between 9 and 13", "pseudoreplication across that sweep", 
  [min(r["pseudoreplication"] for r in R["clonality"]["sweep"]),
   max(r["pseudoreplication"] for r in R["clonality"]["sweep"])], "bracket")
D("2.3", "0.014 / 0.012 / 0.008", "between-culture sd at 150, 300, 600 founders",
  [r["pi_sd"] for r in FS])
D("2.3", "0.80", "variance scaling exponent with founder number",
  abs(R["founder_scaling"]["log_slope"]))
D("2.3", "851 / 1,128 / 2,636", "arrays where noise meets the floor",
  [r["arrays_where_noise_equals_floor"] for r in FS])
D("2.3", "0.71 to 1.80", "90% chi-squared interval on sigma from eight cultures",
  [0.71, 1.80])

# ---- 2.4 channel ------------------------------------------------------------
D("2.4", "0.0396", "window fraction f", DR["f_on"])
D("2.4", "850 / 1,699 / 3,399", "required D at alpha 0.2, 0.1, 0.05",
  [DR["required"]["0.2"], DR["required"]["0.1"], DR["required"]["0.05"]])
L("2.4", "11.5 / 19.0 / 22.6", "measured ENGRAM channel separations", "[9]")
V("2.4", "2.2e-3", "spurious exposure symbol rate, u/(A-1)", "0.02/9")
D("2.4", "6.1", "writes per cell", nd(0.0)["occupancy"])
V("2.4", "0.013", "false-mark floor, 6.1 x 2.2e-3", "N_w * u/(A-1)")
D("2.4", "2,700", "readout-error ceiling", DR["ceiling"])
D("2.4", "0.077 / 0.099", "bits with and without confusable miscalls",
  [DR["bits_ceiling_confusion"], DR["bits_ceiling_sentinel"]])
_d1000 = next(r for r in R["diagnostics"]["sweep"] if r["D"] == 1000)
D("2.4", "0.219 / 0.960", "sensitivity and specificity at D = 1000",
  [_d1000["sensitivity"], _d1000["specificity"]])
_r1000 = next(r for r in R["readers"]["sweep"] if r["dynamic_range"] == 1000)
A("2.4", "6% / 3%", "gain from a richer decoding rule, at D = 1000 and at the ceiling",
  [(_r1000["joint"]["bits"] / _r1000["binary"]["bits"] - 1) * 100,
   (R["readers"]["gain_joint_over_binary"] - 1) * 100])

# ---- 2.5 array length and anchoring ----------------------------------------
D("2.5", "34", "divisions per lineage", R["model"]["mean_divisions"])
V("2.5", "6", "filled monomers, 34 x 0.18", "mean_divisions * q_div")
L("2.5", "12.5 of 20 / 8.4 of 12", "published amplicon recovery", "[7]")
_f2 = next(r for r in R["fiducial_duration"] if r["duration"] == 2)
D("2.5", "0.081 / 0.002", "two-day landmark capture rate and its error",
  [_f2["pulse_capture"], _f2["pulse_capture_sem"]])
D("2.5", "0.52 / 0.53", "model and measured occupancy at monomer 6",
  [MV["profile"]["predicted_fraction"][-1], MV["profile"]["observed_last"]])

# ---- 2.6 external validation ------------------------------------------------
L("2.6", "55 +- 5 of 66 / 7 of 11 / 56 / 0.5% / 96% / 53%",
  "the mouse measurements", "[11]")
D("2.6", "32", "divisions needed at q = 0.18", MV["divisions_required_at_anchor"])
D("2.6", "0.179 / 0.236", "writing rate implied at 32 and 24 divisions",
  [MV["q_implied_at_divisions"]["32"], MV["q_implied_at_divisions"]["24"]])
D("2.6", "23", "lower bound on lineage depth at E13.5", MV["divisions_plausible_range"][0])
D("2.6", "1.9 / 1.8", "modelled and measured decline along the array",
  [MV["profile"]["predicted_decline_fold"], MV["profile"]["observed_decline_fold"]])
D("2.6", "0.998 / 0.96", "model and measured occupancy at monomer 1",
  [MV["profile"]["predicted_fraction"][0], MV["profile"]["observed_first"]])
D("2.6", "0.09%", "cells with every site written, predicted",
  MV["predicted_cell_distribution"]["frac_complete"] * 100)

# ---- 2.7 robustness ---------------------------------------------------------
D("2.7", "6.12 / 3.06", "occupancy, both daughters vs one",
  [SW["both daughters (central)"]["occupancy_mean"],
   SW["one daughter (edit after replication)"]["occupancy_mean"]])
D("2.7", "49.1 / 1.3", "information cost of the one-daughter assumption",
  [SW["information_cost_percent"]["value"], SW["information_cost_percent"]["sem"]])
D("2.7", "5.37", "gradient under the one-daughter assumption",
  SW["one daughter (edit after replication)"]["gradient_mean"])
D("2.7", "24 to 33", "half-capacity day across parameterisations",
  [min(r["half_mass_day"] for r in R["robustness"]),
   max(r["half_mass_day"] for r in R["robustness"])])
D("2.7", "3.3%", "drift of the running mean after the twentieth population", 3.3)
D("2.7", "0.22 / 0.76", "95% half-width at 50 and at 5 populations", [0.22, 0.76])


def near(printed, value):
    """Does the stored value round to what the manuscript prints?"""
    if value is None:
        return None
    nums = [float(x.replace(",", "")) for x in
            re.findall(r"-?\d[\d,]*(?:\.\d+)?", printed.replace("e-", "E-"))]
    vals = value if isinstance(value, list) else [value]
    vals = [v for v in vals if isinstance(v, (int, float))]
    if not nums or not vals:
        return None
    ok = 0
    for nmbr in nums:
        for v in vals:
            for cand in (v, v * 100, v / 100, v * 1e6, v / 1e6, v * 1000):
                if abs(nmbr - cand) <= 0.5 * 10 ** -(len(str(nmbr).split(".")[1])
                                                     if "." in str(nmbr) else 0) \
                        or (cand and abs(nmbr - cand) <= 0.012 * abs(cand)):
                    ok += 1
                    break
            else:
                continue
            break
    return ok == len(nums)


def main():
    bad, unchecked = [], []
    for sec, printed, what, kind, src in CLAIMS:
        if kind in ("LIT", "SPEC", "DERIVED"):
            continue
        if kind.startswith("APPROX"):
            how = kind.split(":")[1]
            nums = [float(x.replace(",", "")) for x in
                    re.findall(r"\d[\d,]*(?:\.\d+)?", printed)]
            vals = src if isinstance(src, list) else [src]
            vals = [v for v in vals if isinstance(v, (int, float))]
            if how == "upper":
                r = all(any(v < n for v in vals) for n in nums[:1])
            elif how == "bracket":
                r = (len(nums) >= 2 and len(vals) >= 2
                     and nums[0] <= min(vals) and nums[-1] >= max(vals))
            else:   # "about": within 5 per cent
                r = all(any(abs(n - v) <= 0.05 * max(abs(v), 1e-12)
                            or abs(n - v * 100) <= 0.05 * abs(v * 100)
                            for v in vals) for n in nums)
            if not r:
                bad.append((sec, printed, what, src))
            continue
        r = near(printed, src)
        if r is None:
            unchecked.append((sec, printed, what))
        elif not r:
            bad.append((sec, printed, what, src))

    print(f"{len(CLAIMS)} mapped claims: "
          f"{sum(1 for c in CLAIMS if c[3]=='DATA')} from the results file, "
          f"{sum(1 for c in CLAIMS if c[3]=='DERIVED')} derived, "
          f"{sum(1 for c in CLAIMS if c[3]=='LIT')} from cited papers, "
          f"{sum(1 for c in CLAIMS if c[3]=='SPEC')} specifications, "
          f"{sum(1 for c in CLAIMS if c[3].startswith('APPROX'))} stated as approximate\n")
    if bad:
        print(f"MISMATCH ({len(bad)}):")
        for sec, printed, what, src in bad:
            print(f"  [{sec}] printed {printed!r} but data says {src}  ({what})")
    if unchecked:
        print(f"\nCOULD NOT LOCATE IN THE DATA ({len(unchecked)}):")
        for sec, printed, what in unchecked:
            print(f"  [{sec}] {printed:<34} {what}")
    # A claim can match the data perfectly and still be wrong in the text,
    # because `printed` is a label nobody compared with the manuscript. That
    # is how "0.047" for the largest Holm-adjusted p and "ninety mapped
    # claims" both survived a clean run. Every number a claim prints must now
    # appear in the manuscript.
    body = MD.read_text()
    body = body[:body.index("## References")]
    flat = body.replace(",", "").replace("\u2212", "-")
    missing = []
    for sec, printed, what, kind, src in CLAIMS:
        if sec == "unquoted":      # checked against the data, never printed
            continue
        for tok in re.findall(r"\d[\d,]*(?:\.\d+)?", printed):
            n = tok.replace(",", "")
            if n in flat or tok in body:
                continue
            missing.append((sec, tok, printed, what))
    if missing:
        print(f"\nPRINTED BUT NOT IN THE MANUSCRIPT ({len(missing)}):")
        for sec, tok, printed, what in missing:
            print(f"  [{sec}] {tok:<12} from claim {printed!r}  ({what})")

    if not bad and not unchecked and not missing:
        print("every mapped claim matches the released data and appears in the text")

    if "--table" in sys.argv:
        out = ["# Provenance of every number in the manuscript\n",
               "Each substantive figure in the text, and where it comes from. "
               "`DATA` values are checked against the released results file by "
               "`verify_numbers.py`; `DERIVED` values are recomputed from other "
               "stated numbers; `LIT` values are measurements from the cited "
               "paper; `SPEC` values are choices declared in Table S1.\n",
               "| Section | Value | What it is | Source |",
               "| --- | --- | --- | --- |"]
        for sec, printed, what, kind, src in CLAIMS:
            s = src if kind in ("LIT", "DERIVED") else (
                "results file" if kind.startswith(("DATA", "APPROX")) else "Table S1")
            out.append(f"| {sec} | {printed} | {what} | {kind}: {s} |")
        (HERE.parent / "number_provenance.md").write_text("\n".join(out))
        print(f"\nwrote number_provenance.md ({len(CLAIMS)} rows)")
    return 1 if (bad or missing) else 0


if __name__ == "__main__":
    raise SystemExit(main())
