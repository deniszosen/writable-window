#!/usr/bin/env python3
"""Quantify the bracket-reader correction behind Supplementary Methods S1.

A conservative reader that bounds an exposure above with the first day
of the following landmark pulse window. That bound is valid only for one-day
pulses: if the pulse lasts d days, the exposure symbol can precede a landmark
written up to d - 1 days later, and the bound has to be the last day of the
window. This script measures the cost of the old bound at a 16-day pulse.

The trick for reproducing the old behaviour exactly: the corrected bound is
fdays[j] + fid_duration - 1, so passing the bracket reader a parameter set whose
fid_duration is 1, while the simulation itself used 16, gives back the loose
bound and nothing else changes.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dataclasses import replace
import numpy as np
from recorder import Params, simulate, readout
import analysis as A

BASE = dict(epochs=(15, 35, 55, 75), dynamic_range=1000.0, seed=1000)
SEEDS = (0, 1, 2)
DUR = 16

old, new, old_all, new_all = [], [], [], []
for sd in SEEDS:
    p = Params(**{**BASE, "fid_duration": DUR, "seed": BASE["seed"] + sd})
    s = simulate(p)
    sym, got, wd, _ = readout(s, p)
    fixed = A.brackets(s, p, sym=sym, got=got, wday=wd, reader="conservative",
                       max_rows=5000)
    buggy = A.brackets(s, replace(p, fid_duration=1), sym=sym, got=got, wday=wd,
                       reader="conservative", max_rows=5000)
    new.append(fixed["acc_two_sided"]); old.append(buggy["acc_two_sided"])
    new_all.append(fixed["acc_all"]); old_all.append(buggy["acc_all"])
    print(f"  seed {sd}: two-sided accuracy {buggy['acc_two_sided']:.3f} -> "
          f"{fixed['acc_two_sided']:.3f}, all records "
          f"{buggy['acc_all']:.3f} -> {fixed['acc_all']:.3f}, "
          f"n = {fixed['n']}, two-sided share {fixed['frac_two_sided']:.3f}")

o, n = float(np.mean(old)), float(np.mean(new))
oa, na = float(np.mean(old_all)), float(np.mean(new_all))
print(f"\n{DUR}-day pulse, {len(SEEDS)} populations. The bound only affects "
      "records with a landmark on both sides,")
print("so the two-sided accuracy is the quantity the correction moves:")
print(f"  uncorrected bound (first day of the window):  two-sided {o:.3f}, all {oa:.3f}")
print(f"  corrected bound (last day of the window):   two-sided {n:.3f}, all {na:.3f}")

claims = [("two-sided accuracy under the uncorrected bound is near 0.61",
           abs(o - 0.61) <= 0.03),
          ("the correction takes it above 0.9", n > 0.9),
          ("records without a right-hand landmark are unaffected", na > oa)]
for label, ok in claims:
    print(("  PASS  " if ok else "  FAIL  ") + label)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")
path = os.path.join(OUT, "results.json")
R = json.load(open(path))
R["bracket_fix"] = dict(pulse_duration=DUR, seeds=len(SEEDS),
                        acc_two_sided_uncorrected=o, acc_two_sided_corrected=n,
                        acc_all_uncorrected=oa, acc_all_corrected=na)
json.dump(R, open(path, "w"), indent=1)
print(f"\nmerged bracket_fix into {os.path.basename(path)}")
sys.exit(0 if all(ok for _, ok in claims) else 1)
