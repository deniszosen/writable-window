# The writable window — simulator, analysis and released results

Code and data behind *The writable window: division-coupled DNA recorders carry
a built-in early bias in human cortical organoids*.

Denis Zosen, Norwegian Centre for Molecular Biosciences and Medicine (NCMBM),
University of Oslo. ORCID [0000-0003-4722-9921](https://orcid.org/0000-0003-4722-9921).

No new experimental data were generated. Every input is a parameter taken from a
published measurement, a specification the experimenter chooses, or a swept
assumption; Table S1 of the supplementary material says which, for each one.

## Contents

    model/src/recorder.py       the generative model: branching population + sequential writing
    model/src/analysis.py       quantities computed from a simulated population
    model/src/power.py          between-culture variance and the noncentral t power solve
    model/src/run_all.py        entry point; writes model/out/results.json
    model/src/run_*.py          the remaining sweep drivers
    model/src/run_glia.py       glial compartment, fate-midpoint spread, and the
                                decomposition of the confound (added in review)
    model/src/check_keys.py     every key a consumer reads must exist in the results file
    model/src/test_recorder.py  unit tests on the model itself, not on its output
    model/src/run_sensitivity.py  interaction surface, writing off division, the no-exit null
    model/out/results.json      the released results file

    code/verify_numbers.py      118 mapped claims in this manuscript against the data
    code/review.py              structure, duplication, citations, abbreviations, length
    code/multiplicity.py        Holm and Benjamini-Hochberg over the fourteen arms
    code/interaction_test.py    is the sensitivity surface additive, and how much
                                interaction could it have seen
    code/mouse_validation.py    the external test against the DNA Typewriter mouse
    code/growth_trajectory.py   the cell-number and cycling-fraction trajectory
    code/figures.py             figures 2 to 5 and S1, S2
    code/fig1_timeline.py       figure 1
    code/md2docx.py             renders the markdown to .docx
    code/results.json           same file as model/out/results.json, read from beside the scripts

    markdown/                   manuscript and supplementary source
    figures/                    every figure as PDF, SVG and PNG
    number_provenance.md        each substantive number and where it comes from
    RESULTS_SCHEMA.md           what is in results.json, in what units, written by which script

## Rebuild everything

Run everything with one command. The order is not obvious and prose
instructions for it kept going stale, so it is a script:

    pip install -r requirements.txt
    ./run_everything.sh --check     # verify every path, runs nothing
    ./run_everything.sh             # the real thing, about two hours

`N_POPULATIONS=3 ./run_everything.sh` gives a fast pass for smoke-testing,
which will **not** reproduce the published numbers. The per-script population
counts are not uniform — 50 for `run_nondiv`, `run_diagnostics`,
`run_mechanisms`, `run_glia`
and `run_sensitivity`, 10 per grid cell for the sensitivity surface, 5 for
`run_landmarks`, 3 for `run_extra`, and 3 and 6 for `run_all` — and changing
them changes the headline results.

Two things make the order non-obvious, and both have caught people, including
while preparing this deposit:

- `run_glia.py` reads `multiplicity.reference_gradient`, so `multiplicity.py`
  has to run partway through rather than with the other post-processing at the
  end. Running every simulation first and leaving all the analysis to the end
  fails with a `KeyError`.
- The simulation scripts live in `model/src` and read `../out/results.json`;
  the analysis scripts live in `code/` and read `results.json` beside
  themselves. The file is copied between the two at three points.

Every check exits non-zero on failure, so the sequence runs unattended.

## What was checked, and when

Against the released `results.json` on 9 October 2026, under Python 3.14.7,
numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2:

| Check | Result |
| --- | --- |
| `verify_numbers.py` | 118 mapped claims, all matching |
| `review.py` | clean |
| `multiplicity.py` | stored block matches a fresh run |
| `interaction_test.py` | stored block matches a fresh run |
| `test_recorder.py` | 12 unit tests on the model, all passing |
| `mouse_validation.py` | reproduces the quantities in section 2.6 |
| `check_keys.py` | every consumed key present |
| all seven figures | regenerate from the released results file |

**The whole pipeline was also re-run from an empty `out/` directory**, in the
order below, under the versions above, on 10 October 2026. Of the 10,072
numeric values shared with the released file, **10,071 reproduce to within
1 × 10⁻⁹ (99.99%)** and 9,794 are bit-identical. The headline early-to-late
gradient comes back as `5.283359245525049`, identical to the last bit, despite
the five-year interpreter gap and the numpy major version. The released
`results.json` **is** that regenerated file, so the data deposited here is the
output of the code deposited beside it.

The single exception is
`power_curves_t.curves[0].points[1].min_detectable_fold_sd_high`, differing
by 1.4 × 10⁻⁹ relative. It is a fold change computed as net/(net − delta) at
two cultures per arm at the pessimistic end of the sigma interval, where
net − delta is about 0.0015; the two solvers agree on delta itself to roughly
1 × 10⁻¹¹ relative and the division by a near-zero denominator amplifies that.
The value is not quoted in the manuscript.

**Two portability fixes were needed**, both the same scipy 1.18 behaviour:
`scipy.stats.nct.cdf` returns NaN for the opposite-tail term at moderate and
large positive noncentrality instead of underflowing to zero, which aborts the
root finder. It broke `figures.py` (the minimum-detectable-effect bars)
and `power.py` (`detectable_delta`), the latter making `run_mechanisms.py` exit 1 and
so blocking `silencing`, `death`, `ipc`, `position_efficiency` and
`power_curves_t` from regenerating at all. A noncentral *t* is (Z + ncp)/sqrt(V/df) with a
positive denominator, so T < 0 requires Z < −ncp and the term is bounded above
by Phi(−ncp); every noncentrality at which scipy returns NaN here puts that
bound below 1e-17, so treating it as zero is exact in double precision. No
reported quantity moves.

## Licence

Code under MIT (`LICENSE`). Manuscript text, supplementary material and figures
under CC BY 4.0.
