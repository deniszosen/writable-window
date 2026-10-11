# `results.json`: what is in it, which script writes it, and in what units

The released results file is one JSON object of about 10,000 leaf values. It
is written by seven scripts in sequence and every `DATA` number in the
manuscript is read from it. This file exists because a referee should not have
to infer the schema from the code.

## Conventions

- Every probability is a fraction in [0, 1], never a percentage.
- Every time is in **days** from the onset of neurogenesis, except cell-cycle
  length, which is in **hours**.
- `*_sem` is a standard error over simulated populations, `*_sd` a standard
  deviation, `*_seeds` the per-population values behind them.
- `gradient` is always the early-to-late recording gradient: the probability
  that a cell recovered at day 150 records a one-day exposure on day 15,
  divided by the same for day 75, formed per population and then averaged.
- `expansion` is endpoint cells divided by founders, for the tracked lineages
  only, not the whole aggregate.
- `occupancy` and `slot` are filled monomers per array.
- A `sweep` is a list of rows, one per setting, each carrying the swept
  parameter under its own name.

## Blocks, and the script that writes each

| Block | Written by | What it holds |
| --- | --- | --- |
| `env` | `run_all.py` | interpreter and library versions of the run |
| `_provenance` | by hand | which blocks came from which run |
| `model` | `run_all.py` | parameters, expansion, cycling fraction by day |
| `writable` | `run_all.py`, then `run_nondiv.py` | per-epoch recording probability, gradient, half-capacity day |
| `central_contrast`, `clock`, `order`, `pulse_rule` | `run_all.py` | mark rates, information, landmark behaviour |
| `dynamic_range` | `run_all.py` | channel sweep under both error models |
| `writing_rate`, `array_length`, `feasibility` | `run_all.py` | one-parameter sweeps |
| `exact_tree`, `independence` | `run_all.py` | full lineage tree, pseudoreplication |
| `variance`, `power` | `run_all.py` | between-culture variance, sample size |
| `robustness`, `switches` | `run_all.py` | alternative parameterisations, structural switches |
| `fiducials`, `fiducial_leak`, `fiducial_duration`, `landmark_alphabet` | `run_all.py`, then `run_landmarks.py` | landmark channel |
| `array_length_sweep`, `array_span`, `anchors`, `clonality`, `founder_scaling` | `run_extra.py` | length, anchoring, clonal structure |
| `nondivision` | `run_nondiv.py` | division-independent writing, 50 populations |
| `persistent_compartment`, `bracket_distribution`, `diagnostics`, `alphabet_ceiling`, `event_size_distribution`, `power_curves_normal` | `run_diagnostics.py` | compartments, decoding, diagnostics |
| `silencing`, `death`, `ipc`, `position_efficiency`, `readers`, `power_curves_t` | `run_mechanisms.py` | mechanism arms, 50 populations |
| `gliogenesis`, `fate_heterogeneity`, `confound_decomposition` | `run_glia.py` | glial compartment, regional heterogeneity, the decomposition behind Table S4 |
| `sensitivity_2d`, `time_coupled`, `no_exit_null` | `run_sensitivity.py` | interaction surface, writing moved off division, the no-exit null |
| `bracket_fix` | `check_bracket_fix.py` | the bracket-reader correction |
| `multiplicity` | `multiplicity.py` | Welch *t*, Holm and Benjamini-Hochberg over the sixteen arms |
| `mouse_validation` | `mouse_validation.py` | the external test against Yu *et al* |

## Three traps

1. **Three different between-culture standard deviations exist.** The ones used
   for every recommendation come from `founder_scaling.sweep[*].pi_sd`, computed
   at a raised population cap. The entry in `variance.between_culture` computed
   at the default cap of 40,000 is inflated roughly threefold by subsampling and
   enters nothing. Table S2 of the supplementary material states which is which.
2. **`writable.gradient` is overwritten.** `run_all.py` writes a three-population
   value; `run_nondiv.py` overwrites it with the 50-population value that the
   manuscript quotes. Running `run_all.py` alone leaves the smaller number in
   place.
3. **Population counts are not uniform across scripts.** `run_nondiv.py`,
   `run_diagnostics.py`, `run_mechanisms.py`, `run_glia.py` and `run_sensitivity.py` default to
   50; `run_landmarks.py` uses 5, `run_extra.py` 3, and `run_all.py` 3 and 6.
   `sensitivity_2d` uses 10 per grid cell. Each block carries its own `n_seeds`
   where it is not the default.

## Which seed produces which figure

Every figure is drawn from this file rather than from a fresh simulation, so no
figure depends on a seed at drawing time. The underlying runs all start from
`seed = 1000` and use offsets 0, 1, 2, … up to the population count for that
block. Two panels show a single realisation rather than a mean, and both use
offset 0, that is `seed = 1000`: figure 3a, the distribution of nuclei per
writing event, and figure 1a and 1b, the growth trajectory, which
`growth_trajectory.py` writes to `growth_trajectory.json`.
