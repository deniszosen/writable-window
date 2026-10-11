# Provenance of every number in the manuscript

Each substantive figure in the text, and where it comes from. `DATA` values are checked against the released results file by `verify_numbers.py`; `DERIVED` values are recomputed from other stated numbers; `LIT` values are measurements from the cited paper; `SPEC` values are choices declared in Table S1.

| Section | Value | What it is | Source |
| --- | --- | --- | --- |
| model | 1,500 | founder progenitors tracked | SPEC: Table S1 |
| model | 50 | monomers per array | SPEC: Table S1 |
| model | 150 | days followed | SPEC: Table S1 |
| model | 0.18 | writing probability per array per division | SPEC: Table S1 |
| model | 0.60 | array recovery probability | SPEC: Table S1 |
| model | 0.02 | per-symbol miscall rate | SPEC: Table S1 |
| model | 0.01 | per-monomer contraction probability | SPEC: Table S1 |
| model | 10 | symbol alphabet A = 2 + k | SPEC: Table S1 |
| model | 15, 35, 55, 75 | exposure days | SPEC: Table S1 |
| model | 40,000 | tracked-population cap | SPEC: Table S1 |
| model | 30 to 110 | cell-cycle length, hours | SPEC: Table S1 |
| model | 45 | fate-switch midpoint, days | SPEC: Table S1 |
| model | 0.06 | fate-switch steepness per day | SPEC: Table S1 |
| model | 0.6 | intermediate-progenitor cycle factor | SPEC: Table S1 |
| model | 8 / 142 | landmark pulses, first and last day | SPEC: Table S1 |
| model | 50 | simulated populations behind headline arms | SPEC: Table S1 |
| 2.1 | 26.8 | half-capacity day | DATA: results file |
| 2.1 | 0.1 | its standard error | DATA: results file |
| 2.1 | 82% | share of writing in the first 60 days | DATA: results file |
| 2.1 | fewer than 5% | share of writing after day 100 | APPROX:upper: results file |
| 2.1 | 0.096 / 0.063 / 0.034 / 0.019 | per-epoch record probability | DATA: results file |
| 2.1 | 5.3 | early-to-late gradient | DATA: results file |
| 2.1 | 0.11 | its standard error | DATA: results file |
| 2.1 | 100% / 20% / 1.5% | cycling fraction at days 0, 27, 149 | APPROX:about: results file |
| 2.1 | 1,598-fold | expansion of the tracked lineages | DATA: results file |
| 2.1 | 2.4 million | cells at day 150 | DATA: results file |
| 2.1 | 2.7 to 6.6 | gradient across seven parameterisations | DATA: results file |
| 2.2 | 42% / 14% | prime editing in neonatal and adult mouse cortex | LIT: [20] |
| 2.2 | 2.6e-2 / 7.2e-3 | those, as per-day rates over 21 days | DERIVED: 1-(1-p)^(1/21) |
| 2.2 | 5.28 | reference gradient | DATA: results file |
| 2.2 | 4.17 | gradient at the adult-cortex rate | DATA: results file |
| 2.2 | 3.53 | gradient at 0.02 per day | DATA: results file |
| 2.2 | about 1,700 / about 2,300 | channel requirement across that sweep | APPROX:about: results file |
| 2.2 | 0.5% / 0.3% | endpoint cells and writes from a 20% persistent compartment | DATA: results file |
| 2.2 | 1.4 | expansion of an expanding late compartment, mantissa of 1.4 x 10^8 | APPROX:about: results file |
| 2.2 | 1,624-fold | calibrated expansion | DATA: results file |
| unquoted | 4,049-fold | unmatched expansion, 50% of daughters, 2 rounds | DATA: results file |
| 2.2 | 5.1 / 5.1 / 5.4 / 5.5 | matched transit-amplifying gradients | DATA: results file |
| 2.2 | 3.86 | gradient at 0.010 per day post-mitotic death | DATA: results file |
| 2.2 | 831-fold | expansion at that death rate | DATA: results file |
| 2.2 | 12.2 / 19.2 | gradient at 50- and 30-day editor half-life | DATA: results file |
| 2.2 | 7.34 | gradient at a six-monomer efficiency half-life | DATA: results file |
| 2.2 | 8% / 29% | resolvable change at 50 and at 5 populations | DERIVED: mde / reference gradient |
| 2.2 | 0.034 | largest Holm-adjusted p among the resolvable arms | APPROX:about: results file |
| 2.2 | 14 / 8 | arms tested, and arms resolvable | DATA: results file |
| 2.2 | 9 / 9 / 8 | arms significant uncorrected, under BH, under Holm | DATA: results file |
| 2.2 | four of the 14 | arms resolvable at 50 populations but not at 5 | DERIVED: 4 |
| model | about 6% | information lost going from 1,500 to 6,000 founders | APPROX:about: results file |
| model | 24,000 | barcoded founder lineages in a cerebral organoid | LIT: Lindenhofer 2024 [18], stored in clonality.measured |
| 2.2 | 5.04 / 5.00 / 4.83 / 4.14 | always-on glia, gradient as p_PP rises 0.06 to 0.30 | DATA: results file |
| 2.2 | 4.50 / 4.79 / 4.96 | gradient with gliogenesis switching on at day 80, 100, 120 | DATA: results file |
| 2.2 | 22 | glia as a per cent of endpoint cells, switch on day 100 | APPROX:about: results file |
| 2.2 | 23 | glia as a per cent of all writes, switch on day 100 | APPROX:about: results file |
| 2.2 | 1 | per cent by which the day-100 switch misses the calibrated expansion | APPROX:upper: results file |
| 2.2 | 1.00 | gradient for a writer coupled to time rather than division | DATA: results file |
| 2.2 | 4.92 | gradient at 40% glial founders, p_PP 0.20 | DATA: results file |
| 2.2 | 30 | glia as a per cent of endpoint cells at p_PP 0.30 | APPROX:about: results file |
| 2.2 | 25 | glia as a per cent of all writes at p_PP 0.30 | APPROX:about: results file |
| unquoted | 1,843 | expansion at p_PP 0.30 | DATA: results file |
| 2.2 | 20 | glial founder fraction, per cent | SPEC: Table S1 |
| 2.2 | 4.67 / 3.65 / 2.35 / 1.48 | gradient at fate-midpoint spread 5 to 20 days | DATA: results file |
| 2.2 | 2,107 / 5,042 / 26,976 / 260,420 | expansion over the same spread | DATA: results file |
| 2.1 | 5.2 | reported early-to-late with no biological window | DATA: results file |
| 2.1 | 15.6 | reported early-to-late given a true three-fold early window | DATA: results file |
| 2.1 | 1.7 | reported early-to-late given a true three-fold late window | DATA: results file |
| 2.1 | 0.99 | reported early-to-late given a true 5.3-fold late window | DATA: results file |
| 2.1 | 1.00 | gradient when nothing leaves the cycle | DATA: results file |
| 2.1 | 0.02 | its standard error | DATA: results file |
| 2.2 | 0.03 | gradient with writing removed from division | APPROX:about: results file |
| 2.7 | 3.16 / 6.32 | gradient range over the grid | DATA: results file |
| 2.7 | 1.49 | interaction sum of squares observed | DATA: results file |
| 2.7 | 1.80 | the same from Monte Carlo noise alone | DATA: results file |
| 2.7 | 0.65 | p for any interaction beyond noise | APPROX:about: results file |
| 2.7 | 0.41 | smallest detectable interaction, RMS | DATA: results file |
| 2.1 | 22 | monomers filled per array in the no-exit null | DATA: results file |
| 2.3 | 27,000 | expansion at a 15-day fate-midpoint spread | APPROX:about: results file |
| 2.3 | about 226,000 | cells from 150 founders, full tree | APPROX:about: results file |
| 2.3 | about 29,000 | recovered nuclei carrying an exposure symbol | APPROX:about: results file |
| 2.3 | about 2,700 | distinct writing events behind them | APPROX:about: results file |
| 2.3 | 11.0 | pseudoreplication factor | DATA: results file |
| 2.3 | 0.2 | its standard error | DATA: results file |
| 2.3 | 81 of 2,908 | events holding half of all marked records | DATA: results file |
| 2.3 | five orders / 80% / 5% | measured clonal spread by day 40 | LIT: [18] |
| 2.3 | 3.7 / 0.23 | clone-size orders and top-5% share, one schedule | DATA: results file |
| 2.3 | 5.2 / 0.70 | the same with a 15-day spread on the fate switch | DATA: results file |
| 2.3 | between 9 and 13 | pseudoreplication across that sweep | APPROX:bracket: results file |
| 2.3 | 0.014 / 0.012 / 0.008 | between-culture sd at 150, 300, 600 founders | DATA: results file |
| 2.3 | 0.80 | variance scaling exponent with founder number | DATA: results file |
| 2.3 | 851 / 1,128 / 2,636 | arrays where noise meets the floor | DATA: results file |
| 2.3 | 0.71 to 1.80 | 90% chi-squared interval on sigma from eight cultures | DATA: results file |
| 2.4 | 0.0396 | window fraction f | DATA: results file |
| 2.4 | 850 / 1,699 / 3,399 | required D at alpha 0.2, 0.1, 0.05 | DATA: results file |
| 2.4 | 11.5 / 19.0 / 22.6 | measured ENGRAM channel separations | LIT: [9] |
| 2.4 | 2.2e-3 | spurious exposure symbol rate, u/(A-1) | DERIVED: 0.02/9 |
| 2.4 | 6.1 | writes per cell | DATA: results file |
| 2.4 | 0.013 | false-mark floor, 6.1 x 2.2e-3 | DERIVED: N_w * u/(A-1) |
| 2.4 | 2,700 | readout-error ceiling | DATA: results file |
| 2.4 | 0.077 / 0.099 | bits with and without confusable miscalls | DATA: results file |
| 2.4 | 0.219 / 0.960 | sensitivity and specificity at D = 1000 | DATA: results file |
| 2.4 | 6% / 3% | gain from a richer decoding rule, at D = 1000 and at the ceiling | APPROX:about: results file |
| 2.5 | 34 | divisions per lineage | DATA: results file |
| 2.5 | 6 | filled monomers, 34 x 0.18 | DERIVED: mean_divisions * q_div |
| 2.5 | 12.5 of 20 / 8.4 of 12 | published amplicon recovery | LIT: [7] |
| 2.5 | 0.081 / 0.002 | two-day landmark capture rate and its error | DATA: results file |
| 2.5 | 0.52 / 0.53 | model and measured occupancy at monomer 6 | DATA: results file |
| 2.6 | 55 +- 5 of 66 / 7 of 11 / 56 / 0.5% / 96% / 53% | the mouse measurements | LIT: [11] |
| 2.6 | 32 | divisions needed at q = 0.18 | DATA: results file |
| 2.6 | 0.179 / 0.236 | writing rate implied at 32 and 24 divisions | DATA: results file |
| 2.6 | 23 | lower bound on lineage depth at E13.5 | DATA: results file |
| 2.6 | 1.9 / 1.8 | modelled and measured decline along the array | DATA: results file |
| 2.6 | 0.998 / 0.96 | model and measured occupancy at monomer 1 | DATA: results file |
| 2.6 | 0.09% | cells with every site written, predicted | DATA: results file |
| 2.7 | 6.12 / 3.06 | occupancy, both daughters vs one | DATA: results file |
| 2.7 | 49.1 / 1.3 | information cost of the one-daughter assumption | DATA: results file |
| 2.7 | 5.37 | gradient under the one-daughter assumption | DATA: results file |
| 2.7 | 24 to 33 | half-capacity day across parameterisations | DATA: results file |
| 2.7 | 3.3% | drift of the running mean after the twentieth population | DATA: results file |
| 2.7 | 0.22 / 0.76 | 95% half-width at 50 and at 5 populations | DATA: results file |