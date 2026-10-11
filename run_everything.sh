#!/bin/bash
# Regenerate everything from an empty directory, in the only order that works.
#
# Two things make the order non-obvious, and both have bitten people:
#   * run_glia.py reads multiplicity.reference_gradient, so the correction has
#     to be computed partway through rather than at the end;
#   * the simulation scripts live in model/src and read ../out/results.json,
#     while the analysis scripts live in code/ and read results.json beside
#     themselves, so the file is copied between them.
#
#   ./run_everything.sh --check   verify every path without running anything
#   ./run_everything.sh           the real thing, about two hours
#
# N_POPULATIONS=3 ./run_everything.sh gives a fast pass that will NOT
# reproduce the published numbers.
set -euo pipefail
cd "$(dirname "$0")"
ROOT="$PWD"
CHECK=0
[ "${1:-}" = "--check" ] && CHECK=1

SIM="run_all.py run_extra.py run_nondiv.py run_landmarks.py run_diagnostics.py run_mechanisms.py"
SIM2="run_glia.py run_sensitivity.py check_bracket_fix.py check_keys.py"
POST="multiplicity.py interaction_test.py mouse_validation.py growth_trajectory.py"
FIGS="fig1_timeline.py figures.py"
VERIFY="verify_numbers.py review.py"

for s in $SIM $SIM2 test_recorder.py; do
  [ -f "model/src/$s" ] || { echo "MISSING model/src/$s"; exit 1; }
done
for s in $POST $FIGS $VERIFY; do
  [ -f "code/$s" ] || { echo "MISSING code/$s"; exit 1; }
done
if [ "$CHECK" = 1 ]; then
  echo "all $(echo $SIM $SIM2 $POST $FIGS $VERIFY test_recorder.py | wc -w) referenced scripts exist"
  exit 0
fi

run() { echo "== $* ($(date +%H:%M:%S))"; "$@"; }

cd "$ROOT/model/src"
for s in $SIM; do run python3 "$s"; done

# the correction has to exist before run_glia.py
cp "$ROOT/model/out/results.json" "$ROOT/code/"
cd "$ROOT/code" && run python3 multiplicity.py --write
cp "$ROOT/code/results.json" "$ROOT/model/out/"

cd "$ROOT/model/src"
for s in $SIM2; do run python3 "$s"; done

# everything after this point is post-processing of the finished file
cp "$ROOT/model/out/results.json" "$ROOT/code/"
cd "$ROOT/code"
run python3 multiplicity.py --write        # the arm family changed
run python3 interaction_test.py --write
run python3 mouse_validation.py --write
run python3 growth_trajectory.py
cp "$ROOT/code/results.json" "$ROOT/model/out/"

for s in $FIGS; do run python3 "$s"; done
for s in $VERIFY; do run python3 "$s"; done
run python3 "$ROOT/model/src/test_recorder.py"
echo "== done ($(date +%H:%M:%S))"
