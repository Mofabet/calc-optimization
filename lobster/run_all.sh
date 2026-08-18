#!/bin/bash
# Run every stage in order. Stops at the first failure.
# Single stage:  ./run_all.sh 4      Range:  ./run_all.sh 2 5
set -e
CONF=lobster.conf
FROM=${1:-1}
TO=${2:-${1:-6}}

run() {
    local n=$1 label=$2
    shift 2
    if [ "$n" -ge "$FROM" ] && [ "$n" -le "$TO" ]; then
        echo "== stage $n: $label"
        "$@"
    fi
}

run 1 "make scf.in"    python3 01_make_scf.py       "$CONF"
run 2 "make lobsterin" python3 02_make_lobsterin.py "$CONF"
run 3 "run SCF"        bash    03_run_scf.sh        "$CONF"
run 4 "run LOBSTER"    bash    04_run_lobster.sh    "$CONF"
run 5 "run LOPOSTER"   bash    05_run_loposter.sh   "$CONF"
run 6 "collect series" python3 06_collect.py        "$CONF"
echo "== done"
