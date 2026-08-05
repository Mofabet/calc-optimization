#!/bin/bash
# Wannier90 stage: -pp -> pw2wannier90 -> main run, for both spin channels.
# Usage:  bin/run_wannier.sh config/GdMnSi.conf [up|dn|both]

set -euo pipefail

CONF="${1:?usage: run_wannier.sh CONF [up|dn|both]}"
WHICH="${2:-both}"

get() { grep -E "^$1=" "$CONF" | head -1 | cut -d= -f2- | sed 's/#.*//' | xargs; }

SYSTEM=$(get SYSTEM)
QE_BIN=$(get QE_BIN)
W90_BIN=$(get W90_BIN)
SRUN_OPTS=$(get SRUN_OPTS)

PW2WAN="${QE_BIN:+$QE_BIN/}pw2wannier90.x"
W90="${W90_BIN:+$W90_BIN/}wannier90.x"

command -v "$W90" >/dev/null 2>&1 || [ -x "$W90" ] || {
    echo "ERROR: wannier90.x not found ($W90). Set W90_BIN in $CONF."; exit 1; }
[ -x "$PW2WAN" ] || command -v "$PW2WAN" >/dev/null 2>&1 || {
    echo "ERROR: pw2wannier90.x not found ($PW2WAN). Set QE_BIN in $CONF."; exit 1; }

need() { [ -s "$1" ] || { echo "FAILED: $1 was not produced. See $2"; exit 1; }; }

run_spin() {
    local spin="$1" seed="${SYSTEM}_$1"
    echo "=== $seed : preprocessing ==="
    "$W90" -pp "$seed" > "pp_${spin}.out" 2>&1
    need "${seed}.nnkp" "pp_${spin}.out"

    echo "=== $seed : pw2wannier90 ==="
    # shellcheck disable=SC2086
    srun $SRUN_OPTS --error="err_pw2wan_${spin}.log" \
         "$PW2WAN" -in "pw2wan_${spin}.in" > "pw2wan_${spin}.out"
    for ext in mmn amn eig; do
        need "${seed}.${ext}" "pw2wan_${spin}.out"
    done
    if grep -qi "error\|%%%%" "pw2wan_${spin}.out"; then
        echo "FAILED: pw2wannier90 reported an error, see pw2wan_${spin}.out"
        exit 1
    fi

    echo "=== $seed : main run ==="
    "$W90" "$seed" > "wannier_${spin}.out" 2>&1
    need "${seed}_hr.dat" "${seed}.wout"

    echo "--- $seed spreads ---"
    grep -A2 "Final State" "${seed}.wout" | tail -2 || true
    grep "Sum of centres and spreads" "${seed}.wout" | tail -1 || true
}

case "$WHICH" in
    up|dn) run_spin "$WHICH" ;;
    both)  run_spin up; run_spin dn ;;
    *)     echo "ERROR: argument must be up, dn or both"; exit 1 ;;
esac

echo
echo "Done:"
ls -lh "${SYSTEM}"_*_hr.dat
