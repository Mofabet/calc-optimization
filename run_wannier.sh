#!/bin/bash
# Wannier90 stage for one spin channel, split into two parts.
#
#   overlaps   wannier90.x -pp  ->  pw2wannier90.x     (produces .mmn .amn .eig)
#   minimise   wannier90.x                             (produces _hr.dat)
#
# They are separate because the disentanglement windows can only be validated
# once the .eig file exists, while .mmn/.amn/.eig do not depend on the windows
# at all -- they depend on the projections and the k-mesh. Running the overlaps
# for both spins first, then setting the windows, then minimising, means the
# expensive pw2wannier90 step runs exactly once per spin.
#
# Usage:  run_wannier.sh CONF SPIN [overlaps|minimise|all]
#           SPIN = up | dn | both

set -euo pipefail

CONF="${1:?usage: run_wannier.sh CONF SPIN [overlaps|minimise|all]}"
WHICH="${2:-both}"
STAGE="${3:-all}"

# Trim with sed, not xargs: xargs with no command runs echo, and echo eats a
# leading -n or -e as its own option, so SRUN_OPTS=-n 16 came back as "16".
get() { sed -n "s/^$1=//p" "$CONF" | head -1 \
        | sed 's/#.*//; s/^[[:space:]]*//; s/[[:space:]]*$//'; }

SYSTEM=$(get SYSTEM)
QE_BIN=$(get QE_BIN)
W90_BIN=$(get W90_BIN)
SRUN_OPTS=$(get SRUN_OPTS)
W90_SRUN_OPTS=$(get W90_SRUN_OPTS)

PW2WAN="${QE_BIN:+$QE_BIN/}pw2wannier90.x"
W90="${W90_BIN:+$W90_BIN/}wannier90.x"

have() { command -v "$1" >/dev/null 2>&1 || [ -x "$1" ]; }
have "$W90"    || { echo "ERROR: wannier90.x not found ($W90). Set W90_BIN."; exit 1; }
have "$PW2WAN" || { echo "ERROR: pw2wannier90.x not found ($PW2WAN). Set QE_BIN."; exit 1; }

# wannier90.x is single-threaded unless the build is MPI-enabled. Left unset it
# runs wherever this script runs, which on a login node means a long serial job
# on a shared machine. Set W90_SRUN_OPTS (e.g. "-n 1 -p <partition>") to push it
# onto a compute node instead.
run_w90() {
    if [ -n "$W90_SRUN_OPTS" ]; then
        echo "    srun $W90_SRUN_OPTS $W90 $*"
        # shellcheck disable=SC2086
        srun $W90_SRUN_OPTS "$W90" "$@"
    else
        "$W90" "$@"
    fi
}

# QE writes fatal errors to CRASH and to the stderr log, not to stdout, so a
# stdout file can look healthy while the run has aborted.
show_qe_error() {
    local out="$1" errlog="${2:-}" f
    for f in CRASH "$errlog" "$out"; do
        [ -n "$f" ] && [ -s "$f" ] || continue
        if grep -qi "%%%%\|Error in routine\|error #\|ERROR AT K-POINT" "$f"; then
            echo "--- $f ---"
            grep -m1 -B3 -A10 -i \
                 "%%%%\|Error in routine\|error #\|ERROR AT K-POINT" "$f"
            return 0
        fi
    done
    echo "--- no error text found; tail of $out ---"
    tail -15 "$out" 2>/dev/null
}

need() {
    [ -s "$1" ] && return 0
    echo "FAILED: $1 was not produced."
    show_qe_error "$2" "${3:-}"
    exit 1
}

do_overlaps() {
    local spin="$1" seed="${SYSTEM}_$1" ext ok=1
    for ext in mmn amn eig; do
        [ -s "${seed}.${ext}" ] || ok=0
    done
    if [ $ok -eq 1 ]; then
        echo "=== $seed : overlaps already present, reusing ==="
        return 0
    fi

    echo "=== $seed : preprocessing ==="
    rm -f CRASH
    run_w90 -pp "$seed" > "pp_${spin}.out" 2>&1
    need "${seed}.nnkp" "pp_${spin}.out"

    echo "=== $seed : pw2wannier90 ==="
    rm -f CRASH
    echo "    srun $SRUN_OPTS $PW2WAN -in pw2wan_${spin}.in"
    # shellcheck disable=SC2086
    srun $SRUN_OPTS --error="err_pw2wan_${spin}.log" \
         "$PW2WAN" -in "pw2wan_${spin}.in" > "pw2wan_${spin}.out"
    for ext in mmn amn eig; do
        need "${seed}.${ext}" "pw2wan_${spin}.out" "err_pw2wan_${spin}.log"
    done
    if [ -s CRASH ] || grep -qi "%%%%\|Error in routine" \
            "pw2wan_${spin}.out" "err_pw2wan_${spin}.log" 2>/dev/null; then
        echo "FAILED: pw2wannier90 reported an error."
        show_qe_error "pw2wan_${spin}.out" "err_pw2wan_${spin}.log"
        exit 1
    fi
}

do_minimise() {
    local spin="$1" seed="${SYSTEM}_$1"
    echo "=== $seed : minimisation ==="
    run_w90 "$seed" > "wannier_${spin}.out" 2>&1 || true
    need "${seed}_hr.dat" "${seed}.wout"

    echo "--- $seed spreads ---"
    grep "Sum of centres and spreads" "${seed}.wout" | tail -1 || true
}

spin_list() {
    case "$1" in
        up|dn) echo "$1" ;;
        both)  echo "up dn" ;;
        *) echo "ERROR: spin must be up, dn or both" >&2; exit 1 ;;
    esac
}

for s in $(spin_list "$WHICH"); do
    case "$STAGE" in
        overlaps) do_overlaps "$s" ;;
        minimise) do_minimise "$s" ;;
        all)      do_overlaps "$s"; do_minimise "$s" ;;
        *) echo "ERROR: stage must be overlaps, minimise or all"; exit 1 ;;
    esac
done
