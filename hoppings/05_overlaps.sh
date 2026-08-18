#!/bin/bash
# 05 -- overlaps: wannier90 -pp, then pw2wannier90.
#
# Produces .nnkp, .mmn, .amn, .eig. These do not depend on the energy windows,
# only on the projections and the k-mesh, which is why they are a separate step
# from the minimisation: the windows can only be validated once .eig exists.
#
#   05_overlaps.sh              both spins
#   05_overlaps.sh up           one channel
#   05_overlaps.sh both --force recompute even if present

set -uo pipefail
BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$BIN/common.sh"

WHICH="both"; FORCE=0
for a in "$@"; do
    case "$a" in
        up|dn|both) WHICH="$a" ;;
        --force) FORCE=1 ;;
        *) die "unknown argument '$a'" ;;
    esac
done

conf_load
W90=$(resolve "${W90_BIN:+$W90_BIN/}wannier90.x") \
    || die "wannier90.x not found; set W90_BIN"
PW2WAN=$(resolve "${QE_BIN:+$QE_BIN/}pw2wannier90.x") \
    || die "pw2wannier90.x not found; set QE_BIN"
[ -n "$OMP_THREADS" ] && export OMP_NUM_THREADS="$OMP_THREADS"

run_w90() {
    if [ -n "$W90_SRUN_OPTS" ]; then
        # shellcheck disable=SC2086
        srun $W90_SRUN_OPTS "$W90" "$@"
    else
        "$W90" "$@"
    fi
}

for spin in $( [ "$WHICH" = both ] && echo "up dn" || echo "$WHICH" ); do
    seed="${SYSTEM}_${spin}"
    [ -f "${seed}.win" ] || die "${seed}.win missing -- run 02_inputs.py"

    # .mmn depends only on the mesh and band count, but .amn depends on the
    # projections. Reusing a stale .amn after changing PROJECTIONS would feed
    # wannier90 the wrong trial orbitals, so the header is checked.
    if [ $FORCE -eq 0 ] && [ -s "${seed}.mmn" ] && [ -s "${seed}.amn" ] \
       && [ -s "${seed}.eig" ]; then
        wantw=$(sed -n 's/^ *num_wann *= *\([0-9]*\).*/\1/p' "${seed}.win" | head -1)
        wantb=$(sed -n 's/^ *num_bands *= *\([0-9]*\).*/\1/p' "${seed}.win" | head -1)
        haveb=$(sed -n '2p' "${seed}.amn" | awk '{print $1}')
        havew=$(sed -n '2p' "${seed}.amn" | awk '{print $3}')
        if [ "$wantw" = "$havew" ] && [ "$wantb" = "$haveb" ]; then
            say "$seed: overlaps present and consistent, reusing"
            continue
        fi
        say "$seed: .amn has num_bands=$haveb num_wann=$havew but .win wants"
        say "        num_bands=$wantb num_wann=$wantw -- recomputing"
    fi

    say "$seed: preprocessing"
    rm -f CRASH
    run_w90 -pp "$seed" > "pp_${spin}.out" 2>&1 || true
    need_file "${seed}.nnkp" "pp_${spin}.out"

    say "$seed: pw2wannier90"
    rm -f CRASH
    say "srun $SRUN_OPTS $PW2WAN -in pw2wan_${spin}.in"
    # shellcheck disable=SC2086
    srun $SRUN_OPTS --error="err_pw2wan_${spin}.log" \
         "$PW2WAN" -in "pw2wan_${spin}.in" > "pw2wan_${spin}.out"
    for e in mmn amn eig; do
        need_file "${seed}.${e}" "pw2wan_${spin}.out" "err_pw2wan_${spin}.log"
    done
    if [ -s CRASH ] || grep -qi "%%%%\|Error in routine" \
            "pw2wan_${spin}.out" "err_pw2wan_${spin}.log" 2>/dev/null; then
        show_error "pw2wan_${spin}.out" "err_pw2wan_${spin}.log"
        die "pw2wannier90 reported an error"
    fi
    say "$seed: done"
done
