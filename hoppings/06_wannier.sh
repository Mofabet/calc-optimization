#!/bin/bash
# 06 -- disentanglement and localisation: wannier90 -> _hr.dat
#
# Reports the two numbers that say whether the result is usable: Omega_I is the
# invariant part fixed by the subspace, Omega minus Omega_I is what localisation
# removes, and their ratio belongs in 1.1-1.5.
#
#   06_wannier.sh              both spins
#   06_wannier.sh up
#   06_wannier.sh both --iter 500    override dis_num_iter for this run

set -uo pipefail
BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$BIN/common.sh"

WHICH="both"; ITER=""
while [ $# -gt 0 ]; do
    case "$1" in
        up|dn|both) WHICH="$1"; shift ;;
        --iter) ITER="$2"; shift 2 ;;
        *) die "unknown argument '$1'" ;;
    esac
done

conf_load
W90=$(resolve "${W90_BIN:+$W90_BIN/}wannier90.x") \
    || die "wannier90.x not found; set W90_BIN"
[ -n "$OMP_THREADS" ] && export OMP_NUM_THREADS="$OMP_THREADS"

for spin in $( [ "$WHICH" = both ] && echo "up dn" || echo "$WHICH" ); do
    seed="${SYSTEM}_${spin}"
    for e in win mmn amn eig; do
        [ -s "${seed}.${e}" ] || die "${seed}.${e} missing -- run 05_overlaps.sh"
    done
    grep -q "^dis_win_max *= *0\.0*$" "${seed}.win" \
        && die "${seed}.win still has placeholder windows -- run 04_windows.py"

    keep=""
    if [ -n "$ITER" ]; then
        keep="${seed}.win.keep-$$"
        cp "${seed}.win" "$keep"
        sed -i "s/^dis_num_iter.*/dis_num_iter    = $ITER/" "${seed}.win"
        trap 'mv -f "$keep" "${seed}.win" 2>/dev/null' EXIT INT TERM
        say "$seed: dis_num_iter overridden to $ITER for this run"
    fi

    say "$seed: minimising"
    if [ -n "$W90_SRUN_OPTS" ]; then
        say "srun $W90_SRUN_OPTS $W90 $seed"
        # shellcheck disable=SC2086
        srun $W90_SRUN_OPTS "$W90" "$seed" > "wannier_${spin}.out" 2>&1 || true
    else
        "$W90" "$seed" > "wannier_${spin}.out" 2>&1 || true
    fi
    [ -n "$keep" ] && { mv -f "$keep" "${seed}.win"; trap - EXIT INT TERM; keep=""; }
    need_file "${seed}_hr.dat" "${seed}.wout"

    awk -v s="$seed" '
        /Final Omega_I/       {oi=$3}
        /Sum of centres and spreads/ {ot=$NF}
        /num_wann *:/         {nw=$NF}
        END{ if(oi>0) printf "  %s: Omega_I %.2f  Omega %.2f  ratio %.3f\n", s, oi, ot, ot/oi }
    ' "${seed}.wout"
    grep -q "Wannierisation convergence criteria satisfied" "${seed}.wout" \
        && say "$seed: localisation converged" \
        || say "$seed: localisation hit the iteration limit"
    grep -q "Maximum number of disentanglement iterations reached" "${seed}.wout" \
        && say "$seed: disentanglement hit the iteration limit"
done
