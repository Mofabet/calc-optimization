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
        [ -s "${seed}.${e}" ] || die "${seed}.${e} missing -- run 04_overlaps.sh"
    done
    grep -q "^dis_win_max *= *0\.0*$" "${seed}.win" \
        && die "${seed}.win still has placeholder windows -- run 05_windows.py"

    # wannier90 answers a header mismatch with "param_read: mismatch in
    # <seed>.eig" and nothing else, so check it here where the two numbers can
    # be named. This bites whenever the .win is regenerated after the overlaps
    # were computed: step 04 checks it, running 06 alone did not.
    wantw=$(sed -n 's/^ *num_wann *= *\([0-9]*\).*/\1/p' "${seed}.win" | head -1)
    wantb=$(sed -n 's/^ *num_bands *= *\([0-9]*\).*/\1/p' "${seed}.win" | head -1)
    haveb=$(tail -1 "${seed}.eig" | awk '{print $1}')
    havek=$(tail -1 "${seed}.eig" | awk '{print $2}')
    amnb=$(sed -n '2p' "${seed}.amn" | awk '{print $1}')
    amnw=$(sed -n '2p' "${seed}.amn" | awk '{print $3}')
    if [ "$wantb" != "$haveb" ] || [ "$wantb" != "$amnb" ] \
       || [ "$wantw" != "$amnw" ]; then
        echo "ERROR: ${seed}.win does not match the overlaps." >&2
        echo "         .win  num_bands=$wantb  num_wann=$wantw" >&2
        echo "         .eig  num_bands=$haveb  ($havek k-points)" >&2
        echo "         .amn  num_bands=$amnb  num_wann=$amnw" >&2
        echo "       The .win was regenerated after the overlaps were made." >&2
        echo "       Either pin NBND_NSCF=$haveb in system.conf and rerun" >&2
        echo "       steps 02 and 05, which keeps the existing overlaps, or" >&2
        echo "       rerun step 04 with --force to recompute them." >&2
        exit 1
    fi

    # The overlaps carry their own band and function counts. wannier90 checks
    # them too, but only after the job has reached a compute node, so a stale
    # .win costs a queue slot to discover. Check here instead.
    wb=$(sed -n 's/^ *num_bands *= *\([0-9]*\).*/\1/p' "${seed}.win" | head -1)
    ww=$(sed -n 's/^ *num_wann *= *\([0-9]*\).*/\1/p' "${seed}.win" | head -1)
    eb=$(tail -1 "${seed}.eig" | awk '{print $1}')
    ab=$(sed -n '2p' "${seed}.amn" | awk '{print $1}')
    aw=$(sed -n '2p' "${seed}.amn" | awk '{print $3}')
    if [ "$wb" != "$eb" ] || [ "$wb" != "$ab" ] || [ "$ww" != "$aw" ]; then
        echo "ERROR: ${seed}.win does not match the overlaps."
        echo "         .win  num_bands $wb   num_wann $ww"
        echo "         .eig  num_bands $eb"
        echo "         .amn  num_bands $ab   num_wann $aw"
        echo "       The .win was regenerated after the overlaps were computed."
        echo "       Either set NBND_NSCF=$eb in system.conf and rerun step 02,"
        echo "       or recompute the overlaps:  04_overlaps.sh $spin --force"
        exit 1
    fi

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
    # Plain if, not `grep && say`: as the last statement of the loop a failing
    # grep becomes the script's exit status, and a converged run -- where the
    # "limit reached" message is absent -- would then look like a failure and
    # stop the driver before the hoppings are extracted.
    if grep -q "Wannierisation convergence criteria satisfied" "${seed}.wout"; then
        say "$seed: localisation converged"
    else
        say "$seed: localisation hit the iteration limit"
    fi
    if grep -q "Maximum number of disentanglement iterations reached" \
              "${seed}.wout"; then
        say "$seed: disentanglement hit the iteration limit"
    else
        say "$seed: disentanglement converged"
    fi
done

exit 0