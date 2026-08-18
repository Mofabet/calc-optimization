#!/bin/bash
# Decide dis_num_iter by measuring t_eff, not by watching the optimiser.
#
# Runs the minimisation at a ladder of iteration counts, extracts the hoppings
# after each, and reports how much they moved. The overlaps are computed once
# and reused, so each rung costs only the minimisation.
#
# Why this and not the convergence delta: the delta describes how fast the
# optimiser is still moving, which is a property of the search, not of the
# answer. t_eff is a Frobenius norm of an orbital block, and such a norm is
# invariant under any unitary mixing of orbitals within each of the two atoms.
# The localisation step is mostly exactly that kind of mixing, so t_eff is far
# more stable than Omega -- it typically stops changing long before Omega_I
# does. Measuring it is the only way to know where.
#
# Usage, from inside a structure folder:
#     ../bin/converge.sh up 250 500 1000 2000
#     ../bin/converge.sh up 250 500 1000 2000 --detach
#
# Leaves conv_<n>/ subfolders and prints a comparison table.

set -uo pipefail

BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$BIN/$(basename "${BASH_SOURCE[0]}")"
CONF="system.conf"

DETACH=0
ARGS=()
for a in "$@"; do
    case "$a" in
        --detach) DETACH=1 ;;
        *) ARGS+=("$a") ;;
    esac
done
SPIN="${ARGS[0]:-up}"
STEPS=("${ARGS[@]:1}")
[ ${#STEPS[@]} -ge 2 ] || { echo "usage: converge.sh SPIN N1 N2 [N3 ...]"; exit 1; }

if [ "$DETACH" = 1 ] && [ -z "${CONV_CHILD:-}" ]; then
    log="converge_${SPIN}_$(date +%Y%m%d-%H%M%S).log"
    CONV_CHILD=1 setsid nohup "$SELF" ${ARGS[@]+"${ARGS[@]}"} >> "$log" 2>&1 &
    echo "running in the background"
    echo "  pid $!   log $PWD/$log"
    echo "  follow  tail -f $PWD/$log"
    exit 0
fi

get() { sed -n "s/^$1=//p" "$CONF" | head -1 \
        | sed 's/#.*//; s/^[[:space:]]*//; s/[[:space:]]*$//'; }

[ -f "$CONF" ] || { echo "ERROR: no $CONF here"; exit 1; }

# Check the helpers before starting: this run takes hours, and discovering a
# missing or renamed script only at the comparison step wastes all of it.
for h in 06_wannier.sh 07_hoppings.py compare.py; do
    [ -f "$BIN/$h" ] || { echo "ERROR: $BIN/$h not found."; \
        echo "       Copy the whole bin/ directory; the names must match."; \
        exit 1; }
done
SYSTEM=$(get SYSTEM); SEED="${SYSTEM}_${SPIN}"
PAIRS=$(get PAIRS); DMAX=$(get DMAX); TMIN=$(get TMIN)

for e in win mmn amn eig; do
    [ -s "${SEED}.${e}" ] || { echo "ERROR: ${SEED}.${e} missing -- run the overlaps first"; exit 1; }
done

echo "structure $PWD"
echo "seed      $SEED"
echo "ladder    ${STEPS[*]}"
echo

# One saved copy plus a trap: a mid-run failure or a Ctrl-C must not leave the
# .win holding a scratch iteration count that a later run would silently use.
WIN_ORIG="${SEED}.win.orig-$$"
cp "${SEED}.win" "$WIN_ORIG"
restore() { [ -f "$WIN_ORIG" ] && mv -f "$WIN_ORIG" "${SEED}.win"; }
trap restore EXIT INT TERM

for n in "${STEPS[@]}"; do
    d="conv_${n}"
    if [ -s "$d/avg_chem.dat" ]; then
        echo "=== $n iterations: already done, reusing $d ==="
        continue
    fi
    echo "=== $n iterations   $(date '+%H:%M:%S') ==="
    mkdir -p "$d"
    sed "s/^dis_num_iter.*/dis_num_iter   = $n/" "$WIN_ORIG" > "${SEED}.win"

    "$BIN/06_wannier.sh" "$SPIN" || {
        echo "  minimisation failed at $n"; continue; }

    grep "Final Omega_I" "${SEED}.wout" | tail -1
    grep "Sum of centres and spreads" "${SEED}.wout" | tail -1
    grep -q "Wannierisation convergence criteria satisfied" "${SEED}.wout" \
        && echo "  localisation converged" \
        || echo "  localisation DID NOT converge -- raise NUM_ITER"

    python3 "$BIN/07_hoppings.py" --spin "$SPIN" >/dev/null
    cp "avg_chem_${SPIN}.dat" "$d/avg_chem.dat"
    cp "hoppings_${SPIN}.dat" "$d/hoppings.dat" 2>/dev/null
    cp "assign_${SPIN}.log" "$d/assign.log" 2>/dev/null
    cp "${SEED}.wout" "${SEED}.win" "$d/" 2>/dev/null
    echo
done

echo "############################################################"
echo "# gauge quality per rung"
echo "############################################################"
echo "  Omega_I is the invariant part, fixed by the choice of subspace and"
echo "  untouched by localisation. Omega minus Omega_I is what localisation"
echo "  removes. Their ratio is the one number to watch: 1.1-1.5 is a healthy"
echo "  gauge, and a ratio that climbs with dis_num_iter means the"
echo "  localisation is losing ground, not that the model is improving."
echo
printf "  %6s %11s %11s %11s %8s  %s\n" iter Omega_I Omega ratio "per WF" localisation
for n in "${STEPS[@]}"; do
    w="conv_${n}/${SEED}.wout"
    [ -s "$w" ] || continue
    oi=$(grep "Final Omega_I" "$w" | tail -1 | awk '{print $3}')
    ot=$(grep "Sum of centres and spreads" "$w" | tail -1 | awk '{print $NF}')
    nw=$(sed -n 's/^ *num_wann *= *\([0-9]*\).*/\1/p' "conv_${n}/${SEED}.win" 2>/dev/null | head -1)
    conv=$(grep -q "Wannierisation convergence criteria satisfied" "$w" \
           && echo yes || echo NO)
    awk -v n="$n" -v oi="$oi" -v ot="$ot" -v nw="${nw:-1}" -v c="$conv" \
        'BEGIN{printf "  %6s %11.2f %11.2f %11.3f %8.3f  %s\n", n, oi, ot, ot/oi, ot/nw, c}'
done
echo
w0=$(head -1 <<< "${STEPS[*]}")

echo "############################################################"
echo "# t_eff against the most converged rung"
echo "############################################################"
last="conv_${STEPS[${#STEPS[@]}-1]}/avg_chem.dat"
for n in "${STEPS[@]}"; do
    f="conv_${n}/avg_chem.dat"
    [ -s "$f" ] || continue
    [ "$f" = "$last" ] && continue
    echo
    echo "--- $n vs ${STEPS[${#STEPS[@]}-1]} ---"
    python3 "$BIN/compare.py" "$f" "$last" --tol 1 | tail -12
done

echo
echo "Pick the smallest rung whose nearest-neighbour bonds already agree with"
echo "the top rung to within the precision you intend to quote. Apply that"
echo "value to every structure in the series."
