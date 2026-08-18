#!/bin/bash
# Run the numbered steps, all of them or a range.
#
#   run.sh              01 through 08
#   run.sh 05           just step 05
#   run.sh 04-06        steps 04 to 06
#   run.sh 05-          step 05 onwards
#   run.sh --detach     in the background, survives logout
#
# Anything after -- is passed to every step that accepts it, e.g.
#   run.sh 05-06 -- up      one spin channel only
#
# Each step is also a standalone command; this only chains them.

set -uo pipefail
BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$BIN/$(basename "${BASH_SOURCE[0]}")"

DETACH=0; RANGE=""; PASS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --detach) DETACH=1; shift ;;
        --) shift; PASS=("$@"); break ;;
        *) RANGE="$1"; shift ;;
    esac
done

# Closing the terminal sends SIGHUP to the driver, which kills it together with
# whatever srun step it is waiting on. setsid gives the run its own session, so
# there is no controlling terminal to lose.
if [ $DETACH -eq 1 ] && [ -z "${RUN_CHILD:-}" ]; then
    log="run_$(date +%Y%m%d-%H%M%S).log"
    RUN_CHILD=1 setsid nohup "$SELF" ${RANGE:+"$RANGE"} \
        ${PASS[@]+-- "${PASS[@]}"} >> "$log" 2>&1 &
    echo "background pid $!   log $PWD/$log"
    echo "follow: tail -f $PWD/$log     stop: kill $!"
    exit 0
fi

ALL=(01_conf.py 02_inputs.py 03_nscf.sh 04_overlaps.sh \
     05_windows.py 06_wannier.sh 07_hoppings.py 08_table.py)

first=1; last=8
if [ -n "$RANGE" ]; then
    case "$RANGE" in
        *-*) first=${RANGE%-*}; last=${RANGE#*-}; [ -n "$last" ] || last=8 ;;
        *)   first=$RANGE; last=$RANGE ;;
    esac
fi
first=$((10#$first)); last=$((10#$last))
{ [ "$first" -ge 1 ] && [ "$last" -le 8 ] && [ "$first" -le "$last" ]; } \
    || { echo "ERROR: bad range '$RANGE' (steps are 01..08)"; exit 1; }

echo "=== $(basename "$PWD")   steps $(printf '%02d' "$first")-$(printf '%02d' "$last")   $(date '+%Y-%m-%d %H:%M:%S')"

for ((i = first; i <= last; i++)); do
    step="${ALL[$((i - 1))]}"
    echo
    echo "--- $step"
    case "$step" in
        01_conf.py)   python3 "$BIN/$step" -o system.conf ;;
        02_inputs.py) python3 "$BIN/$step" ;;
        05_windows.py) python3 "$BIN/$step" system.conf ;;
        07_hoppings.py) python3 "$BIN/$step" ;;
        08_table.py)
            up=$(ls avg_chem_up.dat 2>/dev/null)
            dn=$(ls avg_chem_dn.dat 2>/dev/null)
            if [ -n "$up" ] && [ -n "$dn" ]; then
                python3 "$BIN/$step" "$up" "$dn" \
                        --label "$(basename "$PWD")" --nearest
            else
                echo "  need avg_chem_up.dat and avg_chem_dn.dat"
            fi ;;
        *) "$BIN/$step" ${PASS[@]+"${PASS[@]}"} ;;
    esac
    rc=$?
    [ $rc -eq 0 ] || { echo; echo "STOPPED at $step (exit $rc)"; exit $rc; }
done

echo
echo "=== done   $(date '+%Y-%m-%d %H:%M:%S')"
