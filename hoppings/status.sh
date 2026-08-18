#!/bin/bash
# Status of every structure in a series, in one output.
#
# Answers, per folder: how far the pipeline got, what each stage produced, and
# -- for a Wannier90 minimisation that appears stuck -- whether it is actually
# advancing or has stopped. Reads only; changes nothing.
#
# Usage:
#   status.sh                 # every z* folder below the current directory
#   status.sh z12 z20         # named folders
#   status.sh --watch z12     # sample .wout twice, 60 s apart

set -o pipefail

WATCH=0
DIRS=()
for a in "$@"; do
    case "$a" in
        --watch) WATCH=1 ;;
        *) DIRS+=("$a") ;;
    esac
done
if [ ${#DIRS[@]} -eq 0 ]; then
    for d in z*; do [ -d "$d" ] && DIRS+=("$d"); done
fi
[ ${#DIRS[@]} -gt 0 ] || { echo "no folders found"; exit 1; }

now=$(date +%s)
age() {  # mtime of $1 as "12m ago" / "3h14m ago"
    [ -e "$1" ] || { echo "-"; return; }
    local t=$(( now - $(stat -c %Y "$1") ))
    if   [ $t -lt 3600 ]; then echo "$((t/60))m"
    elif [ $t -lt 86400 ]; then echo "$((t/3600))h$(( (t%3600)/60 ))m"
    else echo "$((t/86400))d"; fi
}
sz() { [ -e "$1" ] && du -h "$1" 2>/dev/null | cut -f1 || echo "-"; }

echo "############################################################"
echo "# host $(hostname)   $(date '+%Y-%m-%d %H:%M:%S')"
echo "############################################################"

echo
echo "== processes owned by ${USER:-$(id -un)} on THIS machine =="
ps -u "$(id -un)" -o pid,etime,pcpu,pmem,rss,comm --sort=-pcpu 2>/dev/null \
    | grep -Ei "wannier|pw\.x|pw2wan|projwfc|PID" | head -15 \
    || echo "  none"
echo "  (anything listed here is running on $(hostname), not on a compute node)"

echo
echo "== load on $(hostname) =="
printf "  cores: %s   load average:%s\n" "$(nproc 2>/dev/null || echo ?)" \
       "$(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null)"
free -h 2>/dev/null | sed 's/^/  /'
echo "  (load far above the core count means the jobs are fighting each other)"

echo
echo "== slurm queue =="
squeue -u "$(id -un)" 2>/dev/null | head -15 || echo "  squeue unavailable"

echo
echo "== disk =="
df -h . 2>/dev/null | tail -1
du -sh "${DIRS[@]}" 2>/dev/null | sed 's/^/  /'

for d in "${DIRS[@]}"; do
    echo
    echo "############################################################"
    echo "# $d"
    echo "############################################################"
    cd "$d" || continue

    sys=$(grep -E "^SYSTEM=" system.conf 2>/dev/null | cut -d= -f2)
    sys=${sys:-GdMnSi}

    if [ -f system.conf ]; then
        echo "-- settings --"
        grep -E "^(PREFIX|PROJECTIONS|MP_GRID|NBND_NSCF|DIS_WIN_LO|DIS_WIN_HI|DIS_FROZ_LO|DIS_FROZ_HI|DIS_NUM_ITER|DIS_CONV_TOL|W90_SRUN_OPTS)=" \
            system.conf | sed 's/^/   /'
    else
        echo "   no system.conf"
    fi

    echo "-- stages --"
    printf "   %-14s %s\n" "scf.out" \
        "$(grep -c 'convergence has been achieved' scf.out 2>/dev/null || echo 0) converged, E_F=$(grep 'the Fermi energy is' scf.out 2>/dev/null | tail -1 | awk '{print $5}')"
    printf "   %-14s %s\n" "nscf.out" \
        "$(grep -c 'End of band structure' nscf.out 2>/dev/null || echo 0) finished  (mtime $(age nscf.out))"

    for f in "${sys}_up.eig" "${sys}_up.amn" "${sys}_up.mmn" \
             "${sys}_dn.eig" "${sys}_dn.amn" "${sys}_dn.mmn"; do
        printf "   %-24s %8s   %s\n" "$f" "$(sz "$f")" "$(age "$f")"
    done

    for s in up dn; do
        w="${sys}_${s}.wout"
        [ -f "$w" ] || continue
        echo "-- ${sys}_${s}.wout   (mtime $(age "$w"), $(sz "$w")) --"
        if grep -q "Final State" "$w"; then
            echo "   minimisation FINISHED"
            grep "Sum of centres and spreads" "$w" | tail -1 | sed 's/^/   /'
        else
            last_dis=$(grep -- "<-- DIS" "$w" | tail -1)
            last_wan=$(grep -- "<-- CONV" "$w" | tail -1)
            if [ -n "$last_wan" ]; then
                echo "   in localisation step, last line:"
                echo "     $last_wan"
            elif [ -n "$last_dis" ]; then
                echo "   in disentanglement, last line:"
                echo "     $last_dis"
                echo "   iterations so far: $(grep -c -- '<-- DIS' "$w")"
            else
                echo "   no iteration lines yet; tail:"
                tail -4 "$w" | sed 's/^/     /'
            fi
            grep -i "Error\|Exiting\|Warning" "$w" | tail -3 | sed 's/^/   /'
        fi
        printf "   %-24s %8s   %s\n" "${sys}_${s}_hr.dat" \
               "$(sz "${sys}_${s}_hr.dat")" "$(age "${sys}_${s}_hr.dat")"
    done

    if [ $WATCH -eq 1 ]; then
        echo "-- advancing? sampling for 60 s --"
        # CPU time is the reliable signal. Fortran buffers its output, so a
        # .wout can sit unchanged for many minutes while the run is working
        # normally; conversely a process can spin without producing anything.
        pids=$(pgrep -u "$(id -un)" wannier90.x 2>/dev/null)
        declare -A t0
        for pp in $pids; do
            t0[$pp]=$(awk '{print $14+$15}' /proc/$pp/stat 2>/dev/null || echo 0)
        done
        for s in up dn; do
            w="${sys}_${s}.wout"
            [ -f "$w" ] || continue
            a=$(stat -c %s "$w")
            sleep 60
            b=$(stat -c %s "$w")
            if [ "$a" -eq "$b" ]; then
                echo "   $w: unchanged over 60 s (output is buffered, see CPU below)"
            else
                echo "   $w: grew by $((b-a)) bytes in 60 s"
            fi
        done
        for pp in $pids; do
            t1=$(awk '{print $14+$15}' /proc/$pp/stat 2>/dev/null || echo 0)
            d=$(( t1 - ${t0[$pp]:-0} ))
            # 100 clock ticks per second per core
            printf "   pid %s: %s s of CPU in 60 s -> %.1f cores busy\n" \
                   "$pp" "$((d/100))" "$(echo "$d/6000" | bc -l 2>/dev/null || echo 0)"
        done
        echo "   (near zero cores busy = truly stuck; many cores busy = running,"
        echo "    and if the machine is oversubscribed, running very slowly)"
    fi

    cd - > /dev/null
done

echo
echo "############################################################"
echo "Done. Nothing was modified."
