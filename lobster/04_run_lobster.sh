#!/bin/bash
# Stage 4. Run LOBSTER for every point. One task, all cores via OpenMP.
set -u
CONF=${1:-lobster.conf}
source "$CONF"

for p in $POINTS; do
    d="$SERIES_DIR/$p"
    [ -f "$d/lobsterin" ] || { echo "  [$p] SKIP, no lobsterin"; continue; }
    grep -q "convergence has been achieved" "$d/scf.out" 2>/dev/null || {
        echo "  [$p] SKIP, SCF not converged"; continue; }
    echo "  [$p] LOBSTER starting"
    ( cd "$d" && OMP_NUM_THREADS="$NCORES" srun -n 1 -c "$NCORES" \
        -p "$PARTITION" --exclusive "$LOBSTER_BIN" > lobster.log 2>&1 )
    if [ -f "$d/lobsterout" ]; then
        echo -n "  [$p] done, spilling:"
        grep -i "spilling" "$d/lobsterout" | head -4 | sed 's/^/    /'
    else
        echo "  [$p] LOBSTER FAILED, see $d/lobster.log"
    fi
done
