#!/bin/bash
# Stage 3. Run the LOBSTER-ready SCF for every point, one after another.
set -u
CONF=${1:-lobster.conf}
source "$CONF"

for p in $POINTS; do
    d="$SERIES_DIR/$p"
    [ -f "$d/scf.in" ] || { echo "  [$p] SKIP, no scf.in"; continue; }
    if grep -q "convergence has been achieved" "$d/scf.out" 2>/dev/null; then
        echo "  [$p] SCF already converged, skipping"
        continue
    fi
    echo "  [$p] SCF starting"
    ( cd "$d" && srun --mpi=pmix -n "$NCORES" -p "$PARTITION" --exclusive \
        "$QE_BIN/pw.x" -i scf.in > scf.out 2>&1 )
    if grep -q "convergence has been achieved" "$d/scf.out"; then
        echo "  [$p] SCF done"
    else
        echo "  [$p] SCF FAILED, see $d/scf.out"
    fi
done
