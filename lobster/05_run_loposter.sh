#!/bin/bash
# Stage 5. Postprocess each point with LOPOSTER (plots + summary files).
# LOPOSTER is called without arguments and reads whatever LOBSTER left behind.
set -u
CONF=${1:-lobster.conf}
source "$CONF"

for p in $POINTS; do
    d="$SERIES_DIR/$p"
    [ -f "$d/lobsterout" ] || { echo "  [$p] SKIP, no lobsterout"; continue; }
    echo "  [$p] LOPOSTER starting"
    ( cd "$d" && srun -n 1 -c "$NCORES" -p "$PARTITION" \
        "$LOPOSTER_BIN" > loposter.log 2>&1 )
    n=$(find "$d" -maxdepth 2 -newer "$d/lobsterout" \
        \( -name '*.png' -o -name '*.pdf' \) 2>/dev/null | wc -l)
    echo "  [$p] done, $n figures"
done
