#!/bin/bash
# Stage 7. Quality gate. Run this BEFORE using any LOBSTER number.
# Prints the spilling, the electron count, and how negative the projected
# DOS actually gets (as a fraction of the peak height).
set -u
CONF=${1:-lobster.conf}
source "$CONF"

for p in $POINTS; do
    d="$SERIES_DIR/$p"
    [ -f "$d/lobsterout" ] || { echo "[$p] no lobsterout"; continue; }
    echo "=== $p"

    echo "-- spilling"
    grep -i -A3 "spilling" "$d/lobsterout" | grep -i "%" | sed 's/^/   /'

    echo "-- electrons / basis"
    grep -i -E "electron|setting up the basis|charge" "$d/lobsterout" \
        | head -8 | sed 's/^/   /'

    echo "-- most negative pDOS value, relative to that file's peak"
    for f in "$d"/pDOSCAR_*.loposter; do
        [ -f "$f" ] || continue
        awk -v name="$(basename "$f")" '
            # Skip the 6-line preamble and every repeated block header.
            # Header lines carry NEDOS as a bare integer in field 3;
            # data lines always have a decimal point there.
            NR <= 6 { next }
            $3 ~ /^[0-9]+$/ { next }
            NF < 3 { next }
            {
                for (i = 2; i <= NF; i++)
                    if ($i + 0 == $i) {
                        if ($i < mn) mn = $i
                        if ($i > mx) mx = $i
                    }
                seen = 1
            }
            END {
                if (seen && mx > 0)
                    printf "   %-34s min %10.4f   %6.2f %% of peak\n",
                           name, mn, -100 * mn / mx
                else
                    printf "   %-34s no data rows parsed\n", name
            }' "$f"
    done
done
