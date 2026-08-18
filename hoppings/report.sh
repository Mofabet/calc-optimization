#!/bin/bash
# Collect everything worth showing into one plain-text file.
#
# Run from the series root -- the directory holding bin/ and the structure
# folders. Every section is guarded: whatever has not been computed yet is
# reported as missing rather than silently omitted, so the file doubles as a
# checklist of what is still outstanding.
#
# Usage:
#   ./bin/report.sh                     # every z* folder, auto-named output
#   ./bin/report.sh -o report.txt z12 z2074
#   ./bin/report.sh --ref z2074         # folder used for the single-structure
#                                       # sections (histogram, bands, ICOHP)

set -uo pipefail

BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT=""
REF=""
DIRS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -o) OUT="$2"; shift 2 ;;
        --ref) REF="$2"; shift 2 ;;
        *) DIRS+=("$1"); shift ;;
    esac
done
if [ ${#DIRS[@]} -eq 0 ]; then
    for d in z*; do [ -d "$d" ] && DIRS+=("$d"); done
fi
[ ${#DIRS[@]} -gt 0 ] || { echo "no structure folders found"; exit 1; }
[ -n "$REF" ] || REF="${DIRS[${#DIRS[@]}-1]}"
[ -n "$OUT" ] || OUT="report_$(date +%Y%m%d-%H%M%S).txt"

py() { python3 "$BIN/$1" "${@:2}" 2>&1; }
sec() { printf '\n\n======================================================================\n%s\n======================================================================\n' "$1"; }
sub() { printf '\n--- %s ---\n' "$1"; }
miss() { printf '  [not computed: %s]\n' "$1"; }

# Everything below writes to stdout; the whole block is redirected once at the
# end, which keeps the section logic free of tee plumbing.
{

printf 'GdMnSi hoppings -- collected results\n'
printf 'generated %s on %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$(hostname)"
printf 'directory %s\n' "$PWD"
printf 'structures %s   (reference for single-structure sections: %s)\n' \
       "${DIRS[*]}" "$REF"

sec "1. SETUP AND GROUND STATE"
for d in "${DIRS[@]}"; do
    sub "$d"
    if [ -f "$d/scf.in" ]; then
        ( cd "$d" && py 01_conf.py -o /dev/null )
    else
        miss "$d/scf.in"
    fi
    if [ -f "$d/scf.out" ]; then
        e=$(grep '^!  *total energy' "$d/scf.out" | tail -1)
        [ -n "$e" ] && printf '  total energy%s\n' "${e#*=}"
    fi
done

sub "total energies side by side (for magnetic comparisons)"
printf '  %-16s %22s %14s\n' folder "total energy, Ry" "vs first, meV"
first=""
for d in "${DIRS[@]}"; do
    e=$(grep '^!  *total energy' "$d/scf.out" 2>/dev/null | tail -1 \
        | awk '{print $(NF-1)}')
    [ -n "$e" ] || { printf '  %-16s %22s\n' "$d" "-"; continue; }
    [ -n "$first" ] || first="$e"
    awk -v d="$d" -v e="$e" -v f="$first" \
        'BEGIN{printf "  %-16s %22.8f %14.1f\n", d, e, (e-f)*13605.693}'
done
printf '  (only meaningful between runs on the same cell and composition)\n'

sec "2. ORBITAL WEIGHT AT E_F  (justifies the projection set)"
found=0
for d in "$REF" "${DIRS[@]}"; do
    if ls "$d"/*.pdos_atm* >/dev/null 2>&1; then
        sub "$d"
        py pdos.py "$d" --bin 1.0
        found=1
        break
    fi
done
[ $found -eq 1 ] || miss "projwfc.x output (*.pdos_atm*) in any folder"

sec "3. EIGENVALUE DENSITY AND WINDOW PLACEMENT"
if [ -f "$REF/scf.out" ]; then
    ( cd "$REF" && py spectrum.py scf.out --histogram )
else
    miss "$REF/scf.out"
fi

sec "4. WANNIER FUNCTION QUALITY"
printf '  Omega_I is the invariant part, set by the subspace; Omega minus\n'
printf '  Omega_I is what localisation removes. Ratio 1.1-1.5 is healthy.\n'
printf '  Spread per function: 1-3 A^2 is what atom-centred d/p looks like.\n\n'
printf '  %-16s %4s %10s %10s %8s %9s %10s %s\n' \
       folder spin Omega_I Omega ratio "per WF" "worst ctr" localisation
for d in "${DIRS[@]}"; do
    for s in up dn; do
        w=$(ls "$d"/*_${s}.wout 2>/dev/null | head -1)
        [ -n "$w" ] || { printf '  %-16s %4s %10s\n' "$d" "$s" "-"; continue; }
        oi=$(grep "Final Omega_I" "$w" | tail -1 | awk '{print $3}')
        ot=$(grep "Sum of centres and spreads" "$w" | tail -1 | awk '{print $NF}')
        win=$(ls "$d"/*_${s}.win 2>/dev/null | head -1)
        nw=$(sed -n 's/^ *num_wann *= *\([0-9]*\).*/\1/p' "$win" 2>/dev/null | head -1)
        wc_=$(sed -n 's/.*worst centre-to-atom distance: *\([0-9.]*\).*/\1/p' \
              "$d/assign_${s}.log" 2>/dev/null | head -1)
        loc=$(grep -q "Wannierisation convergence criteria satisfied" "$w" \
              && echo conv || echo "NOT conv")
        awk -v d="$d" -v s="$s" -v oi="${oi:-0}" -v ot="${ot:-0}" \
            -v nw="${nw:-1}" -v w="${wc_:--}" -v l="$loc" \
            'BEGIN{if(oi>0) printf "  %-16s %4s %10.2f %10.2f %8.3f %9.3f %10s %s\n", d,s,oi,ot,ot/oi,ot/nw,w,l;
                   else printf "  %-16s %4s %10s\n", d,s,"-"}'
    done
done

sub "spread by element, per structure"
for d in "${DIRS[@]}"; do
    for s in up dn; do
        w=$(ls "$d"/*_${s}.wout 2>/dev/null | head -1)
        win=$(ls "$d"/*_${s}.win 2>/dev/null | head -1)
        [ -n "$w" ] && [ -n "$win" ] || continue
        printf '\n  == %s / %s ==\n' "$d" "$s"
        py spreads.py "$w" "$win" | sed 's/^/  /'
    done
done

sub "assignment warnings (empty is good)"
for d in "${DIRS[@]}"; do
    for s in up dn; do
        f="$d/assign_${s}.log"
        [ -f "$f" ] || continue
        n=$(grep -c "WARNING\|do NOT sit" "$f")
        printf '  %-16s %4s  %s\n' "$d" "$s" \
               "$( [ "$n" -eq 0 ] && echo clean || echo "$n warning(s)")"
        [ "$n" -gt 0 ] && sed 's/^/      /' "$f"
    done
done

sec "5. CONVERGENCE LADDER  (dis_num_iter dependence)"
laddered=0
for d in "${DIRS[@]}"; do
    ls -d "$d"/conv_* >/dev/null 2>&1 || continue
    laddered=1
    sub "$d"
    printf '  %8s %11s %11s %9s %9s  %s\n' \
           dis_iter Omega_I Omega ratio "per WF" localisation
    for c in $(ls -d "$d"/conv_* | sort -t_ -k2 -n); do
        n=${c##*conv_}
        w=$(ls "$c"/*_up.wout 2>/dev/null | head -1)
        [ -n "$w" ] || continue
        oi=$(grep "Final Omega_I" "$w" | tail -1 | awk '{print $3}')
        ot=$(grep "Sum of centres and spreads" "$w" | tail -1 | awk '{print $NF}')
        nw=$(sed -n 's/^ *num_wann *= *\([0-9]*\).*/\1/p' "$c"/*_up.win 2>/dev/null | head -1)
        loc=$(grep -q "Wannierisation convergence criteria satisfied" "$w" \
              && echo conv || echo "NOT conv")
        awk -v n="$n" -v oi="${oi:-0}" -v ot="${ot:-0}" -v nw="${nw:-1}" -v l="$loc" \
            'BEGIN{printf "  %8s %11.2f %11.2f %9.3f %9.3f  %s\n", n,oi,ot,ot/oi,ot/nw,l}'
    done
    printf '\n  t_eff per rung, eV  (column 5 of avg_chem.dat: RMS over equivalent bonds)\n'
    # Built in python rather than a shell pipeline: the natural `... | while
    # read` form runs the loop in a subshell, so the accumulated row never
    # makes it back out and the table silently comes out empty.
    python3 - "$d" <<'PYEOF'
import glob, os, re, sys
d = sys.argv[1]
rungs = sorted(glob.glob(os.path.join(d, "conv_*")),
               key=lambda p: int(re.sub(r".*conv_", "", p) or 0))
data, bonds = {}, []
for r in rungs:
    n = re.sub(r".*conv_", "", r)
    f = os.path.join(r, "avg_chem.dat")
    if not os.path.exists(f):
        continue
    data[n] = {}
    for line in open(f):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 5:
            continue
        key = (p[0], round(float(p[1]), 3))
        data[n][key] = float(p[4])
        if key not in bonds:
            bonds.append(key)
if not data:
    print("    [no avg_chem.dat in the conv_* folders]")
else:
    cols = list(data)
    print("    " + f"{'pair':<12s}{'d, A':>9s}" +
          "".join(f"{c:>10s}" for c in cols) + f"{'drift':>9s}")
    for key in sorted(bonds, key=lambda k: (k[0], k[1])):
        vals = [data[c].get(key) for c in cols]
        row = f"    {key[0]:<12s}{key[1]:9.3f}"
        for v in vals:
            row += f"{v:10.4f}" if v is not None else f"{'-':>10s}"
        got = [v for v in vals if v is not None]
        if len(got) > 1 and got[0]:
            row += f"{(got[-1] - got[0]) / got[0] * 100:+8.1f}%"
        print(row)
    print("    drift = change from the first rung to the last, in percent")
PYEOF
done
[ $laddered -eq 1 ] || miss "conv_* folders (run bin/converge.sh)"

sec "6. BAND VALIDATION  (DFT vs Wannier interpolation)"
did=0
for d in "$REF" "${DIRS[@]}"; do
    [ -f "$d/bands.out" ] || continue
    for s in up dn; do
        b=$(ls "$d"/*_${s}_band.dat 2>/dev/null | head -1)
        [ -n "$b" ] || continue
        sub "$d / $s"
        ( cd "$d" && py bands_check.py bands.out "$(basename "$b")" )
        did=1
    done
    [ $did -eq 1 ] && break
done
[ $did -eq 1 ] || miss "bands.out plus *_band.dat (see bin/bands_prep.py)"

sec "7. HOPPINGS"
for d in "${DIRS[@]}"; do
    if [ -s "$d/avg_chem_up.dat" ] && [ -s "$d/avg_chem_dn.dat" ]; then
        sub "$d"
        py 08_table.py "$d/avg_chem_up.dat" "$d/avg_chem_dn.dat" \
           --label "$d" --nearest
        printf '\n  all bonds, spin up\n'
        sed 's/^/    /' "$d/avg_chem_up.dat"
        printf '\n  symmetry check (equivalent bonds should agree), spin up\n'
        head -8 "$d/avg_atom_up.dat" 2>/dev/null | sed 's/^/    /' \
            || printf '    [avg_atom_up.dat missing]\n'
    else
        sub "$d"
        miss "$d/avg_chem_{up,dn}.dat"
    fi
done

sub "series"
py 09_series.py "${DIRS[@]}"

sec "8. LOBSTER  -ICOHP"
did=0
for d in "$REF" "${DIRS[@]}"; do
    f=$(ls "$d"/ICOHPLIST.lobster "$d"/*/ICOHPLIST.lobster 2>/dev/null | head -1)
    [ -n "$f" ] || continue
    sub "$d"
    py icohp.py "$f" --dmax 4.0
    for l in "$d"/lobsterout "$d"/*/lobsterout; do
        [ -f "$l" ] || continue
        printf '\n  charge spilling and electron count\n'
        grep -i "spilling\|charge spilling\|Total number of electrons" "$l" \
            | head -6 | sed 's/^/    /'
        break
    done
    did=1
    break
done
[ $did -eq 1 ] || miss "ICOHPLIST.lobster"

sec "9. OUTSTANDING"
printf '  Sections marked [not computed] above list what is still missing.\n'
printf '  Read section 4 before section 7: a t_eff table is only as good as\n'
printf '  the basis behind it, and the spread and assignment lines are what\n'
printf '  say whether that basis is atom-centred.\n'

printf '\n\nend of report\n'

} > "$OUT" 2>&1

echo "wrote $OUT  ($(wc -l < "$OUT") lines, $(du -h "$OUT" | cut -f1))"
echo
echo "sections:"
grep -n "^[0-9]\+\. " "$OUT" | sed 's/^/  /'
echo
echo "still missing:"
grep -c "not computed" "$OUT" | sed 's/^/  /'
grep "not computed" "$OUT" | sed 's/^/  /' | sort -u
