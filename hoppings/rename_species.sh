#!/bin/bash
# Rename a QE species label in scf.in and in the saved XML.
#
# Why this is needed: QE 7.5 segfaults when restarting a DFT+U calculation
# whose species labels it cannot resolve to a chemical element. The trace is
#   read_file_new -> post_xml_init -> init_hubbard -> offset_atom_wfc
# and it hits every post-processing tool, pw2wannier90 and projwfc alike,
# while the SCF itself runs fine because it reads the HUBBARD card from the
# input rather than from the XML. Labels of the form Mn1 / Mn2 resolve;
# Mn_up / Mn_dn do not.
#
# Only text is edited. The charge density and wavefunctions contain no species
# labels, so nothing has to be recomputed.
#
# Usage:
#   rename_species.sh Mn_up Mn1 Mn_dn Mn2 -- z12 z18 z20
#   rename_species.sh Mn_up Mn1 Mn_dn Mn2 -- z*        # every folder
#   rename_species.sh --check -- z12 z18 z20           # report only

set -euo pipefail

CHECK=0
PAIRS=()
DIRS=()
mode=pairs
for a in "$@"; do
    case "$a" in
        --check) CHECK=1 ;;
        --) mode=dirs ;;
        *) if [ "$mode" = pairs ]; then PAIRS+=("$a"); else DIRS+=("$a"); fi ;;
    esac
done

[ ${#DIRS[@]} -gt 0 ] || { echo "usage: $0 OLD NEW [OLD NEW ...] -- DIR [DIR ...]"; exit 1; }
if [ $CHECK -eq 0 ] && [ $(( ${#PAIRS[@]} % 2 )) -ne 0 ]; then
    echo "ERROR: rename pairs must come in twos (OLD NEW)"; exit 1
fi

for d in "${DIRS[@]}"; do
    [ -d "$d" ] || continue
    echo "=== $d"

    scf="$d/scf.in"
    if [ ! -f "$scf" ]; then echo "    no scf.in, skipped"; continue; fi

    xml=$(ls -1 "$d"/out/*.save/data-file-schema.xml 2>/dev/null | head -1 || true)

    echo -n "    species in scf.in : "
    sed -n '/ATOMIC_SPECIES/,/^[[:space:]]*[A-Z_]\{4,\}/p' "$scf" \
        | awk 'NF>=3 && $1!="ATOMIC_SPECIES" {printf "%s ", $1}'
    echo
    if [ -n "$xml" ]; then
        echo -n "    species in XML    : "
        grep -o '<species name="[^"]*"' "$xml" | sed 's/.*name="//; s/"//' \
            | tr '\n' ' '
        echo
    else
        echo "    no save XML yet (nscf not run)"
    fi

    [ $CHECK -eq 1 ] && continue

    ts=$(date +%Y%m%d-%H%M%S)
    cp "$scf" "$scf.bak-$ts"
    [ -n "$xml" ] && cp "$xml" "$xml.bak-$ts"

    i=0
    while [ $i -lt ${#PAIRS[@]} ]; do
        old="${PAIRS[$i]}"; new="${PAIRS[$((i+1))]}"
        # word boundaries so Mn_up is not matched inside a longer token
        sed -i "s/\\b${old}\\b/${new}/g" "$scf"
        [ -n "$xml" ] && sed -i "s/\\b${old}\\b/${new}/g" "$xml"
        echo "    $old -> $new"
        i=$((i+2))
    done

    echo -n "    now in scf.in     : "
    sed -n '/ATOMIC_SPECIES/,/^[[:space:]]*[A-Z_]\{4,\}/p' "$scf" \
        | awk 'NF>=3 && $1!="ATOMIC_SPECIES" {printf "%s ", $1}'
    echo
    echo "    backups: *.bak-$ts"
done

echo
echo "Species labels in the charge density and wavefunctions do not exist,"
echo "so no SCF has to be repeated. If a folder had already run its nscf,"
echo "resume with:  ../bin/run.sh 04-06"
echo "If it had not, run:  ../bin/run.sh"
