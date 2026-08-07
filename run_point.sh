#!/bin/bash
# One z-point, from a finished SCF to averaged hoppings.
#
# Expects to be run inside a folder that already contains:
#   scf.in   scf.out   out/          (the whole out/, not just the text files)
#
# Usage:
#   ../bin/run_point.sh              # everything
#   ../bin/run_point.sh nscf         # stop after the nscf
#   ../bin/run_point.sh wannier      # from the nscf output onwards
#   ../bin/run_point.sh analyse      # only re-extract hoppings from _hr.dat
#   ../bin/run_point.sh conf         # rebuild system.conf and the inputs only

set -euo pipefail

BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="${1:-all}"
CONF="system.conf"

# Read one key. Deliberately not "| xargs": with no command xargs runs echo,
# and a value starting with -n or -e is then swallowed as an echo flag, which
# silently turns "SRUN_OPTS=-n 16" into "16".
get() {
    local v
    v=$(grep -E "^$1=" "$CONF" | head -1 | cut -d= -f2- | sed 's/[[:space:]]*#.*$//')
    v="${v#"${v%%[![:space:]]*}"}"
    v="${v%"${v##*[![:space:]]}"}"
    printf '%s' "$v"
}

# ---------------------------------------------------------------- checks --
[ -f scf.in ]  || { echo "ERROR: no scf.in here"; exit 1; }
[ -f scf.out ] || { echo "ERROR: no scf.out here"; exit 1; }
grep -q "convergence has been achieved" scf.out || {
    echo "ERROR: scf.out does not report convergence. Fix the SCF first."; exit 1; }
[ -d out ] || { echo "ERROR: no out/ directory. The nscf needs the charge"
                echo "       density from the SCF -- copy the whole out/."; exit 1; }

# ---------------------------------------------------------------- config --
if [ ! -f "$CONF" ] || [ "$STAGE" = "conf" ]; then
    echo "=== deriving $CONF from scf.in ==="
    python3 "$BIN/scf2conf.py" -o "$CONF"
fi

SYSTEM=$(get SYSTEM)
PAIRS=$(get PAIRS)
DMAX=$(get DMAX)
TMIN=$(get TMIN)
QE_BIN=$(get QE_BIN)
SRUN_OPTS=$(get SRUN_OPTS)
PW="${QE_BIN:+$QE_BIN/}pw.x"
if ! command -v "$PW" >/dev/null 2>&1 && [ ! -x "$PW" ]; then
    echo "ERROR: pw.x not found at '$PW'. Set QE_BIN in defaults.conf."
    exit 1
fi

echo "=== generating inputs ==="
python3 "$BIN/make_inputs.py" "$CONF"

if [ "$STAGE" = "conf" ]; then
    echo
    echo "Config and inputs rebuilt. Nothing was run."
    exit 0
fi

# ------------------------------------------------------------------ nscf --
if [ "$STAGE" = "all" ] || [ "$STAGE" = "nscf" ]; then
    echo "=== nscf ==="
    echo "    srun $SRUN_OPTS $PW -in nscf.in"
    # shellcheck disable=SC2086
    srun $SRUN_OPTS "$PW" -in nscf.in > nscf.out
    grep -q "End of band structure calculation" nscf.out || {
        echo "FAILED: nscf did not finish, see nscf.out"; exit 1; }
    grep "the Fermi energy is" nscf.out | tail -1
fi
[ "$STAGE" = "nscf" ] && exit 0

# --------------------------------------------------------------- wannier --
if [ "$STAGE" = "all" ] || [ "$STAGE" = "wannier" ]; then
    # A first, rough set of windows: -pp ignores dis_* entirely, so this only
    # has to be syntactically present for the preprocessing step.
    echo "=== windows (provisional, from E_F) ==="
    python3 "$BIN/set_windows.py" "$CONF"

    # Overlaps first for both spins, then the windows, then the minimisation.
    # The .mmn/.amn/.eig files depend on the projections and the k-mesh, not on
    # the windows, so they are computed once; validating the windows needs the
    # .eig of both channels, and minimising before that validation is what made
    # the first attempt fail on the frozen window.
    echo "=== overlaps (pw2wannier90), both spins ==="
    "$BIN/run_wannier.sh" "$CONF" both overlaps

    for s in up dn; do
        if [ ! -s "${SYSTEM}_${s}.eig" ]; then
            echo "STOPPING: ${SYSTEM}_${s}.eig was not produced."
            exit 1
        fi
    done

    echo "=== windows, validated against both .eig ==="
    python3 "$BIN/set_windows.py" "$CONF"

    echo "=== minimisation, both spins ==="
    "$BIN/run_wannier.sh" "$CONF" both minimise
fi

# --------------------------------------------------------------- analyse --
echo "=== hoppings ==="
for s in up dn; do
    python3 "$BIN/extract_hoppings.py" \
        "${SYSTEM}_${s}_hr.dat" "${SYSTEM}_${s}.win" \
        --pairs "$PAIRS" --dmax "$DMAX" --tmin "$TMIN" --onsite \
        > "hoppings_${s}.dat" 2> "assign_${s}.log"
    python3 "$BIN/average_hoppings.py" "hoppings_${s}.dat" --by chem > "avg_chem_${s}.dat"
    python3 "$BIN/average_hoppings.py" "hoppings_${s}.dat" --by atom > "avg_atom_${s}.dat"
done

echo
echo "--- assignment warnings (empty is good) ---"
grep -h "WARNING\|do NOT sit" assign_up.log assign_dn.log || echo "  none"
echo
echo "--- result ---"
python3 "$BIN/make_table.py" avg_chem_up.dat avg_chem_dn.dat \
        --label "$(basename "$PWD")" --nearest
