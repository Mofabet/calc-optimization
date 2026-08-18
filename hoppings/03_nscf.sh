#!/bin/bash
# 03 -- non-self-consistent run on the uniform mesh.
#
# Reuses the charge density in out/. nosym/noinv are already in nscf.in:
# Wannier90 needs the full mesh, not the irreducible wedge.
#
#   03_nscf.sh              run it
#   03_nscf.sh --check      report status only, run nothing

set -uo pipefail
BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$BIN/common.sh"

CHECK=0
[ "${1:-}" = "--check" ] && CHECK=1

conf_load
[ -f nscf.in ] || die "no nscf.in -- run 02_inputs.py first"
[ -d out ] || die "no out/ -- the nscf needs the density from the SCF"

if grep -q "End of band structure calculation" nscf.out 2>/dev/null; then
    say "already finished: $(grep 'the Fermi energy is' nscf.out | tail -1)"
    exit 0
fi
[ $CHECK -eq 1 ] && { say "not yet run"; exit 0; }

PW=$(resolve "${QE_BIN:+$QE_BIN/}pw.x") || die "pw.x not found; set QE_BIN"
say "srun $SRUN_OPTS $PW -in nscf.in"
# shellcheck disable=SC2086
srun $SRUN_OPTS "$PW" -in nscf.in > nscf.out
grep -q "End of band structure calculation" nscf.out \
    || die "nscf did not finish; see nscf.out"
say "$(grep 'the Fermi energy is' nscf.out | tail -1)"
