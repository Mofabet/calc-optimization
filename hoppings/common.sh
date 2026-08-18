#!/bin/bash
# Shared shell helpers. Sourced by the numbered steps; not run directly.

CONF="${CONF:-system.conf}"

say() { printf '  %s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# Read one key. Deliberately not "| xargs": with no command xargs runs echo,
# and a value starting with -n or -e is swallowed as an echo flag, which
# silently turns "SRUN_OPTS=-n 16" into "16".
get() {
    sed -n "s/^$1=//p" "$CONF" 2>/dev/null | head -1 \
        | sed 's/#.*//; s/^[[:space:]]*//; s/[[:space:]]*$//'
}

conf_load() {
    [ -f "$CONF" ] || die "no $CONF here -- run 01_conf.py first"
    SYSTEM=$(get SYSTEM);      PREFIX=$(get PREFIX)
    QE_BIN=$(get QE_BIN);      W90_BIN=$(get W90_BIN)
    SRUN_OPTS=$(get SRUN_OPTS); W90_SRUN_OPTS=$(get W90_SRUN_OPTS)
    OMP_THREADS=$(get OMP_THREADS)
    PAIRS=$(get PAIRS); DMAX=$(get DMAX); TMIN=$(get TMIN)
}

# Resolve to an absolute path before handing anything to srun: a bare name
# relies on PATH being identical on the compute node, which it often is not,
# and the job then dies with "No such file or directory" there while the check
# passed here.
resolve() {
    if [ -x "$1" ]; then readlink -f "$1"; return 0; fi
    local p
    p=$(command -v "$1" 2>/dev/null) && readlink -f "$p" && return 0
    return 1
}

# QE writes fatal errors to CRASH and to stderr, not to stdout, so a stdout
# file can look healthy while the run has aborted.
show_error() {
    local out="$1" err="${2:-}" f
    for f in CRASH "$err" "$out"; do
        [ -n "$f" ] && [ -s "$f" ] || continue
        if grep -qi "%%%%\|Error in routine\|error #\|ERROR AT K-POINT" "$f"; then
            printf -- '--- %s ---\n' "$f"
            grep -m1 -B3 -A10 -i \
                 "%%%%\|Error in routine\|error #\|ERROR AT K-POINT" "$f"
            return 0
        fi
    done
    if grep -qi "execve()\|No such file or directory" "$out" 2>/dev/null; then
        printf -- '--- %s ---\n' "$out"
        grep -m2 -i "execve()\|No such file or directory" "$out"
        printf '    Found here but not on the compute node. A binary under\n'
        printf '    /usr/bin is installed per node; put it on shared storage\n'
        printf '    and set W90_BIN, or leave W90_SRUN_OPTS empty to run here.\n'
        return 0
    fi
    printf -- '--- tail of %s ---\n' "$out"
    tail -15 "$out" 2>/dev/null
}

need_file() {
    [ -s "$1" ] && return 0
    printf 'FAILED: %s was not produced.\n' "$1"
    show_error "$2" "${3:-}"
    exit 1
}
