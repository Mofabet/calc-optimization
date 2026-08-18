#!/usr/bin/env python3
"""
Summarise LOBSTER's ICOHPLIST by bond type.

ICOHPLIST lists every symmetry-inequivalent contact once per spin channel, so
the raw file has to be aggregated twice over: sum the two spin channels to get
the ICOHP of one bond, then average over crystallographically equivalent bonds
to get a number per contact type.

Sign convention: ICOHP is negative for a bonding interaction, and the
literature quotes -ICOHP, so that a larger number means a stronger bond. That
is what this prints.

Usage:
    icohp.py [ICOHPLIST.lobster] [--dmax 4.0] [--min 0.05]
"""
import argparse
import re
import sys
from collections import defaultdict


def element(label):
    return re.sub(r"[_0-9].*$", "", label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="?", default="ICOHPLIST.lobster")
    ap.add_argument("--dmax", type=float, default=4.0)
    ap.add_argument("--min", type=float, default=0.0,
                    help="drop bond types below this -ICOHP, in eV")
    args = ap.parse_args()

    try:
        lines = open(args.path, errors="ignore").read().splitlines()
    except OSError as e:
        sys.exit(f"cannot read {args.path}: {e}")

    # index label label distance tx ty tz icohp [...]
    pat = re.compile(
        r"^\s*\d+\s+([A-Za-z][A-Za-z_0-9]*)\s+([A-Za-z][A-Za-z_0-9]*)\s+"
        r"([0-9.]+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?[0-9.]+)")

    per_bond = defaultdict(float)     # unique contact -> summed over spins
    bond_type = {}                    # unique contact -> (pair, distance)
    for line in lines:
        m = pat.match(line)
        if not m:
            continue
        l1, l2, d = m.group(1), m.group(2), float(m.group(3))
        tr = (m.group(4), m.group(5), m.group(6))
        val = float(m.group(7))
        if d > args.dmax:
            continue
        key = (l1, l2, tr, round(d, 4))
        per_bond[key] += val
        pair = "-".join(sorted((element(l1), element(l2))))
        bond_type[key] = (pair, round(d, 3))

    if not per_bond:
        sys.exit(f"no data rows recognised in {args.path} -- check the format")

    groups = defaultdict(list)
    for key, v in per_bond.items():
        groups[bond_type[key]].append(v)

    print(f"# {args.path}")
    print(f"# {len(per_bond)} individual contacts within {args.dmax} A, "
          f"summed over spin channels")
    print(f"# -ICOHP in eV; larger means stronger bonding")
    print()
    print(f"  {'bond':<10s} {'d, A':>8s} {'n':>4s} {'-ICOHP':>9s} "
          f"{'min':>9s} {'max':>9s}")
    rows = []
    for (pair, d), v in groups.items():
        mean = -sum(v) / len(v)
        if mean < args.min:
            continue
        rows.append((-mean, pair, d, len(v), mean, -max(v), -min(v)))
    for _, pair, d, n, mean, lo, hi in sorted(rows):
        print(f"  {pair:<10s} {d:8.3f} {n:4d} {mean:9.3f} {lo:9.3f} {hi:9.3f}")

    if not rows:
        print("  (nothing above the --min threshold)")


if __name__ == "__main__":
    main()
