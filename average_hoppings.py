#!/usr/bin/env python3
"""
Collapse the per-bond hopping list into a readable table.

Usage:
    average_hoppings.py hoppings_MnSi_up.dat --by chem
"""
import argparse
import math
import re
from collections import defaultdict


def chem(label):
    return re.sub(r"\d+$", "", label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--by", choices=["atom", "chem"], default="chem")
    ap.add_argument("--round", type=int, default=3,
                    help="decimals used to consider two distances equal")
    args = ap.parse_args()

    groups = defaultdict(list)
    for line in open(args.file):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 8:
            continue
        A, B = (p[0], p[1]) if args.by == "atom" else (chem(p[0]), chem(p[1]))
        pair = "-".join(sorted([A, B]))
        groups[(pair, round(float(p[5]), args.round))].append(float(p[6]))

    print(f"# source: {args.file}   grouping: by {args.by}")
    print("# pair      distance_Ang    n    avg_t      rms_t      min_t      max_t")
    for (pair, d), v in sorted(groups.items(), key=lambda x: (x[0][0], x[0][1])):
        avg = sum(v) / len(v)
        rms = math.sqrt(sum(t * t for t in v) / len(v))
        print(f"{pair:<10s} {d:12.3f} {len(v):4d} {avg:10.5f} {rms:10.5f} "
              f"{min(v):10.5f} {max(v):10.5f}")


if __name__ == "__main__":
    main()
