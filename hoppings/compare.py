#!/usr/bin/env python3
"""
Compare t_eff across two or more runs of the same structure.

Written to settle convergence questions with the observable rather than with
the optimiser's own diagnostics. Omega_I and the per-iteration delta describe
how the minimisation is progressing; they are not the quantity that goes into
a paper. A setting that moves Omega_I by 20% while leaving t_eff within a
percent was already sufficient, and the run time saved across a series is
worth more than the extra digits.

For a series the criterion is stricter in one respect and looser in another:
every structure must use identical settings, but the residual error need only
be small compared with the variation across the series, since a systematic
offset common to all points leaves the trend intact.

Usage:
    compare.py avg_chem_up_1000.dat avg_chem_up_5000.dat
    compare.py it500/avg_chem_up.dat it1000/avg_chem_up.dat \
        it2000/avg_chem_up.dat --labels 500 1000 2000
    compare.py A.dat B.dat --tol 1
"""
import argparse
import os
import sys


def load(path):
    """-> {(pair, distance): rms_t}"""
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) < 5:
            continue
        out[(p[0], round(float(p[1]), 3))] = float(p[4])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--labels", nargs="*")
    ap.add_argument("--tol", type=float, default=1.0,
                    help="relative spread in %% considered negligible")
    args = ap.parse_args()

    if len(args.files) < 2:
        sys.exit("need at least two files")

    data = [load(f) for f in args.files]
    for f, d in zip(args.files, data):
        if not d:
            sys.exit(f"{f} has no data rows")

    labels = args.labels or []
    if len(labels) != len(args.files):
        labels = []
        for f in args.files:
            parent = os.path.basename(os.path.dirname(os.path.abspath(f)))
            stem = os.path.splitext(os.path.basename(f))[0]
            labels.append(parent if parent and parent != "." else stem)
    w = max(8, max(len(l) for l in labels) + 1)

    keys = sorted(set.intersection(*(set(d) for d in data)),
                  key=lambda k: (k[0], k[1]))
    if not keys:
        sys.exit("no bond appears in every file")

    for f, l in zip(args.files, labels):
        print(f"  {l:<{w}s} {f}")
    print()

    head = f"  {'pair':<10s} {'d, A':>7s}" + \
           "".join(f"{l:>{w}s}" for l in labels) + f"{'spread':>9s}"
    print(head)

    worst, worst_key = 0.0, None
    rows = []
    for k in keys:
        vals = [d[k] for d in data]
        lo, hi = min(vals), max(vals)
        rel = 100 * (hi - lo) / lo if lo else float("nan")
        if rel > worst:
            worst, worst_key = rel, k
        rows.append((k, vals, rel))

    for k, vals, rel in rows:
        flag = "" if rel <= args.tol else "  <--"
        print(f"  {k[0]:<10s} {k[1]:7.3f}"
              + "".join(f"{v:>{w}.4f}" for v in vals)
              + f"{rel:8.2f}%{flag}")

    # The nearest bond of each pair is the one the physics rests on; a distant
    # bond with t_eff near the cutoff can swing by a lot in relative terms
    # while being irrelevant in absolute terms.
    print()
    nearest = {}
    for k, vals, rel in rows:
        if k[0] not in nearest or k[1] < nearest[k[0]][0][1]:
            nearest[k[0]] = (k, vals, rel)
    print("  nearest bond of each pair:")
    for pair, (k, vals, rel) in sorted(nearest.items()):
        span = max(vals) - min(vals)
        print(f"    {pair:<10s} {k[1]:.3f} A   spread {span:.4f} eV "
              f"({rel:.2f}%)")

    print()
    if worst <= args.tol:
        print(f"Every bond agrees to within {args.tol:g}% across all runs. "
              f"The cheapest of these settings is sufficient.")
    else:
        print(f"Largest spread {worst:.2f}% on {worst_key[0]} at "
              f"{worst_key[1]:.3f} A, above the {args.tol:g}% threshold.")
        print("Check whether it is a near bond or a distant one: a swing in a "
              "small t_eff near the cutoff rarely changes any conclusion, "
              "while the same swing on the nearest bond does.")

    missing = set()
    for d in data:
        missing |= set(d) - set(keys)
    if missing:
        print()
        print(f"{len(missing)} bond(s) absent from at least one file "
              f"(usually one side of --tmin): "
              + ", ".join(f"{p} {dd:.3f}" for p, dd in sorted(missing)[:6]))


if __name__ == "__main__":
    main()
