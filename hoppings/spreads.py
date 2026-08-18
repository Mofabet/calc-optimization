#!/usr/bin/env python3
"""
Summarise the final Wannier spreads from a .wout, grouped by element.

The spread of a Wannier function is the mean-square spatial extent of the
orbital it represents. Atom-centred d and p functions in a solid land around
1-3 Ang^2; much larger means the function is smeared over several atoms and a
matrix element between two such objects is not a hopping between those atoms.

Grouping by element answers the question a projection set raises directly:
whether the added rare-earth functions are compact and self-contained (in which
case they take the unwanted character out of the d and p functions, which is
what they were added for) or whether they are diffuse (in which case they are
mixing into everything).

Usage:
    spreads.py SEED.wout [SEED.win]
"""
import math
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp


def parse_wout(path):
    """-> [(centre xyz, spread)] from the last 'Final State' block."""
    text = open(path, errors="ignore").read()
    i = text.rfind("Final State")
    if i < 0:
        return [], None
    out, total = [], None
    for line in text[i:].splitlines():
        m = re.search(r"WF centre and spread\s+\d+\s*\(([^)]*)\)\s+"
                      r"([-0-9.]+)", line)
        if m:
            xyz = [float(v) for v in m.group(1).split(",")]
            out.append((xyz, float(m.group(2))))
        m = re.search(r"Sum of centres and spreads.*?\)\s+([-0-9.]+)", line)
        if m:
            total = float(m.group(1))
    return out, total


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__.strip())
    wout = sys.argv[1]
    win = sys.argv[2] if len(sys.argv) > 2 else \
        re.sub(r"\.wout$", ".win", wout)

    wfs, total = parse_wout(wout)
    if not wfs:
        sys.exit(f"no 'Final State' block in {wout} -- the run has not "
                 f"finished localising")
    if not os.path.exists(win):
        sys.exit(f"{win} not found; pass it as the second argument")

    lattice, atoms, blocks = qp.parse_win(win)
    labels, worst = qp.assign_by_centres([c for c, _ in wfs], atoms,
                                         lattice, blocks)

    n = len(wfs)
    print(f"{wout}")
    print(f"  {n} Wannier functions")
    if total is not None:
        print(f"  total spread {total:.3f} A^2   ->  {total / n:.3f} A^2 "
              f"per function")
    print(f"  worst centre-to-atom distance {worst:.3f} A")

    by_el = defaultdict(list)
    for (c, s), lab in zip(wfs, labels):
        by_el[re.sub(r"\d+$", "", lab)].append(s)

    print()
    print("  element   n     mean    min      max     comment")
    for el in sorted(by_el):
        v = by_el[el]
        mean = sum(v) / len(v)
        if mean < 0.5:
            note = "very compact, essentially atomic"
        elif mean < 3.0:
            note = "atom-centred, as expected for d/p"
        elif mean < 6.0:
            note = "loose"
        else:
            note = "smeared over several atoms"
        print(f"  {el:<8s} {len(v):3d}  {mean:7.3f} {min(v):7.3f} "
              f"{max(v):7.3f}   {note}")

    # Name the off-centre functions. "Worst distance 0.85 A" says something is
    # wrong but not what: a function displaced toward a neighbour is a hybrid
    # the projection set cannot represent, and knowing which element it belongs
    # to is what points at the missing projection.
    cart = [qp.frac_to_cart(lattice, f) for _, f in atoms]
    off = []
    for (c, sp), lab in zip(wfs, labels):
        best = 1e30
        for ac in cart:
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    for k in (-1, 0, 1):
                        sh = qp.frac_to_cart(lattice, [i, j, k])
                        dd = math.dist(c, [ac[t] + sh[t] for t in range(3)])
                        best = min(best, dd)
        if best > 0.3:
            off.append((best, lab, sp))
    if off:
        off.sort(reverse=True)
        print(f"\n  {len(off)} functions sit more than 0.3 A from any atom:")
        for dd, lab, sp in off[:10]:
            print(f"    nearest atom {lab:<6s} at {dd:6.3f} A   spread "
                  f"{sp:7.3f} A^2")
        print("    A displacement comparable to half a bond length means the")
        print("    function is bond-centred: the projection set is missing an")
        print("    orbital that the true state hybridises with.")

    hi = sorted(zip(labels, [s for _, s in wfs]), key=lambda x: -x[1])[:5]
    print("\n  widest functions:")
    for lab, s in hi:
        print(f"    {lab:<6s} {s:8.3f} A^2")

    rms = math.sqrt(sum(s for _, s in wfs) / n)
    print(f"\n  RMS radius of an average function: {rms:.2f} A")
    print("  (compare with the nearest bond length: a function whose radius")
    print("   exceeds the bond length overlaps its neighbours substantially)")


if __name__ == "__main__":
    main()
