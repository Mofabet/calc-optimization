#!/usr/bin/env python3
"""
Turn an LMTO symmetry-line path into a KPATH line for the toolkit.

LMTO writes Q vectors in cartesian units of 2*pi/alat, while Wannier90 wants
fractional coordinates of the cell it is given. The two are only the same for a
cubic cell, so the conversion has to go through the actual lattice -- and doing
it explicitly exposes a trap: a path defined for the primitive cell does not
transfer to a supercell unchanged, because the supercell's Brillouin zone is
smaller and points that sat at its boundary land on a reciprocal lattice vector,
i.e. back at Gamma. This script reports that rather than letting it through.

Input is either a CTRL file with a SYML block, or the K_POINTS summary the
lmto_run script produces:

    K_POINTS
    7
      0.000000  0.000000  0.000000   20 !G
      0.000000  0.500000  0.000000   20 !X
      ...

Usage:
    kpath.py CTRL [-s scf.in]
    kpath.py kpoints.txt -s scf.in --units frac
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp


def parse_syml(text):
    """CTRL SYML block -> [(label, (qx, qy, qz))] along a continuous path."""
    m = re.search(r"^SYML(.*?)(?=^[A-Z]{3,}|\Z)", text, re.S | re.M)
    if not m:
        return None
    body = m.group(1)
    segs = []
    for mm in re.finditer(
            r"Q1\s*=\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+LAB1\s*=\s*(\S+)"
            r"\s+Q2\s*=\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+LAB2\s*=\s*(\S+)",
            body, re.S):
        g = mm.groups()
        segs.append(((g[3], tuple(float(x) for x in g[0:3])),
                     (g[7], tuple(float(x) for x in g[4:7]))))
    return segs or None


def parse_summary(text):
    """The lmto_run K_POINTS summary -> [(label, (qx, qy, qz))]."""
    lines = [l for l in text.splitlines() if l.strip()]
    try:
        i = next(i for i, l in enumerate(lines) if l.strip().startswith("K_POINTS"))
    except StopIteration:
        return None
    pts = []
    for l in lines[i + 2:]:
        m = re.match(r"\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+\d+\s*!?\s*(\S*)", l)
        if not m:
            break
        pts.append((m.group(4) or "?",
                    tuple(float(m.group(k)) for k in (1, 2, 3))))
    if len(pts) < 2:
        return None
    return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path_file")
    ap.add_argument("-s", "--scf", default="scf.in")
    ap.add_argument("--units", choices=["cart2pi", "frac"], default="cart2pi",
                    help="units of the input Q vectors (LMTO writes cart2pi)")
    ap.add_argument("--zone", choices=["own", "parent"], default="own",
                    help="own: label the high-symmetry points of the cell in "
                         "scf.in, dividing by its supercell multiplicity. "
                         "parent: keep the literal converted coordinates, "
                         "which trace the parent cell's path through the "
                         "smaller zone")
    args = ap.parse_args()

    text = open(args.path_file, errors="ignore").read()
    segs = parse_syml(text) or parse_summary(text)
    if not segs:
        sys.exit(f"no SYML block or K_POINTS summary found in {args.path_file}")

    if not os.path.exists(args.scf):
        sys.exit(f"{args.scf} not found -- the conversion needs the cell")
    st = qp.read_structure(qp.strip_comments(open(args.scf).read()))
    lat = st["lattice"]
    alat = st["alat_bohr"] * qp.BOHR
    lens = [sum(v * v for v in row) ** 0.5 for row in lat]

    # A fractional coordinate is a fraction of that cell's own b_i, and
    # doubling an axis in real space halves b_i. The same physical k therefore
    # has a fractional coordinate n times larger in an n-fold supercell -- which
    # is why the parent zone boundary comes out at 1.0 and lands on Gamma.
    # Dividing by the multiplicity recovers the number to use for the
    # supercell's own high-symmetry point of the same name.
    mult = qp.supercell_multiplicity(st)

    def conv(q):
        f = list(q) if args.units == "frac" else \
            [sum(q[j] * lat[i][j] / alat for j in range(3)) for i in range(3)]
        if args.zone == "own":
            f = [f[i] / mult[i] for i in range(3)]
        return f

    print(f"cell from {args.scf}: {lens[0]:.4f} x {lens[1]:.4f} x {lens[2]:.4f} A",
          file=sys.stderr)
    print(f"{len(segs)} segments read from {args.path_file}", file=sys.stderr)
    print(f"supercell multiplicity {mult[0]}x{mult[1]}x{mult[2]}"
          + (f", zone = {args.zone}" if max(mult) > 1 else ""), file=sys.stderr)
    if max(mult) > 1 and args.zone == "own":
        print("  components along a multiplied axis are divided by the "
              "multiplicity,", file=sys.stderr)
        print("  so the labels refer to this cell's own zone boundaries.",
              file=sys.stderr)
    print(file=sys.stderr)
    print(f"  {'label':<6s} {'input (2pi/a)':>28s} {'fractional in this cell':>30s}",
          file=sys.stderr)

    seen, folded = [], []
    for (l1, q1), (l2, q2) in segs:
        for lab, q in ((l1, q1), (l2, q2)):
            f = conv(q)
            key = (lab, tuple(round(v, 6) for v in f))
            if key in seen:
                continue
            seen.append(key)
            note = ""
            for v in f:
                if abs(v - round(v)) < 1e-4 and abs(round(v)) >= 1:
                    note = "  <-- equals a reciprocal lattice vector: this is " \
                           "Gamma in this cell"
                    folded.append(lab)
            print(f"  {lab:<6s} {q[0]:9.6f} {q[1]:9.6f} {q[2]:9.6f}"
                  f"   {f[0]:9.6f} {f[1]:9.6f} {f[2]:9.6f}{note}",
                  file=sys.stderr)

    if folded:
        print(file=sys.stderr)
        print(f"WARNING: {', '.join(sorted(set(folded)))} fold onto Gamma.",
              file=sys.stderr)
        print("  The path was defined for a smaller cell than the one in "
              f"{args.scf}.", file=sys.stderr)
        print("  A supercell has a smaller Brillouin zone, so the parent "
              "zone boundary", file=sys.stderr)
        print("  sits outside it. Either halve those components to stay inside "
              "this", file=sys.stderr)
        print("  cell's zone, or keep them and accept that the segment crosses "
              "a zone", file=sys.stderr)
        print("  boundary -- which is fine for a band plot but is not the "
              "high-symmetry", file=sys.stderr)
        print("  point the label claims.", file=sys.stderr)

    parts, prev = [], None
    for (l1, q1), (l2, q2) in segs:
        f1, f2 = conv(q1), conv(q2)
        if prev is None or prev != (l1, tuple(round(v, 6) for v in f1)):
            parts.append((l1, f1))
        parts.append((l2, f2))
        prev = (l2, tuple(round(v, 6) for v in f2))

    kpath = " : ".join(
        f"{lab} {f[0]:.6f} {f[1]:.6f} {f[2]:.6f}" for lab, f in parts)
    print(file=sys.stderr)
    print("Add this line to defaults.conf (or system.conf):", file=sys.stderr)
    print(f"KPATH={kpath}")


if __name__ == "__main__":
    main()
