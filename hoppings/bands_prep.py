#!/usr/bin/env python3
"""
Write a QE 'bands' input on exactly the k-path Wannier90 interpolated along.

Wannier90 saves the path it used to SEED_band.kpt. Feeding that file back to
pw.x as K_POINTS crystal guarantees the two band structures are sampled at the
same points, so they can be compared point by point rather than by eye.

As with nscf.in, the input is produced by editing scf.in rather than rebuilding
it, so non-default settings survive.

Usage:
    bands_prep.py [-c system.conf] [-k SEED_band.kpt] [-o bands.in]
"""
import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp


def set_key(text, nml, key, value):
    m = re.search(rf"(&{nml}\b.*?^\s*/)", text, re.S | re.I | re.M)
    if not m:
        return text
    body = m.group(1)
    if re.search(rf"^\s*{key}\s*=", body, re.I | re.M):
        new = re.sub(rf"^(\s*){key}\s*=[^,\n!]*(,?)\s*$",
                     rf"\g<1>{key} = {value}\g<2>", body, count=1,
                     flags=re.I | re.M)
    else:
        new = re.sub(r"(\n\s*/)$", f",\n    {key} = {value}\\1", body, count=1)
    return text[:m.start(1)] + new + text[m.end(1):]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--conf", default="system.conf")
    ap.add_argument("-k", "--kpt", default=None)
    ap.add_argument("-s", "--scf", default="scf.in")
    ap.add_argument("-o", "--output", default="bands.in")
    args = ap.parse_args()

    kpt = args.kpt
    if kpt is None:
        cands = sorted(glob.glob("*_band.kpt"))
        if not cands:
            sys.exit("no *_band.kpt here -- the Wannier run must have finished "
                     "with bands_plot = true")
        kpt = cands[0]
        if len(cands) > 1:
            print(f"using {kpt} (also found: {', '.join(cands[1:])})")

    lines = [l for l in open(kpt).read().splitlines() if l.strip()]
    nk = int(lines[0].split()[0])
    pts = lines[1:1 + nk]
    if len(pts) != nk:
        sys.exit(f"{kpt} says {nk} points but has {len(pts)}")

    nbnd = None
    if os.path.exists(args.conf):
        for raw in open(args.conf):
            line = raw.split("#", 1)[0].strip()
            if line.startswith("NBND_NSCF="):
                nbnd = line.split("=", 1)[1].strip()

    if not os.path.exists(args.scf):
        sys.exit(f"{args.scf} not found")
    t = open(args.scf).read()
    t = set_key(t, "CONTROL", "calculation", "'bands'")
    t = set_key(t, "CONTROL", "verbosity", "'high'")
    t = set_key(t, "ELECTRONS", "startingpot", "'file'")
    if nbnd:
        t = set_key(t, "SYSTEM", "nbnd", nbnd)

    block = f"K_POINTS crystal\n{nk}\n"
    for p in pts:
        f = p.split()
        w = f[3] if len(f) > 3 else "1.0"
        block += f"{float(f[0]):16.10f} {float(f[1]):16.10f} " \
                 f"{float(f[2]):16.10f}  {w}\n"
    t = re.sub(r"^[ \t]*K_POINTS.*\n(?:(?![ \t]*[A-Z_]{4,}).*\n)*",
               block + "\n", t, count=1, flags=re.M)

    open(args.output, "w").write(t)
    print(f"wrote {args.output}")
    print(f"  path from {kpt}: {nk} k-points")
    print(f"  nbnd {nbnd or 'unchanged from scf.in'}")
    print()
    print("Run it, then compare:")
    print(f"  srun <opts> pw.x -in {args.output} > bands.out")
    print(f"  python3 <bin>/bands_check.py bands.out "
          f"{kpt.replace('_band.kpt', '')}_band.dat")


if __name__ == "__main__":
    main()
