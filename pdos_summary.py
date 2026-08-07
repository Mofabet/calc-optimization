#!/usr/bin/env python3
"""
Summarise projwfc.x output by element and orbital character.

The choice of Wannier projections and of the outer window is a question about
which orbitals carry weight where, and projwfc answers it directly. This reads
every pdos_atm#*(X)_wfc#*(l) file, sums over atoms of the same element and the
same l, and prints the result on a coarse energy grid relative to E_F, per spin
channel.

Read two things off the table: which orbitals have weight at E_F (they must
appear in PROJECTIONS, or the Wannier functions cannot be localised), and where
each manifold ends (that sets DIS_WIN_HI and DIS_WIN_LO).

Usage:
    pdos_summary.py [DIR] [--ef 13.2472] [--bin 1.0] [--lo -12] [--hi 11]
"""
import argparse
import glob
import os
import re
import sys
from collections import defaultdict


def find_ef(folder):
    for name in ("scf.out", "nscf.out"):
        p = os.path.join(folder, name)
        if not os.path.exists(p):
            continue
        hits = re.findall(r"the Fermi energy is\s+([-0-9.]+)\s*ev",
                          open(p, errors="ignore").read())
        if hits:
            return float(hits[-1]), name
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?", default=".")
    ap.add_argument("--ef", type=float)
    ap.add_argument("--bin", type=float, default=1.0)
    ap.add_argument("--lo", type=float, default=-12.0)
    ap.add_argument("--hi", type=float, default=11.0)
    args = ap.parse_args()

    files = glob.glob(os.path.join(args.folder, "*.pdos_atm#*"))
    if not files:
        sys.exit(f"no pdos_atm files in {args.folder} -- run projwfc.x first")

    ef = args.ef
    if ef is None:
        ef, src = find_ef(args.folder)
        if ef is None:
            sys.exit("no Fermi energy found; pass --ef")
        print(f"E_F = {ef:.4f} eV (from {src})")
    else:
        print(f"E_F = {ef:.4f} eV (given)")

    # spin resolution is visible in the header
    head = open(files[0]).readline()
    spin = "ldosup" in head or "dosup" in head
    print(f"spin-resolved: {'yes' if spin else 'no'}   files: {len(files)}")

    nb = int(round((args.hi - args.lo) / args.bin))
    acc = defaultdict(lambda: [[0.0, 0.0] for _ in range(nb)])
    emin_seen, emax_seen = 1e9, -1e9

    for f in files:
        m = re.search(r"atm#\d+\(([A-Za-z_0-9]+)\)_wfc#\d+\(([spdf])\)",
                      os.path.basename(f))
        if not m:
            continue
        el = re.sub(r"[_0-9].*$", "", m.group(1))     # Mn1 / Mn_up -> Mn
        key = f"{el}_{m.group(2)}"
        for line in open(f):
            if line.lstrip().startswith("#"):
                continue
            p = line.split()
            if len(p) < 2:
                continue
            e = float(p[0]) - ef
            emin_seen = min(emin_seen, e)
            emax_seen = max(emax_seen, e)
            i = int((e - args.lo) / args.bin)
            if not (0 <= i < nb):
                continue
            acc[key][i][0] += float(p[1])
            if spin and len(p) > 2:
                acc[key][i][1] += float(p[2])

    print(f"file coverage: E_F{emin_seen:+.1f} .. E_F{emax_seen:+.1f} eV")
    if emax_seen < 3:
        print("\n!! The files stop below E_F+3, so nothing can be said about the")
        print("   unoccupied states. projwfc's Emin/Emax are ABSOLUTE energies,")
        print(f"   not relative to E_F: rerun with Emin={ef - 12:.1f}, "
              f"Emax={ef + 11:.1f}")

    keys = sorted(acc, key=lambda k: (k.split("_")[0], "spdf".index(k[-1])))
    if not keys:
        sys.exit("no recognisable pdos files")

    def table(col, title):
        print(f"\n{title}   (states/eV, summed over atoms of each element)")
        print("  E-E_F  " + "".join(f"{k:>10s}" for k in keys))
        peak = max(max(acc[k][i][col] for k in keys) for i in range(nb)) or 1.0
        for i in range(nb):
            e0 = args.lo + i * args.bin
            vals = [acc[k][i][col] for k in keys]
            bar = "#" * int(round(24 * max(vals) / peak))
            mark = " <- E_F" if e0 <= 0 < e0 + args.bin else ""
            print(f"  {e0:+6.1f} " + "".join(f"{v:10.2f}" for v in vals)
                  + f"  {bar}{mark}")

    table(0, "spin up" if spin else "total")
    if spin:
        table(1, "spin down")

    print("\nweight at E_F (the bin containing it):")
    i_ef = int((0 - args.lo) / args.bin)
    for k in keys:
        u = acc[k][i_ef][0]
        d = acc[k][i_ef][1] if spin else 0.0
        print(f"  {k:<8s} up {u:8.2f}" + (f"   dn {d:8.2f}" if spin else ""))

    print("\nupper edge of each manifold (last bin above 5% of that "
          "orbital's peak):")
    for k in keys:
        for col, lbl in ((0, "up"), (1, "dn")) if spin else ((0, "tot"),):
            series = [acc[k][i][col] for i in range(nb)]
            pk = max(series) or 1.0
            top = max((i for i, v in enumerate(series) if v > 0.05 * pk),
                      default=None)
            bot = min((i for i, v in enumerate(series) if v > 0.05 * pk),
                      default=None)
            if top is None:
                continue
            print(f"  {k:<8s} {lbl}:  E_F{args.lo + bot * args.bin:+.0f} .. "
                  f"E_F{args.lo + (top + 1) * args.bin:+.0f} eV")


if __name__ == "__main__":
    main()
