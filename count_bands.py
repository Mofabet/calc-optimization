#!/usr/bin/env python3
"""
Count how many bands lie below E_F + delta, to size NBND_NSCF.

The .mmn overlap file scales as num_bands^2, so carrying an SCF band count
sized for a density of states into the Wannier stage is expensive. Only bands
reaching a few eV above E_F are needed. This reports the largest count over
all k-points, which is the number that has to fit.

Parsing is line-based rather than regex-based: QE prints eigenvalues eight to
a line, splits them by spin channel, and with verbosity='high' follows each
block with occupation numbers that must not be counted.

Usage:
    count_bands.py [scf.out] [--ef VALUE] [--margin 10]
"""
import argparse
import re
import sys


def numeric_line(line):
    s = line.strip()
    if not s:
        return None
    toks = s.split()
    try:
        return [float(t) for t in toks]
    except ValueError:
        return None


def parse(path):
    ef, nelec, nbnd_scf = None, None, None
    blocks, spin = [], "unpolarised"
    cur, collecting = None, False

    for line in open(path, errors="ignore"):
        m = re.search(r"the Fermi energy is\s+([-0-9.]+)\s*ev", line)
        if m:
            ef = float(m.group(1))
        m = re.search(r"highest occupied.*?:\s+([-0-9.]+)", line)
        if m and ef is None:
            ef = float(m.group(1))
        m = re.search(r"number of electrons\s*=\s*([0-9.]+)", line)
        if m:
            nelec = float(m.group(1))
        m = re.search(r"number of Kohn-Sham states\s*=\s*(\d+)", line)
        if m:
            nbnd_scf = int(m.group(1))
        m = re.search(r"-+\s*SPIN\s+(UP|DOWN)", line, re.I)
        if m:
            spin = m.group(1).lower()

        if "bands (ev)" in line:
            if cur:
                blocks.append((spin, cur))
            cur, collecting = [], True
            continue

        if collecting:
            vals = numeric_line(line)
            if vals is None:
                if line.strip():          # 'occupation numbers', new k, etc.
                    blocks.append((spin, cur))
                    cur, collecting = None, False
                continue                  # blank lines are allowed inside
            cur.extend(vals)

    if cur:
        blocks.append((spin, cur))
    return ef, nelec, nbnd_scf, blocks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log", nargs="?", default="scf.out")
    ap.add_argument("--ef", type=float)
    ap.add_argument("--margin", type=int, default=10,
                    help="spare bands added on top (default 10)")
    ap.add_argument("--histogram", action="store_true",
                    help="print the eigenvalue density around E_F, which "
                         "reveals flat manifolds such as rare-earth 4f")
    args = ap.parse_args()

    ef, nelec, nbnd_scf, blocks = parse(args.log)
    if args.ef is not None:
        ef = args.ef

    if not blocks:
        print(f"No eigenvalue blocks found in {args.log}.")
        print("The SCF must have been run with verbosity = 'high'; otherwise")
        print("eigenvalues are not written. Either rerun the last iteration")
        print("with high verbosity, or skip this check -- set_windows.py")
        print("validates the band count against the .eig file after the NSCF.")
        sys.exit(1)
    if ef is None:
        sys.exit(f"No Fermi energy in {args.log}. Pass --ef.")

    sizes = {len(b) for _, b in blocks}
    spins = sorted({s for s, _ in blocks})
    print(f"{args.log}")
    print(f"  E_F                {ef:.4f} eV")
    if nelec:
        print(f"  electrons          {nelec:.0f}")
    if nbnd_scf:
        print(f"  bands in the scf   {nbnd_scf}")
    print(f"  eigenvalue blocks  {len(blocks)}  "
          f"(spin channels: {', '.join(spins)})")
    print(f"  bands per block    {'/'.join(str(s) for s in sorted(sizes))}")
    if len(sizes) > 1:
        print("  NOTE: blocks differ in length; the file may be truncated.")

    print()
    print("  delta   bands below E_F+delta      suggested NBND_NSCF")
    best = None
    for d in (3, 5, 8, 10):
        counts = [sum(1 for e in b if e < ef + d) for _, b in blocks]
        hi = max(counts)
        sug = hi + args.margin
        print(f"  {d:>2} eV   min {min(counts):4d}   max {hi:4d}"
              f"            {sug:4d}")
        if d == 5:
            best = sug

    print()
    print(f"  With the default window (E_F+5 eV) set NBND_NSCF = {best}.")
    if nbnd_scf:
        ratio = (best / nbnd_scf) ** 2
        print(f"  That makes the .mmn about {1 - ratio:.0%} smaller than "
              f"carrying all {nbnd_scf} bands.")
    print("  Raise it if set_windows.py reports too few bands in the outer "
          "window.")

    if args.histogram:
        histogram(blocks, ef, spins)


def histogram(blocks, ef, spins, lo=-12, hi=12, step=1.0):
    """
    Bands per k-point per eV, relative to E_F.

    A flat manifold -- rare-earth 4f is the usual one -- shows up as a spike.
    Flat bands inside the outer window are what degrade a disentanglement
    based on d and p projections, so this is the picture to consult when
    choosing dis_win_max rather than reusing a range chosen for a DOS plot.
    """
    print()
    print(f"  Eigenvalue density, bands per k-point per {step:g} eV bin")
    print(f"  (E - E_F in eV; a spike is a flat, weakly dispersing manifold)")
    nbins = int((hi - lo) / step)
    for spin in spins:
        sel = [b for s, b in blocks if s == spin]
        if not sel:
            continue
        bins = [0] * nbins
        for b in sel:
            for e in b:
                i = int((e - ef - lo) / step)
                if 0 <= i < nbins:
                    bins[i] += 1
        dens = [c / len(sel) for c in bins]
        peak = max(dens) or 1.0
        print(f"\n  spin {spin}")
        for i, d in enumerate(dens):
            e0 = lo + i * step
            bar = "#" * int(round(40 * d / peak))
            mark = "  <- E_F" if e0 <= 0 < e0 + step else ""
            print(f"   {e0:+6.0f}  {d:6.1f}  {bar}{mark}")


if __name__ == "__main__":
    main()
