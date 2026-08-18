#!/usr/bin/env python3
"""
Compare DFT bands with the Wannier interpolation, point by point.

This is the test that says whether the tight-binding model actually reproduces
the electronic structure it claims to describe. Localisation measures say the
Wannier functions are compact; only this says the Hamiltonian built from them
is right.

What to expect. Inside the frozen window the interpolation should be exact to
a few meV -- those states were constrained. Between the frozen window and the
edge of the outer window, deviations of tens of meV are normal. Outside the
outer window the model says nothing and large deviations are meaningless.

Usage:
    bands_check.py bands.out SEED_band.dat [--win SEED.win] [--ef E]
                    [--plot bands_compare.dat]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp

# A copy of bin/ made at a different time can leave these out of step. Say so
# now rather than failing halfway through with an AttributeError.
for _fn in ("read_structure", "strip_comments"):
    if not hasattr(qp, _fn):
        sys.exit(f"ERROR: {os.path.join(os.path.dirname(qp.__file__), 'qeparse.py')}\n"
                 f"       is missing {_fn}() -- it is older than this script.\n"
                 f"       Copy the whole bin/ directory, not individual files.")


def read_wannier_bands(path):
    """SEED_band.dat -> [[E over k] per band], and the x axis."""
    blocks, cur = [], []
    for line in open(path):
        if not line.strip():
            if cur:
                blocks.append(cur)
                cur = []
            continue
        p = line.split()
        if len(p) >= 2:
            cur.append((float(p[0]), float(p[1])))
    if cur:
        blocks.append(cur)
    if not blocks:
        sys.exit(f"no data in {path}")
    xs = [x for x, _ in blocks[0]]
    bands = [[e for _, e in b] for b in blocks]
    return xs, bands


def read_kblocks(text):
    """
    -> [(spin, (kx, ky, kz), [energies])], one entry per k-point per channel.

    Coordinates and eigenvalues are taken in a single pass so they cannot get
    out of step. Reading them separately looks equivalent and is not: QE prints
    k components in a fixed-width field, so adjacent negative values run
    together as "-0.1250-0.1250" with no separating space, and a whitespace
    regex then silently drops that point while the eigenvalue scan keeps it.
    """
    out = []
    spin, cur_spin, k, e = "unpolarised", None, None, None

    def flush():
        if e is not None and k is not None and e:
            out.append((cur_spin, k, e))

    for line in text.splitlines():
        m = re.search(r"-+\s*SPIN\s+(UP|DOWN)", line, re.I)
        if m:
            spin = m.group(1).lower()
            continue
        if "bands (ev)" in line:
            flush()
            head = line.split("(")[0]
            nums = re.findall(r"-?\d+\.\d+", head)
            k = tuple(float(x) for x in nums[:3]) if len(nums) >= 3 else None
            e, cur_spin = [], spin
            continue
        if e is None:
            continue
        t = line.strip()
        if not t:
            continue
        try:
            e.extend(float(v) for v in t.split())
        except ValueError:          # occupation numbers, a new section, ...
            flush()
            k, e = None, None
    flush()
    return out


def cart2pi_to_crystal(k, lattice, alat_ang):
    """
    Cartesian k in units of 2*pi/alat -> crystal coordinates.

    From b_i . a_j = 2*pi*delta_ij: the component along b_i is the dot product
    of the cartesian k with the real lattice vector a_i, both in units of alat.
    """
    return [sum(k[j] * lattice[i][j] / alat_ang for j in range(3))
            for i in range(3)]


def read_windows(win):
    out = {}
    if win and os.path.exists(win):
        for line in open(win):
            m = re.match(r"\s*(dis_(?:win|froz)_(?:min|max))\s*=\s*([-0-9.eEdD+]+)",
                         line)
            if m:
                out[m.group(1)] = float(m.group(2).replace("d", "e"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("qe_out")
    ap.add_argument("wannier_band")
    ap.add_argument("--win")
    ap.add_argument("--ef", type=float)
    ap.add_argument("--spin", choices=["up", "down", "auto"], default="auto")
    ap.add_argument("--plot", default="bands_compare.dat")
    ap.add_argument("--kpt", help="default: SEED_band.kpt beside the .dat")
    ap.add_argument("--scf", default="scf.in")
    ap.add_argument("--ktol", type=float, default=1e-4)
    args = ap.parse_args()

    win = args.win or re.sub(r"_band\.dat$", ".win", args.wannier_band)
    w = read_windows(win)

    text = open(args.qe_out, errors="ignore").read()
    ef = args.ef
    if ef is None:
        m = re.findall(r"the Fermi energy is\s+([-0-9.]+)\s*ev", text)
        if m:
            ef = float(m[-1])
    if ef is None and w.get("dis_win_min") is not None:
        print("no Fermi energy in the bands output (a 'bands' run does not "
              "compute one); pass --ef or it will be taken as 0", file=sys.stderr)
    ef = ef or 0.0

    kblocks = read_kblocks(text)
    if not kblocks:
        sys.exit(f"no eigenvalue blocks in {args.qe_out} -- needs "
                 f"verbosity = 'high'")

    spin = args.spin
    if spin == "auto":
        spin = "up" if "_up" in os.path.basename(args.wannier_band) else \
               ("down" if "_dn" in os.path.basename(args.wannier_band)
                else "unpolarised")
    chan = [(k, e) for s_, k, e in kblocks if s_ == spin] or \
           [(k, e) for _, k, e in kblocks]
    sel = [e for _, e in chan]
    qk_cart = [k for k, _ in chan]
    nb = {len(e) for e in sel}
    print(f"spin channel: {spin}   k-points in {args.qe_out}: {len(sel)}"
          f"   bands per point: {'/'.join(str(n) for n in sorted(nb))}")
    if len(nb) > 1:
        print("  NOTE: the eigenvalue blocks are not all the same length; "
              "the output may be truncated")

    xs, wb = read_wannier_bands(args.wannier_band)
    print(f"wannier bands: {len(wb)} bands over {len(xs)} k-points")

    # Match by coordinates rather than by position in the file. The two runs
    # are generated by different programs and need not produce identical point
    # counts -- QE expands a crystal_b path with its own rounding, and a
    # segment junction may or may not be duplicated. Comparing index i against
    # index i is only valid if the lists happen to coincide exactly, and
    # silently wrong if they do not.
    kpt = args.kpt or re.sub(r"_band\.dat$", "_band.kpt", args.wannier_band)
    pairs = None
    if os.path.exists(kpt) and os.path.exists(args.scf):
        klines = [l for l in open(kpt).read().splitlines() if l.strip()][1:]
        wk = [[float(v) for v in l.split()[:3]] for l in klines]
        st = qp.read_structure(qp.strip_comments(open(args.scf).read()))
        alat_ang = st["alat_bohr"] * qp.BOHR
        qk = [cart2pi_to_crystal(k, st["lattice"], alat_ang)
              for k in qk_cart]

        def wrapdiff(a, b):
            d = a - b
            return abs(d - round(d))

        pairs, unmatched = [], 0
        for iq, kq in enumerate(qk):
            best, bd = None, 1e9
            for iw, kw in enumerate(wk):
                d = max(wrapdiff(kq[t], kw[t]) for t in range(3))
                if d < bd:
                    bd, best = d, iw
            if bd <= args.ktol:
                pairs.append((iq, best))
            else:
                unmatched += 1
        print(f"matched {len(pairs)} of {len(qk)} DFT k-points against "
              f"{len(wk)} in {os.path.basename(kpt)}"
              + (f"; {unmatched} unmatched" if unmatched else ""))
        if not pairs:
            sys.exit("no k-point matched. The two runs are on different paths; "
                     "rebuild the bands input with bands_prep.py, which "
                     "takes the path from *_band.kpt.")

    if pairs is None:
        if len(sel) != len(xs):
            sys.exit(f"k-point count differs ({len(sel)} vs {len(xs)}) and "
                     f"coordinates could not be read.\n"
                     f"Need {os.path.basename(kpt)} and {args.scf} beside the "
                     f"data, or rebuild with bands_prep.py.")
        pairs = [(i, i) for i in range(len(xs))]

    froz = (w.get("dis_froz_min"), w.get("dis_froz_max"))
    outer = (w.get("dis_win_min"), w.get("dis_win_max"))
    print(f"frozen window {froz[0]} .. {froz[1]} eV")
    print(f"outer  window {outer[0]} .. {outer[1]} eV")

    # For each Wannier band find the nearest DFT band, not the reverse.
    # The model spans a subspace: inside the outer window the DFT bands
    # outnumber the Wannier ones, and the surplus is not represented by
    # construction. Measuring DFT-to-Wannier therefore charges the model for
    # states it was never asked to reproduce and inflates every number.
    def collect(lo, hi):
        devs = []
        for iq, iw in pairs:
            qk = sorted(sel[iq])
            for b in wb:
                e = b[iw]
                if lo is not None and not (lo <= e <= hi):
                    continue
                devs.append(min(abs(e - x) for x in qk))
        return sorted(devs)

    def show(lbl, lo, hi):
        d = collect(lo, hi)
        if not d:
            return
        n = len(d)
        print(f"  {lbl:<26s} {n:6d} {d[n // 2] * 1000:9.2f} "
              f"{sum(d) / n * 1000:9.2f} {d[int(0.95 * n)] * 1000:9.2f} "
              f"{d[-1] * 1000:10.2f}")

    fl, fh = froz
    print()
    print(f"  deviation of each Wannier band from the nearest DFT band, meV")
    print(f"  {'region':<26s} {'n':>6s} {'median':>9s} {'mean':>9s} "
          f"{'95%':>9s} {'max':>10s}")
    if fl is not None:
        show("frozen window", fl, fh)
        show("below frozen", outer[0], fl)
        show("above frozen", fh, outer[1])
    show("E_F -3 .. E_F +3", ef - 3, ef + 3)

    alld = collect(None, None)
    if alld:
        n = len(alld)
        u1 = sum(1 for x in alld if x < 0.001) / n
        u10 = sum(1 for x in alld if x < 0.010) / n
        print(f"\n  {u1:.0%} of interpolated states agree with a DFT band to "
              f"better than 1 meV, {u10:.0%} to better than 10 meV")

    print()
    print("  Read the median, not the mean: a handful of band crossings, where")
    print("  the nearest-band rule picks the wrong partner, dominate the mean")
    print("  while saying nothing about the fit.")
    print("  Inside the frozen window those states were constrained exactly, so")
    print("  the median belongs at or below 1 meV. Above it the disentangled")
    print("  subspace is a smooth mixture rather than a selection of DFT")
    print("  states, and tens of meV are expected there, not a defect.")

    if args.plot:
        with open(args.plot, "w") as f:
            f.write("# x  E-E_F  source(0=DFT 1=Wannier)\n")
            for iq, iw in pairs:
                for e in sel[iq]:
                    f.write(f"{xs[iw]:12.6f} {e - ef:12.6f} 0\n")
            f.write("\n")
            for b in wb:
                for x, e in zip(xs, b):
                    f.write(f"{x:12.6f} {e - ef:12.6f} 1\n")
                f.write("\n")
        print(f"\nplot data: {args.plot}")
        print("  gnuplot> plot '<file>' u 1:($3==0?$2:1/0) w p pt 7 ps 0.3 "
              "t 'DFT', '' u 1:($3==1?$2:1/0) w l t 'Wannier'")


if __name__ == "__main__":
    main()
