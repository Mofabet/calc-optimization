#!/usr/bin/env python3
"""
Set the Wannier90 disentanglement windows from the actual Fermi level, and
check that they are feasible before you burn a run on them.

This is the step that was hand-tuned for GdFeSi (dis_win_min = 8.0 etc. were
written for E_F = 13.78 eV). For Mn and Ru E_F will be somewhere else, so
hardcoding those numbers would silently produce a bad Wannierisation.

Two feasibility conditions, both checked per k-point:
  * bands inside the outer  window  >= num_wann   (otherwise nothing to pick)
  * bands inside the frozen window  <= num_wann   (otherwise W90 aborts)
If either fails the window is walked until it fits, and the change is printed.

Usage:
    set_windows.py CONF                      # after nscf, before the main run
    set_windows.py CONF --ef 13.78           # override the detected E_F
    set_windows.py CONF --dry-run            # report only, do not touch .win
"""
import argparse
import glob
import os
import re
import sys


def read_conf(path):
    conf = {}
    for raw in open(path):
        line = raw.split("#", 1)[0].strip()
        if line and "=" in line:
            k, v = line.split("=", 1)
            conf[k.strip()] = v.strip()
    return conf


def find_fermi(*candidates):
    for path in candidates:
        if not os.path.exists(path):
            continue
        ef = None
        for line in open(path, errors="ignore"):
            m = re.search(r"the Fermi energy is\s+([-0-9.]+)\s*ev", line)
            if m:
                ef = float(m.group(1))
            m = re.search(r"highest occupied level \(ev\):\s+([-0-9.]+)", line)
            if m and ef is None:
                ef = float(m.group(1))
        if ef is not None:
            return ef, path
    return None, None


def read_eig(path):
    """seedname.eig -> {ik: [energies]}"""
    per_k = {}
    for line in open(path):
        p = line.split()
        if len(p) < 3:
            continue
        ib, ik, e = int(p[0]), int(p[1]), float(p[2])
        per_k.setdefault(ik, []).append(e)
    return per_k


def counts(per_k, lo, hi):
    return [sum(1 for e in ens if lo <= e <= hi) for ens in per_k.values()]


def patch_win(path, vals, dry):
    text = open(path).read()
    for key, v in vals.items():
        new = f"{key:<12s} = {v:.4f}"
        text, n = re.subn(rf"^{key}\s*=.*$", new, text, count=1, flags=re.M)
        if n == 0:
            sys.exit(f"ERROR: no '{key}' line in {path}")
    text = text.replace(
        "! WINDOWS ARE PLACEHOLDERS -- run bin/set_windows.py after nscf",
        "! windows set by bin/set_windows.py")
    if not dry:
        open(path, "w").write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("conf")
    ap.add_argument("--ef", type=float)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    c = read_conf(args.conf)
    sysname = c["SYSTEM"]

    ef = args.ef
    if ef is None:
        # Prefer the SCF Fermi level. The NSCF mesh is often deliberately
        # coarser than the SCF one -- it only has to resolve the Wannier
        # Fourier transform, not integrate the density -- and a Fermi level
        # from a coarse mesh is the less accurate of the two.
        ef, src = find_fermi("scf.out", "nscf.out")
        if ef is None:
            sys.exit("ERROR: no Fermi energy in scf.out/nscf.out. "
                     "Run the scf first, or pass --ef.")
        print(f"E_F = {ef:.4f} eV   (from {src})")
        other, osrc = find_fermi("nscf.out" if src.startswith("scf")
                                 else "scf.out")
        if other is not None and abs(other - ef) > 0.05:
            print(f"      {osrc} gives {other:.4f} eV, a difference of "
                  f"{abs(other - ef):.3f} eV -- expected when the two meshes "
                  f"differ; override with --ef if you prefer the other.")
    else:
        print(f"E_F = {ef:.4f} eV   (given)")

    win_lo = ef - float(c["DIS_WIN_LO"])
    win_hi = ef + float(c["DIS_WIN_HI"])
    froz_lo = ef - float(c["DIS_FROZ_LO"])
    froz_hi = ef + float(c["DIS_FROZ_HI"])

    for spin in ("up", "dn"):
        seed = f"{sysname}_{spin}"
        win = f"{seed}.win"
        if not os.path.exists(win):
            print(f"  {win}: not found, skipped")
            continue

        m = re.search(r"^num_wann\s*=\s*(\d+)", open(win).read(), re.M)
        num_wann = int(m.group(1))

        lo, hi, flo, fhi = win_lo, win_hi, froz_lo, froz_hi
        eig = f"{seed}.eig"
        notes = []

        if os.path.exists(eig):
            per_k = read_eig(eig)
            nk = len(per_k)

            # outer window must hold at least num_wann bands everywhere
            guard = 0
            while min(counts(per_k, lo, hi)) < num_wann and guard < 200:
                hi += 0.1
                lo -= 0.05
                guard += 1
            if guard:
                notes.append(f"outer widened by {guard * 0.1:.1f} eV upward "
                             f"and {guard * 0.05:.2f} eV downward to fit "
                             f"num_wann")
                if guard * 0.1 > 0.5:
                    notes.append("!! that is a large widening. It usually "
                                 "means num_wann is too big for the clean "
                                 "part of the spectrum, and the window has "
                                 "been pushed into a neighbouring manifold "
                                 "(rare-earth 4f is the usual culprit). "
                                 "Removing a projection is normally better "
                                 "than swallowing those bands -- check "
                                 "count_bands.py --histogram")

            # frozen window must hold at most num_wann bands everywhere
            guard = 0
            while max(counts(per_k, flo, fhi)) > num_wann and guard < 400:
                fhi -= 0.05
                guard += 1
            if guard:
                notes.append(f"frozen top lowered by {guard * 0.05:.2f} eV")
            if fhi < ef:
                notes.append("!! frozen window no longer reaches E_F -- "
                             "num_wann is too small for this energy range, "
                             "add another projection (e.g. Gd:d)")

            cnt_out = counts(per_k, lo, hi)
            cnt_frz = counts(per_k, flo, fhi)
            stat = (f"    bands in outer  window: {min(cnt_out)}..{max(cnt_out)}"
                    f"  (need >= {num_wann})\n"
                    f"    bands in frozen window: {min(cnt_frz)}..{max(cnt_frz)}"
                    f"  (need <= {num_wann})\n"
                    f"    checked over {nk} k-points")
        else:
            stat = (f"    {eig} not found -- windows set from E_F only.\n"
                    f"    Re-run this script after pw2wannier90.x to validate them.")

        patch_win(win, {"dis_win_min": lo, "dis_win_max": hi,
                        "dis_froz_min": flo, "dis_froz_max": fhi},
                  args.dry_run)

        print(f"\n  {win}   num_wann = {num_wann}")
        print(f"    outer  [{lo:8.3f} , {hi:8.3f}]")
        print(f"    frozen [{flo:8.3f} , {fhi:8.3f}]")
        print(stat)
        for n in notes:
            print(f"    NOTE: {n}")

    if args.dry_run:
        print("\n(dry run -- no file was modified)")


if __name__ == "__main__":
    main()
