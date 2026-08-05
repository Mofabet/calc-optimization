#!/usr/bin/env python3
"""
Read an LMTO CTRL file and emit the structural part of a system .conf.

This replaces the step you used to do by hand: taking STRUC ALAT / PLAT and
the SITE ATOM=/POS= block out of CTRL and turning it into QE coordinates.

LMTO POS is in units of alat (like QE's ATOMIC_POSITIONS alat), so the
fractional z used by Wannier90 is  z_frac = |POS_z| / (c/a).

Usage:
    ctrl2conf.py CTRL --tm Mn
    ctrl2conf.py CTRL --tm Ru --template config/GdRuSi.conf > config/GdRuSi.conf
"""
import argparse
import re
import sys


def parse_ctrl(path):
    text = open(path).read()

    m = re.search(r"ALAT\s*=\s*([0-9.eEdD+-]+)", text)
    if not m:
        sys.exit("ERROR: no 'ALAT=' found in CTRL")
    alat = float(m.group(1).replace("D", "E").replace("d", "e"))

    # PLAT: three vectors of three numbers, possibly spread over lines
    m = re.search(r"PLAT\s*=\s*((?:[-.0-9eEdD+\s]+))", text)
    if not m:
        sys.exit("ERROR: no 'PLAT=' found in CTRL")
    nums = re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eEdD][-+]?[0-9]+)?", m.group(1))
    if len(nums) < 9:
        sys.exit("ERROR: could not read 9 PLAT components")
    plat = [float(x.replace("D", "E").replace("d", "e")) for x in nums[:9]]
    coa = plat[8]

    # SITE block: ATOM=Xx POS=x y z, repeated
    sites = []
    for m in re.finditer(
        r"ATOM\s*=\s*([A-Za-z]{1,2})\s+POS\s*=\s*"
        r"([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)", text):
        label = m.group(1)
        # LMTO writes '-.25' and '.25' which float() handles fine
        pos = [float(m.group(i)) for i in (2, 3, 4)]
        sites.append((label, pos))
    if not sites:
        sys.exit("ERROR: no 'ATOM= ... POS=' lines found in CTRL")

    return alat, plat, coa, sites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ctrl")
    ap.add_argument("--tm", required=True, help="transition metal symbol, e.g. Mn")
    ap.add_argument("--template", help="existing .conf to patch instead of a bare dump")
    args = ap.parse_args()

    alat, plat, coa, sites = parse_ctrl(args.ctrl)

    z = {}
    for label, pos in sites:
        z.setdefault(label, abs(pos[2]))

    if "Gd" not in z:
        sys.exit("ERROR: no Gd site in CTRL -- check the file")
    if "Si" not in z:
        sys.exit("ERROR: no Si site in CTRL -- check the file")

    z_gd = z["Gd"] / coa
    z_si = z["Si"] / coa

    vals = {
        "ALAT_BOHR": f"{alat:.8f}",
        "COA": f"{coa:.8f}",
        "Z_GD": f"{z_gd:.6f}",
        "Z_SI": f"{z_si:.6f}",
    }

    # sanity: PLAT must be diagonal-tetragonal for the generator below
    off = [plat[i] for i in (1, 2, 3, 5, 6, 7)]
    if any(abs(v) > 1e-8 for v in off):
        print("# WARNING: PLAT is not diagonal. The input generator assumes a",
              file=sys.stderr)
        print("#          tetragonal primitive cell -- edit CELL_PARAMETERS by hand.",
              file=sys.stderr)

    if args.template:
        out = []
        for line in open(args.template):
            key = line.split("=", 1)[0].strip()
            if key in vals:
                comment = ""
                if "#" in line:
                    comment = "  #" + line.split("#", 1)[1].rstrip()
                out.append(f"{key}={vals[key]}{comment}\n")
            else:
                out.append(line)
        sys.stdout.write("".join(out))
    else:
        print(f"# generated from {args.ctrl}")
        for k, v in vals.items():
            print(f"{k}={v}")

    print(f"\n# --- check ---", file=sys.stderr)
    print(f"# alat  = {alat:.8f} Bohr = {alat * 0.529177210903:.4f} Ang",
          file=sys.stderr)
    print(f"# c/a   = {coa:.8f}", file=sys.stderr)
    print(f"# sites in CTRL: {', '.join(s[0] for s in sites)}", file=sys.stderr)
    a_ang = alat * 0.529177210903
    d = a_ang * ((0.5 ** 2) + (z_si * coa) ** 2) ** 0.5
    print(f"# nearest {args.tm}-Si distance = {d:.4f} Ang", file=sys.stderr)


if __name__ == "__main__":
    main()
