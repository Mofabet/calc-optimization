#!/usr/bin/env python3
"""
Build a .conf from an SCF calculation that has already been run.
Usage:
    scf2conf.py                          # scf.in / scf.out in the current dir
    scf2conf.py -i scf.in -o out.conf
    scf2conf.py --nscf-mesh 8 8 8 --projections 'Gd:d;Mn:d;Si:p'
"""
import argparse
import os
import re
import sys

BOHR = 0.529177210903


# ------------------------------------------------------------ parsing ----
def strip_comments(text):
    return re.sub(r"!.*", "", text)


def namelist(text, name):
    m = re.search(rf"&{name}(.*?)^\s*/", text, re.S | re.I | re.M)
    if not m:
        return {}
    out = {}
    for k, v in re.findall(r"([A-Za-z_]+(?:\(\d+\))?)\s*=\s*([^,\n]+)",
                           m.group(1)):
        out[k.strip().lower()] = v.strip().strip("',\" ")
    return out


def card(text, name):
    """Return (unit, [lines]) for a QE card."""
    pat = (rf"^[ \t]*{name}[ \t]*(\(?[A-Za-z_-]*\)?)[ \t]*$\n"
           r"((?:(?![ \t]*[A-Z_]{4,}).*\n)*)")
    m = re.search(pat, text, re.M)
    if not m:
        return None, []
    unit = m.group(1).strip("() ").lower()
    lines = [l.strip() for l in m.group(2).splitlines() if l.strip()]
    return unit, lines


def fnum(s):
    return float(s.replace("d", "e").replace("D", "E"))


# --------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="scf.in")
    ap.add_argument("-l", "--log", default="scf.out")
    ap.add_argument("-o", "--output", default="-")
    ap.add_argument("--nscf-mesh", nargs=3, type=int, default=[8, 8, 8])
    ap.add_argument("--projections", default=None,
                    help="default: Gd:d;<TM>:d;Si:p")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: {args.input} not found")
    text = strip_comments(open(args.input).read())

    ctrl = namelist(text, "CONTROL")
    syst = namelist(text, "SYSTEM")

    prefix = ctrl.get("prefix")
    if not prefix:
        sys.exit("ERROR: no prefix in &CONTROL")
    outdir = ctrl.get("outdir", "./out")

    ibrav = int(fnum(syst.get("ibrav", "0")))
    if "celldm(1)" in syst:
        alat = fnum(syst["celldm(1)"])
    elif "a" in syst:
        alat = fnum(syst["a"]) / BOHR
    else:
        sys.exit("ERROR: neither celldm(1) nor A found in &SYSTEM")

    # ---- c/a ----------------------------------------------------------
    if ibrav == 0:
        unit, lines = card(text, "CELL_PARAMETERS")
        if not lines:
            sys.exit("ERROR: ibrav=0 but no CELL_PARAMETERS card")
        vec = [[fnum(x) for x in l.split()[:3]] for l in lines[:3]]
        if unit in ("angstrom", "ang"):
            coa = vec[2][2] / vec[0][0]
        elif unit == "bohr":
            coa = vec[2][2] / vec[0][0]
        else:                                   # alat
            coa = vec[2][2] / vec[0][0]
        off = [vec[0][1], vec[0][2], vec[1][0], vec[1][2], vec[2][0], vec[2][1]]
        if any(abs(v) > 1e-6 for v in off):
            sys.exit("ERROR: the cell is not diagonal. This toolkit assumes a "
                     "tetragonal primitive cell; edit CELL_PARAMETERS by hand.")
    elif ibrav in (6, 7):
        coa = fnum(syst.get("celldm(3)", "0"))
        if coa == 0:
            sys.exit("ERROR: ibrav=6 but celldm(3) missing")
    else:
        sys.exit(f"ERROR: ibrav={ibrav} not supported (use 0 or 6)")

    # ---- species ------------------------------------------------------
    _, sp_lines = card(text, "ATOMIC_SPECIES")
    pseudo = {}
    for l in sp_lines:
        p = l.split()
        if len(p) >= 3:
            pseudo[p[0]] = p[2]
    if not pseudo:
        sys.exit("ERROR: ATOMIC_SPECIES card not found")

    tm = next((e for e in pseudo if e not in ("Gd", "Si")), None)
    if tm is None:
        sys.exit(f"ERROR: cannot identify the transition metal among "
                 f"{list(pseudo)}")

    # ---- positions -> fractional z ------------------------------------
    unit, pos_lines = card(text, "ATOMIC_POSITIONS")
    if not pos_lines:
        sys.exit("ERROR: ATOMIC_POSITIONS card not found")
    zfrac = {}
    for l in pos_lines:
        p = l.split()
        if len(p) < 4:
            continue
        el, z = p[0], fnum(p[3])
        if unit in ("crystal", "crystal_sg"):
            zf = abs(z) % 1.0
        elif unit == "angstrom":
            zf = abs(z) / (alat * BOHR * coa)
        elif unit == "bohr":
            zf = abs(z) / (alat * coa)
        else:                                   # alat (also the QE default)
            zf = abs(z) / coa
        zf = min(zf, 1.0 - zf) if zf > 0 else 0.0
        zfrac.setdefault(el, zf)

    for need_el in ("Gd", "Si"):
        if need_el not in zfrac:
            sys.exit(f"ERROR: no {need_el} in ATOMIC_POSITIONS")

    # ---- magnetisation, Hubbard, cutoffs -------------------------------
    order = [l.split()[0] for l in sp_lines if l.split()]
    mag = {}
    for k, v in syst.items():
        m = re.match(r"starting_magnetization\((\d+)\)", k)
        if m:
            i = int(m.group(1)) - 1
            if i < len(order):
                mag[order[i]] = fnum(v)

    hub_u, hub_j = "6.7", "0.7"
    for l in card(text, "HUBBARD")[1]:
        p = l.split()
        if len(p) >= 3 and p[1].lower().startswith("gd"):
            if p[0].upper() == "U":
                hub_u = p[2]
            elif p[0].upper().startswith("J"):
                hub_j = p[2]

    kunit, klines = card(text, "K_POINTS")
    kscf = klines[0].split()[:3] if klines and kunit == "automatic" else \
        ["12", "12", "12"]

    # ---- scf.out -------------------------------------------------------
    ef, converged, moments, tot_mag = None, None, [], None
    if os.path.exists(args.log):
        log = open(args.log, errors="ignore").read()
        converged = "convergence has been achieved" in log
        m = re.findall(r"the Fermi energy is\s+([-0-9.]+)\s*ev", log)
        if m:
            ef = float(m[-1])
        m = re.findall(r"total magnetization\s*=\s*([-0-9.]+)", log)
        if m:
            tot_mag = float(m[-1])
        blocks = re.findall(
            r"Magnetic moment per site.*?\n((?:\s*atom.*\n)+)", log)
        if blocks:
            for l in blocks[-1].splitlines():
                mm = re.search(r"atom\s+(\d+).*?magn=\s*([-0-9.]+)", l)
                if mm:
                    moments.append((int(mm.group(1)), float(mm.group(2))))

    # ---- emit -----------------------------------------------------------
    sysname = f"Gd{tm}Si"
    projections = args.projections or f"Gd:d;{tm}:d;Si:p"
    mesh = " ".join(str(n) for n in args.nscf_mesh)
    a_ang = alat * BOHR
    d_tm_si = a_ang * (0.25 + (zfrac["Si"] * coa) ** 2) ** 0.5

    conf = f"""# generated by scf2conf.py from {args.input}
# z(Si) = {zfrac['Si']:.6f}   d({tm}-Si) = {d_tm_si:.4f} Ang

SYSTEM={sysname}
PREFIX={prefix}
TM={tm}

ALAT_BOHR={alat:.8f}
COA={coa:.8f}
Z_GD={zfrac['Gd']:.6f}
Z_SI={zfrac['Si']:.6f}

PSEUDO_DIR={ctrl.get('pseudo_dir', './pseudo')}
PP_GD={pseudo.get('Gd', 'FIXME')}
PP_TM={pseudo.get(tm, 'FIXME')}
PP_SI={pseudo.get('Si', 'FIXME')}

ECUTWFC={fnum(syst.get('ecutwfc', '60')):.0f}
ECUTRHO={fnum(syst.get('ecutrho', '600')):.0f}
NBND={int(fnum(syst.get('nbnd', '80')))}
DEGAUSS={syst.get('degauss', '0.01')}
KMESH_SCF={' '.join(kscf)}
KMESH_NSCF={mesh}

START_MAG_GD={mag.get('Gd', 0.6)}
START_MAG_TM={mag.get(tm, 0.0)}
START_MAG_SI={mag.get('Si', 0.0)}

HUBBARD_U_GD={hub_u}
HUBBARD_J_GD={hub_j}

MP_GRID={mesh}
PROJECTIONS={projections}
DIS_WIN_LO=6.0
DIS_WIN_HI=5.0
DIS_FROZ_LO=5.0
DIS_FROZ_HI=1.5
DIS_NUM_ITER=1000
NUM_ITER=500

PAIRS={tm}-Si,{tm}-{tm},Gd-{tm}
DMAX=5.0
TMIN=0.02

QE_BIN=/mnt/wizardry/qe/build/bin
W90_BIN=
SRUN_OPTS=--mpi=pmix -n 48 -p trollen --exclusive
"""
    if args.output == "-":
        sys.stdout.write(conf)
    else:
        open(args.output, "w").write(conf)
        print(f"wrote {args.output}")

    # ---- report to stderr -------------------------------------------------
    e = sys.stderr
    print(f"\n--- {args.input} ---", file=e)
    print(f"  prefix   {prefix}   outdir {outdir}", file=e)
    print(f"  {tm} system, a = {a_ang:.5f} A, c/a = {coa:.6f}", file=e)
    print(f"  z(Gd) = {zfrac['Gd']:.5f}   z(Si) = {zfrac['Si']:.5f}", file=e)
    print(f"  nearest {tm}-Si = {d_tm_si:.4f} A", file=e)
    if converged is None:
        print(f"  {args.log} not found -- SCF state unknown", file=e)
    else:
        print(f"  SCF converged: {'yes' if converged else 'NO'}", file=e)
        if ef is not None:
            print(f"  E_F = {ef:.4f} eV", file=e)
        if tot_mag is not None:
            print(f"  total magnetization = {tot_mag:.3f} uB/cell", file=e)
        if moments:
            txt = "  ".join(f"{i}:{v:+.3f}" for i, v in moments)
            print(f"  site moments (uB): {txt}", file=e)

    d = os.path.join(os.path.dirname(args.input) or ".",
                     outdir.lstrip("./"))
    if not os.path.isdir(d):
        print(f"\n  WARNING: {outdir} not found next to {args.input}.", file=e)
        print(f"  The nscf needs the charge density from the SCF run.",
              file=e)
        print(f"  Copy the whole {outdir} directory, not just the .in/.out.",
              file=e)


if __name__ == "__main__":
    main()
