#!/usr/bin/env python3
"""
Generate every input file for the hopping pipeline from one .conf.
Produces: scf.in, nscf.in, SEED_up.win, SEED_dn.win, pw2wan_up.in, pw2wan_dn.in
Usage:
    make_inputs.py config/GdMnSi.conf [-o workdir]
"""
import argparse
import os
import re
import sys

BOHR = 0.529177210903

# valence orbital counts for Wannier90 projection specifiers
ORB_COUNT = {"s": 1, "p": 3, "d": 5, "f": 7,
             "sp": 2, "sp2": 3, "sp3": 4, "sp3d": 5, "sp3d2": 6}


def read_conf(path):
    conf = {}
    for raw in open(path):
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        conf[k.strip()] = v.strip()
    return conf


def need(conf, key):
    v = conf.get(key)
    if v is None:
        sys.exit(f"ERROR: {key} missing from config")
    if v.upper().startswith("FIXME"):
        sys.exit(f"ERROR: {key} is still FIXME -- fill it in "
                 f"(bin/ctrl2conf.py can do the structural ones)")
    return v


def orbitals(spec):
    """'d' -> 5, 'dxy;dyz' -> 2, 'l=2' -> 5"""
    spec = spec.strip()
    if spec.startswith("l="):
        return 2 * int(spec[2:]) + 1
    if ";" in spec:
        return len(spec.split(";"))
    return ORB_COUNT.get(spec, None) or sys.exit(
        f"ERROR: cannot count orbitals in projection '{spec}'")


def kmesh(nk):
    n1, n2, n3 = nk
    return [(i / n1, j / n2, k / n3)
            for i in range(n1) for j in range(n2) for k in range(n3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("conf")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--force-scf", action="store_true",
                    help="overwrite an existing scf.in")
    args = ap.parse_args()

    c = read_conf(args.conf)
    os.makedirs(args.outdir, exist_ok=True)

    sysname = need(c, "SYSTEM")
    prefix = need(c, "PREFIX")
    tm = need(c, "TM")

    alat = float(need(c, "ALAT_BOHR"))
    coa = float(need(c, "COA"))
    z_gd = float(need(c, "Z_GD"))
    z_si = float(need(c, "Z_SI"))

    a_ang = alat * BOHR
    c_ang = a_ang * coa

    nk_nscf = tuple(int(x) for x in need(c, "KMESH_NSCF").split())
    mp_grid = tuple(int(x) for x in need(c, "MP_GRID").split())
    if nk_nscf != mp_grid:
        sys.exit(f"ERROR: KMESH_NSCF {nk_nscf} != MP_GRID {mp_grid}. "
                 "pw2wannier90 requires the identical mesh.")

    # ---- geometry -----------------------------------------------------
    # QE: ATOMIC_POSITIONS alat (z in units of alat, i.e. z_frac * c/a)
    qe_pos = [
        ("Gd", 0.25, 0.25, -z_gd * coa),
        ("Gd", -0.25, -0.25, z_gd * coa),
        (tm, 0.25, -0.25, 0.0),
        (tm, -0.25, 0.25, 0.0),
        ("Si", 0.25, 0.25, z_si * coa),
        ("Si", -0.25, -0.25, -z_si * coa),
    ]
    # Wannier90: same sites, wrapped into [0,1)
    w90_pos = [
        ("Gd", 0.25, 0.25, 1.0 - z_gd),
        ("Gd", 0.75, 0.75, z_gd),
        (tm, 0.25, 0.75, 0.0),
        (tm, 0.75, 0.25, 0.0),
        ("Si", 0.25, 0.25, z_si),
        ("Si", 0.75, 0.75, 1.0 - z_si),
    ]

    # ---- projections and num_wann ---------------------------------------
    projs = [p.strip() for p in need(c, "PROJECTIONS").split(";") if p.strip()]
    # a projection like 'Gd:d' may itself contain ';' for orbital lists;
    # rejoin those by looking for the ':' separator
    fixed, buf = [], ""
    for p in projs:
        if ":" in p:
            if buf:
                fixed.append(buf)
            buf = p
        else:
            buf += ";" + p
    if buf:
        fixed.append(buf)
    projs = fixed

    num_wann = 0
    proj_blocks = []          # (element, n_orb, n_atoms)
    for p in projs:
        elem, spec = p.split(":", 1)
        elem = elem.strip()
        n_orb = orbitals(spec)
        n_at = sum(1 for s in w90_pos if s[0] == elem)
        if n_at == 0:
            sys.exit(f"ERROR: projection '{p}' but no {elem} atom in the cell")
        num_wann += n_orb * n_at
        proj_blocks.append((elem, n_orb, n_at))

    nbnd = int(need(c, "NBND"))
    if num_wann > nbnd:
        sys.exit(f"ERROR: num_wann={num_wann} > NBND={nbnd}")

    # ---- shared blocks ---------------------------------------------------
    species = (
        f"ATOMIC_SPECIES\n"
        f"Gd   1.0 {need(c, 'PP_GD')}\n"
        f"{tm:<4s} 1.0 {need(c, 'PP_TM')}\n"
        f"Si   1.0 {need(c, 'PP_SI')}\n"
    )
    cell = ("CELL_PARAMETERS alat\n"
            "  1.00000000  0.00000000  0.00000000\n"
            "  0.00000000  1.00000000  0.00000000\n"
            f"  0.00000000  0.00000000 {coa:12.8f}\n")
    positions = "ATOMIC_POSITIONS alat\n" + "".join(
        f"{el:<4s} {x:13.8f} {y:13.8f} {z:13.8f}\n" for el, x, y, z in qe_pos)
    hubbard = (f"HUBBARD (ortho-atomic)\n"
               f"U Gd-4f {need(c, 'HUBBARD_U_GD')}\n"
               f"J0 Gd-4f {need(c, 'HUBBARD_J_GD')}\n")

    def system_block(extra=""):
        return (
            f"&SYSTEM\n"
            f"    ibrav = 0,\n"
            f"    celldm(1) = {alat:.8f},\n"
            f"    nat = 6,\n"
            f"    ntyp = 3,\n"
            f"    nspin = 2,\n"
            f"    nbnd = {nbnd},\n"
            f"    ecutwfc = {need(c, 'ECUTWFC')}.,\n"
            f"    ecutrho = {need(c, 'ECUTRHO')}.,\n"
            f"    starting_magnetization(1) = {need(c, 'START_MAG_GD')},\n"
            f"    starting_magnetization(2) = {need(c, 'START_MAG_TM')},\n"
            f"    starting_magnetization(3) = {need(c, 'START_MAG_SI')},\n"
            f"    occupations = 'smearing',\n"
            f"    smearing = 'marzari-vanderbilt',\n"
            f"    degauss = {need(c, 'DEGAUSS')}{extra}\n"
            f"/\n")

    electrons = ("&ELECTRONS\n"
                 "    startingpot = '{}',\n"
                 "    conv_thr = 1.0d-8,\n"
                 "    mixing_beta = 0.3,\n"
                 "    mixing_mode = 'plain',\n"
                 "    diagonalization = 'david',\n"
                 "    electron_maxstep = 200\n"
                 "/\n")

    def control(calc, title):
        return (f"&CONTROL\n"
                f"    prefix = '{prefix}',\n"
                f"    calculation = '{calc}',\n"
                f"    title = '{title}',\n"
                f"    pseudo_dir = '{need(c, 'PSEUDO_DIR')}',\n"
                f"    disk_io = 'low',\n"
                f"    outdir = './out',\n"
                f"    verbosity = 'high'\n"
                f"/\n")

    def write(name, text):
        path = os.path.join(args.outdir, name)
        with open(path, "w") as f:
            f.write(text)
        print(f"  {name}")

    # ---- scf.in ----------------------------------------------------------
    scf_text = (control("scf", f"{sysname}_scf")
                + system_block()
                + electrons.format("atomic")
                + species + "\n" + cell + "\n" + positions + "\n"
                + f"K_POINTS automatic\n{need(c, 'KMESH_SCF')} 0 0 0\n\n"
                + hubbard)
    scf_path = os.path.join(args.outdir, "scf.in")

    if os.path.exists(scf_path) and not args.force_scf:
        # Do not touch a completed run. Instead verify that what we would have
        # generated agrees with it -- if the geometry differs, the nscf would
        # read a charge density computed for a different structure and the
        # whole chain is silently wrong.
        old = open(scf_path).read()
        def poslines(t):
            m = re.search(r"ATOMIC_POSITIONS[^\n]*\n((?:\s*[A-Z][a-z]?\s+"
                          r"[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+\s*\n|\s*\n)+)", t)
            if not m:
                return None
            out = []
            for l in m.group(1).splitlines():
                p = l.split()
                if len(p) == 4:
                    out.append((p[0], *(round(float(x), 6) for x in p[1:])))
            return out
        a_old, a_new = poslines(old), poslines(scf_text)
        m_old = re.search(r"prefix\s*=\s*'([^']+)'", old)
        problems = []
        if a_old != a_new:
            problems.append("ATOMIC_POSITIONS differ from the config")
        if m_old and m_old.group(1) != prefix:
            problems.append(f"prefix in scf.in is '{m_old.group(1)}', "
                            f"config says '{prefix}'")
        if problems:
            print("\nERROR: the existing scf.in does not match this config:")
            for p in problems:
                print(f"  - {p}")
            print("Regenerate the config with bin/scf2conf.py, or pass "
                  "--force-scf to overwrite scf.in.")
            sys.exit(1)
        print("  scf.in (already present, kept -- geometry and prefix match)")
    else:
        write("scf.in", scf_text)

    # ---- nscf.in with an explicit k list ---------------------------------
    kpts = kmesh(nk_nscf)
    klist = (f"K_POINTS crystal\n{len(kpts)}\n"
             + "".join(f"{x:16.10f} {y:16.10f} {z:16.10f}  1.0\n"
                       for x, y, z in kpts))

    if os.path.exists(scf_path):
        # Derive nscf.in from the real scf.in by targeted edits, not from the
        # template. An scf.in for an AFM start (or any hand-tuned run) carries
        # keys this generator does not know about -- tot_magnetization,
        # constrained_magnetization, input_dft, per-type occupations. Rebuilding
        # from the template would drop them and the nscf would then describe a
        # different system than the charge density it reads.
        t = open(scf_path).read()

        def set_key(text, namelist, key, value):
            pat = rf"(&{namelist}\b.*?^\s*/)"
            m = re.search(pat, text, re.S | re.I | re.M)
            if not m:
                return text
            body = m.group(1)
            if re.search(rf"^\s*{key}\s*=", body, re.I | re.M):
                new = re.sub(rf"^(\s*){key}\s*=[^,\n]*(,?)\s*$",
                             rf"\g<1>{key} = {value}\g<2>",
                             body, count=1, flags=re.I | re.M)
            else:
                new = re.sub(r"(\n\s*/)$", f",\n    {key} = {value}\\1", body,
                             count=1)
            return text[:m.start(1)] + new + text[m.end(1):]

        t = set_key(t, "CONTROL", "calculation", "'nscf'")
        t = set_key(t, "SYSTEM", "nosym", ".TRUE.")
        t = set_key(t, "SYSTEM", "noinv", ".TRUE.")
        t = set_key(t, "ELECTRONS", "startingpot", "'file'")

        # replace the whole K_POINTS card
        t = re.sub(r"^[ \t]*K_POINTS.*\n(?:(?![ \t]*[A-Z_]{4,}).*\n)*",
                   klist + "\n", t, count=1, flags=re.M)

        nscf_text = t
        note = " (derived from scf.in, extra settings preserved)"
    else:
        nscf_text = (control("nscf", f"{sysname}_nscf_wannier")
                     + system_block(",\n    nosym = .TRUE.,\n    noinv = .TRUE.")
                     + electrons.format("file")
                     + species + "\n" + cell + "\n" + positions + "\n"
                     + klist + "\n" + hubbard)
        note = ""

    write("nscf.in", nscf_text)
    if note:
        print(f"        {note.strip()}")

    # ---- .win ------------------------------------------------------------
    kblock = ("begin kpoints\n"
              + "".join(f"{x:16.10f} {y:16.10f} {z:16.10f}\n"
                        for x, y, z in kpts)
              + "end kpoints\n")

    for spin in ("up", "dn"):
        seed = f"{sysname}_{spin}"
        win = (
            f"! generated from {os.path.basename(args.conf)} -- do not edit by hand\n"
            f"! num_wann = {' + '.join(f'{n_at}x{el}({n_orb})' for el, n_orb, n_at in proj_blocks)}\n\n"
            f"num_wann  = {num_wann}\n"
            f"num_bands = {nbnd}\n"
            f"spinors   = false\n\n"
            f"! WINDOWS ARE PLACEHOLDERS -- run bin/set_windows.py after nscf\n"
            f"dis_win_min  = 0.0\n"
            f"dis_win_max  = 0.0\n"
            f"dis_froz_min = 0.0\n"
            f"dis_froz_max = 0.0\n\n"
            f"dis_num_iter = {need(c, 'DIS_NUM_ITER')}\n"
            f"num_iter     = {need(c, 'NUM_ITER')}\n"
            f"conv_window  = 5\n"
            f"conv_tol     = 1.0d-10\n\n"
            f"write_hr  = true\n"
            f"write_tb  = true\n"
            f"write_xyz = true\n\n"
            f"bands_plot = true\n"
            f"bands_num_points = 100\n\n"
            f"mp_grid = {' '.join(str(n) for n in mp_grid)}\n\n"
            f"begin unit_cell_cart\n"
            f"Ang\n"
            f"{a_ang:11.6f} {0.0:11.6f} {0.0:11.6f}\n"
            f"{0.0:11.6f} {a_ang:11.6f} {0.0:11.6f}\n"
            f"{0.0:11.6f} {0.0:11.6f} {c_ang:11.6f}\n"
            f"end unit_cell_cart\n\n"
            f"begin atoms_frac\n"
            + "".join(f"{el:<3s} {x:12.8f} {y:12.8f} {z:12.8f}\n"
                      for el, x, y, z in w90_pos)
            + f"end atoms_frac\n\n"
            f"begin projections\n"
            + "".join(f"{p}\n" for p in projs)
            + f"end projections\n\n"
            f"begin kpoint_path\n"
            f"G  0.0000 0.0000 0.0000   X  0.5000 0.0000 0.0000\n"
            f"X  0.5000 0.0000 0.0000   M  0.5000 0.5000 0.0000\n"
            f"M  0.5000 0.5000 0.0000   G  0.0000 0.0000 0.0000\n"
            f"G  0.0000 0.0000 0.0000   Z  0.0000 0.0000 0.5000\n"
            f"end kpoint_path\n\n"
            + kblock)
        write(f"{seed}.win", win)

        pw2wan = (f"&inputpp\n"
                  f"    outdir = './out',\n"
                  f"    prefix = '{prefix}',\n"
                  f"    seedname = '{seed}',\n"
                  f"    spin_component = '{'up' if spin == 'up' else 'down'}',\n"
                  f"    write_mmn = .true.,\n"
                  f"    write_amn = .true.,\n"
                  f"    write_eig = .true.\n"
                  f"/\n")
        write(f"pw2wan_{spin}.in", pw2wan)

    # ---- report -----------------------------------------------------------
    d_tm_si = a_ang * (0.25 + (z_si * coa) ** 2) ** 0.5
    print(f"\n{sysname}:  a = {a_ang:.5f} A   c = {c_ang:.5f} A   c/a = {coa:.6f}")
    print(f"  num_wann = {num_wann}   num_bands = {nbnd}   k-mesh = "
          f"{'x'.join(str(n) for n in mp_grid)} ({len(kpts)} points)")
    print(f"  nearest {tm}-Si = {d_tm_si:.4f} A")
    print(f"\nNext: run scf, then nscf, then bin/set_windows.py")


if __name__ == "__main__":
    main()
