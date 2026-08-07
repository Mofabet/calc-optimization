#!/usr/bin/env python3
"""
Generate the downstream inputs for the hopping pipeline.

Produces: nscf.in, <SYSTEM>_up.win, <SYSTEM>_dn.win, pw2wan_up.in, pw2wan_dn.in
(and scf.in only when one is not already present).

The cell and the atomic positions are read from scf.in, not rebuilt from
idealised coordinates, so magnetic supercells and relaxed structures pass
through unchanged. QE species that split a magnetic sublattice (Mn_up, Mn_dn)
are merged back into the chemical element for the Wannier projections: the
projections act on the element, and the spin channel is chosen by
pw2wannier90 through spin_component.

Usage:
    make_inputs.py system.conf [-o workdir]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp

ORB_COUNT = {"s": 1, "p": 3, "d": 5, "f": 7,
             "sp": 2, "sp2": 3, "sp3": 4, "sp3d": 5, "sp3d2": 6}


def read_conf(path):
    conf = {}
    for raw in open(path):
        line = raw.split("#", 1)[0].strip()
        if line and "=" in line:
            k, v = line.split("=", 1)
            conf[k.strip()] = v.strip()
    return conf


def need(conf, key):
    v = conf.get(key)
    if v is None:
        sys.exit(f"ERROR: {key} missing from config")
    if v.upper().startswith("FIXME"):
        sys.exit(f"ERROR: {key} is still FIXME")
    return v


def orbitals(spec):
    spec = spec.strip()
    if spec.startswith("l="):
        return 2 * int(spec[2:]) + 1
    if ";" in spec:
        return len(spec.split(";"))
    n = ORB_COUNT.get(spec)
    if n is None:
        sys.exit(f"ERROR: cannot count orbitals in projection '{spec}'")
    return n


def split_projections(raw):
    """'Gd:d;Mn:d;Si:p' -> ['Gd:d', 'Mn:d', 'Si:p'], tolerating 'Fe:dxy;dyz'."""
    out, buf = [], ""
    for piece in raw.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        if ":" in piece:
            if buf:
                out.append(buf)
            buf = piece
        else:
            buf += ";" + piece
    if buf:
        out.append(buf)
    return out


def kmesh(nk):
    n1, n2, n3 = nk
    return [(i / n1, j / n2, k / n3)
            for i in range(n1) for j in range(n2) for k in range(n3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("conf")
    ap.add_argument("-o", "--outdir", default=None,
                    help="where to write; defaults to the directory holding "
                         "the config, not the current directory")
    args = ap.parse_args()

    c = read_conf(args.conf)
    # The inputs belong next to the SCF they are derived from, which is where
    # the config lives. Defaulting to the current directory made
    # "make_inputs.py z12/system.conf" write into the parent and then fail on
    # a missing key, because no scf.in was found there.
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.conf))
    args.outdir = outdir
    os.makedirs(outdir, exist_ok=True)
    scf_path = os.path.join(outdir, "scf.in")

    sysname = need(c, "SYSTEM")
    prefix = need(c, "PREFIX")
    tm = need(c, "TM")

    nk = tuple(int(x) for x in need(c, "KMESH_NSCF").split())
    mp = tuple(int(x) for x in need(c, "MP_GRID").split())
    if nk != mp:
        sys.exit(f"ERROR: KMESH_NSCF {nk} != MP_GRID {mp}. "
                 "pw2wannier90 requires the identical mesh.")

    # ---- geometry --------------------------------------------------------
    # scf.in is the authority on the structure. Reading it rather than
    # rebuilding from idealised coordinates is what lets relaxed structures
    # and magnetic supercells pass through untouched.
    if not os.path.exists(scf_path):
        sys.exit(f"ERROR: no scf.in in {os.path.abspath(outdir)}.\n"
                 f"       This toolkit starts from a completed SCF. Run one "
                 f"first, or point -o at the directory that holds it.")
    text = qp.strip_comments(open(scf_path).read())
    got = qp.namelist(text, "CONTROL").get("prefix")
    if got and got != prefix:
        sys.exit(f"ERROR: prefix in scf.in is '{got}' but the config says "
                 f"'{prefix}'. The config belongs to a different run.")
    st = qp.read_structure(text)
    lattice = st["lattice"]
    w90_atoms = [(el, f[0], f[1], f[2]) for _, el, f in st["atoms"]]
    nbnd = int(c.get("NBND_NSCF")
               or qp.fnum(qp.namelist(text, "SYSTEM").get("nbnd",
                                                          need(c, "NBND"))))

    counts = {}
    for el, *_ in w90_atoms:
        counts[el] = counts.get(el, 0) + 1

    # ---- projections ------------------------------------------------------
    projs = split_projections(need(c, "PROJECTIONS"))
    num_wann, blocks, unprojected = 0, [], []
    for p in projs:
        el, spec = p.split(":", 1)
        el = el.strip()
        if el not in counts:
            sys.exit(f"ERROR: projection '{p}' but no {el} atom in the cell. "
                     f"Present: {', '.join(sorted(counts))}")
        n = orbitals(spec)
        num_wann += n * counts[el]
        blocks.append((el, n, counts[el]))
    projected = {p.split(":", 1)[0].strip() for p in projs}
    unprojected = [e for e in counts if e not in projected]
    if num_wann > nbnd:
        sys.exit(f"ERROR: num_wann={num_wann} > num_bands={nbnd}")

    def write(name, text):
        with open(os.path.join(args.outdir, name), "w") as f:
            f.write(text)
        print(f"  {name}")

    kpts = kmesh(nk)
    klist = (f"K_POINTS crystal\n{len(kpts)}\n"
             + "".join(f"{x:16.10f} {y:16.10f} {z:16.10f}  1.0\n"
                       for x, y, z in kpts))

    # ---- nscf.in -----------------------------------------------------------
    # Four targeted edits to the real scf.in. Rebuilding from a template would
    # drop non-default keys -- tot_magnetization, constrained_magnetization,
    # input_dft, per-species settings -- and the nscf would then describe a
    # different system than the charge density it reads.
    t = open(scf_path).read()

    def set_key(text, nml, key, value):
        m = re.search(rf"(&{nml}\b.*?^\s*/)", text, re.S | re.I | re.M)
        if not m:
            return text
        body = m.group(1)
        if re.search(rf"^\s*{key}\s*=", body, re.I | re.M):
            new = re.sub(rf"^(\s*){key}\s*=[^,\n!]*(,?)\s*$",
                         rf"\g<1>{key} = {value}\g<2>", body,
                         count=1, flags=re.I | re.M)
        else:
            new = re.sub(r"(\n\s*/)$", f",\n    {key} = {value}\\1", body,
                         count=1)
        return text[:m.start(1)] + new + text[m.end(1):]

    t = set_key(t, "CONTROL", "calculation", "'nscf'")
    t = set_key(t, "SYSTEM", "nbnd", str(nbnd))
    t = set_key(t, "SYSTEM", "nosym", ".TRUE.")
    t = set_key(t, "SYSTEM", "noinv", ".TRUE.")
    t = set_key(t, "ELECTRONS", "startingpot", "'file'")
    t = re.sub(r"^[ \t]*K_POINTS.*\n(?:(?![ \t]*[A-Z_]{4,}).*\n)*",
               klist + "\n", t, count=1, flags=re.M)
    write("nscf.in", t)
    print("        (derived from scf.in, other settings preserved)")

    # ---- .win and pw2wan ----------------------------------------------------
    for spin in ("up", "dn"):
        seed = f"{sysname}_{spin}"
        win = (
            f"! generated by make_inputs.py from {os.path.basename(args.conf)}\n"
            f"! geometry source: scf.in\n"
            f"! num_wann = "
            + " + ".join(f"{na}x{el}({no})" for el, no, na in blocks) + "\n\n"
            f"num_wann  = {num_wann}\n"
            f"num_bands = {nbnd}\n"
            f"spinors   = false\n\n"
            f"! placeholders -- run bin/set_windows.py\n"
            f"dis_win_min  = 0.0\n"
            f"dis_win_max  = 0.0\n"
            f"dis_froz_min = 0.0\n"
            f"dis_froz_max = 0.0\n\n"
            f"dis_num_iter   = {need(c, 'DIS_NUM_ITER')}\n"
            f"dis_conv_tol    = {c.get('DIS_CONV_TOL', '1.0d-9')}\n"
            f"dis_conv_window = {c.get('DIS_CONV_WINDOW', '3')}\n"
            f"num_iter     = {need(c, 'NUM_ITER')}\n"
            f"conv_tol     = {c.get('CONV_TOL', '1.0d-10')}\n"
            f"conv_window  = {c.get('CONV_WINDOW', '5')}\n\n"
            f"write_hr  = true\n"
            f"write_tb  = true\n"
            f"write_xyz = true\n\n"
            f"bands_plot = true\n"
            f"bands_num_points = 100\n\n"
            f"mp_grid = {' '.join(str(n) for n in mp)}\n\n"
            f"begin unit_cell_cart\n"
            f"Ang\n"
            + "".join(f"{v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}\n"
                      for v in lattice)
            + "end unit_cell_cart\n\n"
            f"begin atoms_frac\n"
            + "".join(f"{el:<3s} {x:12.8f} {y:12.8f} {z:12.8f}\n"
                      for el, x, y, z in w90_atoms)
            + "end atoms_frac\n\n"
            f"begin projections\n"
            + "".join(f"{p}\n" for p in projs)
            + "end projections\n\n"
            f"begin kpoint_path\n"
            f"G  0.0000 0.0000 0.0000   X  0.5000 0.0000 0.0000\n"
            f"X  0.5000 0.0000 0.0000   M  0.5000 0.5000 0.0000\n"
            f"M  0.5000 0.5000 0.0000   G  0.0000 0.0000 0.0000\n"
            f"G  0.0000 0.0000 0.0000   Z  0.0000 0.0000 0.5000\n"
            f"end kpoint_path\n\n"
            f"begin kpoints\n"
            + "".join(f"{x:16.10f} {y:16.10f} {z:16.10f}\n"
                      for x, y, z in kpts)
            + "end kpoints\n")
        write(f"{seed}.win", win)
        write(f"pw2wan_{spin}.in",
              f"&inputpp\n"
              f"    outdir = './out',\n"
              f"    prefix = '{prefix}',\n"
              f"    seedname = '{seed}',\n"
              f"    spin_component = '{'up' if spin == 'up' else 'down'}',\n"
              f"    write_mmn = .true.,\n"
              f"    write_amn = .true.,\n"
              f"    write_eig = .true.\n"
              f"/\n")

    # ---- report --------------------------------------------------------------
    lens = [sum(v * v for v in row) ** 0.5 for row in lattice]
    print(f"\n{sysname}:  cell {lens[0]:.4f} x {lens[1]:.4f} x {lens[2]:.4f} A"
          f"   {len(w90_atoms)} atoms")
    print(f"  composition  {', '.join(f'{n} {e}' for e, n in counts.items())}")
    print(f"  num_wann {num_wann}   num_bands {nbnd}   k-mesh "
          f"{'x'.join(str(n) for n in mp)} ({len(kpts)} points)")
    if unprojected:
        print(f"  WARNING: no projections on {', '.join(unprojected)} -- "
              f"those atoms contribute no Wannier functions")


if __name__ == "__main__":
    main()
