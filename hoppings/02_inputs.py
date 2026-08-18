#!/usr/bin/env python3
"""
02 -- system.conf -> nscf.in, <SYSTEM>_{up,dn}.win, pw2wan_{up,dn}.in

The cell and the atomic positions come from scf.in, not from idealised
coordinates, so relaxed structures and magnetic supercells pass through
unchanged. QE species that split a magnetic sublattice (Mn1/Mn2) are merged
into the chemical element for the projections: projections act on the element,
and the spin channel is chosen by pw2wannier90.

    02_inputs.py                       everything
    02_inputs.py --only win            just the .win files
    02_inputs.py --only nscf           just nscf.in
    02_inputs.py --only win --spin up  one channel
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp


def orbitals(spec):
    spec = spec.strip()
    if spec.startswith("l="):
        return 2 * int(spec[2:]) + 1
    if ";" in spec:
        return len(spec.split(";"))
    n = qp.ORB_COUNT.get(spec)
    if n is None:
        sys.exit(f"ERROR: cannot count orbitals in projection '{spec}'")
    return n


def split_projections(raw):
    """'Gd:d;Mn:d' -> ['Gd:d', 'Mn:d'], tolerating 'Fe:dxy;dyz'."""
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
    ap.add_argument("conf", nargs="?", default="system.conf")
    ap.add_argument("-o", "--outdir", default=None)
    ap.add_argument("--only", choices=["all", "nscf", "win", "pw2wan"],
                    default="all")
    ap.add_argument("--spin", choices=["up", "dn", "both"], default="both")
    args = ap.parse_args()

    c = qp.read_conf(args.conf)
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.conf))
    os.makedirs(outdir, exist_ok=True)
    scf = os.path.join(outdir, "scf.in")
    spins = ["up", "dn"] if args.spin == "both" else [args.spin]

    sysname, prefix, tm = (qp.need(c, k) for k in ("SYSTEM", "PREFIX", "TM"))
    nk = tuple(int(x) for x in qp.need(c, "KMESH_NSCF").split())
    mp = tuple(int(x) for x in qp.need(c, "MP_GRID").split())
    if nk != mp:
        sys.exit(f"ERROR: KMESH_NSCF {nk} != MP_GRID {mp}; pw2wannier90 needs "
                 f"the identical mesh")

    if not os.path.exists(scf):
        sys.exit(f"ERROR: no scf.in in {os.path.abspath(outdir)}. This toolkit "
                 f"starts from a completed SCF.")
    text = qp.strip_comments(open(scf).read())
    got = qp.namelist(text, "CONTROL").get("prefix")
    if got and got != prefix:
        sys.exit(f"ERROR: prefix in scf.in is '{got}', config says '{prefix}'")

    st = qp.read_structure(text)
    lattice = st["lattice"]
    atoms = [(el, f) for _, el, f in st["atoms"]]
    nbnd = int(c.get("NBND_NSCF")
               or qp.fnum(qp.namelist(text, "SYSTEM").get("nbnd",
                                                          qp.need(c, "NBND"))))

    counts = {}
    for el, _ in atoms:
        counts[el] = counts.get(el, 0) + 1

    projs = split_projections(qp.need(c, "PROJECTIONS"))
    num_wann, blocks = 0, []
    for p in projs:
        el, spec = p.split(":", 1)
        el = el.strip()
        if el not in counts:
            sys.exit(f"ERROR: projection '{p}' but no {el} in the cell "
                     f"({', '.join(sorted(counts))})")
        n = orbitals(spec)
        num_wann += n * counts[el]
        blocks.append((el, n, counts[el]))
    unprojected = [e for e in counts
                   if e not in {p.split(":", 1)[0].strip() for p in projs}]
    if num_wann > nbnd:
        sys.exit(f"ERROR: num_wann={num_wann} > num_bands={nbnd}")

    kpts = [(i / nk[0], j / nk[1], k / nk[2])
            for i in range(nk[0]) for j in range(nk[1]) for k in range(nk[2])]
    klist = (f"K_POINTS crystal\n{len(kpts)}\n"
             + "".join(f"{x:16.10f} {y:16.10f} {z:16.10f}  1.0\n"
                       for x, y, z in kpts))

    written = []

    def write(name, body):
        open(os.path.join(outdir, name), "w").write(body)
        written.append(name)

    # ---- nscf.in: four targeted edits to the real scf.in ------------------
    # Rebuilding from a template would drop non-default keys and the nscf would
    # then describe a different system than the density it reads.
    if args.only in ("all", "nscf"):
        t = open(scf).read()
        t = set_key(t, "CONTROL", "calculation", "'nscf'")
        t = set_key(t, "SYSTEM", "nbnd", str(nbnd))
        t = set_key(t, "SYSTEM", "nosym", ".TRUE.")
        t = set_key(t, "SYSTEM", "noinv", ".TRUE.")
        t = set_key(t, "ELECTRONS", "startingpot", "'file'")
        t = re.sub(r"^[ \t]*K_POINTS.*\n(?:(?![ \t]*[A-Z_]{4,}).*\n)*",
                   klist + "\n", t, count=1, flags=re.M)
        write("nscf.in", t)

    # ---- band path ---------------------------------------------------------
    # Never assumed: high-symmetry points depend on the lattice, and
    # coordinates written for a primitive cell mean something else in a
    # supercell. tools/kpath.py converts an LMTO SYML block.
    kblock, knote = "", "none set (no band plot)"
    raw = c.get("KPATH", "").strip()
    if raw:
        pts = []
        for chunk in raw.split(":"):
            f = chunk.split()
            if len(f) != 4:
                sys.exit(f"ERROR: bad KPATH entry '{chunk.strip()}'; expected "
                         f"'LABEL x y z' separated by ':'")
            pts.append((f[0], [float(v) for v in f[1:]]))
        rows = [f"{pts[i][0]:<3s}{pts[i][1][0]:9.4f}{pts[i][1][1]:9.4f}"
                f"{pts[i][1][2]:9.4f}   {pts[i+1][0]:<3s}"
                f"{pts[i+1][1][0]:9.4f}{pts[i+1][1][1]:9.4f}"
                f"{pts[i+1][1][2]:9.4f}" for i in range(len(pts) - 1)]
        kblock = ("bands_plot = true\nbands_num_points = "
                  f"{c.get('BANDS_NUM_POINTS', '100')}\n\n"
                  "begin kpoint_path\n" + "\n".join(rows)
                  + "\nend kpoint_path\n\n")
        knote = f"{len(pts)} points from KPATH"

    # ---- .win and pw2wan ---------------------------------------------------
    for spin in spins:
        seed = f"{sysname}_{spin}"
        if args.only in ("all", "win"):
            write(f"{seed}.win",
                  f"! generated by 02_inputs.py from "
                  f"{os.path.basename(args.conf)}\n"
                  f"! geometry from scf.in\n"
                  f"! num_wann = "
                  + " + ".join(f"{na}x{el}({no})" for el, no, na in blocks)
                  + f"\n\nnum_wann  = {num_wann}\n"
                  f"num_bands = {nbnd}\nspinors   = false\n\n"
                  f"! placeholder windows -- 04_windows.py sets these\n"
                  f"dis_win_min  = 0.0\ndis_win_max  = 0.0\n"
                  f"dis_froz_min = 0.0\ndis_froz_max = 0.0\n\n"
                  f"dis_num_iter    = {qp.need(c, 'DIS_NUM_ITER')}\n"
                  f"dis_conv_tol    = {c.get('DIS_CONV_TOL', '1.0d-9')}\n"
                  f"dis_conv_window = {c.get('DIS_CONV_WINDOW', '3')}\n"
                  f"num_iter     = {qp.need(c, 'NUM_ITER')}\n"
                  f"conv_tol     = {c.get('CONV_TOL', '1.0d-10')}\n"
                  f"conv_window  = {c.get('CONV_WINDOW', '5')}\n\n"
                  f"write_hr  = true\nwrite_tb  = true\nwrite_xyz = true\n\n"
                  + kblock
                  + f"mp_grid = {' '.join(str(n) for n in mp)}\n\n"
                  f"begin unit_cell_cart\nAng\n"
                  + "".join(f"{v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}\n"
                            for v in lattice)
                  + "end unit_cell_cart\n\nbegin atoms_frac\n"
                  + "".join(f"{el:<3s} {f[0]:12.8f} {f[1]:12.8f} {f[2]:12.8f}\n"
                            for el, f in atoms)
                  + "end atoms_frac\n\nbegin projections\n"
                  + "".join(f"{p}\n" for p in projs)
                  + "end projections\n\nbegin kpoints\n"
                  + "".join(f"{x:16.10f} {y:16.10f} {z:16.10f}\n"
                            for x, y, z in kpts)
                  + "end kpoints\n")
        if args.only in ("all", "pw2wan"):
            write(f"pw2wan_{spin}.in",
                  f"&inputpp\n    outdir = './out',\n    prefix = '{prefix}',\n"
                  f"    seedname = '{seed}',\n"
                  f"    spin_component = '{'up' if spin == 'up' else 'down'}',\n"
                  f"    write_mmn = .true.,\n    write_amn = .true.,\n"
                  f"    write_eig = .true.\n/\n")

    lens = [sum(v * v for v in row) ** 0.5 for row in lattice]
    print(f"{sysname}  {lens[0]:.4f} x {lens[1]:.4f} x {lens[2]:.4f} A, "
          f"{len(atoms)} atoms: "
          + ", ".join(f"{n} {e}" for e, n in counts.items()))
    print(f"  num_wann {num_wann}   num_bands {nbnd}   mesh "
          f"{'x'.join(str(n) for n in mp)} ({len(kpts)} k-points)")
    print(f"  band path: {knote}")
    if unprojected:
        print(f"  WARNING: no projections on {', '.join(unprojected)}")
    print("  wrote " + ", ".join(written))


if __name__ == "__main__":
    main()
