#!/usr/bin/env python3
"""
Build a .conf from an SCF calculation that has already been run.

Works for any cell: primitive, magnetic supercell, or a structure whose
magnetic sublattices are split into separate QE species (Mn_up / Mn_dn). The
geometry is never re-derived from idealised Wyckoff positions -- it is read
from ATOMIC_POSITIONS as written.

Usage:
    scf2conf.py                          # scf.in / scf.out in the current dir
    scf2conf.py -o system.conf
    scf2conf.py -o /dev/null             # inspect only, write nothing
    scf2conf.py --nscf-mesh 8 8 4 --projections 'Mn:d;Si:p'
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="scf.in")
    ap.add_argument("-l", "--log", default="scf.out")
    ap.add_argument("-o", "--output", default="-")
    ap.add_argument("--nscf-mesh", nargs=3, type=int, default=None,
                    help="explicit mesh; overrides --mesh-base")
    ap.add_argument("--mesh-base", type=int, default=8,
                    help="points along the shortest axis; the others are "
                         "scaled from the cell shape (default 8)")
    ap.add_argument("--projections", default=None)
    ap.add_argument("--defaults", default=None,
                    help="shared policy file applied to every structure "
                         "(default: ../defaults.conf if it exists)")
    ap.add_argument("--kmesh-primitive", nargs=3, type=int, default=None,
                    help="mesh for the primitive cell; divided by the "
                         "detected supercell multiplicity")
    args = ap.parse_args()

    # ---- shared policy ---------------------------------------------------
    POLICY = ("PROJECTIONS", "DIS_WIN_LO", "DIS_WIN_HI", "DIS_FROZ_LO",
              "DIS_FROZ_HI", "DIS_NUM_ITER", "NUM_ITER", "PAIRS", "DMAX",
              "TMIN", "QE_BIN", "W90_BIN", "SRUN_OPTS", "KMESH_PRIMITIVE",
              "NBND_WINDOW", "NBND_MARGIN", "MESH_BASE", "W90_SRUN_OPTS", "DIS_CONV_TOL",
              "DIS_CONV_WINDOW", "CONV_TOL", "CONV_WINDOW")
    dflt = {}
    dpath = args.defaults
    if dpath is None:
        for cand in ("defaults.conf", os.path.join("..", "defaults.conf")):
            if os.path.exists(cand):
                dpath = cand
                break
    if dpath:
        if not os.path.exists(dpath):
            sys.exit(f"ERROR: {dpath} not found")
        for raw in open(dpath):
            line = raw.split("#", 1)[0].strip()
            if line and "=" in line:
                k, v = line.split("=", 1)
                k = k.strip()
                if k in POLICY:
                    dflt[k] = v.strip()
                else:
                    print(f"WARNING: {dpath} sets {k}, which is taken from "
                          f"scf.in and cannot be overridden here.",
                          file=sys.stderr)

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: {args.input} not found")
    text = qp.strip_comments(open(args.input).read())

    ctrl = qp.namelist(text, "CONTROL")
    syst = qp.namelist(text, "SYSTEM")
    prefix = ctrl.get("prefix")
    if not prefix:
        sys.exit("ERROR: no prefix in &CONTROL")
    outdir = ctrl.get("outdir", "./out")

    st = qp.read_structure(text)
    lat = st["lattice"]
    elements = []
    for _, el, _ in st["atoms"]:
        if el not in elements:
            elements.append(el)

    # ---- identify the transition metal --------------------------------
    tm = next((e for e in elements if e not in ("Gd", "Si")), None)
    if tm is None:
        sys.exit(f"ERROR: cannot identify a transition metal among "
                 f"{elements}. Expected something other than Gd and Si.")
    others = [e for e in elements if e not in ("Gd", "Si", tm)]
    if others:
        print(f"WARNING: extra species {others} will get no Wannier "
              f"projections unless you edit PROJECTIONS.", file=sys.stderr)

    counts = {e: sum(1 for a in st["atoms"] if a[1] == e) for e in elements}
    sysname = f"Gd{tm}Si"

    # ---- pseudopotentials, one per chemical element --------------------
    pp = {}
    for label, f in st["species"].items():
        pp.setdefault(qp.element(label), f)

    # ---- cell shape -> default k-mesh ----------------------------------
    lens = [sum(v * v for v in row) ** 0.5 for row in lat]
    mult = qp.supercell_multiplicity(st)
    kprim = args.kmesh_primitive or (
        [int(x) for x in dflt["KMESH_PRIMITIVE"].split()]
        if "KMESH_PRIMITIVE" in dflt else None)
    if args.nscf_mesh:
        mesh = list(args.nscf_mesh)
        mesh_why = "given"
    elif kprim:
        mesh = [max(1, k // m) for k, m in zip(kprim, mult)]
        mesh_why = (f"{'x'.join(str(k) for k in kprim)} on the primitive cell, "
                    f"divided by the {'x'.join(str(m) for m in mult)} "
                    f"supercell multiplicity")
    else:
        # keep the sampling density roughly isotropic in reciprocal space
        base = int(dflt.get("MESH_BASE", args.mesh_base))
        ref = min(lens)
        mesh = [max(2, 2 * -(-int(base * ref / L * 100) // 200)) for L in lens]
        mesh_why = f"base {base}, scaled from the cell shape"
    mesh_s = " ".join(str(n) for n in mesh)

    # ---- magnetisation --------------------------------------------------
    order = list(st["species"].keys())
    mag = {}
    for k, v in syst.items():
        m = re.match(r"starting_magnetization\((\d+)\)", k)
        if m:
            i = int(m.group(1)) - 1
            if i < len(order):
                mag[qp.element(order[i])] = qp.fnum(v)

    hub_u, hub_j = "6.7", "0.7"
    for l in qp.card(text, "HUBBARD")[1]:
        p = l.split()
        if len(p) >= 3 and p[1].lower().startswith("gd"):
            if p[0].upper() == "U":
                hub_u = p[2]
            elif p[0].upper().startswith("J"):
                hub_j = p[2]

    kunit, klines = qp.card(text, "K_POINTS")
    kscf = klines[0].split()[:3] if klines and kunit == "automatic" \
        else ["12", "12", "12"]

    # ---- distances -------------------------------------------------------
    d_tm_si = qp.min_distance(st, tm, "Si")
    d_tm_tm = qp.min_distance(st, tm, tm)

    # ---- projections and num_wann ---------------------------------------
    orb = {"Gd": "d", tm: "d", "Si": "p"}
    n_orb = {"s": 1, "p": 3, "d": 5, "f": 7}
    pairs = (dflt.get("PAIRS", "").replace("{TM}", tm)
             or f"{tm}-Si,{tm}-{tm},Gd-{tm}")
    projections = (args.projections
                   or dflt.get("PROJECTIONS", "").replace("{TM}", tm)
                   or ";".join(f"{e}:{orb[e]}" for e in ("Gd", tm, "Si")
                               if e in counts))
    num_wann = 0
    for spec in projections.split(";"):
        if ":" not in spec:
            continue
        e, o = spec.split(":", 1)
        num_wann += counts.get(e.strip(), 0) * n_orb.get(o.strip(), 0)
    nbnd = int(qp.fnum(syst.get("nbnd", "0")) or 0)

    # ---- scf.out ---------------------------------------------------------
    ef, converged, moments, tot_mag, nelec = None, None, [], None, None
    eig_blocks = []
    if os.path.exists(args.log):
        log = open(args.log, errors="ignore").read()
        converged = "convergence has been achieved" in log
        m = re.search(r"number of electrons\s*=\s*([0-9.]+)", log)
        if m:
            nelec = float(m.group(1))
        m = re.findall(r"the Fermi energy is\s+([-0-9.]+)\s*ev", log)
        if m:
            ef = float(m[-1])
        m = re.findall(r"total magnetization\s*=\s*([-0-9.]+)", log)
        if m:
            tot_mag = float(m[-1])
        eig_blocks = qp.read_eigenvalue_blocks(log)
        blocks = re.findall(r"Magnetic moment per site.*?\n((?:\s*atom.*\n)+)",
                            log)
        if blocks:
            for l in blocks[-1].splitlines():
                mm = re.search(r"atom\s+(\d+).*?magn=\s*([-0-9.]+)", l)
                if mm:
                    moments.append(float(mm.group(2)))

    # ---- bands for the nscf -------------------------------------------------
    # The .mmn overlap file scales as num_bands^2 * nkpts * nntot, so carrying
    # an SCF band count sized for a DOS calculation into the Wannier stage
    # costs tens of gigabytes and hours of pw2wannier90 for no benefit. Only
    # bands reaching a few eV above E_F are needed.
    win_ev = float(dflt.get("NBND_WINDOW", 10))
    margin = int(dflt.get("NBND_MARGIN", 10))
    if eig_blocks and ef is not None:
        hi = max(sum(1 for e in b if e < ef + win_ev) for _, b in eig_blocks)
        nbnd_nscf = hi + margin
        nbnd_why = (f"bands below E_F+{win_ev:g} eV ({hi}) plus {margin}")
    elif nelec:
        nbnd_nscf = max(num_wann + 20, int(nelec / 2 * 1.35) + 10)
        nbnd_why = "estimated from the electron count (no eigenvalues in the log)"
    else:
        nbnd_nscf = nbnd
        nbnd_why = "unchanged from the scf"
    if nbnd:
        nbnd_nscf = min(nbnd_nscf, nbnd)
    nbnd_nscf = max(nbnd_nscf, num_wann + 5)

    # ---- emit -------------------------------------------------------------
    def mag_of(e):
        return mag.get(e, 0.0)

    conf = f"""# generated by scf2conf.py from {args.input}
# {st['nat']} atoms: {', '.join(f'{n} {e}' for e, n in counts.items())}
# nearest {tm}-Si = {d_tm_si:.4f} A   nearest {tm}-{tm} = {d_tm_tm:.4f} A
#
# The geometry below is informational. While scf.in is present, the cell and
# the atomic positions used to build the .win files are read from it directly,
# so these keys do not have to describe the structure exactly.

SYSTEM={sysname}
PREFIX={prefix}
TM={tm}
NAT={st['nat']}
D_TM_SI={d_tm_si:.6f}

ALAT_BOHR={st['alat_bohr']:.8f}
COA={lens[2] / lens[0]:.8f}

PSEUDO_DIR={ctrl.get('pseudo_dir', './pseudo')}
PP_GD={pp.get('Gd', 'FIXME')}
PP_TM={pp.get(tm, 'FIXME')}
PP_SI={pp.get('Si', 'FIXME')}

ECUTWFC={qp.fnum(syst.get('ecutwfc', '60')):.0f}
ECUTRHO={qp.fnum(syst.get('ecutrho', '600')):.0f}
NBND={nbnd}
NBND_NSCF={nbnd_nscf}      # bands kept for the nscf and the Wannier stage
DEGAUSS={syst.get('degauss', '0.01')}
KMESH_SCF={' '.join(kscf)}
KMESH_NSCF={mesh_s}

START_MAG_GD={mag_of('Gd')}
START_MAG_TM={mag_of(tm)}
START_MAG_SI={mag_of('Si')}

HUBBARD_U_GD={hub_u}
HUBBARD_J_GD={hub_j}

MP_GRID={mesh_s}
PROJECTIONS={projections}
DIS_WIN_LO={dflt.get('DIS_WIN_LO', '6.0')}
DIS_WIN_HI={dflt.get('DIS_WIN_HI', '5.0')}
DIS_FROZ_LO={dflt.get('DIS_FROZ_LO', '5.0')}
DIS_FROZ_HI={dflt.get('DIS_FROZ_HI', '1.5')}
DIS_NUM_ITER={dflt.get('DIS_NUM_ITER', '1000')}
DIS_CONV_TOL={dflt.get('DIS_CONV_TOL', '1.0d-9')}
DIS_CONV_WINDOW={dflt.get('DIS_CONV_WINDOW', '3')}
NUM_ITER={dflt.get('NUM_ITER', '500')}

PAIRS={pairs}
DMAX={dflt.get('DMAX', '5.0')}
TMIN={dflt.get('TMIN', '0.02')}

QE_BIN={dflt.get('QE_BIN', '')}
W90_BIN={dflt.get('W90_BIN', '')}
SRUN_OPTS={dflt.get('SRUN_OPTS', '-n 16')}
W90_SRUN_OPTS={dflt.get('W90_SRUN_OPTS', '')}
"""
    if args.output == "-":
        sys.stdout.write(conf)
    else:
        open(args.output, "w").write(conf)
        if args.output != "/dev/null":
            print(f"wrote {args.output}")

    # ---- report ------------------------------------------------------------
    e = sys.stderr
    print(f"\n--- {args.input} ---", file=e)
    print(f"  prefix     {prefix}     outdir {outdir}", file=e)
    print(f"  cell       {lens[0]:.4f} x {lens[1]:.4f} x {lens[2]:.4f} A"
          f"   ({st['nat']} atoms)", file=e)
    labels = sorted({a[0] for a in st['atoms']})
    print(f"  species    {', '.join(labels)}"
          f"  ->  {', '.join(f'{n} {el}' for el, n in counts.items())}", file=e)
    if len(labels) > len(counts):
        print(f"  NOTE: magnetic sublattices are split into separate QE "
              f"species; they are merged into chemical elements for the "
              f"Wannier projections.", file=e)
    print(f"  nearest    {tm}-Si {d_tm_si:.4f} A   {tm}-{tm} {d_tm_tm:.4f} A",
          file=e)
    print(f"  projections {projections}  ->  num_wann = {num_wann}", file=e)
    if nbnd and num_wann > nbnd:
        print(f"  ERROR: num_wann {num_wann} exceeds nbnd {nbnd}", file=e)
    if nelec:
        print(f"  electrons  {nelec:.0f}", file=e)
    if nbnd_nscf != nbnd:
        saving = 1 - (nbnd_nscf / nbnd) ** 2
        print(f"  nbnd       {nbnd} in the scf -> {nbnd_nscf} for the nscf"
              f"   ({nbnd_why}; .mmn about {saving:.0%} smaller)", file=e)
    if dpath:
        print(f"  policy     {dpath}", file=e)
    nk = mesh[0] * mesh[1] * mesh[2]
    print(f"  k-mesh     nscf {mesh_s}  ({nk} points, {mesh_why})", file=e)
    # The mesh sets the real-space range of H(R): with N points along an axis
    # the Wigner-Seitz supercell reaches about N/2 cells, so this is the
    # quantity to match between calculations you intend to compare, not the
    # reciprocal-space density.
    cut = [n * L / 2 for n, L in zip(mesh, lens)]
    print(f"             H(R) reaches ~{cut[0]:.0f} x {cut[1]:.0f} x "
          f"{cut[2]:.0f} A. Match this across a series rather than the "
          f"mesh itself.", file=e)

    if converged is None:
        print(f"  {args.log} not found -- SCF state unknown", file=e)
    else:
        print(f"  SCF        {'converged' if converged else 'NOT CONVERGED'}",
              file=e)
        if ef is not None:
            print(f"  E_F        {ef:.4f} eV", file=e)
        if tot_mag is not None:
            print(f"  total mag  {tot_mag:.3f} uB/cell", file=e)
        if moments and len(moments) == st["nat"]:
            by_el = {}
            for (_, el, _), mv in zip(st["atoms"], moments):
                by_el.setdefault(el, []).append(mv)
            for el, vals in by_el.items():
                lo, hi = min(vals), max(vals)
                rng = f"{lo:+.3f}" if hi - lo < 5e-3 else f"{lo:+.3f}..{hi:+.3f}"
                afm = "  <- AFM" if lo < -0.05 < 0.05 < hi else ""
                print(f"    m({el:<3s}) {rng}{afm}", file=e)

    d = os.path.join(os.path.dirname(args.input) or ".", outdir.lstrip("./"))
    if not os.path.isdir(d):
        print(f"\n  WARNING: {outdir} not found. The nscf needs the charge "
              f"density from the SCF -- copy the whole directory.", file=e)


if __name__ == "__main__":
    main()
