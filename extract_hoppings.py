#!/usr/bin/env python3
"""
Extract effective hoppings t_eff(A,B,R) from a Wannier90 _hr.dat.

    t_eff(A,B,R) = sqrt( < |H_ij(R)|^2 > )   over i in block A, j in block B

Usage:
    extract_hoppings.py GdMnSi_up_hr.dat GdMnSi_up.win \
        --pairs Mn-Si,Mn-Mn --dmax 5.0 --tmin 0.02 > hoppings_MnSi_up.dat
"""
import argparse
import math
import os
import re
import sys
from collections import defaultdict

ORB_COUNT = {"s": 1, "p": 3, "d": 5, "f": 7,
             "sp": 2, "sp2": 3, "sp3": 4, "sp3d": 5, "sp3d2": 6}
BOHR = 0.529177210903


# ---------------------------------------------------------------- .win ----
def parse_win(path):
    text = open(path).read()

    def block(name):
        m = re.search(rf"begin\s+{name}(.*?)end\s+{name}", text,
                      re.S | re.I)
        if not m:
            return None
        return [l.strip() for l in m.group(1).splitlines() if l.strip()
                and not l.strip().startswith("!")]

    # lattice
    lines = block("unit_cell_cart")
    if lines is None:
        sys.exit(f"ERROR: no 'begin unit_cell_cart' in {path}")
    scale = 1.0
    if lines[0].lower() in ("ang", "angstrom"):
        lines = lines[1:]
    elif lines[0].lower() == "bohr":
        scale, lines = BOHR, lines[1:]
    lattice = [[float(x) * scale for x in l.split()[:3]] for l in lines[:3]]

    # atoms
    frac = block("atoms_frac")
    cart = block("atoms_cart")
    atoms = []                      # (element, fractional xyz)
    if frac:
        for l in frac:
            p = l.split()
            atoms.append((p[0], [float(x) for x in p[1:4]]))
    elif cart:
        inv = invert3(lattice)
        start = 1 if cart[0].lower() in ("ang", "angstrom", "bohr") else 0
        for l in cart[start:]:
            p = l.split()
            c = [float(x) for x in p[1:4]]
            atoms.append((p[0], matvec(inv, c)))
    else:
        sys.exit(f"ERROR: no atoms block in {path}")

    # projections -> ordered list of (label, n_orb)
    proj = block("projections")
    if proj is None:
        sys.exit(f"ERROR: no 'begin projections' in {path}")
    blocks = []
    for spec in proj:
        if spec.lower() in ("ang", "bohr", "random"):
            continue
        # The atom counter restarts for every projection line. An element may
        # legitimately appear more than once -- 'Gd:f' and 'Gd:d' are two
        # blocks on the same four atoms -- and a shared counter would invent
        # Gd5..Gd8 that no atom corresponds to.
        counter = defaultdict(int)
        elem, orb = spec.split(":", 1)
        elem = elem.strip()
        orb = orb.strip()
        if orb.startswith("l="):
            n_orb = 2 * int(re.match(r"l=(\d+)", orb).group(1)) + 1
        elif ";" in orb:
            n_orb = len(orb.split(";"))
        else:
            n_orb = ORB_COUNT.get(orb)
        if n_orb is None:
            sys.exit(f"ERROR: cannot count orbitals in '{spec}'")
        for el, _ in atoms:
            if el == elem:
                counter[elem] += 1
                blocks.append((f"{elem}{counter[elem]}", n_orb))
    return lattice, atoms, blocks


def matvec(m, v):
    return [sum(m[i][j] * v[j] for j in range(3)) for i in range(3)]


def invert3(m):
    a, b, c = m[0], m[1], m[2]
    det = (a[0] * (b[1] * c[2] - b[2] * c[1])
           - a[1] * (b[0] * c[2] - b[2] * c[0])
           + a[2] * (b[0] * c[1] - b[1] * c[0]))
    if abs(det) < 1e-12:
        sys.exit("ERROR: singular lattice matrix")
    co = [[(b[1] * c[2] - b[2] * c[1]), -(a[1] * c[2] - a[2] * c[1]),
           (a[1] * b[2] - a[2] * b[1])],
          [-(b[0] * c[2] - b[2] * c[0]), (a[0] * c[2] - a[2] * c[0]),
           -(a[0] * b[2] - a[2] * b[0])],
          [(b[0] * c[1] - b[1] * c[0]), -(a[0] * c[1] - a[1] * c[0]),
           (a[0] * b[1] - a[1] * b[0])]]
    return [[co[i][j] / det for j in range(3)] for i in range(3)]


def frac_to_cart(lattice, f):
    return [sum(f[i] * lattice[i][j] for i in range(3)) for j in range(3)]


# ------------------------------------------------------------ centres ----
def parse_centres(path, num_wann):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    n_tot = int(lines[0].split()[0])
    body = lines[2:2 + n_tot]
    centres = [[float(x) for x in l.split()[1:4]]
               for l in body if l.split()[0].upper() == "X"]
    if len(centres) != num_wann:
        return None
    return centres


def assign_by_centres(centres, atoms, lattice, labels):
    """Map each Wannier centre to its nearest atomic image."""
    cart = [frac_to_cart(lattice, f) for _, f in atoms]
    out, worst = [], 0.0
    for wc in centres:
        best, best_d = None, 1e30
        for ia, ac in enumerate(cart):
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    for k in (-1, 0, 1):
                        sh = frac_to_cart(lattice, [i, j, k])
                        d = math.dist(wc, [ac[t] + sh[t] for t in range(3)])
                        if d < best_d:
                            best_d, best = d, ia
        out.append(best)
        worst = max(worst, best_d)
    # atom index -> label like Mn1
    counter = defaultdict(int)
    names = []
    for el, _ in atoms:
        counter[el] += 1
        names.append(f"{el}{counter[el]}")
    return [names[i] for i in out], worst


# ------------------------------------------------------------- hr.dat ----
def read_hr(path):
    with open(path) as f:
        lines = f.readlines()
    num_wann = int(lines[1].split()[0])
    nrpts = int(lines[2].split()[0])
    n_deg = math.ceil(nrpts / 15)
    start = 3 + n_deg
    # tolerate a degeneracy block written with a different line width
    while start < len(lines) and len(lines[start].split()) != 7:
        start += 1
    return num_wann, nrpts, lines[start:]


# ---------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hr")
    ap.add_argument("win")
    ap.add_argument("--centres", help="default: <seed>_centres.xyz")
    ap.add_argument("--assign", choices=["centres", "projections"],
                    default="centres")
    ap.add_argument("--pairs", default="",
                    help="e.g. Mn-Si,Mn-Mn  (default: every pair)")
    ap.add_argument("--dmax", type=float, default=5.0)
    ap.add_argument("--tmin", type=float, default=0.02)
    ap.add_argument("--onsite", action="store_true",
                    help="also print on-site energies (for eps_d - eps_p)")
    args = ap.parse_args()

    lattice, atoms, proj_blocks = parse_win(args.win)
    num_wann, nrpts, data = read_hr(args.hr)

    n_proj = sum(n for _, n in proj_blocks)
    if n_proj != num_wann:
        sys.exit(f"ERROR: .win projections give {n_proj} Wannier functions but "
                 f"{args.hr} has {num_wann}. The two files are from different "
                 f"runs.")

    # --- Wannier function -> atom -----------------------------------------
    seed = re.sub(r"_hr\.dat$", "", os.path.basename(args.hr))
    cfile = args.centres or os.path.join(os.path.dirname(args.hr) or ".",
                                         f"{seed}_centres.xyz")
    wf_atom, method, worst = None, "projections", None

    if args.assign == "centres" and os.path.exists(cfile):
        centres = parse_centres(cfile, num_wann)
        if centres:
            wf_atom, worst = assign_by_centres(centres, atoms, lattice,
                                               proj_blocks)
            method = "centres"

    if wf_atom is None:
        wf_atom = []
        for label, n_orb in proj_blocks:
            wf_atom += [label] * n_orb
        if args.assign == "centres":
            print(f"# NOTE: {os.path.basename(cfile)} not usable, "
                  f"fell back to projection order", file=sys.stderr)

    # index (1-based) -> atom label
    idx2atom = {i + 1: wf_atom[i] for i in range(num_wann)}

    # --- diagnostics -------------------------------------------------------
    order = []
    for label, n_orb in proj_blocks:
        order += [label] * n_orb
    mism = [(i + 1, order[i], wf_atom[i]) for i in range(num_wann)
            if order[i] != wf_atom[i]]
    print(f"# assignment: {method}", file=sys.stderr)
    if worst is not None:
        print(f"# worst centre-to-atom distance: {worst:.3f} Ang",
              file=sys.stderr)
        if worst > 1.0:
            print("# WARNING: a Wannier centre sits far from every atom. "
                  "Check the spreads in the .wout before trusting these "
                  "numbers.", file=sys.stderr)
    if mism:
        print(f"# {len(mism)} Wannier functions do NOT sit where the "
              f"projection order says:", file=sys.stderr)
        for i, exp, got in mism[:12]:
            print(f"#   WF {i:3d}: expected {exp:>5s}, found on {got}",
                  file=sys.stderr)

    # --- atom positions and pair filter ------------------------------------
    counter = defaultdict(int)
    pos = {}
    for el, f in atoms:
        counter[el] += 1
        pos[f"{el}{counter[el]}"] = f

    want = None
    if args.pairs.strip():
        want = set()
        for p in args.pairs.split(","):
            a, b = p.strip().split("-")
            want.add(tuple(sorted((a, b))))

    def elem(lbl):
        return re.sub(r"\d+$", "", lbl)

    # --- accumulate ---------------------------------------------------------
    # H(R) is Hermitian: H_ij(R) = H_ji(-R)*, so (A,B,R) and (B,A,-R) are the
    # same bond written twice. Keeping both would double every count and, for
    # same-species pairs like Mn-Mn, list each bond twice in the final table.
    # Only the canonical half is accumulated; the other half is used as a
    # consistency check on the file.
    def canon(A, B, R):
        return min((A, B, R), (B, A, (-R[0], -R[1], -R[2])))

    hop = defaultdict(list)
    mirror = defaultdict(list)
    onsite = defaultdict(list)
    for line in data:
        p = line.split()
        if len(p) < 7:
            continue
        R = (int(p[0]), int(p[1]), int(p[2]))
        i, j = int(p[3]), int(p[4])
        val2 = float(p[5]) ** 2 + float(p[6]) ** 2
        A, B = idx2atom[i], idx2atom[j]
        if R == (0, 0, 0) and i == j:
            onsite[A].append(float(p[5]))
            continue
        if R == (0, 0, 0) and A == B:
            continue                      # intra-atomic crystal field
        if want and tuple(sorted((elem(A), elem(B)))) not in want:
            continue
        key = canon(A, B, R)
        (hop if (A, B, R) == key else mirror)[key].append(val2)

    worst_herm = 0.0
    for key, vals in hop.items():
        mv = mirror.get(key)
        if not mv or len(mv) != len(vals):
            continue
        t1 = math.sqrt(sum(vals) / len(vals))
        t2 = math.sqrt(sum(mv) / len(mv))
        if t1 > 1e-6:
            worst_herm = max(worst_herm, abs(t1 - t2) / t1)
    if worst_herm > 1e-3:
        print(f"# WARNING: H(R) is not Hermitian to better than "
              f"{worst_herm:.1%}. The _hr.dat may be truncated or the run "
              f"did not converge.", file=sys.stderr)

    # --- distances and output ------------------------------------------------
    rows = []
    for (A, B, R), vals in hop.items():
        df = [pos[B][t] + R[t] - pos[A][t] for t in range(3)]
        dc = frac_to_cart(lattice, df)
        d = math.sqrt(sum(x * x for x in dc))
        if d > args.dmax:
            continue
        t = math.sqrt(sum(vals) / len(vals))
        if t < args.tmin:
            continue
        rows.append((d, -t, A, B, R, t, len(vals)))

    print(f"# file: {args.hr}")
    print(f"# win:  {args.win}   assignment: {method}")
    print(f"# dmax = {args.dmax} Ang   tmin = {args.tmin} eV")
    print("# A B R1 R2 R3 distance_Ang t_eff_eV n_matrix_elements")
    for d, _, A, B, R, t, n in sorted(rows):
        print(f"{A:4s} {B:4s} {R[0]:3d} {R[1]:3d} {R[2]:3d} "
              f"{d:9.4f} {t:10.5f} {n:4d}")

    if args.onsite:
        print("\n# on-site energies (eV), mean over the orbitals of each atom")
        print("# atom n_orb eps_mean eps_min eps_max")
        by_elem = defaultdict(list)
        for A, vals in sorted(onsite.items()):
            m = sum(vals) / len(vals)
            print(f"# {A:5s} {len(vals):3d} {m:10.4f} "
                  f"{min(vals):10.4f} {max(vals):10.4f}")
            by_elem[elem(A)].append(m)
        print("# ---- charge-transfer energies ----")
        keys = sorted(by_elem)
        for a in keys:
            for b in keys:
                if a < b:
                    ea = sum(by_elem[a]) / len(by_elem[a])
                    eb = sum(by_elem[b]) / len(by_elem[b])
                    print(f"# eps({a}) - eps({b}) = {ea - eb:8.4f} eV")

    if not rows:
        print("# (nothing passed the filters -- loosen --tmin/--dmax "
              "or check --pairs)", file=sys.stderr)


if __name__ == "__main__":
    main()
