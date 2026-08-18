#!/usr/bin/env python3
"""
Shared library: Quantum ESPRESSO and Wannier90 file readers.

Every script imports from here so that a file format is parsed in exactly one
place. Nothing here assumes a particular cell, atom count or composition.

  read_structure          scf.in    -> cell, atoms, species
  read_eigenvalue_blocks  scf.out   -> eigenvalues per k, per spin
  supercell_multiplicity  structure -> repeats along each axis
  parse_win               .win      -> lattice, atoms, projection blocks
  parse_centres           _centres.xyz
  assign_by_centres       centres   -> nearest atom label
  read_hr                 _hr.dat   -> matrix element lines
  read_conf / need        KEY=value config files
"""
import math
import re
from collections import defaultdict
import sys

BOHR = 0.529177210903

CARDS = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS",
         "HUBBARD", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_VELOCITIES",
         "ATOMIC_FORCES", "ADDITIONAL_K_POINTS", "SOLVENTS")


def strip_comments(text):
    return re.sub(r"!.*", "", text)


def namelist(text, name):
    m = re.search(rf"&{name}\b(.*?)^\s*/", text, re.S | re.I | re.M)
    if not m:
        return {}
    out = {}
    for k, v in re.findall(r"([A-Za-z_]+(?:\(\d+\))?)\s*=\s*([^,\n]+)",
                           m.group(1)):
        out[k.strip().lower()] = v.strip().strip("',\" ")
    return out


def card(text, name):
    """-> (unit, [body lines]). Stops at the next card keyword."""
    lines = text.splitlines()
    unit, body, grabbing = None, [], False
    for line in lines:
        m = re.match(r"\s*([A-Z_]+)\b(.*)$", line)
        if m and m.group(1) in CARDS:
            if grabbing:
                break
            if m.group(1) == name:
                unit = m.group(2).strip().strip("(){} ").lower() or None
                grabbing = True
            continue
        if grabbing and line.strip():
            body.append(line.strip())
    return unit, body


def fnum(s):
    return float(s.replace("d", "e").replace("D", "E"))


def element(label):
    """
    Chemical element behind a QE species label.

    QE distinguishes magnetic sublattices by inventing species: Mn_up and
    Mn_dn are both manganese with the same pseudopotential. Wannier90
    projections act on the chemical species, and the spin channel is selected
    by pw2wannier90, so the labels must be collapsed before building the .win
    or half the Mn atoms would be left without projections.
    """
    lab = label.split("_")[0]
    lab = re.sub(r"[0-9]+$", "", lab)
    return lab[0].upper() + lab[1:].lower() if lab else lab


def invert3(m):
    a, b, c = m
    det = (a[0] * (b[1] * c[2] - b[2] * c[1])
           - a[1] * (b[0] * c[2] - b[2] * c[0])
           + a[2] * (b[0] * c[1] - b[1] * c[0]))
    if abs(det) < 1e-12:
        sys.exit("ERROR: singular cell matrix")
    co = [[(b[1] * c[2] - b[2] * c[1]), -(a[1] * c[2] - a[2] * c[1]),
           (a[1] * b[2] - a[2] * b[1])],
          [-(b[0] * c[2] - b[2] * c[0]), (a[0] * c[2] - a[2] * c[0]),
           -(a[0] * b[2] - a[2] * b[0])],
          [(b[0] * c[1] - b[1] * c[0]), -(a[0] * c[1] - a[1] * c[0]),
           (a[0] * b[1] - a[1] * b[0])]]
    return [[co[i][j] / det for j in range(3)] for i in range(3)]


def frac_to_cart(lattice, f):
    return [sum(f[i] * lattice[i][j] for i in range(3)) for j in range(3)]


def cart_to_frac(lattice, c):
    inv = invert3(lattice)
    return [sum(c[k] * inv[k][j] for k in range(3)) for j in range(3)]


def read_structure(text):
    """
    -> dict(alat_bohr, lattice (3x3 Angstrom), atoms [(label, element, frac)],
            species {label: pseudo}, nat)
    """
    syst = namelist(text, "SYSTEM")

    if "celldm(1)" in syst:
        alat_b = fnum(syst["celldm(1)"])
    elif "a" in syst:
        alat_b = fnum(syst["a"]) / BOHR
    else:
        sys.exit("ERROR: neither celldm(1) nor A in &SYSTEM")
    a_ang = alat_b * BOHR

    ibrav = int(fnum(syst.get("ibrav", "0")))
    if ibrav == 0:
        unit, lines = card(text, "CELL_PARAMETERS")
        if len(lines) < 3:
            sys.exit("ERROR: ibrav=0 but no usable CELL_PARAMETERS card")
        vec = [[fnum(x) for x in l.split()[:3]] for l in lines[:3]]
        scale = {"angstrom": 1.0, "ang": 1.0, "bohr": BOHR}.get(unit, a_ang)
        lattice = [[v * scale for v in row] for row in vec]
    elif ibrav == 6:
        coa = fnum(syst.get("celldm(3)", "0"))
        if coa == 0:
            sys.exit("ERROR: ibrav=6 but celldm(3) missing")
        lattice = [[a_ang, 0, 0], [0, a_ang, 0], [0, 0, a_ang * coa]]
    else:
        sys.exit(f"ERROR: ibrav={ibrav} is not supported (use 0 or 6)")

    _, sp_lines = card(text, "ATOMIC_SPECIES")
    species = {}
    for l in sp_lines:
        p = l.split()
        if len(p) >= 3:
            species[p[0]] = p[2]
    if not species:
        sys.exit("ERROR: ATOMIC_SPECIES card not found")

    unit, pos_lines = card(text, "ATOMIC_POSITIONS")
    if not pos_lines:
        sys.exit("ERROR: ATOMIC_POSITIONS card not found")
    atoms = []
    for l in pos_lines:
        p = l.split()
        if len(p) < 4:
            continue
        label = p[0]
        v = [fnum(x) for x in p[1:4]]
        if unit in ("crystal", "crystal_sg"):
            f = v
        else:
            mult = {"angstrom": 1.0, "bohr": BOHR}.get(unit, a_ang)
            f = cart_to_frac(lattice, [x * mult for x in v])
        f = [x - round(x) + (1.0 if x - round(x) < 0 else 0.0) for x in f]
        f = [x % 1.0 for x in f]
        atoms.append((label, element(label), f))

    nat = int(fnum(syst.get("nat", str(len(atoms)))))
    if nat != len(atoms):
        sys.exit(f"ERROR: nat={nat} but {len(atoms)} ATOMIC_POSITIONS lines")

    return dict(alat_bohr=alat_b, lattice=lattice, atoms=atoms,
                species=species, nat=nat)


def min_distance(struct, el_a, el_b):
    """Shortest distance between two chemical species, over image shifts."""
    lat = struct["lattice"]
    A = [a for a in struct["atoms"] if a[1] == el_a]
    B = [a for a in struct["atoms"] if a[1] == el_b]
    best = None
    for _, _, fa in A:
        for _, _, fb in B:
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    for k in (-1, 0, 1):
                        df = [fb[0] + i - fa[0], fb[1] + j - fa[1],
                              fb[2] + k - fa[2]]
                        dc = frac_to_cart(lat, df)
                        d = sum(x * x for x in dc) ** 0.5
                        if d > 1e-6 and (best is None or d < best):
                            best = d
    return best


def read_eigenvalue_blocks(text):
    """
    -> [(spin, [energies])], one entry per k-point per spin channel.

    Line-based rather than regex-based: QE prints eigenvalues eight to a line,
    splits them by spin channel, and with verbosity='high' follows each block
    with occupation numbers that must not be counted.
    Requires verbosity='high'; returns [] otherwise.
    """
    blocks, spin, cur, collecting = [], "unpolarised", None, False
    for line in text.splitlines():
        m = re.search(r"-+\s*SPIN\s+(UP|DOWN)", line, re.I)
        if m:
            spin = m.group(1).lower()
        if "bands (ev)" in line:
            if cur:
                blocks.append((spin, cur))
            cur, collecting = [], True
            continue
        if not collecting:
            continue
        s = line.strip()
        if not s:
            continue
        try:
            cur.extend(float(t) for t in s.split())
        except ValueError:
            blocks.append((spin, cur))
            cur, collecting = None, False
    if cur:
        blocks.append((spin, cur))
    return blocks


def supercell_multiplicity(struct, tol=1e-4):
    """
    -> (n1, n2, n3): how many primitive repeats fit along each axis.

    Found by testing whether translating every atom by 1/n along an axis maps
    the structure onto itself, comparing chemical elements. A magnetic cell
    built by splitting one sublattice into Mn_up and Mn_dn is a genuine
    structural doubling, and the Brillouin zone shrinks accordingly, so the
    k-mesh along that axis needs proportionally fewer points to reach the same
    real-space cutoff for H(R).
    """
    atoms = struct["atoms"]
    out = []
    for axis in range(3):
        found = 1
        for n in (4, 3, 2):
            shift = 1.0 / n
            ok = True
            for _, el, f in atoms:
                target = list(f)
                target[axis] = (target[axis] + shift) % 1.0
                if not any(
                        el2 == el and all(
                            min(abs(target[k] - f2[k]),
                                1 - abs(target[k] - f2[k])) < tol
                            for k in range(3))
                        for _, el2, f2 in atoms):
                    ok = False
                    break
            if ok:
                found = n
                break
        out.append(found)
    return tuple(out)


def read_kpoints_from_out(text):
    """
    -> [(kx, ky, kz)] in cartesian units of 2*pi/alat, in output order.

    QE prints the k-point above each eigenvalue block. This is the only place
    the coordinates appear in the output, so it is what a comparison against an
    independently generated path has to be matched on.
    """
    out = []
    for m in re.finditer(r"^\s*k\s*=\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)"
                         r".*?bands \(ev\)", text, re.M):
        out.append(tuple(float(m.group(i)) for i in (1, 2, 3)))
    return out


def cart2pi_to_crystal(k, lattice, alat_ang):
    """
    Cartesian k in units of 2*pi/alat -> crystal (fractional) coordinates.

    Uses b_i . a_j = 2*pi*delta_ij: the crystal component along b_i is the dot
    product of the cartesian k with the real lattice vector a_i, both expressed
    in units of alat.
    """
    return [sum(k[j] * lattice[i][j] / alat_ang for j in range(3))
            for i in range(3)]


# ---- Wannier90 files ----

ORB_COUNT = {"s": 1, "p": 3, "d": 5, "f": 7,
             "sp": 2, "sp2": 3, "sp3": 4, "sp3d": 5, "sp3d2": 6}


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
            atoms.append((p[0], cart_to_frac(lattice, c)))
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


# ---- config files ----

def read_conf(path):
    """KEY=value, '#' comments, whitespace-tolerant."""
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
