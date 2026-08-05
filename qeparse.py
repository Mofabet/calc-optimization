#!/usr/bin/env python3
"""
Shared parsing of Quantum ESPRESSO input files.

Used by scf2conf.py and make_inputs.py so that the structure is read from
scf.in in exactly one place. Nothing here assumes a particular cell, atom
count or composition.
"""
import re
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
