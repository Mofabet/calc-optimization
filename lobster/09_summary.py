#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_summary.py -- quality control and post-processing summary for a LOBSTER series.

Walks a set of calculation folders, extracts everything needed to judge the
projection quality and to fill the thesis tables, prints a readable report and
writes a flat CSV.

Reads (anything missing is skipped silently):
    lobsterout              charge spilling, electrons recovered by projection
    lobsterin               basis set, keywords actually requested
    ICOHPLIST.lobster       -ICOHP per bond, per spin channel
    COHPCAR.lobster         -COHP at E_F, occupied antibonding share
    pDOSCAR_*.loposter      size of the negative excursions in the projected DOS
    scf.out                 E_F, magnetic moments, total energy
    scf.in                  atom order, so Mn1 / Mn2 sublattices can be told apart

Examples:
    ./09_summary.py z*  --csv lobster_qc.csv --out lobster_qc.txt
    ./09_summary.py                          # folders taken from lobster.conf
    ./09_summary.py z12 z2074 --species      # split Mn1-Mn1 from Mn1-Mn2
    ./09_summary.py z* --no-cohp             # fast: ICOHP and moments only
"""

import argparse
import glob
import math
import os
import re
import sys
from collections import defaultdict

# ----------------------------------------------------------------------------
# Warning thresholds -- edit here
# ----------------------------------------------------------------------------
SPILLING_WARN = 3.0          # %, abs. charge spilling above this is suspect
SPILLING_BAD = 5.0           # %, above this nothing downstream is reliable
ELECTRON_TOL = 0.05          # electrons missing from the projection
PDOS_NEG_WARN = 5.0          # %, negative pDOS excursion relative to the peak
FIXED_TOL = 0.01             # A, spread below which a distance counts as fixed
CTRL_SPREAD_WARN = 0.15      # -ICOHP spread at a fixed distance, fraction
QUENCHED_MOMENT = 0.5        # muB, below this the moment counts as quenched
SLOW_SCF_WARN = 120          # iterations

PAIR_CUTOFF = 4.0            # A, pair radius
MAX_SHELLS = 3               # how many shells per pair to report
RY_TO_MEV = 13605.693


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
def fnum(x, nd=3, dash="-"):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return dash
    return f"{x:.{nd}f}"


def table(headers, rows, aligns=None, indent="  "):
    if not rows:
        return indent + "(no data)"
    n = len(headers)
    aligns = list(aligns or []) + [">"] * (n - len(aligns or []))
    cells = [[str(c) for c in r] + [""] * (n - len(r)) for r in rows]
    w = [max(len(headers[i]), max(len(r[i]) for r in cells)) for i in range(n)]
    out = [indent + "  ".join(f"{headers[i]:{aligns[i]}{w[i]}}" for i in range(n)),
           indent + "  ".join("-" * w[i] for i in range(n))]
    for r in cells:
        out.append(indent + "  ".join(f"{r[i]:{aligns[i]}{w[i]}}" for i in range(n)))
    return "\n".join(out)


def section(title):
    return "\n" + "=" * 78 + f"\n{title}\n" + "=" * 78 + "\n"


def loglinfit(x, y):
    """y = A exp(-b x); -> A, b, RMS of the residuals in the logarithm."""
    n = len(x)
    ly = [math.log(v) for v in y]
    mx = sum(x) / n
    my = sum(ly) / n
    den = sum((a - mx) ** 2 for a in x)
    if den == 0:
        return None, None, None
    b = sum((a - mx) * (c - my) for a, c in zip(x, ly)) / den
    a0 = my - b * mx
    rms = math.sqrt(sum((c - (a0 + b * a)) ** 2 for a, c in zip(x, ly)) / n)
    return math.exp(a0), -b, rms


def sort_key(path):
    m = re.match(r"^z(\d+)", os.path.basename(path.rstrip("/")))
    return (float("0." + m.group(1)) if m else 9.9, path)


def z_of(name):
    m = re.match(r"^z(\d+)", name)
    return float("0." + m.group(1)) if m else None


def config_of(name):
    m = re.search(r"_(FM|AFM|FIM)$", name, re.I)
    return m.group(1).upper() if m else "gs"


def pattern(vals, thr=0.05):
    return "".join("+" if v > thr else "-" if v < -thr else "0" for v in vals)


def tm_moments(scf):
    """Moments of the transition metal, whichever of Mn / Fe is present."""
    mom = (scf or {}).get("moments", {})
    for el in ("Mn", "Fe", "Co", "Ni", "Ru"):
        if mom.get(el):
            return mom[el]
    return []


# ----------------------------------------------------------------------------
# Atom labels.  LOBSTER writes element+index and loses the QE species names,
# so the sublattice split has to come back from ATOMIC_POSITIONS.
# ----------------------------------------------------------------------------
def split_label(lab):
    m = re.match(r"^([A-Za-z][A-Za-z]?\d*)_(\d+)$", lab)
    if m:
        return m.group(1), int(m.group(2))
    m = re.match(r"^([A-Za-z]+)(\d+)$", lab)
    if m:
        return m.group(1), int(m.group(2))
    return lab, None


def species_by_index(scf_in):
    if not os.path.isfile(scf_in):
        return {}
    txt = open(scf_in, errors="replace").read()
    m = re.search(r"ATOMIC_POSITIONS[^\n]*\n(.*?)(?=\n\s*[A-Z_]{4,}|\Z)", txt, re.S)
    out, i = {}, 0
    if m:
        for line in m.group(1).strip().splitlines():
            f = line.split()
            if len(f) >= 4:
                i += 1
                out[i] = f[0]
    return out


def make_atom_id(by_species, smap):
    def atom_id(lab):
        name, idx = split_label(lab)
        if by_species and idx is not None and idx in smap:
            return smap[idx]
        return re.sub(r"\d+$", "", name)
    return atom_id


def pair_of(a, b, atom_id):
    x, y = sorted([atom_id(a), atom_id(b)])
    return f"{x}-{y}"


# ----------------------------------------------------------------------------
# lobsterout / lobsterin / pDOS
# ----------------------------------------------------------------------------
def parse_lobsterout(path):
    if not os.path.isfile(path):
        return None
    d = {"spilling": [], "e_rec": None, "e_tot": None}
    for line in open(path, errors="replace"):
        if "spilling" in line.lower():
            d["spilling"] += [float(v) for v in re.findall(r"([-\d.]+)\s*%", line)]
        m = re.search(r"electrons recovered by projection:\s*([\d.]+)\s*of\s*([\d.]+)",
                      line)
        if m:
            d["e_rec"], d["e_tot"] = float(m.group(1)), float(m.group(2))
    return d


def parse_lobsterin(path):
    d = {"basis": None, "rsh": False, "gen_from": None}
    if not os.path.isfile(path):
        return d
    for line in open(path, errors="replace"):
        s = line.strip()
        m = re.match(r"basisSet\s+(\S+)", s, re.I)
        if m:
            d["basis"] = m.group(1)
        if re.match(r"realspaceHamiltonian", s, re.I):
            d["rsh"] = True
        m = re.match(r"cohpGenerator\s+from\s+([\d.]+)\s+to\s+([\d.]+)", s, re.I)
        if m:
            d["gen_from"] = (float(m.group(1)), float(m.group(2)))
    return d


def pdos_negativity(folder):
    """Deepest negative pDOS excursion, as a percentage of that file's peak."""
    worst = None
    for f in sorted(glob.glob(os.path.join(folder, "pDOSCAR_*.loposter"))):
        mn, mx, seen = 0.0, 0.0, False
        for i, line in enumerate(open(f, errors="replace")):
            if i < 6:                      # DOSCAR-style preamble
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            if re.fullmatch(r"\d+", fields[2]):   # repeated block header
                continue
            for v in fields[1:]:
                try:
                    x = float(v)
                except ValueError:
                    continue
                mn, mx, seen = min(mn, x), max(mx, x), True
        if seen and mx > 0:
            pct = -100.0 * mn / mx
            if worst is None or pct > worst[1]:
                worst = (os.path.basename(f), pct)
    return worst


# ----------------------------------------------------------------------------
# ICOHPLIST / COHPCAR
# ----------------------------------------------------------------------------
def parse_icohplist(path, atom_id, cutoff):
    """{(pair, d_round): {'d','n','up','dn','tot'}}, eV per bond.

    The spin channels are written as consecutive blocks with the bond index
    restarting, which is how they are told apart.
    """
    if not os.path.isfile(path):
        return None
    blocks, cur, last = [], [], None
    for line in open(path, errors="replace"):
        f = line.split()
        if len(f) < 8 or not f[0].rstrip(".").isdigit():
            continue
        try:
            i = int(f[0].rstrip("."))
            a, b, dist, val = f[1], f[2], float(f[3]), float(f[7])
        except (ValueError, IndexError):
            continue
        if last is not None and i <= last:
            blocks.append(cur)
            cur = []
        cur.append((a, b, dist, val))
        last = i
    if cur:
        blocks.append(cur)
    if not blocks:
        return None

    nblocks = min(len(blocks), 2)
    per = defaultdict(lambda: {"up": None, "dn": None, "d": 0.0, "n": 0,
                               "tot": 0.0, "split": "-"})
    for s, block in enumerate(blocks[:2]):
        chan = "up" if s == 0 else "dn"
        seen = defaultdict(list)
        for a, b, dist, val in block:
            if dist < 0.1 or dist > cutoff:        # on-site terms are not bonds
                continue
            seen[(pair_of(a, b, atom_id), round(dist / 0.05) * 0.05)].append(
                (dist, val))
        for key, vals in seen.items():
            per[key][chan] = -sum(v for _, v in vals) / len(vals)
            per[key]["d"] = sum(x for x, _ in vals) / len(vals)
            per[key]["n"] = len(vals)
    for e in per.values():
        if nblocks == 2:
            # two blocks: the file really is spin resolved
            e["tot"] = (e["up"] or 0.0) + (e["dn"] or 0.0)
            e["split"] = "file"
        else:
            # one block: the value is already summed over spin. The channels
            # are unknown here and get filled in from COHPCAR later.
            e["tot"] = e["up"] or 0.0
            e["up"] = e["dn"] = None
        e["nblocks"] = nblocks
    return dict(per)


def parse_cohpcar(path):
    """(energies, {i: {spin: [COHP]}}, labels, nspin) or None.

    The numeric block starts after the last "No.N:" label line -- the header
    lines above it also parse as numbers and are not data.
    """
    if not os.path.isfile(path):
        return None
    lines = open(path, errors="replace").read().splitlines()
    lab_re = re.compile(r"^No\.\d+:(\S+?)->(\S+?)\((-?[\d.]+)\)")
    labels, last = [], -1
    for i, line in enumerate(lines):
        m = lab_re.match(line.strip())
        if m:
            labels.append((m.group(1), m.group(2), float(m.group(3))))
            last = i
    if last < 0:
        return None
    rows = []
    for line in lines[last + 1:]:
        f = line.split()
        if len(f) < 4:
            continue
        try:
            rows.append([float(x) for x in f])
        except ValueError:
            continue
    if not rows:
        return None
    widths = defaultdict(int)
    for r in rows:
        widths[len(r)] += 1
    ncol = max(widths, key=widths.get)
    data = [r for r in rows if len(r) == ncol]

    layout = None
    for nspin in (1, 2):
        for extra in (1, 0):
            if ncol == 1 + 2 * (len(labels) + extra) * nspin:
                layout = (nspin, extra)
                break
        if layout:
            break
    if not layout:
        return None
    nspin, extra = layout
    energies = [r[0] for r in data]
    cohp = {}
    for i in range(len(labels)):
        cohp[i] = {s: [r[1 + s * 2 * (len(labels) + extra) + 2 * (extra + i)]
                       for r in data] for s in range(nspin)}
    return energies, cohp, labels, nspin


def cohp_at_fermi(cc, atom_id, cutoff):
    """{(pair, d): {'ef': [-COHP per spin], 'anti': % of occupied weight}}"""
    if not cc:
        return {}
    energies, cohp, labels, nspin = cc
    if len(energies) < 2:
        return {}
    groups = defaultdict(list)
    for i, (a, b, dist) in enumerate(labels):
        if dist < 0.1 or dist > cutoff:
            continue
        groups[(pair_of(a, b, atom_id), round(dist / 0.05) * 0.05)].append(i)
    dE = energies[1] - energies[0]
    j0 = min(range(len(energies)), key=lambda k: abs(energies[k]))
    out = {}
    for key, idx in groups.items():
        ef, occ, anti, bond = [], [], 0.0, 0.0
        for s in range(nspin):
            y = [-sum(cohp[i][s][j] for i in idx) / len(idx)
                 for j in range(len(energies))]
            ef.append(y[j0])
            occ.append(sum(v for e, v in zip(energies, y) if e <= 0) * dE)
            for e, v in zip(energies, y):
                if e <= 0:
                    if v < 0:
                        anti += -v * dE
                    else:
                        bond += v * dE
        out[key] = {"ef": ef, "occ": occ, "nspin": nspin,
                    "anti": 100 * anti / (anti + bond) if (anti + bond) else None}
    return out


# ----------------------------------------------------------------------------
# scf.out
# ----------------------------------------------------------------------------
def parse_scf_out(path, smap):
    if not os.path.isfile(path):
        return None
    d, moments = {}, {}
    for line in open(path, errors="replace"):
        m = re.search(r"^!\s+total energy\s*=\s*(-?[\d.]+)", line)
        if m:
            d["E_Ry"] = float(m.group(1))
        m = re.search(r"the Fermi energy is\s+(-?[\d.]+)", line)
        if m:
            d["Ef"] = float(m.group(1))
        m = re.search(r"total magnetization\s*=\s*(-?[\d.]+)", line)
        if m:
            d["M_tot"] = float(m.group(1))
        m = re.search(r"absolute magnetization\s*=\s*(-?[\d.]+)", line)
        if m:
            d["M_abs"] = float(m.group(1))
        m = re.search(r"convergence has been achieved in\s+(\d+)", line)
        if m:
            d["n_iter"] = int(m.group(1))
        m = re.search(r"^\s*atom\s+(\d+).*magn=\s*(-?[\d.]+)", line)
        if m:
            moments[int(m.group(1))] = float(m.group(2))
    by_el = defaultdict(list)
    for i, mu in sorted(moments.items()):
        by_el[re.sub(r"\d+$", "", smap.get(i, "?"))].append(mu)
    d["moments"] = dict(by_el)
    return d


# ----------------------------------------------------------------------------
# Collection
# ----------------------------------------------------------------------------
def collect(folder, args):
    name = os.path.basename(folder.rstrip("/"))
    r = {"name": name, "path": folder, "notes": []}
    smap = species_by_index(os.path.join(folder, "scf.in"))
    atom_id = make_atom_id(args.species, smap)
    if args.species and not smap:
        r["notes"].append("--species asked for but scf.in has no "
                          "ATOMIC_POSITIONS; fell back to elements")

    r["lob"] = parse_lobsterout(os.path.join(folder, "lobsterout"))
    if r["lob"] is None:
        r["notes"].append("no lobsterout")
    r["lin"] = parse_lobsterin(os.path.join(folder, "lobsterin"))

    r["scf"] = parse_scf_out(os.path.join(folder, "scf.out"), smap)
    if r["scf"] is None:
        r["notes"].append("no scf.out")

    ic = parse_icohplist(os.path.join(folder, "ICOHPLIST.lobster"),
                         atom_id, args.cutoff)
    if ic is None:
        r["notes"].append("no ICOHPLIST.lobster")
    r["icohp"] = ic or {}

    if args.no_cohp:
        r["cohp"] = {}
    else:
        cc = parse_cohpcar(os.path.join(folder, "COHPCAR.lobster"))
        if cc is None:
            r["notes"].append("COHPCAR.lobster missing or unreadable")
        r["cohp"] = cohp_at_fermi(cc, atom_id, args.cutoff)

    fill_spin_from_cohp(r)
    r["pdos"] = pdos_negativity(folder)
    return r


def fill_spin_from_cohp(r):
    """Split -ICOHP into channels when ICOHPLIST holds only the sum.

    Some LOBSTER versions write ICOHPLIST already summed over spin. The exact
    total then comes from the list, and the ratio between the channels from
    integrating -COHP over the occupied states -- which is truncated at
    COHPstartEnergy, so only the ratio is used, never the absolute value.
    """
    for key, e in r["icohp"].items():
        if e.get("up") is not None:
            continue
        c = r["cohp"].get(key)
        occ = (c or {}).get("occ") or []
        if len(occ) < 2 or sum(occ) <= 0:
            e["split"] = "-"
            continue
        f_up = occ[0] / sum(occ)
        e["up"] = e["tot"] * f_up
        e["dn"] = e["tot"] * (1 - f_up)
        e["split"] = "COHP"


def shell_index(results):
    """Number the shells of each pair WITHIN a folder, then describe them.

    Shell 1 is the nearest neighbour of that pair in that folder. Numbering
    per folder rather than globally is what lets a bond whose length varies
    across the series stay "shell 1" everywhere.

    Returns (folder, pair, d_round) -> shell, and (pair, shell) -> geometry,
    where the geometry is taken over the ground-state folders only: a forced
    magnetic configuration sits at the same geometry as its parent and would
    otherwise be counted twice.
    """
    idx = {}
    for r in results:
        by_pair = defaultdict(set)
        for pair, d in r["icohp"]:
            by_pair[pair].add(round(d, 2))
        for pair, dists in by_pair.items():
            for sh, d in enumerate(sorted(dists), 1):
                idx[(r["name"], pair, d)] = sh

    meta = defaultdict(list)
    for r in results:
        if config_of(r["name"]) != "gs":
            continue
        for (pair, d), e in r["icohp"].items():
            meta[(pair, idx[(r["name"], pair, round(d, 2))])].append(e["d"])
    out = {}
    for key, ds in meta.items():
        out[key] = {"mean": sum(ds) / len(ds),
                    "spread": max(ds) - min(ds),
                    "fixed": (max(ds) - min(ds)) < FIXED_TOL,
                    "nfold": len(ds)}
    return idx, out


def shell_of(idx, folder, pair, d):
    return idx.get((folder, pair, round(d, 2)))


def series_values(results, idx, pair, sh, gs_only=True):
    """[(folder, entry)] for one pair/shell across the series."""
    out = []
    for r in results:
        if gs_only and config_of(r["name"]) != "gs":
            continue
        for (p, d), e in r["icohp"].items():
            if p == pair and shell_of(idx, r["name"], p, d) == sh:
                out.append((r["name"], e))
    return out


def reference_energy(results):
    ref = None
    for r in results:
        if config_of(r["name"]) == "gs" and r["scf"] and "E_Ry" in r["scf"]:
            e = r["scf"]["E_Ry"]
            ref = e if ref is None else min(ref, e)
    return ref


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def render(results, args):
    idx, meta = shell_index(results)
    L = ["LOBSTER SERIES QUALITY REPORT",
         f"folders: {len(results)}   pairs up to {args.cutoff:.1f} A   "
         f"shells per pair: {args.shells}   "
         f"grouping: {'species (Mn1/Mn2 split)' if args.species else 'element'}"]

    # 1 --------------------------------------------------------------------
    L.append(section("1. Projection quality"))
    rows = []
    for r in results:
        lo, lin = r["lob"] or {}, r["lin"]
        sp = lo.get("spilling") or []
        miss = (lo["e_tot"] - lo["e_rec"]
                if lo.get("e_rec") is not None and lo.get("e_tot") else None)
        rows.append([r["name"], lin.get("basis") or "-",
                     fnum(sp[0], 2) if len(sp) > 0 else "-",
                     fnum(sp[1], 2) if len(sp) > 1 else "-",
                     fnum(lo.get("e_rec"), 4), fnum(miss, 4),
                     fnum(r["pdos"][1], 2) if r["pdos"] else "-",
                     "yes" if lin.get("rsh") else "no",
                     len(r["icohp"])])
    L.append(table(["folder", "basis", "spill up", "spill dn", "e recov",
                    "e lost", "pDOS neg %", "H(R)", "pairs"], rows, ["<", "<"]))
    L.append("""
  spill        abs. charge spilling, %. Under 3 is good; over 5-8 the local
               basis does not reproduce the plane-wave wavefunction and
               nothing further down is reliable. It measures relocation of
               charge, not its absence
  e recov      electrons recovered by projection; e lost is the shortfall
  pDOS neg %   deepest negative excursion of any projected DOS file, as a
               percentage of that file's own peak. A few percent is the
               normal price of a Mulliken partition, not a defect
  H(R)         realspaceHamiltonian requested, so the tight-binding matrix
               is on disk next to the COHP data""")

    # 2 --------------------------------------------------------------------
    L.append(section("2. Magnetic state and energies"))
    ref = reference_energy(results)
    rows = []
    for r in results:
        s = r["scf"] or {}
        mn, gd = tm_moments(s), s.get("moments", {}).get("Gd") or []
        dE = ((s["E_Ry"] - ref) * RY_TO_MEV
              if ref is not None and "E_Ry" in s else None)
        rows.append([r["name"], fnum(z_of(r["name"]), 4), config_of(r["name"]),
                     fnum(s.get("E_Ry"), 6), fnum(dE, 1),
                     fnum(sum(abs(x) for x in mn) / len(mn), 3) if mn else "-",
                     pattern(mn) or "-",
                     fnum(sum(gd) / len(gd), 3) if gd else "-",
                     fnum(s.get("M_tot"), 2), fnum(s.get("M_abs"), 2),
                     fnum(s.get("Ef"), 3), s.get("n_iter", "-")])
    L.append(table(["folder", "z", "conf", "E, Ry", "dE, meV", "|m| TM",
                    "pattern", "m Gd", "M tot", "M abs", "E_F, eV", "iter"],
                   rows, ["<", ">", "<"]))
    L.append("""
  dE           relative to the lowest ground-state folder. Folders whose name
               ends in _FM / _AFM are forced configurations and are left out
               of the choice of that reference
  pattern      sign of the moment on each TM atom, in ATOMIC_POSITIONS order
  M tot/M abs  cell magnetisation. M abs rising while M tot falls is the
               signature of compensated moments forming""")

    # 3 --------------------------------------------------------------------
    L.append(section("3. -ICOHP by pair and shell (eV per bond)"))
    rows = []
    for r in results:
        for (pair, d), e in sorted(r["icohp"].items(),
                                   key=lambda kv: (kv[0][0], kv[0][1])):
            sh = shell_of(idx, r["name"], pair, d)
            if sh is None or sh > args.shells:
                continue
            c = r["cohp"].get((pair, d), {})
            ef = c.get("ef") or []
            rows.append([r["name"], pair, sh, fnum(e["d"], 3), e["n"],
                         fnum(e["tot"], 3),
                         fnum(e.get("up"), 3), fnum(e.get("dn"), 3),
                         e.get("split", "-"),
                         fnum(ef[0], 3) if len(ef) > 0 else "-",
                         fnum(ef[1], 3) if len(ef) > 1 else "-",
                         fnum(c.get("anti"), 1)])
    L.append(table(["folder", "pair", "sh", "d, A", "N", "total", "up", "dn",
                    "src", "-COHP@EF up", "dn", "anti %"], rows,
                   ["<", "<", ">", ">", ">", ">", ">", "<"]))
    L.append("""
  total        -ICOHP integrated to E_F, per bond, summed over spin. Positive
               means bonding. This is the share of the band energy belonging to
               one pair, not a dissociation energy: compare bonds with each
               other, never add them up into an energy difference
  up/dn, src   the channels. src=file means ICOHPLIST was written spin
               resolved; src=COHP means the file held only the sum and the
               ratio between the channels was taken from COHPCAR, with the
               exact total kept. src=- means no split was available
  -COHP@EF     the value at the Fermi level. Negative means the states at E_F
               are antibonding for that pair, i.e. E_F sits above the
               bonding/antibonding crossover
  anti %       share of the occupied weight that is antibonding
  Equal up and dn on a TM-TM pair says the two atoms are antiparallel:
  swapping the spin index then maps the bond onto itself. Compare the channels
  within one row, not across folders -- E_F moves between them.""")

    # 4 --------------------------------------------------------------------
    L.append(section("4. Internal control: bonds at a fixed distance"))
    rows = []
    for (pair, sh), m in sorted(meta.items()):
        if not m["fixed"] or m["nfold"] < 3 or sh > args.shells:
            continue
        vals = [e["tot"] for _, e in series_values(results, idx, pair, sh)]
        if len(vals) < 3:
            continue
        mean = sum(vals) / len(vals)
        rel = (max(vals) - min(vals)) / mean if mean else float("nan")
        rows.append([pair, sh, fnum(m["mean"], 3), len(vals), fnum(min(vals), 3),
                     fnum(max(vals), 3), fnum(100 * rel, 1)])
    L.append(table(["pair", "sh", "d, A", "points", "min", "max", "spread %"],
                   rows, ["<"]))
    L.append("""
  Ground-state folders only. These distances do not change across the
  series, so whatever -ICOHP does here is purely electronic. A small spread means the bond is a spectator and
  cannot be what drives any change in the series; a large one means it
  responds to something other than its own geometry.""")

    # 5 --------------------------------------------------------------------
    L.append(section("5. Bond-length dependence"))
    rows = []
    for (pair, sh), m in sorted(meta.items()):
        if m["fixed"] or sh > 1:
            continue
        pts = sorted((e["d"], e["tot"])
                     for _, e in series_values(results, idx, pair, sh)
                     if e["tot"] > 0)
        if len(pts) < 3:
            continue
        A, kappa, rms = loglinfit([x for x, _ in pts], [y for _, y in pts])
        if kappa is None:
            continue
        rows.append([pair, len(pts), f"{pts[0][0]:.3f}-{pts[-1][0]:.3f}",
                     fnum(A, 1), fnum(kappa, 2), fnum(rms, 4),
                     fnum(100 * (pts[-1][1] / pts[0][1] - 1), 1)])
    L.append(table(["pair", "points", "d range, A", "A", "kappa, 1/A",
                    "RMS log", "change %"], rows, ["<"]))
    L.append("""
  Ground-state folders only. Fit -ICOHP = A exp(-kappa d) over the folders
  where this distance varies.
  kappa around 1.5-2 1/A is what d-p overlap normally gives. RMS log is the
  scatter in the logarithm: below 0.01 the points sit on the curve.""")

    # 6 --------------------------------------------------------------------
    L.append(section("6. Moment against hybridisation"))
    rows = []
    for r in results:
        mn = tm_moments(r["scf"])
        best = None
        for (pair, d), e in r["icohp"].items():
            a, b = pair.split("-")
            if (a.startswith(("Mn", "Fe")) and b.startswith("Si")) or \
               (b.startswith(("Mn", "Fe")) and a.startswith("Si")):
                if best is None or e["d"] < best[0]:
                    best = (e["d"], e["tot"], pair)
        if not (mn and best):
            continue
        avg = sum(abs(x) for x in mn) / len(mn)
        rows.append([r["name"], fnum(z_of(r["name"]), 4), best[2],
                     fnum(best[0], 3), fnum(best[1], 3), fnum(avg, 3),
                     "quenched" if avg < QUENCHED_MOMENT else pattern(mn)])
    L.append(table(["folder", "z", "pair", "d, A", "-ICOHP", "|m| TM", "state"],
                   rows, ["<", ">", "<"]))
    L.append("""
  The shortest TM-Si bond against the local moment. Strong hybridisation
  broadens the d band and suppresses the moment; as the bond weakens the band
  narrows and the moment appears. Read together with section 4: if the TM-TM
  distance is fixed there, direct overlap cannot be what drives the change.""")

    # 7 --------------------------------------------------------------------
    L.append(section("7. Warnings"))
    fl = flags(results, idx, meta, args)
    L.append("\n".join("  " + f for f in fl) if fl else "  Nothing triggered.")

    notes = [f"  {r['name']}: {n}" for r in results for n in r["notes"]]
    if notes:
        L.append(section("8. Could not be read"))
        L.append("\n".join(notes))
    return "\n".join(L)


def flags(results, idx, meta, args):
    out = []
    for r in results:
        lo = r["lob"] or {}
        for i, v in enumerate((lo.get("spilling") or [])[:2]):
            ch = "up" if i == 0 else "dn"
            if v >= SPILLING_BAD:
                out.append(f"{r['name']}: charge spilling {v:.2f}% ({ch}) -- "
                           f"above {SPILLING_BAD}%, results not usable")
            elif v >= SPILLING_WARN:
                out.append(f"{r['name']}: charge spilling {v:.2f}% ({ch}) -- "
                           f"above {SPILLING_WARN}%, check the basis")
        if lo.get("e_rec") is not None and lo.get("e_tot"):
            miss = lo["e_tot"] - lo["e_rec"]
            if abs(miss) > ELECTRON_TOL:
                out.append(f"{r['name']}: {miss:.3f} electrons missing from "
                           f"the projection")
        if r["pdos"] and r["pdos"][1] > PDOS_NEG_WARN:
            out.append(f"{r['name']}: pDOS negative excursion "
                       f"{r['pdos'][1]:.1f}% of peak in {r['pdos'][0]}")
        s = r["scf"] or {}
        if s.get("n_iter") and s["n_iter"] > SLOW_SCF_WARN:
            out.append(f"{r['name']}: SCF took {s['n_iter']} iterations -- "
                       f"check starting_magnetization and the final state")

    onsite = [r["name"] for r in results
              if r["lin"].get("gen_from") and r["lin"]["gen_from"][0] < 0.05]
    if onsite:
        out.append(f"cohpGenerator starts at 0 A in {len(onsite)} folder(s), "
                   f"so on-site terms are in the files -- this report excludes "
                   f"them, but 08_plot_cohp.py needs its own switch")

    derived = [r["name"] for r in results
               if any(e.get("split") == "COHP" for e in r["icohp"].values())]
    if derived:
        out.append(f"ICOHPLIST holds one spin-summed block in {len(derived)} "
                   f"folder(s); the up/dn split shown was derived from COHPCAR "
                   f"ratios, the totals are exact")

    pats = defaultdict(list)
    for r in results:
        mn = tm_moments(r["scf"])
        if mn and sum(abs(x) for x in mn) / len(mn) >= QUENCHED_MOMENT:
            pats[pattern(mn)].append(r["name"])
    if len(pats) > 1:
        for p, names in sorted(pats.items(), key=lambda kv: -len(kv[1])):
            out.append(f"moment pattern {p}: {', '.join(names)} -- folders "
                       f"with different patterns are not directly comparable")

    for (pair, sh), m in sorted(meta.items()):
        if not m["fixed"] or m["nfold"] < 3:
            continue
        vals = [e["tot"] for _, e in series_values(results, idx, pair, sh)]
        mean = sum(vals) / len(vals) if vals else 0
        if len(vals) >= 3 and mean:
            rel = (max(vals) - min(vals)) / mean
            if rel > CTRL_SPREAD_WARN:
                out.append(f"{pair} shell {sh} sits at a fixed distance yet "
                           f"-ICOHP varies by {100 * rel:.0f}% -- an electronic "
                           f"effect, worth explaining in the text")
    return out


# ----------------------------------------------------------------------------
def write_csv(results, path, idx, args):
    import csv
    cols = ["folder", "z", "conf", "pair", "shell", "d_A", "N",
            "icohp_up", "icohp_dn", "icohp_tot", "spin_src",
            "cohp_ef_up", "cohp_ef_dn",
            "anti_pct", "spill_up", "spill_dn", "e_lost", "E_Ry", "dE_meV",
            "m_TM", "pattern", "m_Gd", "M_tot", "M_abs", "Ef_eV", "n_iter"]
    ref = reference_energy(results)
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in results:
            s, lo = r["scf"] or {}, r["lob"] or {}
            mn, gd = tm_moments(s), s.get("moments", {}).get("Gd") or []
            sp = lo.get("spilling") or []
            dE = ((s["E_Ry"] - ref) * RY_TO_MEV
                  if ref is not None and "E_Ry" in s else None)
            tail = [fnum(sp[0], 2) if len(sp) > 0 else "",
                    fnum(sp[1], 2) if len(sp) > 1 else "",
                    fnum(lo["e_tot"] - lo["e_rec"], 4)
                    if lo.get("e_rec") is not None and lo.get("e_tot") else "",
                    fnum(s.get("E_Ry"), 6), fnum(dE, 1),
                    fnum(sum(abs(x) for x in mn) / len(mn), 3) if mn else "",
                    pattern(mn),
                    fnum(sum(gd) / len(gd), 3) if gd else "",
                    fnum(s.get("M_tot"), 2), fnum(s.get("M_abs"), 2),
                    fnum(s.get("Ef"), 3), s.get("n_iter", "")]
            for (pair, d), e in sorted(r["icohp"].items()):
                c = r["cohp"].get((pair, d), {})
                ef = c.get("ef") or []
                w.writerow([r["name"], fnum(z_of(r["name"]), 4),
                            config_of(r["name"]), pair,
                            shell_of(idx, r["name"], pair, d) or "",
                            fnum(e["d"], 3),
                            e["n"], fnum(e.get("up"), 4), fnum(e.get("dn"), 4),
                            fnum(e["tot"], 4), e.get("split", "-"),
                            fnum(ef[0], 4) if len(ef) > 0 else "",
                            fnum(ef[1], 4) if len(ef) > 1 else "",
                            fnum(c.get("anti"), 2)] + tail)
                n += 1
    return n, len(cols)


def folders_from_conf(path):
    if not os.path.isfile(path):
        return []
    c = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            c[k.strip()] = v.strip().strip("\"'")
    root = c.get("SERIES_DIR")
    if not root or not os.path.isdir(root):
        return []
    return [os.path.join(root, p) for p in sorted(os.listdir(root))
            if os.path.isfile(os.path.join(root, p, "ICOHPLIST.lobster"))]


def main():
    ap = argparse.ArgumentParser(
        description="Quality control summary for a series of LOBSTER calculations.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*",
                    help="calculation folders (shell globs are fine); with "
                         "none given, they are taken from lobster.conf")
    ap.add_argument("--conf", default="lobster.conf",
                    help="config used when no folders are given")
    ap.add_argument("--csv", help="write the flat table to a CSV file")
    ap.add_argument("--out", help="write the report to a file (default: stdout)")
    ap.add_argument("--species", action="store_true",
                    help="group by QE species, splitting Mn1 from Mn2")
    ap.add_argument("--no-cohp", action="store_true",
                    help="skip COHPCAR: no -COHP at E_F, no antibonding share")
    ap.add_argument("--cutoff", type=float, default=PAIR_CUTOFF,
                    help=f"pair radius in A (default {PAIR_CUTOFF})")
    ap.add_argument("--shells", type=int, default=MAX_SHELLS,
                    help=f"shells per pair to report (default {MAX_SHELLS})")
    args = ap.parse_args()

    folders = []
    for d in args.dirs:
        for h in (sorted(glob.glob(d)) or [d]):
            if os.path.isdir(h):
                folders.append(h)
    if not folders:
        folders = folders_from_conf(args.conf)
    if not folders:
        sys.exit("No folders found.")
    folders.sort(key=sort_key)

    results = []
    for f in folders:
        sys.stderr.write(f"  reading {f} ...\n")
        sys.stderr.flush()
        results.append(collect(f, args))

    report = render(results, args)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(report + "\n")
        sys.stderr.write(f"report: {args.out}\n")
    else:
        print(report)
    if args.csv:
        idx, _ = shell_index(results)
        n, k = write_csv(results, args.csv, idx, args)
        sys.stderr.write(f"CSV: {args.csv} -- {n} rows, {k} columns\n")


if __name__ == "__main__":
    main()
