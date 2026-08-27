#!/usr/bin/env python3
"""Stage 8. Publication figure + LaTeX table from COHPCAR/ICOHPLIST.

    python3 08_plot_cohp.py lobster.conf z2074 Mn-Si Mn-Mn Gd-Mn
    python3 08_plot_cohp.py lobster.conf z2074 --species Mn1-Mn1 Mn1-Mn2

On-site terms (d = 0) are dropped unless --onsite is given. With --species
the atom index is mapped back to the QE species via ATOMIC_POSITIONS, which
separates the Mn1 / Mn2 magnetic sublattices.

Bonds are grouped by element pair and distance (0.05 A tolerance), then
averaged over the symmetry-equivalent members of each group. Plots -COHP,
so that positive means bonding, which is the convention used in the COHP
literature. Writes cohp_<point>.pdf, cohp_<point>.dat and icohp_<point>.tex.

The column layout of COHPCAR differs between LOBSTER versions, so the
script infers it from the column count and PRINTS what it inferred --
check that line once against your file before trusting the figure.
"""
import os
import re
import sys
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def load_conf(path="lobster.conf"):
    c = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            c[k.strip()] = v.strip().strip("\"'")
    return c


def split_label(lab):
    """LOBSTER atom label -> (name, index). Handles 'Mn3' and 'Mn1_3'."""
    m = re.match(r"^([A-Za-z][A-Za-z]?\d*)_(\d+)$", lab)
    if m:
        return m.group(1), int(m.group(2))
    m = re.match(r"^([A-Za-z]+)(\d+)$", lab)
    if m:
        return m.group(1), int(m.group(2))
    return lab, None


def element(label):
    """Mn1_5 -> Mn ; Si9 -> Si"""
    name, _ = split_label(label)
    return re.sub(r"\d+$", "", name)


def species_map(scf_in):
    """1-based atom index -> species label, from ATOMIC_POSITIONS order.

    LOBSTER writes labels as element+index and loses the QE species names,
    so this map is the only way back to the Mn1 / Mn2 sublattices.
    """
    txt = open(scf_in).read()
    m = re.search(r"ATOMIC_POSITIONS[^\n]*\n(.*?)(?=\n\s*[A-Z_]{4,}|\Z)",
                  txt, re.S)
    out = {}
    if not m:
        return out
    i = 0
    for line in m.group(1).strip().splitlines():
        f = line.split()
        if len(f) >= 4:
            i += 1
            out[i] = f[0]
    return out


def make_atom_id(by_species, smap):
    """Return the function that maps a LOBSTER label to a group name."""
    def atom_id(lab):
        name, idx = split_label(lab)
        if by_species:
            if idx is not None and idx in smap:
                return smap[idx]
            return name
        return re.sub(r"\d+$", "", name)
    return atom_id


def bond_key(a, b, dist, atom_id):
    ea, eb = sorted([atom_id(a), atom_id(b)])
    return "%s-%s" % (ea, eb), round(dist / 0.05) * 0.05


def parse_cohpcar(path):
    """Return (energies, {label_index: {spin: [COHP]}}, labels).

    The numeric block starts after the last "No.N:" label line -- the two
    header lines above it also parse as numbers and must not be mistaken
    for data.
    """
    lines = open(path, errors="ignore").read().splitlines()
    lab_re = re.compile(r"^No\.\d+:(\S+?)->(\S+?)\((-?[\d.]+)\)")

    labels, last_label = [], -1
    for i, line in enumerate(lines):
        m = lab_re.match(line.strip())
        if m:
            labels.append((m.group(1), m.group(2), float(m.group(3))))
            last_label = i
    if last_label < 0:
        raise SystemExit("no 'No.N:' interaction labels found in %s" % path)

    rows = []
    for line in lines[last_label + 1:]:
        f = line.split()
        if len(f) < 4:
            continue
        try:
            rows.append([float(x) for x in f])
        except ValueError:
            continue
    if not rows:
        raise SystemExit("no numeric data found in %s" % path)

    # keep only the dominant row width, in case of stray lines
    widths = {}
    for r in rows:
        widths[len(r)] = widths.get(len(r), 0) + 1
    ncol = max(widths, key=widths.get)
    data = [r for r in rows if len(r) == ncol]
    layout = None
    for nspin in (1, 2):
        for extra in (1, 0):          # leading "average" pair, present or not
            if ncol == 1 + 2 * (len(labels) + extra) * nspin:
                layout = (nspin, extra)
                break
        if layout:
            break
    if not layout:
        raise SystemExit("cannot match %d columns to %d interactions in %s"
                         % (ncol, len(labels), path))
    nspin, extra = layout
    print("  layout: %d interactions, %d columns -> nspin=%d, average block=%s"
          % (len(labels), ncol, nspin, bool(extra)))

    energies = [row[0] for row in data]
    cohp = {}
    for i in range(len(labels)):
        cohp[i] = {}
        for s in range(nspin):
            base = 1 + s * 2 * (len(labels) + extra)
            col = base + 2 * (extra + i)
            cohp[i][s] = [row[col] for row in data]
    return energies, cohp, labels


def parse_icohplist(path, atom_id):
    """{(pair, dist): [ICOHP summed over spin, ...]} keyed per bond instance."""
    per_bond = defaultdict(float)
    meta = {}
    for line in open(path, errors="ignore"):
        f = line.split()
        if len(f) < 8 or not f[0].rstrip(".").isdigit():
            continue
        try:
            a, b, dist, val = f[1], f[2], float(f[3]), float(f[7])
        except (ValueError, IndexError):
            continue
        ident = (a, b, round(dist, 4), tuple(f[4:7]))
        per_bond[ident] += val
        meta[ident] = bond_key(a, b, dist, atom_id)

    groups = defaultdict(list)
    for ident, val in per_bond.items():
        groups[meta[ident]].append((ident[2], val))
    return groups


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = {a for a in sys.argv[1:] if a.startswith("--")}
    by_species = "--species" in flags
    keep_onsite = "--onsite" in flags

    conf = args[0] if args else "lobster.conf"
    c = load_conf(conf)
    point = args[1] if len(args) > 1 else c["POINTS"].split()[0]
    wanted = args[2:] or None

    d = os.path.join(c["SERIES_DIR"], point)
    smap = species_map(os.path.join(d, "scf.in"))
    atom_id = make_atom_id(by_species, smap)
    if by_species and not smap:
        print("  WARNING: no ATOMIC_POSITIONS found, falling back to labels")
    if by_species:
        print("  species mode: %d atoms mapped, species %s"
              % (len(smap), sorted(set(smap.values()))))

    energies, cohp, labels = parse_cohpcar(os.path.join(d, "COHPCAR.lobster"))

    # average -COHP over the members of each bond group
    groups = defaultdict(list)
    for i, (a, b, dist) in enumerate(labels):
        groups[bond_key(a, b, dist, atom_id)].append(i)

    keys = sorted(groups, key=lambda k: k[1])
    if not keep_onsite:
        keys = [k for k in keys if k[1] > 0.1]   # drop on-site terms
    if wanted:
        keys = [k for k in keys if k[0] in wanted]
    if not keys:
        raise SystemExit("no bond groups selected; available: %s"
                         % sorted({k[0] for k in groups}))

    nspin = len(cohp[0])
    curves = {}
    for key in keys:
        idx = groups[key]
        for s in range(nspin):
            avg = [-sum(cohp[i][s][j] for i in idx) / len(idx)
                   for j in range(len(energies))]
            curves[(key, s)] = avg
        print("  %s at %.2f A: %d equivalent bonds" % (key[0], key[1], len(idx)))

    # ---- data file ----
    dat = os.path.join(d, "cohp_%s.dat" % point)
    with open(dat, "w") as fh:
        head = ["E-Ef"] + ["-COHP_%s_%.2fA_spin%d" % (k[0], k[1], s)
                           for k in keys for s in range(nspin)]
        fh.write("# " + "  ".join(head) + "\n")
        for j, e in enumerate(energies):
            row = [e] + [curves[(k, s)][j] for k in keys for s in range(nspin)]
            fh.write("  ".join("%12.6f" % v for v in row) + "\n")

    # ---- figure ----
    if plt:
        fig, axes = plt.subplots(1, len(keys), figsize=(3.0 * len(keys), 4.2),
                                 sharey=True, squeeze=False)
        for ax, key in zip(axes[0], keys):
            for s in range(nspin):
                sign = 1 if s == 0 else -1
                y = [sign * v for v in curves[(key, s)]]
                ax.plot(y, energies, lw=1.1,
                        label="spin up" if s == 0 else "spin down")
            ax.axhline(0, color="k", lw=0.6, ls="--")
            ax.axvline(0, color="k", lw=0.6)
            ax.set_title("%s, %.2f $\\AA$" % (key[0], key[1]), fontsize=10)
            ax.set_xlabel(r"$-$COHP")
        axes[0][0].set_ylabel(r"$E - E_\mathrm{F}$ (eV)")
        axes[0][0].set_ylim(float(c["COHP_EMIN"]), float(c["COHP_EMAX"]))
        axes[0][-1].legend(fontsize=8, frameon=False)
        fig.tight_layout()
        pdf = os.path.join(d, "cohp_%s.pdf" % point)
        fig.savefig(pdf)
        print("  wrote %s" % pdf)
    else:
        print("  matplotlib not available, only the .dat was written")

    # ---- LaTeX table ----
    ic = parse_icohplist(os.path.join(d, "ICOHPLIST.lobster"), atom_id)
    tex = os.path.join(d, "icohp_%s.tex" % point)
    with open(tex, "w") as fh:
        fh.write("\\begin{tabular}{lccc}\n\\hline\n")
        fh.write("Bond & $d$ (\\AA) & $N$ & "
                 "$-$ICOHP (eV/bond) \\\\\n\\hline\n")
        for key in sorted(ic, key=lambda k: k[1]):
            if not keep_onsite and key[1] <= 0.1:
                continue
            if wanted and key[0] not in wanted:
                continue
            vals = [v for _, v in ic[key]]
            dists = [x for x, _ in ic[key]]
            fh.write("%s & %.3f & %d & %.3f \\\\\n"
                     % (key[0], sum(dists) / len(dists), len(vals),
                        -sum(vals) / len(vals)))
        fh.write("\\hline\n\\end{tabular}\n")
    print("  wrote %s and %s" % (dat, tex))


if __name__ == "__main__":
    main()
