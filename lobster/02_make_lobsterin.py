#!/usr/bin/env python3
"""Stage 2. Write a lobsterin next to every scf.in.

Species labels are read from ATOMIC_SPECIES; the element is the label with
trailing digits and any _suffix stripped (Mn1 -> Mn, Mn_up -> Mn).
"""
import os
import re
import sys

# Valence shells offered to LOBSTER, per element. Must match the valence of
# the pseudopotential; add elements here as needed.
BASIS = {
    "H": "1s",
    "C": "2s 2p",
    "N": "2s 2p",
    "O": "2s 2p",
    "Si": "3s 3p",
    "P": "3s 3p",
    "Mn": "3s 3p 3d 4s",
    "Fe": "3s 3p 3d 4s",
    "Co": "3s 3p 3d 4s",
    "Ni": "3s 3p 3d 4s",
    "Ru": "4s 4p 4d 5s",
    "Gd": "5s 5p 4f 5d 6s",
}

# Rough count of basis functions per shell, for the nbnd sanity check.
NORB = {"s": 1, "p": 3, "d": 5, "f": 7}


def load_conf(path="lobster.conf"):
    c = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            c[k.strip()] = v.strip().strip("\"'")
    return c


def element_of(label):
    return re.sub(r"[_\d].*$", "", label).capitalize()


def read_species(scf_in):
    """Return [(label, element), ...] and the atom count per label."""
    txt = open(scf_in).read()
    m = re.search(r"ATOMIC_SPECIES(.*?)(?=\n\s*[A-Z_]{4,}|\Z)", txt, re.S)
    species = []
    if m:
        for line in m.group(1).strip().splitlines():
            f = line.split()
            if len(f) >= 3:
                species.append((f[0], element_of(f[0])))

    counts = {}
    m = re.search(r"ATOMIC_POSITIONS[^\n]*\n(.*?)(?=\n\s*[A-Z_]{4,}|\Z)",
                  txt, re.S)
    if m:
        for line in m.group(1).strip().splitlines():
            f = line.split()
            if len(f) >= 4:
                counts[f[0]] = counts.get(f[0], 0) + 1
    return species, counts


def write_lobsterin(path, species, conf):
    L = []
    L.append("COHPstartEnergy %s" % conf["COHP_EMIN"])
    L.append("COHPendEnergy   %s" % conf["COHP_EMAX"])
    L.append("gaussianSmearingWidth %s" % conf["SMEARING"])
    L.append("")
    L.append("basisSet %s" % conf["BASIS_SET"])
    for label, elem in species:
        if elem not in BASIS:
            raise SystemExit("no basis entry for element %s (label %s); "
                             "add it to BASIS in 02_make_lobsterin.py"
                             % (elem, label))
        L.append("basisFunctions %s %s" % (label, BASIS[elem]))
    L.append("")
    L.append("cohpGenerator from 0.0 to %s" % conf["COHP_RMAX"])
    L.append("orbitalwise")
    L.append("saveProjectionToFile")
    if conf.get("REAL_SPACE_H", "no").lower().startswith("y"):
        L.append("")
        L.append("realspaceHamiltonian")
        L.append("realspaceOverlap")
    extra = conf.get("EXTRA_LINES", "").strip()
    if extra:
        L.append("")
        L.extend(x.strip() for x in extra.split(";") if x.strip())
    L.append("")
    L.append("!xrange %s %s" % (conf["COHP_EMIN"], conf["COHP_EMAX"]))
    open(path, "w").write("\n".join(L) + "\n")


def nbasis(species, counts):
    tot = 0
    for label, elem in species:
        per = sum(NORB[s[-1]] for s in BASIS[elem].split())
        tot += per * counts.get(label, 0)
    return tot


def main():
    c = load_conf(sys.argv[1] if len(sys.argv) > 1 else "lobster.conf")
    for p in c["POINTS"].split():
        d = os.path.join(c["SERIES_DIR"], p)
        scf_in = os.path.join(d, "scf.in")
        if not os.path.isfile(scf_in):
            print("  [%s] SKIP, run stage 1 first" % p)
            continue
        species, counts = read_species(scf_in)
        write_lobsterin(os.path.join(d, "lobsterin"), species, c)
        nb = nbasis(species, counts)
        flag = "" if int(c["NBND"]) >= nb else "  <-- RAISE NBND"
        print("  [%s] lobsterin written, %d basis functions vs nbnd=%s%s"
              % (p, nb, c["NBND"], flag))


if __name__ == "__main__":
    main()
