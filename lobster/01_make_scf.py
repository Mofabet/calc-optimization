#!/usr/bin/env python3
"""Stage 1. Build a LOBSTER-ready scf.in for every point.

Takes SRC_DIR/<point>/scf.in, forces nosym/noinv/nbnd, replaces the
automatic k-mesh by the full explicit list, writes SERIES_DIR/<point>/scf.in.
Also checks that the pseudopotentials are not ultrasoft.
"""
import os
import re
import sys


def load_conf(path="lobster.conf"):
    c = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            c[k.strip()] = v.strip().strip("\"'")
    return c


def kpoint_block(n1, n2, n3):
    lines = ["K_POINTS crystal", "%d" % (n1 * n2 * n3)]
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                lines.append("  %12.8f %12.8f %12.8f  1.0"
                             % (i / n1, j / n2, k / n3))
    return lines


def check_pseudos(src_in, tag):
    """Warn if any UPF is ultrasoft -- LOBSTER wants PAW or norm-conserving."""
    txt = open(src_in).read()
    pdir = re.search(r"pseudo_dir\s*=\s*['\"]([^'\"]+)", txt)
    if not pdir:
        return
    pdir = pdir.group(1)
    m = re.search(r"ATOMIC_SPECIES(.*?)(?=\n\s*[A-Z_]{4,}|\Z)", txt, re.S)
    if not m:
        return
    for line in m.group(1).strip().splitlines():
        f = line.split()
        if len(f) < 3:
            continue
        path = os.path.join(pdir, f[2])
        if not os.path.isfile(path):
            print("  [%s] pseudo not found: %s" % (tag, f[2]))
            continue
        head = open(path, errors="ignore").read(4000)
        t = re.search(r'pseudo_type\s*=\s*"?([A-Za-z]+)', head)
        t = t.group(1).upper() if t else "?"
        if t.startswith("US") or t.startswith("USPP"):
            print("  [%s] WARNING ultrasoft pseudo: %s" % (tag, f[2]))


def rewrite(src_in, dst_in, nbnd, mesh, tag):
    lines = open(src_in).read().splitlines()
    out, i, n = [], 0, len(lines)
    drop = re.compile(r"^\s*(nosym|noinv|nbnd|wf_collect)\s*=", re.I)

    while i < n:
        line = lines[i]
        low = line.strip().lower()

        if drop.match(line):
            i += 1
            continue

        if low.startswith("&control"):
            # LOBSTER greps scf.in for this keyword; QE 7 accepts it and
            # ignores it (wavefunctions are always collected).
            out.append(line)
            out.append("  wf_collect = .true.")
            i += 1
            continue

        if low.startswith("&system"):
            out.append(line)
            out.append("  nosym = .true.")
            out.append("  noinv = .true.")
            out.append("  nbnd = %d" % nbnd)
            i += 1
            continue

        if low.startswith("k_points"):
            if "automatic" in low:
                i += 2                      # header + mesh line
            elif "gamma" in low:
                i += 1
            else:                           # explicit list: header + count + N
                cnt = int(lines[i + 1].split()[0])
                i += 2 + cnt
            out.extend(kpoint_block(*mesh))
            continue

        out.append(line)
        i += 1

    open(dst_in, "w").write("\n".join(out) + "\n")
    print("  [%s] scf.in written, %d k-points, nbnd=%d"
          % (tag, mesh[0] * mesh[1] * mesh[2], nbnd))


def main():
    c = load_conf(sys.argv[1] if len(sys.argv) > 1 else "lobster.conf")
    mesh = tuple(int(x) for x in c["KMESH"].split())
    nbnd = int(c["NBND"])

    for p in c["POINTS"].split():
        src = os.path.join(c["SRC_DIR"], p, "scf.in")
        if not os.path.isfile(src):
            print("  [%s] SKIP, no %s" % (p, src))
            continue
        dst_dir = os.path.join(c["SERIES_DIR"], p)
        os.makedirs(dst_dir, exist_ok=True)
        check_pseudos(src, p)
        rewrite(src, os.path.join(dst_dir, "scf.in"), nbnd, mesh, p)


if __name__ == "__main__":
    main()
