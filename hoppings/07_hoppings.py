#!/usr/bin/env python3
"""
07 -- _hr.dat -> hoppings_<spin>.dat, avg_chem_<spin>.dat, avg_atom_<spin>.dat

    t_eff(A,B,R) = sqrt( < |H_ij(R)|^2 > ),  i in A, j in B

Three things this does not do naively:

Wannier functions are matched to atoms by their final centres, read from
<seed>_centres.xyz, not by projection index order. After disentanglement a
function can move onto a bond or swap with a neighbour, and trusting the index
order then mislabels the output while leaving the numbers plausible.

H(R) is Hermitian, so (A,B,R) and (B,A,-R) are one bond. Only the canonical
half is accumulated; the other half is used as a consistency check.

The averages report both the mean and the RMS over equivalent bonds. Quote the
RMS -- t_eff is itself an RMS -- and read a gap between the two as a warning
that supposedly equivalent bonds are not.

    07_hoppings.py                      both spins, everything
    07_hoppings.py --spin up
    07_hoppings.py --pairs Mn-Si --tmin 0.05
    07_hoppings.py --onsite             also print on-site energies
"""
import argparse
import math
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp


def elem(label):
    return re.sub(r"\d+$", "", label)


def extract(hr, win, centres, pairs, dmax, tmin, assign_mode, log):
    lattice, atoms, blocks = qp.parse_win(win)
    num_wann, nrpts, data = qp.read_hr(hr)

    n_proj = sum(n for _, n in blocks)
    if n_proj != num_wann:
        sys.exit(f"ERROR: {win} implies {n_proj} Wannier functions, {hr} has "
                 f"{num_wann}; the two are from different runs")

    order = []
    for label, n in blocks:
        order += [label] * n

    wf_atom, method, worst = None, "projections", None
    if assign_mode == "centres" and os.path.exists(centres):
        c = qp.parse_centres(centres, num_wann)
        if c:
            wf_atom, worst = qp.assign_by_centres(c, atoms, lattice, blocks)
            method = "centres"
    if wf_atom is None:
        wf_atom = order

    print(f"# assignment: {method}", file=log)
    if worst is not None:
        print(f"# worst centre-to-atom distance: {worst:.3f} Ang", file=log)
        if worst > 1.0:
            print("# WARNING: a Wannier centre sits far from every atom; check "
                  "the spreads before trusting these numbers", file=log)
    mism = [(i + 1, order[i], wf_atom[i]) for i in range(num_wann)
            if order[i] != wf_atom[i]]
    if mism:
        print(f"# {len(mism)} functions are not where the projection order "
              f"says:", file=log)
        for i, exp, got in mism[:12]:
            print(f"#   WF {i:3d}: expected {exp:>5s}, found on {got}", file=log)

    idx2atom = {i + 1: wf_atom[i] for i in range(num_wann)}
    counter = defaultdict(int)
    pos = {}
    for el, f in atoms:
        counter[el] += 1
        pos[f"{el}{counter[el]}"] = f

    want = None
    if pairs.strip():
        want = {tuple(sorted(p.strip().split("-"))) for p in pairs.split(",")}

    def canon(a, b, r):
        return min((a, b, r), (b, a, (-r[0], -r[1], -r[2])))

    hop, mirror, onsite = defaultdict(list), defaultdict(list), defaultdict(list)
    for line in data:
        p = line.split()
        if len(p) < 7:
            continue
        R = (int(p[0]), int(p[1]), int(p[2]))
        i, j = int(p[3]), int(p[4])
        v2 = float(p[5]) ** 2 + float(p[6]) ** 2
        A, B = idx2atom[i], idx2atom[j]
        if R == (0, 0, 0) and i == j:
            onsite[A].append(float(p[5]))
            continue
        if R == (0, 0, 0) and A == B:
            continue
        if want and tuple(sorted((elem(A), elem(B)))) not in want:
            continue
        key = canon(A, B, R)
        (hop if (A, B, R) == key else mirror)[key].append(v2)

    worst_h = 0.0
    for key, v in hop.items():
        mv = mirror.get(key)
        if mv and len(mv) == len(v):
            t1 = math.sqrt(sum(v) / len(v))
            t2 = math.sqrt(sum(mv) / len(mv))
            if t1 > 1e-6:
                worst_h = max(worst_h, abs(t1 - t2) / t1)
    if worst_h > 1e-3:
        print(f"# WARNING: H(R) is not Hermitian to better than {worst_h:.1%}; "
              f"the _hr.dat may be truncated", file=log)

    rows = []
    for (A, B, R), v in hop.items():
        df = [pos[B][t] + R[t] - pos[A][t] for t in range(3)]
        d = math.sqrt(sum(x * x for x in qp.frac_to_cart(lattice, df)))
        if d > dmax:
            continue
        t = math.sqrt(sum(v) / len(v))
        if t >= tmin:
            rows.append((d, A, B, R, t, len(v)))
    return sorted(rows), onsite


def average(rows, by, ndec=3):
    groups = defaultdict(list)
    for d, A, B, R, t, n in rows:
        a, b = (A, B) if by == "atom" else (elem(A), elem(B))
        groups[("-".join(sorted([a, b])), round(d, ndec))].append(t)
    out = []
    for (pair, d), v in sorted(groups.items()):
        out.append((pair, d, len(v), sum(v) / len(v),
                    math.sqrt(sum(x * x for x in v) / len(v)), min(v), max(v)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--conf", default="system.conf")
    ap.add_argument("--spin", choices=["up", "dn", "both"], default="both")
    ap.add_argument("--pairs")
    ap.add_argument("--dmax", type=float)
    ap.add_argument("--tmin", type=float)
    ap.add_argument("--assign", choices=["centres", "projections"],
                    default="centres")
    ap.add_argument("--onsite", action="store_true")
    args = ap.parse_args()

    c = qp.read_conf(args.conf) if os.path.exists(args.conf) else {}
    sysname = c.get("SYSTEM", "")
    pairs = args.pairs if args.pairs is not None else c.get("PAIRS", "")
    dmax = args.dmax if args.dmax is not None else float(c.get("DMAX", 5.0))
    tmin = args.tmin if args.tmin is not None else float(c.get("TMIN", 0.02))

    for spin in (["up", "dn"] if args.spin == "both" else [args.spin]):
        seed = f"{sysname}_{spin}" if sysname else None
        hr = f"{seed}_hr.dat" if seed else None
        if not hr or not os.path.exists(hr):
            cand = [f for f in os.listdir(".") if f.endswith(f"_{spin}_hr.dat")]
            if not cand:
                print(f"  {spin}: no *_{spin}_hr.dat here, skipped")
                continue
            hr = cand[0]
            seed = hr[:-len("_hr.dat")]
        win, cen = f"{seed}.win", f"{seed}_centres.xyz"

        with open(f"assign_{spin}.log", "w") as log:
            rows, onsite = extract(hr, win, cen, pairs, dmax, tmin,
                                   args.assign, log)

        with open(f"hoppings_{spin}.dat", "w") as f:
            f.write(f"# {hr}   pairs={pairs or 'all'} dmax={dmax} tmin={tmin}\n")
            f.write("# A B R1 R2 R3 distance_Ang t_eff_eV n_elements\n")
            for d, A, B, R, t, n in rows:
                f.write(f"{A:4s} {B:4s} {R[0]:3d} {R[1]:3d} {R[2]:3d} "
                        f"{d:9.4f} {t:10.5f} {n:4d}\n")

        for by in ("chem", "atom"):
            with open(f"avg_{by}_{spin}.dat", "w") as f:
                f.write(f"# from {hr}, grouped by {by}\n")
                f.write("# pair distance_Ang n avg_t rms_t min_t max_t\n")
                for pair, d, n, avg, rms, lo, hi in average(rows, by):
                    f.write(f"{pair:<10s} {d:12.3f} {n:4d} {avg:10.5f} "
                            f"{rms:10.5f} {lo:10.5f} {hi:10.5f}\n")

        warn = sum(1 for l in open(f"assign_{spin}.log") if "WARNING" in l)
        print(f"  {spin}: {len(rows)} bonds -> hoppings_{spin}.dat, "
              f"avg_chem_{spin}.dat, avg_atom_{spin}.dat"
              + (f"   [{warn} assignment warning(s)]" if warn else ""))

        if args.onsite and onsite:
            by_el = defaultdict(list)
            print(f"    on-site energies, {spin}")
            for A in sorted(onsite):
                m = sum(onsite[A]) / len(onsite[A])
                by_el[elem(A)].append(m)
                print(f"      {A:5s} {len(onsite[A]):3d} orbitals  {m:9.4f} eV")
            keys = sorted(by_el)
            for i, a in enumerate(keys):
                for b in keys[i + 1:]:
                    ea = sum(by_el[a]) / len(by_el[a])
                    eb = sum(by_el[b]) / len(by_el[b])
                    print(f"      eps({a}) - eps({b}) = {ea - eb:8.4f} eV")


if __name__ == "__main__":
    main()
