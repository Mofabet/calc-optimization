#!/usr/bin/env python3
"""
Merge the up and down averaged files into the final table.

Usage:
    08_table.py avg_chem_up.dat avg_chem_dn.dat --label GdMnSi
    08_table.py ... --pair Mn-Si --nearest      # only the shortest bond
"""
import argparse
from collections import OrderedDict


def load(path):
    out = OrderedDict()
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        out[(p[0], float(p[1]))] = float(p[4])      # rms_t
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("up")
    ap.add_argument("dn")
    ap.add_argument("--label", default="")
    ap.add_argument("--pair", help="restrict to one chemical pair")
    ap.add_argument("--nearest", action="store_true",
                    help="only the shortest distance of each pair")
    args = ap.parse_args()

    up, dn = load(args.up), load(args.dn)
    keys = [k for k in up if k in dn]
    if args.pair:
        keys = [k for k in keys if k[0] == args.pair]
    if args.nearest:
        seen, keep = {}, []
        for pair, d in sorted(keys, key=lambda x: (x[0], x[1])):
            if pair not in seen:
                seen[pair] = d
                keep.append((pair, d))
        keys = keep

    head = "| " + (" system |" if args.label else "")
    print(head + " pair | d, Å | t_eff↑, eV | t_eff↓, eV | ⟨t_eff⟩, eV |")
    print("|" + ("---:|" if args.label else "") +
          "---:|---:|---:|---:|---:|")
    for pair, d in sorted(keys, key=lambda x: (x[0], x[1])):
        tu, td = up[(pair, d)], dn[(pair, d)]
        row = f"| {args.label} " if args.label else "| "
        print(f"{row}| {pair} | {d:.3f} | {tu:.3f} | {td:.3f} | "
              f"{0.5 * (tu + td):.3f} |")

    missing = [k for k in up if k not in dn] + [k for k in dn if k not in up]
    if missing:
        print(f"\n<!-- {len(missing)} entries present in only one spin "
              f"channel were dropped -->")


if __name__ == "__main__":
    main()
