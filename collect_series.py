#!/usr/bin/env python3
"""
Collect a series of z-points into one table.

Walks the given folders, reads each system.conf for the geometry, each
avg_chem_*.dat for the hoppings and each scf.out for the moments, and prints
the table in the same shape you used for GdFeSi:

    z(Si) | d(M-Si) | t_eff up | t_eff dn | <t_eff> | m(M)

Usage:
    collect_series.py z0.15 z0.19 z0.20 z0.22
    collect_series.py z*                       --pair Mn-Si
    collect_series.py z* --csv results.csv
"""
import argparse
import os
import re
import sys

BOHR = 0.529177210903


def read_conf(path):
    c = {}
    for raw in open(path):
        line = raw.split("#", 1)[0].strip()
        if line and "=" in line:
            k, v = line.split("=", 1)
            c[k.strip()] = v.strip()
    return c


def read_avg(path, pair):
    """-> {distance: rms_t}, restricted to one chemical pair"""
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if len(p) >= 5 and p[0] == pair:
            out[float(p[1])] = float(p[4])
    return out


def read_moments(path, tm, n_gd=2):
    """site moments from scf.out, in the order the atoms appear in the input"""
    if not os.path.exists(path):
        return None, None
    log = open(path, errors="ignore").read()
    tot = None
    m = re.findall(r"total magnetization\s*=\s*([-0-9.]+)", log)
    if m:
        tot = float(m[-1])
    blocks = re.findall(r"Magnetic moment per site.*?\n((?:\s*atom.*\n)+)", log)
    if not blocks:
        return None, tot
    vals = []
    for line in blocks[-1].splitlines():
        mm = re.search(r"atom\s+(\d+).*?magn=\s*([-0-9.]+)", line)
        if mm:
            vals.append(float(mm.group(2)))
    # atom order is Gd, Gd, TM, TM, Si, Si
    tm_moments = vals[n_gd:n_gd + 2] if len(vals) >= n_gd + 2 else []
    return tm_moments, tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--pair", help="default: <TM>-Si")
    ap.add_argument("--csv")
    args = ap.parse_args()

    rows, skipped = [], []
    for d in sorted(args.dirs):
        conf_p = os.path.join(d, "system.conf")
        if not os.path.exists(conf_p):
            skipped.append((d, "no system.conf"))
            continue
        c = read_conf(conf_p)
        tm = c.get("TM", "?")
        pair = args.pair or "-".join(sorted([tm, "Si"]))

        alat, coa = float(c["ALAT_BOHR"]), float(c["COA"])
        z_si = float(c["Z_SI"])
        a_ang = alat * BOHR
        dist = a_ang * (0.25 + (z_si * coa) ** 2) ** 0.5

        up_p = os.path.join(d, "avg_chem_up.dat")
        dn_p = os.path.join(d, "avg_chem_dn.dat")
        if not (os.path.exists(up_p) and os.path.exists(dn_p)):
            skipped.append((d, "no avg_chem_{up,dn}.dat -- not analysed yet"))
            continue
        up, dn = read_avg(up_p, pair), read_avg(dn_p, pair)
        common = sorted(set(up) & set(dn))
        if not common:
            skipped.append((d, f"no '{pair}' entry in both spin channels"))
            continue
        # the bond nearest to the geometric M-Si distance
        dsel = min(common, key=lambda x: abs(x - dist))
        tu, td = up[dsel], dn[dsel]

        mom, tot = read_moments(os.path.join(d, "scf.out"), tm)
        rows.append(dict(dir=d, z=z_si, d=dsel, tu=tu, td=td,
                         tm=(tu + td) / 2, mom=mom, tot=tot, pair=pair))

    if not rows:
        sys.exit("nothing to collect -- " +
                 "; ".join(f"{d}: {why}" for d, why in skipped))

    has_mom = any(r["mom"] for r in rows)
    tmlab = rows[0]["pair"].replace("-Si", "").replace("Si-", "")

    hdr = "| z(Si) | d, Å | t_eff↑, eV | t_eff↓, eV | ⟨t_eff⟩, eV |"
    sep = "|---:|---:|---:|---:|---:|"
    if has_mom:
        hdr += f" m({tmlab}), μB |"
        sep += "---:|"
    print(f"Пара {rows[0]['pair']}\n")
    print(hdr)
    print(sep)
    for r in sorted(rows, key=lambda x: x["z"]):
        line = (f"| {r['z']:.3f} | {r['d']:.3f} | {r['tu']:.3f} | "
                f"{r['td']:.3f} | {r['tm']:.3f} |")
        if has_mom:
            mom = r["mom"]
            if not mom:
                line += " – |"
            elif max(mom) - min(mom) < 0.01:
                line += f" {mom[0]:+.3f} |"
            else:
                line += f" {min(mom):+.3f}…{max(mom):+.3f} |"
        print(line)

    if len(rows) > 1:
        s = sorted(rows, key=lambda x: x["d"])
        dt = (s[0]["tm"] - s[-1]["tm"]) / s[0]["tm"] * 100
        print(f"\nОт d = {s[0]['d']:.3f} Å до {s[-1]['d']:.3f} Å "
              f"⟨t_eff⟩ меняется на {dt:+.1f} %.")

    if skipped:
        print("\nПропущено:")
        for d, why in skipped:
            print(f"  {d}: {why}")

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("dir,z_si,d_ang,t_up,t_dn,t_mean,m_min,m_max,total_mag\n")
            for r in sorted(rows, key=lambda x: x["z"]):
                mom = r["mom"] or []
                f.write(f"{r['dir']},{r['z']:.6f},{r['d']:.4f},{r['tu']:.5f},"
                        f"{r['td']:.5f},{r['tm']:.5f},"
                        f"{min(mom) if mom else ''},{max(mom) if mom else ''},"
                        f"{r['tot'] if r['tot'] is not None else ''}\n")
        print(f"\nCSV: {args.csv}")


if __name__ == "__main__":
    main()
