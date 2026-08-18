#!/usr/bin/env python3
"""
Collect a series of z-points into one table.

Walks the given folders, reads each system.conf for the geometry, each
avg_chem_*.dat for the hoppings and each scf.out for the moments, and prints
the table in the same shape you used for GdFeSi:

    z(Si) | d(M-Si) | t_eff up | t_eff dn | <t_eff> | m(M)

Usage:
    09_series.py z0.15 z0.19 z0.20 z0.22
    09_series.py z*                       --pair Mn-Si
    09_series.py z* --csv results.csv
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qeparse as qp



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


def read_moments(folder, tm):
    """
    Site moments of the transition metal, in Bohr magnetons.

    Which sites are the transition metal is taken from scf.in rather than
    assumed from a fixed atom ordering, because magnetic supercells put the
    species in whatever order the structure was written in.
    """
    log_p = os.path.join(folder, "scf.out")
    scf_p = os.path.join(folder, "scf.in")
    if not os.path.exists(log_p):
        return None, None
    log = open(log_p, errors="ignore").read()

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

    if os.path.exists(scf_p):
        try:
            st = qp.read_structure(qp.strip_comments(open(scf_p).read()))
            if len(st["atoms"]) == len(vals):
                sel = [v for (_, el, _), v in zip(st["atoms"], vals)
                       if el == tm]
                return sel, tot
        except SystemExit:
            pass
    return None, tot


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

        # D_TM_SI is measured from the real ATOMIC_POSITIONS by 01_conf.py.
        # Older configs stored only the idealised Wyckoff parameters, so fall
        # back to the 2-formula-unit formula when it is absent.
        if "D_TM_SI" not in c:
            skipped.append((d, "no D_TM_SI in system.conf -- regenerate it"))
            continue
        dist = float(c["D_TM_SI"])

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

        mom, tot = read_moments(d, tm)
        rows.append(dict(dir=d, d=dsel, tu=tu, td=td,
                         tm=(tu + td) / 2, mom=mom, tot=tot, pair=pair, nat=c.get("NAT","")))

    if not rows:
        sys.exit("nothing to collect -- " +
                 "; ".join(f"{d}: {why}" for d, why in skipped))

    has_mom = any(r["mom"] for r in rows)
    tmlab = rows[0]["pair"].replace("-Si", "").replace("Si-", "")

    hdr = "| structure | d(%s–Si), Å | t_eff↑, eV | t_eff↓, eV | ⟨t_eff⟩, eV |" % tmlab
    sep = "|:---|---:|---:|---:|---:|"
    if has_mom:
        hdr += f" m({tmlab}), μB |"
        sep += "---:|"
    print(f"Пара {rows[0]['pair']}\n")
    print(hdr)
    print(sep)
    for r in sorted(rows, key=lambda x: x["d"]):
        line = (f"| {r['dir']} | {r['d']:.3f} | {r['tu']:.3f} | "
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
            f.write("dir,nat,d_ang,t_up,t_dn,t_mean,m_min,m_max,total_mag\n")
            for r in sorted(rows, key=lambda x: x["d"]):
                mom = r["mom"] or []
                f.write(f"{r['dir']},{r['nat']},{r['d']:.4f},{r['tu']:.5f},"
                        f"{r['td']:.5f},{r['tm']:.5f},"
                        f"{min(mom) if mom else ''},{max(mom) if mom else ''},"
                        f"{r['tot'] if r['tot'] is not None else ''}\n")
        print(f"\nCSV: {args.csv}")


if __name__ == "__main__":
    main()
