#!/usr/bin/env python3
"""Stage 6. Collect the series into one table.

Writes SERIES_DIR/icohp_series.csv with, per point and per atom pair:
distance and ICOHP summed over spin channels. Also prints the charge
spilling of each point -- read that column first, the rest is only
meaningful when it is small.
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


def spilling(lobsterout):
    """Return the spilling percentages found in lobsterout."""
    vals = []
    if os.path.isfile(lobsterout):
        for line in open(lobsterout, errors="ignore"):
            if "spilling" in line.lower():
                vals += re.findall(r"([-\d.]+)\s*%", line)
    return vals


def read_icohp(path):
    """Sum ICOHP over spin channels. Key: (labelA, labelB, distance, transl)."""
    data = {}
    if not os.path.isfile(path):
        return data
    for line in open(path, errors="ignore"):
        f = line.split()
        if len(f) < 8 or not f[0].rstrip(".").isdigit():
            continue
        try:
            a, b = f[1], f[2]
            dist = float(f[3])
            transl = tuple(f[4:7])
            icohp = float(f[7])
        except (ValueError, IndexError):
            continue
        key = (a, b, round(dist, 4), transl)
        data[key] = data.get(key, 0.0) + icohp
    return data


def main():
    c = load_conf(sys.argv[1] if len(sys.argv) > 1 else "lobster.conf")
    rows, header_shown = [], False
    for p in c["POINTS"].split():
        d = os.path.join(c["SERIES_DIR"], p)
        lst = os.path.join(d, "ICOHPLIST.lobster")
        if not header_shown and os.path.isfile(lst):
            # Column order differs between LOBSTER versions -- check it here.
            print("  ICOHPLIST header: %s"
                  % open(lst, errors="ignore").readline().strip())
            header_shown = True
        sp = spilling(os.path.join(d, "lobsterout"))
        sp_str = "/".join(sp) if sp else "n/a"
        icohp = read_icohp(os.path.join(d, "ICOHPLIST.lobster"))
        print("  [%s] spilling %s %%, %d pairs" % (p, sp_str, len(icohp)))
        for (a, b, dist, transl), v in sorted(icohp.items(),
                                              key=lambda x: x[0][2]):
            rows.append([p, sp_str, a, b, "%.4f" % dist,
                         " ".join(transl), "%.5f" % v])

    out = os.path.join(c["SERIES_DIR"], "icohp_series.csv")
    with open(out, "w") as fh:
        fh.write("point,spilling_pct,atomA,atomB,dist_A,translation,ICOHP_eV\n")
        for r in rows:
            fh.write(",".join(r) + "\n")
    print("  wrote %s (%d rows)" % (out, len(rows)))


if __name__ == "__main__":
    main()
