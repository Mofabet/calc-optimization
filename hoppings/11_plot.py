#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
11_plot.py -- the main series figure, built from the CSV written by 10_qc.py.

Top panel: t_eff(d) with an exponential fit, and -ICOHP(d) on a second axis.
Bottom panel: the local moment and the on-site exchange splitting.

    ./11_plot.py qc.csv --pair Mn-Si --out fig_hyb.pdf
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load(path, pair, shell=1):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        sys.exit("Empty CSV.")
    dcol, tcol = f"d_{pair}_s{shell}", f"teff_{pair}_s{shell}"
    if dcol not in rows[0]:
        avail = sorted({re.sub(r"^teff_", "", k) for k in rows[0] if k.startswith("teff_")})
        sys.exit(f"No pair {pair} shell {shell}. Available: {', '.join(avail) or 'none'}")

    per = defaultdict(dict)
    for r in rows:
        e = per[r["folder"]]
        e.setdefault("d", f(r.get(dcol)))
        e.setdefault("t", {})[r["spin"]] = f(r.get(tcol))
        m = [abs(f(v)) for k, v in r.items()
             if k.startswith("magn_") and f(v) is not None]
        if m:
            e["m"] = max(m)
        tm = pair.split("-")[0]
        tm = tm if tm != "Si" else pair.split("-")[1]
        on = f(r.get(f"onsite_{tm}"))
        if on is not None:
            e.setdefault("onsite", {})[r["spin"]] = on
        for k, v in r.items():
            if k.startswith(f"icohp_{pair}_") and f(v) is not None:
                d_i = f(k.rsplit("_", 1)[1])
                cur = e.get("icohp")
                if cur is None or (e.get("d") and abs(d_i - e["d"]) < abs(cur[0] - e["d"])):
                    e["icohp"] = (d_i, f(v))

    pts = []
    for folder, e in per.items():
        if e.get("d") is None:
            continue
        tv = [x for x in e.get("t", {}).values() if x]
        on = e.get("onsite", {})
        pts.append({
            "folder": folder, "d": e["d"],
            "t": float(np.mean(tv)) if tv else None,
            "t_up": e.get("t", {}).get("up"), "t_dn": e.get("t", {}).get("dn"),
            "icohp": e["icohp"][1] if e.get("icohp") else None,
            "m": e.get("m"),
            "dex": (on["dn"] - on["up"]) if ("up" in on and "dn" in on) else None,
        })
    pts.sort(key=lambda p: p["d"])
    return pts


def expfit(d, y):
    b, a = np.polyfit(d, np.log(y), 1)
    return math.exp(a), -b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--pair", default=None, help="e.g. Mn-Si (guessed if omitted)")
    ap.add_argument("--shell", type=int, default=1)
    ap.add_argument("--out", default="fig_hyb.pdf")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    if args.pair is None:
        head = next(csv.reader(open(args.csv)))
        cands = [re.sub(r"^teff_|_s\d+$", "", k) for k in head if k.startswith("teff_")]
        cands = [c for c in cands if "Si" in c and c != "Si-Si"]
        cands = [c for c in cands if "Gd" not in c] or cands
        if not cands:
            sys.exit("Could not guess the pair, use --pair")
        args.pair = cands[0]
        sys.stderr.write(f"pair: {args.pair}\n")

    pts = load(args.csv, args.pair, args.shell)
    if len(pts) < 2:
        sys.exit("Too few points to plot.")
    d = np.array([p["d"] for p in pts])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 6.6), sharex=True,
                                   gridspec_kw={"height_ratios": [1.35, 1]})

    t = np.array([p["t"] if p["t"] else np.nan for p in pts])
    ok = np.isfinite(t) & (t > 0)
    ax1.plot(d[ok], t[ok], "o", ms=7, color="#1f4e79", label=f"$t_{{eff}}$ {args.pair}")
    for key, mk in (("t_up", "^"), ("t_dn", "v")):
        v = np.array([p[key] if p[key] else np.nan for p in pts])
        if np.isfinite(v).sum() > 1:
            ax1.plot(d, v, mk, ms=4, mfc="none", color="#1f4e79", alpha=0.7,
                     label="up" if key == "t_up" else "dn")
    if ok.sum() >= 3:
        A, b = expfit(d[ok], t[ok])
        xs = np.linspace(d[ok].min() - 0.03, d[ok].max() + 0.03, 100)
        ax1.plot(xs, A * np.exp(-b * xs), "-", lw=1.2, color="#1f4e79", alpha=0.6,
                 label=rf"$t_0e^{{-\beta d}}$, $\beta$ = {b:.2f} Å$^{{-1}}$")
    ax1.set_ylabel(r"$t_{\rm eff}$, eV", color="#1f4e79")
    ax1.tick_params(axis="y", labelcolor="#1f4e79")

    ic = np.array([p["icohp"] if p["icohp"] else np.nan for p in pts])
    if np.isfinite(ic).sum() >= 2:
        axb = ax1.twinx()
        oki = np.isfinite(ic) & (ic > 0)
        axb.plot(d[oki], ic[oki], "s", ms=6, color="#a33", label="$-$ICOHP")
        if oki.sum() >= 3:
            A2, b2 = expfit(d[oki], ic[oki])
            xs = np.linspace(d[oki].min() - 0.03, d[oki].max() + 0.03, 100)
            axb.plot(xs, A2 * np.exp(-b2 * xs), "--", lw=1.2, color="#a33", alpha=0.6,
                     label=rf"$\beta$ = {b2:.2f} Å$^{{-1}}$")
        axb.set_ylabel(r"$-$ICOHP, eV", color="#a33")
        axb.tick_params(axis="y", labelcolor="#a33")
        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = axb.get_legend_handles_labels()
        ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right", framealpha=0.9)
    else:
        ax1.legend(fontsize=8, loc="upper right")
    if args.title:
        ax1.set_title(args.title)

    m = np.array([p["m"] if p["m"] is not None else np.nan for p in pts])
    if np.isfinite(m).sum():
        ax2.plot(d, m, "o-", ms=6, color="#2a7", label=r"$|m|$, $\mu_B$")
    dex = np.array([p["dex"] if p["dex"] is not None else np.nan for p in pts])
    if np.isfinite(dex).sum() >= 2:
        axc = ax2.twinx()
        axc.plot(d, np.abs(dex), "d--", ms=5, color="#846", label=r"$\Delta_{ex}$, eV")
        axc.set_ylabel(r"$\Delta_{\rm ex}$, eV", color="#846")
        axc.tick_params(axis="y", labelcolor="#846")
        h1, l1 = ax2.get_legend_handles_labels()
        h2, l2 = axc.get_legend_handles_labels()
        ax2.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower right")
    else:
        ax2.legend(fontsize=8, loc="lower right")

    # threshold: between the last quenched point and the first magnetic one
    if np.isfinite(m).sum() >= 2:
        q = np.where(m < 0.5)[0]
        g = np.where(m >= 0.5)[0]
        if q.size and g.size and q.max() < g.min():
            xc = 0.5 * (d[q.max()] + d[g.min()])
            for ax in (ax1, ax2):
                ax.axvline(xc, color="k", lw=0.8, ls=":", alpha=0.6)
            ax2.annotate("threshold", xy=(xc, ax2.get_ylim()[1] * 0.85),
                         xytext=(4, 0), textcoords="offset points", fontsize=8)

    ax2.set_ylabel(r"$|m|$, $\mu_B$", color="#2a7")
    ax2.tick_params(axis="y", labelcolor="#2a7")
    ax2.set_xlabel(rf"$d$({args.pair.replace('-', '–')}), Å")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    sys.stderr.write(f"figure: {args.out}\n")

    print(f"{'folder':>8}  {'d, A':>7}  {'t_eff':>8}  {'-ICOHP':>8}  "
          f"{'|m|':>6}  {'D_ex':>7}")
    for p in pts:
        print(f"{p['folder']:>8}  {p['d']:7.3f}  "
              f"{p['t'] if p['t'] else float('nan'):8.4f}  "
              f"{p['icohp'] if p['icohp'] else float('nan'):8.3f}  "
              f"{p['m'] if p['m'] is not None else float('nan'):6.3f}  "
              f"{p['dex'] if p['dex'] is not None else float('nan'):7.3f}")


if __name__ == "__main__":
    main()