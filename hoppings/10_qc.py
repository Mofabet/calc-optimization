#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
10_qc.py -- quality control and post-processing summary for a wannier series.

Walks a set of calculation folders, extracts everything needed to judge model
quality and to fill the thesis tables, prints a readable report and writes a
flat CSV.

Reads (anything missing is skipped silently):
    *.wout              Omega_I/D/OD, WF centres and spreads, windows, atoms, lattice
    *.win               k-mesh (for the band interpolation check)
    *.eig               DFT eigenvalues: band counts in the windows, interpolation reference
    *_hr.dat            t_eff per atom pair, on-site block, hermiticity, tail decay
    scf.out / nscf.out  E_F, magnetic moments, total energy
    ICOHPLIST.lobster   -ICOHP per bond

Everything is computed from the primary files -- the output of 07_hoppings.py is
not used -- so this doubles as an independent check of the pipeline.

Examples:
    ./10_qc.py GdMnSi/restart/z*  --csv qc.csv --out qc.txt
    ./10_qc.py z12 z2074 --no-hr            # fast: Omega and windows only
    ./10_qc.py z*  --bands-nk 300           # stricter interpolation check
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from dataclasses import dataclass, field
from itertools import islice

import numpy as np

# ----------------------------------------------------------------------------
# Warning thresholds -- edit here
# ----------------------------------------------------------------------------
RATIO_OK = (1.10, 1.50)      # healthy corridor for Omega/Omega_I
RATIO_SPREAD_WARN = 0.08     # deviation of ratio from the series median, fraction
CENTRE_OFFSET_WARN = 0.90    # A, WF centre displacement from its own atom
BAND_COUNT_SPREAD_WARN = 4   # spread of the frozen-window band count across the series
UPDN_TEFF_WARN = 0.05        # up/dn mismatch when the moment is quenched
QUENCHED_MOMENT = 0.5        # muB, below this the moment counts as quenched
HERMIT_WARN = 1e-6           # eV, max |H(R) - H(-R)^+|
TAIL_WARN = 0.30             # far H(R) as a fraction of the near ones
BANDS_MEDIAN_WARN = 5.0      # meV, median |E_wann - E_dft| inside the frozen window
FIXED_TOL = 0.01             # A, spread below which a distance counts as fixed
CTRL_SPREAD_WARN = 0.05      # t_eff spread at a fixed distance

PAIR_CUTOFF = 4.0            # A, pair radius for t_eff
MAX_SHELLS = 3               # how many shells per pair to report
BANDS_NK = 200               # k-points sampled for the interpolation check

RY_TO_EV = 13.605693122990
FLOAT = r"[-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?"


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


def strip_index(label):
    m = re.match(r"^([A-Za-z]{1,2})", label)
    return m.group(1) if m else label


def _isfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# ----------------------------------------------------------------------------
# .wout
# ----------------------------------------------------------------------------
@dataclass
class Wout:
    path: str = ""
    num_wann: int | None = None
    num_bands: int | None = None
    omega_i: float | None = None
    omega_d: float | None = None
    omega_od: float | None = None
    omega_tot: float | None = None
    win_outer: tuple | None = None
    win_inner: tuple | None = None
    lattice: np.ndarray | None = None
    atoms: list = field(default_factory=list)   # (label, x, y, z) in A
    centres: np.ndarray | None = None
    spreads: np.ndarray | None = None
    converged: bool | None = None

    @property
    def ratio(self):
        return self.omega_tot / self.omega_i if (self.omega_tot and self.omega_i) else None

    @property
    def omega_i_per_wf(self):
        return self.omega_i / self.num_wann if (self.omega_i and self.num_wann) else None


def parse_wout(path):
    try:
        text = open(path, "r", errors="replace").read()
    except OSError:
        return None
    w = Wout(path=path)
    lines = text.splitlines()

    def gi(pat):
        m = re.search(pat, text)
        return int(m.group(1)) if m else None

    w.num_wann = gi(r"Number of Wannier Functions\s*:\s*(\d+)")
    w.num_bands = (gi(r"Number of input Bloch states\s*:\s*(\d+)")
                   or gi(r"Number of Bands\s*:\s*(\d+)"))

    def win(tag, kmin, kmax):
        m = re.search(tag + r":?\s*(" + FLOAT + r")\s*(?:to)?\s+(" + FLOAT + r")", text)
        if m:
            return float(m.group(1)), float(m.group(2))
        a = re.search(kmin + r"\s*[:=]\s*(" + FLOAT + ")", text)
        b = re.search(kmax + r"\s*[:=]\s*(" + FLOAT + ")", text)
        return (float(a.group(1)), float(b.group(1))) if a and b else None

    w.win_outer = win("Outer", "dis_win_min", "dis_win_max")
    w.win_inner = win("(?:Inner|Frozen)", "dis_froz_min", "dis_froz_max")

    lat = []
    for i, ln in enumerate(lines):
        if "Lattice Vectors" in ln:
            for j in range(i + 1, min(i + 5, len(lines))):
                if re.search(r"a_\d", lines[j]):
                    nums = re.findall(FLOAT, lines[j])
                    if len(nums) >= 4:
                        lat.append([float(x) for x in nums[-3:]])
            break
    if len(lat) == 3:
        w.lattice = np.array(lat)

    atom_re = re.compile(
        r"^\|\s*([A-Za-z]{1,2}[0-9]*)\s+(\d+)\s+" + r"(" + FLOAT + r")\s+(" + FLOAT +
        r")\s+(" + FLOAT + r")\s*\|\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" +
        FLOAT + r")\s*\|")
    seen = set()
    for ln in lines:
        m = atom_re.match(ln.strip())
        if m and (m.group(1), m.group(2)) not in seen:
            seen.add((m.group(1), m.group(2)))
            w.atoms.append((m.group(1), float(m.group(6)),
                            float(m.group(7)), float(m.group(8))))

    fs = text.rfind("Final State")
    if fs != -1:
        wf_re = re.compile(
            r"WF cent(?:re|er) and spread\s+\d+\s*\(\s*(" + FLOAT + r")\s*,\s*(" +
            FLOAT + r")\s*,\s*(" + FLOAT + r")\s*\)\s*(" + FLOAT + r")")
        c, s = [], []
        for m in wf_re.finditer(text[fs:]):
            c.append([float(m.group(i)) for i in (1, 2, 3)])
            s.append(float(m.group(4)))
        if c:
            w.centres, w.spreads = np.array(c), np.array(s)

    def last(pat):
        ms = re.findall(pat + r"\s*=\s*(" + FLOAT + ")", text)
        return float(ms[-1].replace("D", "E").replace("d", "e")) if ms else None

    w.omega_i, w.omega_d = last(r"Omega\s*I"), last(r"Omega\s*D")
    w.omega_od, w.omega_tot = last(r"Omega\s*OD"), last(r"Omega\s*Total")
    if w.omega_tot is None and None not in (w.omega_i, w.omega_d, w.omega_od):
        w.omega_tot = w.omega_i + w.omega_d + w.omega_od
    w.converged = ("convergence criteria satisfied" in text.lower()
                   or "all done" in text.lower())
    return w


def assign_wf(w: Wout):
    """Nearest atom for every WF (periodic images included) and the centre offset."""
    if w.centres is None or not w.atoms or w.lattice is None:
        return None, None
    pos = np.array([[a[1], a[2], a[3]] for a in w.atoms])
    sh = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1)
                   for k in (-1, 0, 1)]) @ w.lattice
    img = pos[:, None, :] + sh[None, :, :]
    d = np.linalg.norm(w.centres[:, None, None, :] - img[None, :, :, :], axis=-1)
    dmin = d.min(axis=2)
    return dmin.argmin(axis=1), dmin.min(axis=1)


def parse_win_kpoints(path):
    try:
        text = open(path, "r", errors="replace").read()
    except OSError:
        return None
    m = re.search(r"begin\s+kpoints(.*?)end\s+kpoints", text, re.S | re.I)
    if not m:
        return None
    k = []
    for ln in m.group(1).splitlines():
        p = ln.split()
        if len(p) >= 3 and all(_isfloat(x) for x in p[:3]):
            k.append([float(x) for x in p[:3]])
    return np.array(k) if k else None


# ----------------------------------------------------------------------------
# .eig
# ----------------------------------------------------------------------------
def parse_eig(path):
    try:
        raw = np.loadtxt(path)
    except Exception:
        return None
    if raw.ndim != 2 or raw.shape[1] < 3:
        return None
    nb, nk = int(raw[:, 0].max()), int(raw[:, 1].max())
    if nb * nk != raw.shape[0]:
        return None
    e = np.zeros((nb, nk))
    e[raw[:, 0].astype(int) - 1, raw[:, 1].astype(int) - 1] = raw[:, 2]
    return e


def bands_in_window(e, win):
    if e is None or not win:
        return None
    c = ((e >= win[0]) & (e <= win[1])).sum(axis=0)
    return int(c.min()), int(c.max()), float(np.median(c))


# ----------------------------------------------------------------------------
# _hr.dat
# ----------------------------------------------------------------------------
@dataclass
class HR:
    nwann: int = 0
    rvecs: np.ndarray | None = None
    ham: np.ndarray | None = None
    deg: np.ndarray | None = None


def parse_hr(path):
    """Block-wise read: peak memory is the Hamiltonian itself, not the text."""
    try:
        fh = open(path, "r", errors="replace")
    except OSError:
        return None
    with fh:
        fh.readline()
        try:
            nw = int(fh.readline().split()[0])
            nr = int(fh.readline().split()[0])
        except (ValueError, IndexError):
            return None
        deg = []
        while len(deg) < nr:
            ln = fh.readline()
            if not ln:
                return None
            deg += [int(x) for x in ln.split()]
        deg = np.array(deg[:nr], dtype=float)
        nlin = nw * nw
        ham = np.empty((nr, nw, nw), dtype=complex)
        rvec = np.empty((nr, 3), dtype=int)
        got = 0
        for ir in range(nr):
            buf = list(islice(fh, nlin))
            if len(buf) < nlin:
                break
            arr = np.array("".join(buf).split(), dtype=float)
            if arr.size != nlin * 7:
                break
            arr = arr.reshape(nlin, 7)
            rvec[ir] = arr[0, :3].astype(int)
            # in the file the m index runs fastest -> [n, m], so transpose
            ham[ir] = (arr[:, 5] + 1j * arr[:, 6]).reshape(nw, nw).T
            got += 1
    if got == 0:
        return None
    return HR(nwann=nw, rvecs=rvec[:got], ham=ham[:got], deg=deg[:got])


def hermiticity_error(h: HR):
    idx = {tuple(r): i for i, r in enumerate(h.rvecs)}
    err = 0.0
    for r, i in idx.items():
        j = idx.get((-r[0], -r[1], -r[2]))
        if j is not None:
            err = max(err, float(np.abs(h.ham[i] - h.ham[j].conj().T).max()))
    return err


def tail_ratio(h: HR, lattice):
    if lattice is None:
        return None
    rc = np.linalg.norm(h.rvecs @ lattice, axis=1)
    a = np.abs(h.ham).max(axis=(1, 2))
    if rc.max() <= 0:
        return None
    near = a[(rc > 0.1) & (rc < 0.35 * rc.max())]
    far = a[rc > 0.80 * rc.max()]
    if near.size == 0 or far.size == 0:
        return None
    return float(far.max() / near.max())


def compute_teff(h: HR, w: Wout, owner, cutoff):
    """t_eff(A,B,d) = sqrt(<|H_ij|^2>) over the atom-atom block, averaged over
    equivalent pairs at the same distance."""
    pos = np.array([[a[1], a[2], a[3]] for a in w.atoms])
    lab = [strip_index(a[0]) for a in w.atoms]
    groups = [np.where(owner == i)[0] for i in range(len(w.atoms))]
    rc = h.rvecs @ w.lattice
    keep = np.where(np.linalg.norm(rc, axis=1) <= cutoff + 1.0 +
                    np.linalg.norm(pos - pos.mean(0), axis=1).max() * 2)[0]
    acc = {}
    for ir in keep:
        H = h.ham[ir]
        R = rc[ir]
        for ia in range(len(pos)):
            gi = groups[ia]
            if gi.size == 0:
                continue
            dist = np.linalg.norm(pos + R - pos[ia], axis=1)
            for ib in np.where((dist > 0.1) & (dist <= cutoff))[0]:
                gj = groups[ib]
                if gj.size == 0:
                    continue
                blk = H[np.ix_(gi, gj)]
                pair = "-".join(sorted((lab[ia], lab[ib])))
                key = (pair, round(float(dist[ib]), 3))
                a = acc.setdefault(key, [0.0, 0])
                a[0] += float(np.mean(np.abs(blk) ** 2))
                a[1] += 1
    # merge near-degenerate distances inside one folder
    recs = [{"pair": pair, "dist": d, "s2": s2, "n": n}
            for (pair, d), (s2, n) in acc.items()]
    recs.sort(key=lambda r: (r["pair"], r["dist"]))
    merged = []
    for r in recs:
        if merged and merged[-1]["pair"] == r["pair"] and \
           abs(merged[-1]["dist"] - r["dist"]) < 0.005:
            m = merged[-1]
            m["s2"] += r["s2"]
            m["dist"] = (m["dist"] * m["n"] + r["dist"] * r["n"]) / (m["n"] + r["n"])
            m["n"] += r["n"]
        else:
            merged.append(dict(r))
    out = [{"pair": r["pair"], "dist": r["dist"],
            "teff": math.sqrt(r["s2"] / r["n"]), "n": r["n"]} for r in merged]
    # number the shells within each pair
    by_pair = {}
    for r in out:
        by_pair.setdefault(r["pair"], []).append(r)
    for lst in by_pair.values():
        lst.sort(key=lambda x: x["dist"])
        for i, r in enumerate(lst):
            r["shell"] = i + 1
    return out


def onsite_levels(h: HR, w: Wout, owner):
    i0 = np.where((h.rvecs == 0).all(axis=1))[0]
    if i0.size == 0:
        return {}
    H0 = h.ham[i0[0]]
    per = {}
    for ia, a in enumerate(w.atoms):
        g = np.where(owner == ia)[0]
        if g.size == 0:
            continue
        ev = np.linalg.eigvalsh(H0[np.ix_(g, g)])
        per.setdefault(strip_index(a[0]), []).append(
            (float(ev.mean()), float(ev.max() - ev.min()), int(g.size)))
    return {k: (float(np.mean([x[0] for x in v])),
                float(np.mean([x[1] for x in v])),
                v[0][2], len(v)) for k, v in per.items()}


def bands_check(h: HR, kpts, eig, win_inner, nk_max=BANDS_NK, seed=0):
    """
    Wannier interpolation on the original k-mesh against DFT.
    -> (median, 90th percentile, max |dE| in meV, fraction within 10 meV, n)
    """
    if kpts is None or eig is None or win_inner is None:
        return None
    nk = min(len(kpts), eig.shape[1])
    if nk == 0:
        return None
    rng = np.random.default_rng(seed)
    sel = np.arange(nk) if nk <= nk_max else np.sort(rng.choice(nk, nk_max, replace=False))
    R = h.rvecs.astype(float)
    inv_deg = 1.0 / h.deg
    lo, hi = win_inner
    resid = []
    for ik in sel:
        ph = np.exp(2j * np.pi * (kpts[ik] @ R.T)) * inv_deg
        Hk = np.tensordot(ph, h.ham, axes=(0, 0))
        Hk = 0.5 * (Hk + Hk.conj().T)
        ew = np.linalg.eigvalsh(Hk)
        ed = eig[:, ik]
        ed = ed[(ed >= lo) & (ed <= hi)]
        if ed.size == 0:
            continue
        resid.append(np.abs(ed[:, None] - ew[None, :]).min(axis=1))
    if not resid:
        return None
    r = np.concatenate(resid) * 1000.0     # meV
    return (float(np.median(r)), float(np.percentile(r, 90)), float(r.max()),
            float((r < 10.0).mean()), int(r.size))


# ----------------------------------------------------------------------------
# QE and LOBSTER
# ----------------------------------------------------------------------------
def parse_qe(path):
    out = {"efermi": None, "moments": {}, "etot": None, "nelec": None,
           "abs_magn": None, "tot_magn": None}
    try:
        text = open(path, "r", errors="replace").read()
    except OSError:
        return out
    ms = re.findall(r"the Fermi energy is\s+(" + FLOAT + r")\s*ev", text)
    if ms:
        out["efermi"] = float(ms[-1])
    else:
        ms = re.findall(r"highest occupied[^:]*:\s*(" + FLOAT + ")", text)
        if ms:
            out["efermi"] = float(ms[-1])
    ms = re.findall(r"!\s+total energy\s*=\s*(" + FLOAT + r")\s*Ry", text)
    if ms:
        out["etot"] = float(ms[-1]) * RY_TO_EV
    ms = re.findall(r"number of electrons\s*=\s*(" + FLOAT + ")", text)
    if ms:
        out["nelec"] = float(ms[-1])
    ms = re.findall(r"total magnetization\s*=\s*(" + FLOAT + ")", text)
    if ms:
        out["tot_magn"] = float(ms[-1])
    ms = re.findall(r"absolute magnetization\s*=\s*(" + FLOAT + ")", text)
    if ms:
        out["abs_magn"] = float(ms[-1])
    blocks = re.findall(r"atom\s+(\d+)\s*\(R=[^)]*\)\s*charge=\s*(" + FLOAT +
                        r")\s*magn=\s*(" + FLOAT + ")", text)
    if blocks:
        n = max(int(b[0]) for b in blocks)
        out["moments"] = {int(b[0]): float(b[2]) for b in blocks[-n:]}
    return out


def parse_icohp(path):
    try:
        lines = open(path, "r", errors="replace").read().splitlines()
    except OSError:
        return []
    acc = {}
    for ln in lines:
        p = ln.split()
        if len(p) < 6 or not p[0].isdigit():
            continue
        if not (_isfloat(p[3]) and not _isfloat(p[1]) and not _isfloat(p[2])):
            continue
        dist = float(p[3])
        val = None
        for tok in p[7:]:
            if _isfloat(tok):
                val = float(tok)
                break
        if val is None and _isfloat(p[-1]):
            val = float(p[-1])
        if val is None:
            continue
        pair = "-".join(sorted((strip_index(p[1]), strip_index(p[2]))))
        key = (pair, round(dist, 2))
        a = acc.setdefault(key, {"pair": pair, "d": [], "v": []})
        a["d"].append(dist)
        a["v"].append(val)
    out = [{"pair": v["pair"], "dist": float(np.mean(v["d"])),
            "icohp": float(np.sum(v["v"]) / max(1, len(v["v"]) // 2)) if len(v["v"]) > 1
            else v["v"][0], "n": len(v["v"])} for v in acc.values()]
    out.sort(key=lambda r: (r["pair"], r["dist"]))
    return out


# ----------------------------------------------------------------------------
# Collecting one folder
# ----------------------------------------------------------------------------
def spin_of(name):
    n = os.path.basename(name).lower()
    if re.search(r"(^|[._-])(up|spn1|spin1)([._-]|$)", n):
        return "up"
    if re.search(r"(^|[._-])(dn|down|spn2|spin2)([._-]|$)", n):
        return "dn"
    return "-"


def pick(cands, sp):
    same = [c for c in cands if spin_of(c) == sp]
    return same[0] if same else (cands[0] if cands else None)


def collect(folder, args):
    res = {"folder": folder, "name": os.path.basename(os.path.normpath(folder)),
           "spins": {}, "qe": {}, "icohp": [], "notes": []}

    qe = {}
    for pat in ("scf.out", "*.scf.out", "*scf.out", "nscf.out", "*nscf.out",
                "*/scf.out", "*/nscf.out"):
        for hpath in sorted(glob.glob(os.path.join(folder, pat))):
            q = parse_qe(hpath)
            for k, v in q.items():
                if v not in (None, {}, []):
                    qe.setdefault(k, v)
    res["qe"] = qe
    if not qe:
        res["notes"].append("no scf/nscf.out -- no E_F and no moments")

    for pat in ("ICOHPLIST.lobster", "*/ICOHPLIST.lobster", "*/*/ICOHPLIST.lobster"):
        hits = sorted(glob.glob(os.path.join(folder, pat)))
        if hits:
            res["icohp"] = parse_icohp(hits[0])
            break

    wouts = sorted(glob.glob(os.path.join(folder, "**", "*.wout"), recursive=True))
    if not wouts:
        res["notes"].append("no .wout -- the run never reached Wannier90")
        return res

    for wp in wouts:
        w = parse_wout(wp)
        if w is None or w.num_wann is None:
            res["notes"].append(f"{os.path.basename(wp)}: could not be parsed")
            continue
        sp = spin_of(wp)
        owner, offset = assign_wf(w)
        d = os.path.dirname(wp)
        stem = wp[:-5]
        e = {"wout": w, "owner": owner, "offset": offset, "bands_outer": None,
             "bands_inner": None, "teff": [], "onsite": {}, "herm": None,
             "tail": None, "bands": None}

        eigp = (stem + ".eig") if os.path.exists(stem + ".eig") else \
            pick(sorted(glob.glob(os.path.join(d, "*.eig"))), sp)
        eig = parse_eig(eigp) if eigp else None
        if eig is not None:
            e["bands_outer"] = bands_in_window(eig, w.win_outer)
            e["bands_inner"] = bands_in_window(eig, w.win_inner)

        if not args.no_hr:
            hrp = (stem + "_hr.dat") if os.path.exists(stem + "_hr.dat") else \
                pick(sorted(glob.glob(os.path.join(d, "*_hr.dat"))), sp)
            if hrp:
                h = parse_hr(hrp)
                if h is None:
                    res["notes"].append(f"{os.path.basename(hrp)}: could not be read")
                elif h.nwann != w.num_wann:
                    res["notes"].append(
                        f"{os.path.basename(hrp)}: num_wann {h.nwann} != {w.num_wann}")
                elif owner is None:
                    res["notes"].append(f"{os.path.basename(wp)}: no atoms or lattice")
                else:
                    e["herm"] = hermiticity_error(h)
                    e["tail"] = tail_ratio(h, w.lattice)
                    e["teff"] = compute_teff(h, w, owner, args.cutoff)
                    e["onsite"] = onsite_levels(h, w, owner)
                    if not args.no_bands:
                        winp = (stem + ".win") if os.path.exists(stem + ".win") else \
                            pick(sorted(glob.glob(os.path.join(d, "*.win"))), sp)
                        kp = parse_win_kpoints(winp) if winp else None
                        e["bands"] = bands_check(h, kp, eig, w.win_inner,
                                                 nk_max=args.bands_nk)
                    del h
        res["spins"][sp] = e
    return res


# ----------------------------------------------------------------------------
# Series assembly: shells
# ----------------------------------------------------------------------------
def shell_index(results):
    """{(pair, shell): {folder: {spin: rec}}} plus whether the distance is fixed."""
    idx = {}
    for r in results:
        for sp, e in r["spins"].items():
            for rec in e["teff"]:
                key = (rec["pair"], rec["shell"])
                idx.setdefault(key, {}).setdefault(r["name"], {})[sp] = rec
    meta = {}
    for key, per_folder in idx.items():
        ds = [rec["dist"] for f in per_folder.values() for rec in f.values()]
        meta[key] = {"dmin": min(ds), "dmax": max(ds),
                     "fixed": (max(ds) - min(ds)) < FIXED_TOL,
                     "nfold": len(per_folder)}
    return idx, meta


def elem_spreads(w: Wout, owner):
    if owner is None or w.spreads is None:
        return {}
    per = {}
    for ia, a in enumerate(w.atoms):
        g = np.where(owner == ia)[0]
        if g.size:
            per.setdefault(strip_index(a[0]), []).extend(w.spreads[g].tolist())
    return {k: float(np.mean(v)) for k, v in per.items()}


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def render(results, args):
    idx, meta = shell_index(results)
    L = ["WANNIER SERIES QUALITY REPORT",
         f"folders: {len(results)}   pairs up to {args.cutoff:.1f} A   "
         f"shells per pair: {args.shells}"]

    # 1 --------------------------------------------------------------------
    L.append(section("1. Model quality"))
    rows = []
    for r in results:
        for sp, e in sorted(r["spins"].items()):
            w, off, bi, bo, bc = e["wout"], e["offset"], e["bands_inner"], \
                e["bands_outer"], e["bands"]
            rows.append([
                r["name"], sp, w.num_wann or "-", w.num_bands or "-",
                fnum(w.omega_i, 2), fnum(w.omega_tot, 2), fnum(w.ratio, 3),
                fnum(w.omega_i_per_wf, 2),
                f"{bi[0]}-{bi[1]}" if bi else "-",
                f"{bo[0]}-{bo[1]}" if bo else "-",
                fnum(float(off.max()) if off is not None else None, 2),
                fnum(bc[0], 2) if bc else "-",
                f"{bc[3]*100:.0f}" if bc else "-",
                "yes" if w.converged else "NO"])
    L.append(table(["folder", "spin", "N_wf", "N_bnd", "Omega_I", "Omega", "ratio",
                    "Om_I/N", "bnd frz", "bnd out", "offs", "dE meV", "<10meV %",
                    "conv"], rows, aligns=["<", "<"]))
    L.append("""
  ratio     Omega/Omega_I -- healthy corridor 1.10-1.50; constancy across the
            series matters more than the value itself
  Om_I/N    invariant part per function, A^2
  bnd frz   number of DFT bands inside the frozen window (min-max over k)
  offs      largest distance from a WF centre to its own atom, A
  dE meV    median |E_Wannier - E_DFT| inside the frozen window on the original mesh
  <10meV %  fraction of states within 10 meV""")

    # 2 --------------------------------------------------------------------
    L.append(section("2. Wannier function spreads"))
    elems = []
    for r in results:
        for e in r["spins"].values():
            for k in elem_spreads(e["wout"], e["owner"]):
                if k not in elems:
                    elems.append(k)
    rows = []
    for r in results:
        for sp, e in sorted(r["spins"].items()):
            es = elem_spreads(e["wout"], e["owner"])
            w = e["wout"]
            worst = ""
            if w.spreads is not None and e["owner"] is not None:
                k = np.argsort(w.spreads)[-4:][::-1]
                worst = ", ".join(
                    f"{strip_index(w.atoms[e['owner'][i]][0])} {w.spreads[i]:.2f}"
                    for i in k)
            rows.append([r["name"], sp] + [fnum(es.get(el), 2) for el in elems] + [worst])
    L.append(table(["folder", "spin"] + [f"<{el}>" for el in elems] +
                   ["four widest, A^2"], rows,
                   aligns=["<", "<"] + [">"] * len(elems) + ["<"]))

    # 3 --------------------------------------------------------------------
    L.append(section("3. t_eff by pair and shell (eV)"))
    rows = []
    for r in results:
        pairs_seen = {}
        for sp, e in sorted(r["spins"].items()):
            for rec in e["teff"]:
                if rec["shell"] > args.shells:
                    continue
                pairs_seen.setdefault((rec["pair"], rec["shell"]), {})[sp] = rec
        for (pair, sh), d in sorted(pairs_seen.items(),
                                    key=lambda x: (x[0][0], x[0][1])):
            up, dn, one = d.get("up"), d.get("dn"), d.get("-")
            dist = (up or dn or one)["dist"]
            tu = up["teff"] if up else (one["teff"] if one else None)
            td = dn["teff"] if dn else None
            avg = np.mean([x for x in (tu, td) if x is not None]) if (tu or td) else None
            dev = (abs(td - tu) / tu * 100) if (tu and td and tu > 0) else None
            rows.append([r["name"], pair, sh, f"{dist:.3f}", fnum(tu, 4),
                         fnum(td, 4), fnum(avg, 4), fnum(dev, 1)])
    L.append(table(["folder", "pair", "sh", "d, A", "t(up)", "t(dn)", "t(mean)",
                    "up/dn, %"], rows, aligns=["<", "<", ">"]))
    L.append("\n  'sh' is the coordination shell index within the pair.")

    # 4 --------------------------------------------------------------------
    L.append(section("4. Internal control: shells at a fixed distance"))
    L.append("""  These distances do not change along the series, so t_eff should not
  change either. The spread is the systematic error of the method -- round
  the t_eff values in the text to it.
""")
    rows = []
    ctrl = []
    for (pair, sh), m in sorted(meta.items()):
        if not m["fixed"] or m["nfold"] < 2 or sh > args.shells:
            continue
        t = np.array([rec["teff"] for f in idx[(pair, sh)].values()
                      for rec in f.values()])
        rel = float(t.std(ddof=1) / t.mean() * 100) if t.size > 1 and t.mean() else None
        rows.append([pair, sh, f"{m['dmin']:.3f}", m["nfold"], len(t),
                     fnum(float(t.mean()), 4),
                     fnum(float(t.std(ddof=1)) if t.size > 1 else None, 4),
                     fnum(rel, 1), fnum(float(t.min()), 4), fnum(float(t.max()), 4)])
        ctrl.append((pair, sh, rel))
    L.append(table(["pair", "sh", "d, A", "folders", "n", "mean", "sigma",
                    "sigma/t, %", "min", "max"], rows, aligns=["<", ">"]))
    cand = [c for c in ctrl if c[2] is not None and "Gd" in c[0] and
            not any(tm in c[0] for tm in ("Mn", "Fe"))]
    best = min(cand or [c for c in ctrl if c[2] is not None],
               key=lambda c: c[2], default=None)
    if best:
        L.append(f"""
  Reference: {best[0]}, shell {best[1]} -- spread {best[2]:.1f} %.
  Use a pair whose d is fixed geometrically AND whose electronic structure does
  not change along the series (usually Gd-Gd). Pairs involving the transition
  metal at a fixed d carry a physical signal (the moment changes), not an error.""")

    # 5 --------------------------------------------------------------------
    L.append(section("5. On-site levels (R = 0) and exchange splitting"))
    rows = []
    for r in results:
        els = set()
        for e in r["spins"].values():
            els |= set(e["onsite"])
        for el in sorted(els):
            up = r["spins"].get("up", {}).get("onsite", {}).get(el)
            dn = r["spins"].get("dn", {}).get("onsite", {}).get(el)
            one = r["spins"].get("-", {}).get("onsite", {}).get(el)
            src = up or dn or one
            ex = (dn[0] - up[0]) if (up and dn) else None
            rows.append([r["name"], el, src[2] if src else "-",
                         fnum((up or one)[0] if (up or one) else None, 3),
                         fnum(dn[0] if dn else None, 3), fnum(ex, 3),
                         fnum((up or one)[1] if (up or one) else None, 3)])
    L.append(table(["folder", "elem", "N_wf", "eps(up)", "eps(dn)", "D_ex",
                    "block width"], rows, aligns=["<", "<"]))
    L.append("""
  D_ex = <eps_dn> - <eps_up> over the on-site block: a direct indicator of the
  moment that costs no extra calculation. For an atom with N_wf = 5 this is a
  clean d block; 'block width' is the crystal-field splitting of the levels.""")

    # 6 --------------------------------------------------------------------
    L.append(section("6. Magnetic moments, energies, -ICOHP"))
    rows = []
    for r in results:
        mom = r["qe"].get("moments") or {}
        mx = max((abs(v) for v in mom.values()), default=None)
        s = ", ".join(f"{i}:{v:+.2f}" for i, v in sorted(mom.items()))
        rows.append([r["name"], fnum(r["qe"].get("efermi"), 4),
                     fnum(r["qe"].get("etot"), 3),
                     fnum(r["qe"].get("abs_magn"), 3), fnum(mx, 3),
                     s[:52] + ("..." if len(s) > 52 else "")])
    L.append(table(["folder", "E_F, eV", "E_tot, eV", "|M| cell", "max |m|",
                    "per atom"], rows, aligns=["<", ">", ">", ">", ">", "<"]))
    if any(r["icohp"] for r in results):
        rows = []
        for r in results:
            for rec in r["icohp"]:
                rows.append([r["name"], rec["pair"], f"{rec['dist']:.3f}",
                             fnum(rec["icohp"], 3), rec["n"]])
        L.append("")
        L.append(table(["folder", "bond", "d, A", "-ICOHP, eV", "lines"], rows,
                       aligns=["<", "<"]))

    # 7 --------------------------------------------------------------------
    L.append(section("7. Bond-length dependence"))
    L.append(fit_report(results, idx, meta, args))

    # 8 --------------------------------------------------------------------
    L.append(section("8. Warnings"))
    fl = flags(results, idx, meta)
    L.append("\n".join("  " + f for f in fl) if fl else "  Nothing triggered.")

    notes = [f"  {r['name']}: {n}" for r in results for n in r["notes"]]
    if notes:
        L.append(section("9. Could not be read"))
        L.append("\n".join(notes))
    return "\n".join(L)


def loglinfit(x, y):
    """y = A exp(-b x); -> A, b, RMS of the residuals in the logarithm."""
    b, a = np.polyfit(x, np.log(y), 1)
    rms = float(np.sqrt(np.mean((np.log(y) - (a + b * x)) ** 2)))
    return math.exp(a), -b, rms


def fit_report(results, idx, meta, args):
    out = []
    rows = []
    for (pair, sh), m in sorted(meta.items()):
        if m["fixed"] or m["nfold"] < 3 or sh > 1:
            continue
        pts = [(rec["dist"], rec["teff"])
               for per_sp in idx[(pair, sh)].values() for rec in per_sp.values()]
        d = np.array([p[0] for p in pts])
        t = np.array([p[1] for p in pts])
        ok = t > 0
        if ok.sum() < 3 or d[ok].max() - d[ok].min() < 0.05:
            continue
        A, b, rms = loglinfit(d[ok], t[ok])
        p_eff = b * float(d[ok].mean())
        row = [pair, int(ok.sum()), f"{d[ok].min():.3f}-{d[ok].max():.3f}",
               fnum(A, 2), fnum(b, 3), fnum(p_eff, 2), fnum(rms * 100, 1)]
        # -ICOHP for the same pair
        ipts = []
        for r in results:
            same = [rec for rec in r["icohp"] if rec["pair"] == pair]
            if same:
                ipts.append(min(same, key=lambda x: x["dist"]))
        if len(ipts) >= 3:
            di = np.array([x["dist"] for x in ipts])
            vi = np.array([x["icohp"] for x in ipts])
            msk = vi > 0
            if msk.sum() >= 3 and di[msk].max() - di[msk].min() > 0.05:
                _, bi, _ = loglinfit(di[msk], vi[msk])
                row += [fnum(bi, 3), fnum(bi / (2 * b) if b else None, 2)]
            else:
                row += ["-", "-"]
        else:
            row += ["-", "-"]
        rows.append(row)
    if not rows:
        out.append("  Not enough points with a varying distance yet -- the fit will")
        out.append("  appear once the remaining folders are done.")
        return "\n".join(out)
    out.append(table(["pair", "n", "d range, A", "t_0, eV", "beta, 1/A", "p_eff",
                      "RMS, %", "beta(ICOHP)", "ratio"], rows, aligns=["<", ">"]))
    out.append("""
  t_eff = t_0 * exp(-beta * d), fitted on the nearest shell, both spins pooled.
  p_eff = beta * <d> -- the equivalent power-law exponent for t ~ d^(-p);
          Harrison's rule for p-d gives p = 3.5.
  ratio = beta(ICOHP) / (2 * beta(t_eff)). In the two-level limit -ICOHP is
          proportional to t^2 / dE, so a value near 1 means the two measures
          agree and the 't_eff or -ICOHP' question is settled.""")
    return "\n".join(out)


def flags(results, idx, meta):
    out = []
    ratios = []
    for r in results:
        for sp, e in sorted(r["spins"].items()):
            w, tag = e["wout"], f"{r['name']}/{sp}"
            if w.ratio:
                ratios.append((tag, w.ratio))
                if not (RATIO_OK[0] <= w.ratio <= RATIO_OK[1]):
                    out.append(f"[ratio] {tag}: Omega/Omega_I = {w.ratio:.3f} outside "
                               f"the {RATIO_OK[0]}-{RATIO_OK[1]} corridor")
            if e["offset"] is not None and e["offset"].max() > CENTRE_OFFSET_WARN:
                i = int(e["offset"].argmax())
                out.append(f"[centre] {tag}: function {i+1} sits {e['offset'][i]:.2f} A "
                           f"from its nearest atom -- the atom-block split behind "
                           f"t_eff stops being meaningful")
            if w.converged is False:
                out.append(f"[convergence] {tag}: the minimisation did not report "
                           f"convergence")
            if e["herm"] is not None and e["herm"] > HERMIT_WARN:
                out.append(f"[hermiticity] {tag}: max|H(R)-H(-R)^+| = {e['herm']:.2e} eV")
            if e["tail"] is not None and e["tail"] > TAIL_WARN:
                out.append(f"[tails] {tag}: distant H(R) reach {e['tail']*100:.0f} % of "
                           f"the near ones -- the functions are not localised enough")
            if e["bands"] and e["bands"][0] > BANDS_MEDIAN_WARN:
                out.append(f"[bands] {tag}: median |E_W-E_DFT| = {e['bands'][0]:.1f} meV "
                           f"in the frozen window, only {e['bands'][3]*100:.0f} % of "
                           f"states within 10 meV")
    if len(ratios) > 1:
        v = np.array([x[1] for x in ratios])
        med = float(np.median(v))
        for tag, x in ratios:
            if med and abs(x - med) / med > RATIO_SPREAD_WARN:
                out.append(f"[series] {tag}: ratio {x:.3f} departs from the median "
                           f"{med:.3f} by more than {RATIO_SPREAD_WARN*100:.0f} % -- "
                           f"this point is not comparable with the rest")
    counts = [(f"{r['name']}/{sp}", e["bands_inner"][2])
              for r in results for sp, e in sorted(r["spins"].items())
              if e["bands_inner"]]
    if len(counts) > 1:
        v = np.array([c[1] for c in counts])
        if v.max() - v.min() > BAND_COUNT_SPREAD_WARN:
            out.append(f"[window] the frozen-window band count drifts across the series "
                       f"({v.min():.0f}-{v.max():.0f}): part of the t_eff difference may "
                       f"be a window artefact rather than physics")
    for r in results:
        mom = r["qe"].get("moments") or {}
        if mom and max(abs(x) for x in mom.values()) < QUENCHED_MOMENT:
            up, dn = r["spins"].get("up"), r["spins"].get("dn")
            if up and dn:
                for rec in up["teff"]:
                    if rec["shell"] != 1:
                        continue
                    o = [x for x in dn["teff"] if x["pair"] == rec["pair"]
                         and x["shell"] == 1]
                    if o and rec["teff"] > 1e-6:
                        rel = abs(o[0]["teff"] - rec["teff"]) / rec["teff"]
                        if rel > UPDN_TEFF_WARN:
                            out.append(
                                f"[up/dn] {r['name']}: the moment is quenched, yet "
                                f"t_eff {rec['pair']} differs by {rel*100:.0f} % -- "
                                f"that is a minimisation error, not physics")
                            break
    for (pair, sh), m in sorted(meta.items()):
        if not m["fixed"] or m["nfold"] < 3 or sh > 1:
            continue
        t = np.array([rec["teff"] for f in idx[(pair, sh)].values()
                      for rec in f.values()])
        if t.size > 2 and t.mean() and t.std(ddof=1) / t.mean() > CTRL_SPREAD_WARN:
            out.append(f"[reference] {pair} at the fixed d = {m['dmin']:.3f} A scatters "
                       f"by {t.std(ddof=1)/t.mean()*100:.1f} % -- either physics (the "
                       f"moment changes) or the real uncertainty of t_eff")
    return out


# ----------------------------------------------------------------------------
# CSV
# ----------------------------------------------------------------------------
def write_csv(results, path, args):
    rows = []
    for r in results:
        for sp, e in sorted(r["spins"].items()):
            w = e["wout"]
            row = {
                "folder": r["name"], "spin": sp, "num_wann": w.num_wann,
                "num_bands": w.num_bands, "omega_i": w.omega_i, "omega_d": w.omega_d,
                "omega_od": w.omega_od, "omega_tot": w.omega_tot, "ratio": w.ratio,
                "omega_i_per_wf": w.omega_i_per_wf,
                "win_outer_lo": w.win_outer[0] if w.win_outer else None,
                "win_outer_hi": w.win_outer[1] if w.win_outer else None,
                "win_froz_lo": w.win_inner[0] if w.win_inner else None,
                "win_froz_hi": w.win_inner[1] if w.win_inner else None,
                "bands_froz_min": e["bands_inner"][0] if e["bands_inner"] else None,
                "bands_froz_max": e["bands_inner"][1] if e["bands_inner"] else None,
                "bands_out_min": e["bands_outer"][0] if e["bands_outer"] else None,
                "bands_out_max": e["bands_outer"][1] if e["bands_outer"] else None,
                "max_centre_offset": float(e["offset"].max()) if e["offset"] is not None else None,
                "hermiticity_ev": e["herm"], "tail_ratio": e["tail"],
                "converged": w.converged,
                "bands_med_mev": e["bands"][0] if e["bands"] else None,
                "bands_p90_mev": e["bands"][1] if e["bands"] else None,
                "bands_max_mev": e["bands"][2] if e["bands"] else None,
                "bands_frac_10mev": e["bands"][3] if e["bands"] else None,
                "efermi": r["qe"].get("efermi"), "etot_ev": r["qe"].get("etot"),
                "abs_magn": r["qe"].get("abs_magn"),
                "tot_magn": r["qe"].get("tot_magn"),
            }
            for el, v in elem_spreads(w, e["owner"]).items():
                row[f"spread_{el}"] = v
            for el, v in e["onsite"].items():
                row[f"onsite_{el}"] = v[0]
                row[f"cfwidth_{el}"] = v[1]
            for rec in e["teff"]:
                if rec["shell"] <= args.shells:
                    row[f"d_{rec['pair']}_s{rec['shell']}"] = rec["dist"]
                    row[f"teff_{rec['pair']}_s{rec['shell']}"] = rec["teff"]
            for rec in r["icohp"]:
                row[f"icohp_{rec['pair']}_{rec['dist']:.2f}"] = rec["icohp"]
            for i, v in sorted((r["qe"].get("moments") or {}).items()):
                row[f"magn_{i}"] = v
            rows.append(row)
    keys = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    return len(rows), len(keys)


# ----------------------------------------------------------------------------
def main():
    global PAIR_CUTOFF
    ap = argparse.ArgumentParser(
        description="Quality control summary for a series of wannier calculations.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="calculation folders (shell globs are fine)")
    ap.add_argument("--csv", help="write the flat table to a CSV file")
    ap.add_argument("--out", help="write the report to a file (default: stdout)")
    ap.add_argument("--no-hr", action="store_true",
                    help="skip _hr.dat: fast, but no t_eff, on-site or band check")
    ap.add_argument("--no-bands", action="store_true",
                    help="skip the band interpolation check (the slowest part)")
    ap.add_argument("--cutoff", type=float, default=PAIR_CUTOFF,
                    help=f"pair radius for t_eff in A (default {PAIR_CUTOFF})")
    ap.add_argument("--shells", type=int, default=MAX_SHELLS,
                    help=f"shells per pair to report (default {MAX_SHELLS})")
    ap.add_argument("--bands-nk", type=int, default=BANDS_NK,
                    help=f"k-points for the band check (default {BANDS_NK}, 0 = all)")
    args = ap.parse_args()
    if args.bands_nk <= 0:
        args.bands_nk = 10 ** 9

    folders = []
    for d in args.dirs:
        for h in (sorted(glob.glob(d)) or [d]):
            if os.path.isdir(h):
                folders.append(h)
    if not folders:
        sys.exit("No folders found.")

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
        n, k = write_csv(results, args.csv, args)
        sys.stderr.write(f"CSV: {args.csv} -- {n} rows, {k} columns\n")


if __name__ == "__main__":
    main()