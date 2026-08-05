# Wannier hopping toolkit — quick start

Extracts effective hopping integrals `t_eff` from a completed Quantum ESPRESSO
SCF calculation, via `pw2wannier90.x` and Wannier90, and tabulates them across
a series of structures.

---

## 1. Install

```bash
cp -r hoppings/bin /path/to/project/
cd /path/to/project
chmod +x bin/*.py bin/*.sh
which pw.x pw2wannier90.x wannier90.x
```

Set the site defaults once, in `bin/scf2conf.py` (lines 242–244), so that every
generated configuration inherits them:

```
QE_BIN=/path/to/qe/bin        # directory with pw.x and pw2wannier90.x
W90_BIN=                      # directory with wannier90.x; empty = $PATH
SRUN_OPTS=-n 48               # options passed to srun
```

## 2. Lay out the directories

One directory per structure, `bin/` shared:

```
project/
├── bin/
├── z0.15/
├── z0.19/
└── z0.20/
    ├── scf.in
    ├── scf.out
    └── out/          ← the QE outdir, copied whole
```

`out/` holds the self-consistent charge density that the NSCF reads. Copy the
directory, not just the text files:

```bash
cp -r source/{scf.in,scf.out,out} z0.20/
```

## 3. Check the inputs before spending compute

```bash
for d in z*; do (cd "$d" && python3 ../bin/scf2conf.py -o /dev/null); done
```

Writes nothing. Reports per directory: prefix, lattice parameters, fractional
`z` of each species, nearest metal–metalloid distance, SCF convergence, Fermi
energy and site moments. Verify that the distances match the intended series
and that `out/` is present everywhere.

## 4. Run one structure

```bash
cd z0.20
../bin/run_point.sh
```

Generates the inputs, runs the NSCF, sets and validates the disentanglement
windows, runs Wannier90 for both spin channels, extracts and averages the
hoppings, prints a table.

Individual stages:

| | |
|---|---|
| `../bin/run_point.sh nscf` | NSCF only |
| `../bin/run_point.sh wannier` | windows and Wannier90 |
| `../bin/run_point.sh analyse` | re-extract from an existing `_hr.dat` |

`analyse` is cheap; use it when changing `PAIRS`, `DMAX` or `TMIN` in
`system.conf`.

## 5. Check the result

```bash
cat assign_up.log                                        # expect no warnings
grep "Sum of centres and spreads" *_up.wout | tail -1    # expect 1-3 Å² per WF
head -6 avg_atom_up.dat                                  # equivalent bonds should agree
```

A warning in `assign_up.log` means Wannier functions were not found on the
expected atoms and the labels may be wrong. Spreads of tens of Å² mean the
minimisation failed — widen the outer window or add projections.

## 6. Collect the series

```bash
cd ..
python3 bin/collect_series.py z* --csv results.csv
```

Prints a Markdown table of fractional `z`, bond length, `t_eff` per spin
channel, their mean, and the transition-metal moment. Unfinished directories
are listed separately.

---

## Notes

- `KMESH_NSCF` and `MP_GRID` in `system.conf` must be equal.
- `t_eff` depends on the projections and windows. Keep `PROJECTIONS`,
  `MP_GRID` and the `DIS_*` offsets identical across a series, or the
  comparison is not meaningful.
- An existing `scf.in` is never overwritten; it is checked against
  `system.conf` and the run stops if the geometry or prefix disagree.
- Collinear spin-polarised calculations only.
