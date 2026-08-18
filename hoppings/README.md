# Wannier hopping toolkit

Effective hopping integrals from a completed Quantum ESPRESSO SCF, via
`pw2wannier90` and Wannier90.

    t_eff(A,B,R) = sqrt( < |H_ij(R)|^2 > ),   i in A,  j in B

Steps are numbered in run order. Each is a standalone command; `run.sh` only
chains them.

---

## Layout

```
project/
├── bin/              the scripts
├── defaults.conf     settings shared by the whole series
└── z18/              one folder per structure
    ├── scf.in
    ├── scf.out
    └── out/          the QE outdir, copied whole
```

`out/` holds the charge density the NSCF reads. Copy the directory, not just
the text files.

---

## Pipeline

| | | |
|---|---|---|
| `01_conf.py` | `scf.in`, `scf.out` → `system.conf` | `-o /dev/null` to inspect only |
| `02_inputs.py` | config → `nscf.in`, `*.win`, `pw2wan_*.in` | `--only win\|nscf\|pw2wan`, `--spin` |
| `03_nscf.sh` | run `pw.x` | `--check` |
| `04_overlaps.sh` | `-pp` → `pw2wannier90` → `.mmn .amn .eig` | `up\|dn\|both`, `--force` |
| `05_windows.py` | set `dis_*` in the `.win`, validated against `.eig` | `--dry-run`, `--ef` |
| `06_wannier.sh` | minimisation → `_hr.dat` | `up\|dn\|both`, `--iter N` |
| `07_hoppings.py` | `_hr.dat` → `hoppings_*`, `avg_chem_*`, `avg_atom_*` | `--spin`, `--pairs`, `--tmin`, `--onsite` |
| `08_table.py` | up + dn → one table | `--nearest`, `--pair` |
| `09_series.py` | all folders → series table | `--csv` |

```bash
cd z18
../bin/run.sh              # 01 through 08
../bin/run.sh 06           # one step
../bin/run.sh 04-06        # a range
../bin/run.sh 04-          # from there on
../bin/run.sh 04-06 -- up  # arguments after -- go to the steps
../bin/run.sh --detach     # background, survives logout
```

`run.sh` stops at the first failure. After changing `defaults.conf`, rerun
`01-02` in each folder.

The overlaps come before the windows because the window check needs `.eig`,
which only `pw2wannier90` produces. `wannier90.x -pp` does not read the window
settings, so the placeholders `02_inputs.py` writes are enough for it.

---

## Tools

| | |
|---|---|
| `spectrum.py` | bands below `E_F+Δ`; `--histogram` shows flat manifolds |
| `pdos.py` | `projwfc.x` output by element and orbital — decides the projections |
| `kpath.py` | LMTO `SYML` → `KPATH`, with unit and zone-folding conversion |
| `spreads.py` | final spreads by element, names off-centre functions |
| `bands_prep.py` / `bands_check.py` | QE bands on the Wannier path, then compare |
| `compare.py` | `t_eff` difference between two runs |
| `converge.sh` | `dis_num_iter` ladder, decided on `t_eff` |
| `icohp.py` | LOBSTER `ICOHPLIST` by bond type |
| `status.sh` | what is running, what is done; `--watch` |
| `report.sh` | every result and check in one text file |
| `rename_species.sh` | `Mn_up` → `Mn1` (QE 7.5 cannot restart DFT+U otherwise) |
| `qeparse.py` | shared library, not run directly |

---

## Settings

`defaults.conf` sits beside `bin/` and holds what must be identical across the
series; `01_conf.py` stamps it into every `system.conf`. `{TM}` expands to the
transition metal found in each structure.

```
PROJECTIONS=Gd:f;Gd:d;{TM}:d;Si:p
DIS_WIN_LO=11.0   DIS_WIN_HI=8.0      offsets from E_F, in eV
DIS_FROZ_LO=4.0   DIS_FROZ_HI=1.5
KMESH_PRIMITIVE=12 12 12              divided by the supercell multiplicity
NBND_WINDOW=10  NBND_MARGIN=10        sets NBND_NSCF per structure
KPATH=                                empty: no band plot
QE_BIN=  W90_BIN=  SRUN_OPTS=  W90_SRUN_OPTS=  OMP_THREADS=
```

`system.conf` is generated. Editing it in one folder is how a series stops
being comparable — change `defaults.conf` and rerun `01-02`.

---

## Checking a result

```bash
cat assign_up.log                          # expect no warnings
python3 ../bin/spreads.py *_up.wout        # 1-3 A^2 per function
head -6 avg_atom_up.dat                    # equivalent bonds should agree
```

`Omega_I` is the invariant part of the spread, fixed by the subspace; `Omega`
minus `Omega_I` is what localisation removes. Their ratio belongs in 1.1-1.5,
and a ratio that climbs with `dis_num_iter` means the model is degrading, not
improving. `06_wannier.sh` prints both.

`t_eff` depends on the projections and the windows, so numbers are comparable
only between runs that used the same ones.

Collinear spin-polarised calculations only. Split magnetic sublattices must be
named `Mn1`/`Mn2`, not `Mn_up`/`Mn_dn`.
