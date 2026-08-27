# Wannier hopping toolkit

Effective hopping integrals from a completed Quantum ESPRESSO SCF, via
`pw2wannier90` and Wannier90.

    t_eff(A,B,R) = sqrt( < |H_ij(R)|^2 > ),   i in A,  j in B

Steps are numbered in run order. Each is a standalone command; `run.sh` only
chains them. Steps 01-09 produce one structure or one series; 10-11 read the
finished folders back and turn them into a quality report and the series figure.

---

## Layout

```
project/
├── bin/              the scripts
├── defaults.conf     settings shared by the whole series
├── lobster/          LOBSTER runs, one folder per structure (optional)
└── z18/              one folder per structure
    ├── scf.in
    ├── scf.out
    └── out/          the QE outdir, copied whole
```

`out/` holds the charge density the NSCF reads. Copy the directory, not just
the text files.

Do not leave a second copy of a `.wout` in a subfolder of a structure folder.
`10_qc.py` searches recursively and will tell you which file it picked, but the
only way to be sure is not to have the ambiguity.

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
| `10_qc.py` | all folders → quality report + flat CSV | `--icohp-glob`, `--bands-nk`, `--no-hr` |
| `11_plot.py` | CSV → the series figure | `--pair`, `--shell`, `--title` |

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

`10_qc.py` and `11_plot.py` are not part of `run.sh`. They run over the whole
series at once, after every structure has finished:

```bash
./bin/10_qc.py z12 z15 z18 z20 z2074 z23 \
    --icohp-glob 'lobster/{name}/ICOHPLIST.lobster' \
    --bands-nk 0 \
    --csv qc_$(date +%F).csv --out qc_$(date +%F).txt

./bin/11_plot.py qc_$(date +%F).csv --pair Mn-Si --out fig_hyb.pdf
```

List the folders explicitly rather than using a glob: `z2074` sorts between
`z20` and `z23`, which reads badly in the report. Collect the whole series in a
single run — rows produced by different versions of the script at different
times cannot be compared afterwards.

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
PROJECTIONS=Gd:f;Gd:d;{TM}:d;Si:p;Si:s
DIS_WIN_LO=11.0   DIS_WIN_HI=8.0      offsets from E_F, in eV
DIS_FROZ_LO=4.0   DIS_FROZ_HI=1.5
KMESH_PRIMITIVE=12 12 12              divided by the supercell multiplicity
NBND_WINDOW=10  NBND_MARGIN=10        sets NBND_NSCF per structure
KPATH=                                empty: no band plot
QE_BIN=  W90_BIN=  SRUN_OPTS=  W90_SRUN_OPTS=  OMP_THREADS=
```

`Si:s` is not optional: Si 3s sits inside the outer window around `E_F-9…-8`,
and without it the `s-p` hybridisation pushes the p functions onto the bonds.
With it, a 12-atom cell gives `num_wann = 4*(7+5+5+3+1) = 84`; without it, 80.
`Gd:d` carries the second largest weight at `E_F`, six times that of Si p, and
is likewise required.

`system.conf` is generated. Editing it in one folder is how a series stops
being comparable — change `defaults.conf` and rerun `01-02`.

---

## Checking one result

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

---

## Checking the whole series: `10_qc.py`

Pure numpy, plus matplotlib for `11_plot.py`. Everything is computed from the
primary files — the output of `07_hoppings.py` is never used — so the first pass
also cross-checks two independent `t_eff` implementations against each other.

### What is read

| file | what is taken |
|---|---|
| `*.wout` | Omega_I/D/OD, WF centres and spreads, windows, lattice, atoms |
| `*.win` | k-mesh for the interpolation check |
| `*.eig` | DFT eigenvalues: band counts in the windows, interpolation reference |
| `*_hr.dat` | `t_eff` per pair, on-site block, hermiticity, tail decay |
| `scf.out` / `nscf.out` | `E_F`, magnetic moments, total energy |
| `ICOHPLIST.lobster` | `-ICOHP` per bond |

Missing files are reported in the last section rather than silently skipped.
Spin comes from the file name (`up`/`dn`/`down`/`spn1`/`spn2`); if neither
matches, the run is treated as non-polarised.

Options: `--cutoff` (pair radius, default 4.0 Å), `--shells` (coordination
shells per pair, 3), `--bands-nk` (k-points for the interpolation check, 200;
`0` = all, which is what you want for final numbers), `--no-bands`, `--no-hr`,
`--icohp-glob` (point at LOBSTER when it lives outside the structure folder;
`{name}` expands to the folder name).

### Report sections

1. **Model quality** — `Omega_I`, `Omega`, their ratio, `Omega_I/N`, DFT band
   counts inside the frozen and outer windows, the largest WF-centre offset from
   its own atom, and the median `|E_Wannier - E_DFT|` in three regions: inside
   the frozen window (exact by construction, so a large value there means
   something is broken), below it, and above it.
2. **Spreads** — per-element means plus the four widest functions.
3. **t_eff** — folder, pair, shell, `d`, up, dn, mean, up/dn mismatch in percent.
4. **Internal control** — shells whose `d` does not change along the series. The
   `t_eff` scatter over them bounds the systematic error of the method. Read it
   with care: a fixed distance does not guarantee a fixed environment, and a
   monotonic drift there is physics, not error.
5. **On-site** — `<eps>` over the `H(R=0)` block of each element, referenced to
   `E_F`; `D_ex = eps_dn - eps_up` computed **per atom** and then averaged, since
   in an AFM cell the two sublattices carry opposite splittings and averaging
   `eps` first cancels them. `mean|D_ex|` is the exchange splitting; `<D_ex>`
   near zero beside it is the AFM signature. `N_wf` per element is a hard check
   on the WF-to-atom assignment: 12 for Gd, 5 for the TM, 4 for Si.
6. **Moments, energies, -ICOHP** — moments resolved per element through the
   `.wout` atom order, so the rare-earth moment does not mask the 3d one.
7. **Fit** — `t_eff = t0*exp(-beta*d)` on the nearest shell, the equivalent
   power-law exponent `p = beta*<d>` for comparison with Harrison (3.5), and
   `beta(ICOHP)/(2*beta)` as a consistency check between the two measures. In
   the two-level limit `-ICOHP ~ t^2/dE`, so a value near 1 means the basis-
   dependent and basis-independent measures agree.
8. **Warnings** — the flags below.
9. **Could not be read** — missing files, ambiguous `.wout`, parse failures.

### Flags

| tag | condition |
|---|---|
| `[ratio]` | `Omega/Omega_I` outside 1.10-1.50 |
| `[series]` | ratio departs from the series median by more than 8 % — the point is not comparable |
| `[centre]` | a function sits further than 0.90 Å from its nearest atom — the atom-block split behind `t_eff` stops being meaningful |
| `[window]` | the frozen-window band count drifts by more than 4 across the series |
| `[bands]` | median `\|E_W - E_DFT\|` inside the frozen window above 5 meV |
| `[hermiticity]` | `max \|H(R) - H(-R)+\|` above 1e-6 eV |
| `[tails]` | distant `H(R)` above 30 % of the near ones |
| `[up/dn]` | the moment is quenched, yet `t_eff(up)` and `t_eff(dn)` differ by more than 5 % |
| `[reference]` | `t_eff` at a fixed `d` scatters by more than 5 % |
| `[convergence]` | no convergence marker in `.wout` |

All thresholds are constants at the top of `10_qc.py`, one per line.

### The figure

```bash
./bin/11_plot.py qc.csv --pair Mn-Si --shell 1 --out fig_hyb.pdf
```

Top panel: `t_eff(d)` with the exponential fit, up and dn shown separately, and
`-ICOHP(d)` on a second axis when LOBSTER data is present. Bottom panel: the
local moment of the metal in that pair and the on-site exchange splitting. A
dotted vertical line marks the moment threshold, drawn only when the series
actually contains a point below 0.5 muB.

The pair is guessed from the CSV if `--pair` is omitted; the printed table under
the figure is the data behind it.

---

## Conventions that have already cost time

**QE 7.5 and underscores in species names.** `Mn_up`/`Mn_dn` break the DFT+U
restart from XML: `pw2wannier90` and `projwfc` segfault in `offset_atom_wfc`
while the SCF itself passes. Use `Mn1`/`Mn2`. `rename_species.sh` repairs a
finished calculation textually, with no recomputation.

**LOBSTER 5.1.1** refuses to run unless `scf.in` contains the literal line
`wf_collect`. Adding the line is enough.

**Wannier90 on the head node.** Two earlier attempts died there: six processes
gave 240-380 s per iteration against 13.9 s alone, with OpenMP grabbing about 15
threads each. Always a compute node, one structure at a time.

**The minimal basis `M:d;Si:p`** (32 functions) gives `Omega_I/N = 5.73 Å^2` and
`16.3` in total — the functions smear over several atoms, and the numbers look
plausible while meaning nothing. Early GdFeSi work used that basis, so its
`t_eff` values are not comparable with current ones.

**Off-centre Si functions are not a defect.** `t_eff` is the Frobenius norm of a
block and is invariant under mixing of orbitals within one atom; `s-p`
hybridisation on Si is exactly such a mixing. Directed hybrids therefore do not
affect `t_eff`, as long as every function stays nearest to its own atom. Section
5 of the report checks that through `N_wf` per element.

**`.win` and the overlaps must agree.** `num_bands` and `num_wann` in the `.win`
must match the headers of `.eig` and `.amn`. Rebuilding the `.win` after the
overlaps gives `param_read: mismatch`. Steps 04 and 06 check this and explain it.

**`ICOHPLIST.lobster` sign and multiplicity.** The file column is `ICOHP`,
negative for a bonding contact; `10_qc.py` flips the sign to report `-ICOHP`
positive. One row per symmetry-equivalent bond, so equivalent bonds are
averaged, not summed. On-site rows at `d = 0` are dropped. Only the first ICOHP
column is read, matching `06_collect.py`; if the file turns out to be
spin-resolved with two value columns, absolute values would need doubling, but
the slope `beta` is unaffected either way. Worth checking one line by eye once:

```bash
head -2 lobster/z2074/ICOHPLIST.lobster
awk 'NR==2{print NF, $8, $9}' lobster/z2074/ICOHPLIST.lobster
```

**`t_eff` converges before the spreads do.** Lengthening the minimisation for one
structure moved `Omega` by 1.1 % and redistributed the individual spreads
substantially (the widest Gd functions went from 5.58 to 3.51 Å^2) while `t_eff`
moved by 0.5 %. Convergence of `Omega` is not the criterion; convergence of the
observable is.

**Atom order.** Per-element moments assume `ATOMIC_POSITIONS` and the `.win`
atoms block list atoms in the same order. That is what `02_inputs.py` produces,
but confirm it once against `scf.out`.

---

## Speed

Reading a `_hr.dat` takes roughly 1.5 s per 100 MB, with memory on the order of
the Hamiltonian itself (84 functions × 250 R ≈ 30 MB; 1183 R ≈ 130 MB). The band
interpolation check at 200 k-points takes a fraction of a second, and at the full
864 k-points a few seconds. A six-folder series in two spin channels runs in a
minute or two.