# LOBSTER toolkit

Stages, each runnable on its own:

| stage | script | what it does |
|---|---|---|
| 1 | `01_make_scf.py` | copies geometry from `SRC_DIR/<point>/scf.in`, forces `nosym`/`noinv`/`nbnd`, replaces the automatic mesh by the full explicit k-list, writes `SERIES_DIR/<point>/scf.in`; warns about ultrasoft pseudos |
| 2 | `02_make_lobsterin.py` | writes `lobsterin` per point from the species in `scf.in`; reports basis-function count against `NBND` |
| 3 | `03_run_scf.sh` | `pw.x` per point, skips points already converged |
| 4 | `04_run_lobster.sh` | LOBSTER per point, prints the spilling lines |
| 5 | `05_run_loposter.sh` | LOPOSTER per point |
| 6 | `06_collect.py` | `icohp_series.csv` across the series |
| 7 | `07_check_quality.sh` | spilling, electron count, size of the negative pDOS excursions |
| 8 | `08_plot_cohp.py` | `cohp_<point>.pdf/.dat` + `icohp_<point>.tex` for ONE point |
| 9 | `09_summary.py` | scans every folder with an ICOHPLIST and writes the series tables |

```
chmod +x *.sh *.py
./run_all.sh          # everything
./run_all.sh 4        # only stage 4
./run_all.sh 2 5      # stages 2..5
```

Every script takes the config as an optional first argument, so a second
config file (different mesh, different points) needs no edits to the scripts.

## Before the first run

* **Pseudopotentials.** Stage 1 greps `pseudo_type` in each UPF. An
  ultrasoft one means the SCF must be redone with PAW or norm-conserving
  before LOBSTER is worth trying.
* **`EXTRA_LINES`.** The LOBSTER distribution ships an example `lobsterin`
  for QE. If it contains keywords that point LOBSTER at the QE input/output
  files, put them in `EXTRA_LINES` (semicolon-separated) — the generated
  `lobsterin` does not guess them.
* **Basis functions.** `BASIS` at the top of `02_make_lobsterin.py` lists the
  valence shells per element. They must match the valence of the
  pseudopotential; edit there, not in `lobsterin`.
* **`NBND`.** For the 12-atom AFM-G cell the basis is ~124 functions, so 170
  is comfortable. Stage 2 prints the count and flags it if `NBND` is short.
* **Disk.** Explicit k-list plus `nosym` means the `.save` holds every
  k-point: expect tens of GB per point. Check quota before stage 3.

Stages 7 and 8 are not part of `run_all.sh` -- run them by hand:

```
bash 07_check_quality.sh
python3 08_plot_cohp.py lobster.conf z2074 Mn-Si Mn-Mn Gd-Mn
for p in $POINTS; do python3 08_plot_cohp.py lobster.conf $p --species; done
python3 09_summary.py
```

Stage 8 handles one point per call by design (it draws one multi-panel
figure); loop over the points as above. Stage 9 needs no POINTS list -- it
discovers every folder that has an ICOHPLIST.lobster, so new points appear
automatically.

Stage 8 infers the COHPCAR column layout from the column count and prints
what it inferred. Check that line once against your file.

## After the run

Read `lobsterout` first. Charge spilling under ~3 % is good, over ~5–8 %
means the local basis did not reproduce the plane-wave wavefunction and
nothing downstream is trustworthy. With Gd 4f expect this to be the
sensitive number.

Then:

* `ICOHPLIST.lobster` — pair, distance, ICOHP per spin (stage 6 sums the
  spins into the CSV)
* `COHPCAR.lobster`, `DOSCAR.lobster` — curves; the projected DOS is a good
  cross-check against `projwfc.x`
* `CHARGE.lobster`, `GROSSPOP.lobster` — charges and orbital populations
* `RealSpaceHamiltonians.lobster` — the tight-binding Hamiltonian, readable
  the same way as `wannier90_hr.dat`; note the LOBSTER basis is
  non-orthogonal, so use the Löwdin-orthogonalised variant before comparing
  hoppings with Wannier values

## Suggested order of work

Run one cheap point first (small `POINTS`, `KMESH=4 4 2`) end to end. If the
spilling is fine and the DOS matches `projwfc`, raise the mesh to `6 6 3`,
then `8 8 4`, and add the rest of the series.
