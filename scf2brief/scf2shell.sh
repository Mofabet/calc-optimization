#!/bin/bash
pwd
grep -E "!    total energy|magn=|Fermi|unit-cell|lattice parameter" scf.out
grep -B3 "convergence has been achieved" scf.out