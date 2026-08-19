#!/bin/bash
pwd >> brief
grep -E "!    total energy|magn:|Fermi|unit-cell|lattice parameter" scf.out >> brief
grep -B3 "convergence has been achieved" scf.out >> brief