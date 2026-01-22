#!/bin/bash

# Run macros to compute systematics in batch mode (-b).

root -b -q GeneratorCov.cc
root -b -q EnergyRecoCov.cc
root -b -q RecoClassify3Cat.cc
root -b -q Unfold.cc