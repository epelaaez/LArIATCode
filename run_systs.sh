#!/bin/bash

# Run macros to compute systematics in batch mode (-b).

root -l -b -q MCStatCov.cc
root -l -b -q GeneratorCov.cc
root -l -b -q EnergyRecoCov.cc
root -l -b -q BeamlineMuCov.cc
root -l -b -q BeamlineElCov.cc
root -l -b -q RecoClassify3Cat.cc
root -l -b -q Unfold.cc
root -l -b -q UnfoldFD.cc