#!/bin/bash

# Run each (active) ROOT macro in batch mode (-b).

root -b -q RecoNNShowers.cc
root -b -q RecoAllAnalysis.cc
root -b -q RecoClassifyAll.cc
root -b -q RecoClassifyAllSimplified.cc
root -b -q RecoDataAnalysis.cc