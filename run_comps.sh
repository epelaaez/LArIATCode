#!/bin/bash

# Run macros to get data-MC comparisons in batch mode (-b).

root -l -b -q RecoClassify3Cat.cc
root -l -b -q RecoDataAnalysis.cc