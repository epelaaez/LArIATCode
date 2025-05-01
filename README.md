# Pion Absorption in LArIAT

## LArSoft scripts

The following scripts are available in a separate [GitHub repo](https://github.com/epelaaez/lariatsoft/tree/epelaez/pion_absorption). These scripts run on the generated events to create ROOT histograms that are analyzed with ROOT scripts later on.

Reconstruction modules:
- `RecoEval_module`: rejects non pion events.
- `RecoAllEval_module`: does not reject any event, but stores signal and background information.
- `ShowerRecoEval_module`: analyzes events with electromagnetic showers.

Filter modules:
- `PionAbsorptionSelection_module`: performs reconstruction-level selection and saves truth-level data about events that pass selection.
- `SignalPionabsorption_module`: accepts event passing truth-level signal definition.

## ROOT scripts
- `RecoAnalysis`
- `SelectionAnalysis`
- `ShowerAnalysis`
- `SignalAnalysis`