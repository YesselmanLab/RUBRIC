# Validation of RUBRIC

This summarizes the evidence that RUBRIC behaves as intended on real structures, and
its known limitations. The numbers come from the RUBRIC paper; the data and analysis
code that produce them are in https://github.com/YesselmanLabPublications/2026_Rubric_paper
(engine version `rubric-rna==1.0.0`).

## Agreement with established annotation programs

RUBRIC was compared with DSSR, FR3D, RNAView, barnaba, RNApolis Annotator and BPNet
on the BGSU non-redundant set.

- **Canonical pairs:** RUBRIC's detection agreed with each program with a Jaccard
  index of 0.96-0.98.
- **Non-canonical pairs:** RUBRIC's Jaccard index with each program was 0.54-0.77,
  within or above the range among the established programs themselves (0.42-0.70).
  For class assignment, RUBRIC's Cohen's kappa with each program was 0.76-0.98
  (established programs among themselves: 0.71-0.97).
- **Against DSSR, class by class:** cWW recall 99.7% and precision 99.8%; the other
  eleven classes 85-98%, except recall of tWS (55%), recall of cHH, cSS, cHS and cWS
  (61-76%) and precision of tSS (75%). Where RUBRIC and DSSR assigned different classes
  to the same pair (13,286 pairs), four of the five other programs supported RUBRIC's
  class significantly more often than DSSR's, by 13-41 percentage points.

## Ambiguous pairs

Pairs that RUBRIC flags as ambiguous are the pairs other programs disagree on: among
pairs classified by at least two of the six other programs, they assigned different
classes to 61% of flagged pairs, but to 23% of unflagged pairs of the same classes.

## Flagged problems are correctable

Problems RUBRIC flags were followed into independent, improved models of the same
molecule:

- **Same molecule at higher resolution (162 molecules):** the mean score of flagged
  base pairs increased in 113 of the 115 molecules that had priority-review pairs;
  backbone suiteness increased in 73% of molecules.
- **PDB-REDO re-refinement of the same X-ray data (1,521 structures):** flagged pairs
  improved in 89% of the 364 structures that contained them.

## Backbone and riboses

RUBRIC's structure suiteness reproduces the wwPDB validation reports (R² = 0.992 over
2,604 structures). Suites, conformers and the delta pucker bands follow Suitename.

## Speed

RUBRIC analyzes a complete human 80S ribosome (about 5,500-6,000 nucleotides) in
0.58-0.76 s. Its advantage over DSSR, RNAView and BPNet grows with RNA size, reaching
11- to 24-fold above 5,000 nucleotides. It was faster than FR3D, barnaba and RNApolis
Annotator on every structure they completed.

## Continuous checks

- The test suite (`tests/README.md`) runs on Linux, macOS and Windows on every change,
  including a regression panel of ten real structures chosen for hard cases.
- RUBRIC-PDB-live (https://yesselmanlab.github.io/RUBRIC-PDB-live/) runs every
  RNA-containing PDB entry each week; the latest full run scored all 10,571 entries
  without a failure.

## Known limitations

- **Pairs need a base-base hydrogen bond.** RUBRIC classifies only pairs with at least
  one hydrogen bond between base atoms; a pair held only by a base-to-ribose hydrogen
  bond (counted by DSSR) is not reported.
- **Not every pair is scored.** Ambiguous pairs, pairs involving DNA, and pairs
  involving nucleotides that the chemical component dictionary types as DNA (for
  example 2'-fluoro nucleotides) are reported without a score.
- **Suites need a 2'-O.** A suite is assigned only when the nucleotide and the one
  before it have a C2' bonded to an O2'; DNA, 2'-deoxy and 2'-fluoro nucleotides have
  no suite.
- **First model only.** In a file with several models (NMR ensembles), only the first
  is analyzed.
- **Geometry only.** RUBRIC assesses model geometry, not the fit to the experimental
  data.
- **Calibration.** Scores are calibrated on curated X-ray and cryo-EM RNA structures;
  NMR structures are scored on the same scale but were not used in calibration.
- **Class boundaries.** Near the tSS/tWS boundary RUBRIC assigns some pairs to tSS
  that DSSR assigns to tWS, a boundary on which annotation programs are known to
  disagree.
- **Command line.** Give the input before the options (`rubric <input> [options]`).
  When an unrecognized option comes first, the error message names the wrong argument.
