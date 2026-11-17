# Regression panel

Ten structures from the Protein Data Bank (CC0), each chosen for a case that is hard
to get right. `1GID.cif.gz` sits one folder up, in `tests/data/`. The stored results
are in `expected/`; `tests/regression/test_panel.cpp` compares against them.

| Entry | Method | What it covers |
|---|---|---|
| 1GID | X-ray | Group I intron domain: ambiguous pairs; Fixable, Review and Flag-only suites; pucker outliers and unclassified puckers |
| 1EHZ | X-ray | tRNA-Phe with many modified nucleotides |
| 1AQO | NMR | Several models (only the first is analyzed); priority-review pairs |
| 301D | X-ray | Insertion codes |
| 3IBK | X-ray | RNA G-quadruplex |
| 3P4D | X-ray | 2'-fluoro RNA: no suites (no 2'-O), and its pairs are typed as DNA and not scored |
| 406D | X-ray | A file with four models |
| 8UPT | Cryo-EM | Insertion codes; many priority-review pairs |
| 8YDC | X-ray | Hammerhead ribozyme with a DNA strand and negative residue numbers |
| 9C6K | Cryo-EM | Larger cryo-EM model (497 nucleotides) |

After an intended change to RUBRIC's results, rewrite the stored results with

    RUBRIC_UPDATE_PANEL=1 ctest --test-dir build -R test_panel

and review the change in `git diff` before committing it.
