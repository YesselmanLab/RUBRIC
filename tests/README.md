# RUBRIC tests

The tests check the engine at three levels: its building blocks (unit), the `rubric`
program as a user runs it (cli), and its results on real structures (regression).
They need nothing beyond what RUBRIC itself needs, and run offline in a few seconds.

## Running

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure        # everything
ctest --test-dir build -L unit                    # one level: unit, cli or regression
./build/tests/test_suites                         # one program, with one line per test
```

A program prints `pass` or `FAIL` for each test, with the values involved when a
check fails, and `OK` when all pass. The tests read the bundled `resources/` (ctest
sets `RUBRIC_RESOURCES`) and the files in `tests/data/` and `tests/fragments/`.

## What is tested

| Program | Level | What it checks |
|---|---|---|
| `unit/test_geometry.cpp` | unit | Vector arithmetic, rigid superposition (exact recovery of a known rotation), base reference frames, residue-ID parsing, angular gaps |
| `unit/test_suites.cpp` | unit | All 54 Richardson conformer centres, suiteness falling with distance, outliers, the delta pucker bands |
| `unit/test_backbone.cpp` | unit | Suites, the Fixable / Review / Flag-only recommendation rules and their damped targets, pucker outliers and unclassified puckers, the C2'-O2' suite rule |
| `unit/test_basepairs.cpp` | unit | H-bonds of Watson-Crick pairs and distance limits, classification of every idealized template, scores of ideal and stretched pairs, score bookkeeping, tier boundaries |
| `unit/test_chemistry.cpp` | unit | Parent bases of modified nucleotides, purine/pyrimidine/DNA typing, the modified-nucleotide table, which pairs are scored (not DNA, not ambiguous), resource loading |
| `cli/test_cli.cpp` | cli | Every option (`--out`, `--details`, `--no-score`, `--pml`, `--version`, `--help`), the JSON fields, error exits, the `dump-*` commands, and that rebuilding the resource caches reproduces the bundled ones |
| `cli/test_inputs.cpp` | cli | Negative residue numbers, insertion codes, several models, alternate locations, protein/water/ions next to RNA, files with no RNA |
| `regression/test_panel.cpp` | regression | Ten real structures (`data/panel/README.md`) against stored results |

`check_package.py` checks an installed package (the `rubric` command) against
`data/1GID.expected.json`. The release workflow runs it on every wheel, installed on
each supported Python version (3.8 to 3.15) and built from source, and the conda
build runs it on the conda package.

## What is not tested, and why

| Not tested | Reason |
|---|---|
| Every possible structure | The panel covers the hard cases; RUBRIC-PDB-live runs every RNA-containing PDB entry each week |
| Whether the science is right | Tests check that RUBRIC behaves as designed; `VALIDATION.md` summarizes the evidence that the design is right |
| The PyMOL plug-in display | Needs PyMOL; the `--pml` script it shares logic with is tested |
| Downloading by PDB ID | The tests run offline; the error path without download is tested |
| Speed | Depends on the machine; see the runtime benchmark in the paper |

## Adding a test

Add a `TEST(name) { ... }` to the program that covers the area, or a new
`tests/<level>/test_<area>.cpp` ending with `RUBRIC_TEST_MAIN()` and registered in
`tests/CMakeLists.txt`. Start each test with a one-line comment stating what it
checks and what it expects. Keep inputs small (a fragment or an idealized template,
edited at run time if needed), keep expected values away from RUBRIC's exact
cutoffs so rounding cannot flip them, and use `CHECK_NEAR` for floating-point values.

## Updating the stored panel results

After an intended change to RUBRIC's results:

```bash
RUBRIC_UPDATE_PANEL=1 ctest --test-dir build -R test_panel
git diff tests/data/panel/expected     # review every change before committing
```

## Coverage and sanitizers

```bash
cmake -S . -B build-cov -DCMAKE_BUILD_TYPE=Debug -DRUBRIC_COVERAGE=ON     # line coverage
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"                 # memory errors, undefined behavior
```

The CI runs both on every push and shows the coverage summary on the run page; the
tests cover about 87% of the engine's lines.
