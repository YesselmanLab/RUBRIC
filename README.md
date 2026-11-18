# RUBRIC (Rapid Unified Base-pair, Ribose and backbone Integrity Check)

[![PyPI](https://img.shields.io/pypi/v/rubric-rna?label=pypi%20package)](https://pypi.org/project/rubric-rna/)
[![CI](https://github.com/YesselmanLab/RUBRIC/actions/workflows/ci.yml/badge.svg)](https://github.com/YesselmanLab/RUBRIC/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24324565.svg)](https://doi.org/10.5281/zenodo.24324565)

RUBRIC annotates RNA 3D structures and assesses their local quality. It classifies
base pairs in the Leontis-Westhof scheme, flags ambiguous pairs, scores each pair
from 0 to 100 against high-quality structures, and flags problem backbone suites
and riboses. Results are written to JSON, and a PyMOL plug-in shows them on the
structure. A complete human ribosome takes under a second.

The data and code behind the analyses in the RUBRIC paper are at
https://github.com/YesselmanLabPublications/2026_Rubric_paper.

---

## 1. Install

| | Needed |
|---|---|
| Computer | macOS 11+, Linux x86_64 or ARM, Windows 10+ 64-bit |
| Python | 3.8 or newer |
| PyMOL | only for `rubric-view` |

### With pip (recommended)

```bash
pip install rubric-rna
```

This installs two commands, `rubric` and the optional PyMOL viewer
`rubric-view`. Check the installation with

```bash
rubric --version
```

Upgrade with `pip install --upgrade rubric-rna` and remove with
`pip uninstall rubric-rna`.

### Latest version from GitHub

```bash
pip install git+https://github.com/YesselmanLab/RUBRIC.git
```

This compiles RUBRIC on your computer, so it needs a C++20 compiler and git (see
the requirements in section 2). pip downloads CMake by itself.

### With conda (macOS and Linux)

The conda recipe in `conda/` builds the same `rubric` and `rubric-view` commands
with [rattler-build](https://rattler.build), with no Python dependency:

```bash
git clone https://github.com/YesselmanLab/RUBRIC.git
cd RUBRIC
conda install -c conda-forge rattler-build
rattler-build build --recipe conda/recipe.yaml
conda install -c ./output rubric-rna
```

### Windows

curl, which RUBRIC uses to download PDB entries, is part of Windows 10 and 11.
`rubric-view` needs PyMOL and the bash shell from
[Git for Windows](https://gitforwindows.org); the `rubric` command needs
neither.

Installed with pip or conda? Skip section 2 and go to [section 3](#3-run).

---

## 2. Build from source (for development)

Building from source is only needed to change the code. To use RUBRIC, install
it with pip (section 1).

### Requirements

You need a C++ toolchain, CMake, git, and an internet connection for the first
build. No databases or other tools to install.

| Requirement | Notes |
|---|---|
| **C++20 compiler** | GCC 11 or newer, Clang (Xcode 16 or newer on macOS), or Visual Studio 2022 |
| **CMake ≥ 3.28** | check with `cmake --version` |
| **git** | used to fetch a dependency during the first build |
| **curl** | used at runtime to download PDB entries from RCSB (preinstalled on macOS and most Linux) |

All other dependencies are bundled in the repository (`nlohmann/json` and the
reference scoring data) or downloaded by CMake during the first build (GEMMI, the
structure parser), so nothing else needs to be installed.

### Build

```bash
# 1. clone
git clone https://github.com/YesselmanLab/RUBRIC.git
cd RUBRIC

# 2. configure (Release = optimized)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# 3. build (uses all cores)
cmake --build build -j
```

The first configure/build takes a few minutes because it fetches and compiles
GEMMI; later builds are fast. When it finishes, the executable is at
`build/rubric`, and the reference scoring data is staged next to it
automatically - the tool finds it with no extra setup.

### Install the build (run `rubric` from anywhere)

So you don't have to type `./build/rubric` each time:

```bash
cmake --install build --prefix ~/.local
```

This installs `rubric` (and the optional `rubric-view` viewer) into
`~/.local/bin`, and the reference data into `~/.local/share/rubric`. If
`~/.local/bin` is on your `PATH` (it usually is), you can now run **`rubric
1GID`** and **`rubric-view 1GID`** from any directory.

If they're not found afterwards, add `~/.local/bin` to your PATH:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc   # or ~/.bashrc
```

(You can also install system-wide with `sudo cmake --install build --prefix
/usr/local`.) The examples below use the installed `rubric`; if you skip the
install, just prefix them with `./build/`.

---

## 3. Run

The basic form is `rubric <input> [options]`, where `<input>` is either a
4-character **PDB ID** (downloaded from RCSB) or a path to a local `.cif` / `.pdb`
(`.gz` accepted).

```bash
rubric 1EHZ                          # a PDB ID
rubric 1EHZ --out                    # save the result to ./1EHZ.json
rubric path/to/mymodel.cif --out     # a local file
```

### Options

| Option | Effect |
|---|---|
| `--out [FILE]` | write JSON to a file (default: stdout). Bare `--out` → `<PDB_id>.json` in the current directory; any name or path you give is used exactly as-is |
| `--pml [FILE]` | also write a PyMOL script that shows the flagged pairs and pucker outliers (open with `pymol FILE`; no plugin needed). Bare `--pml` → `<PDB_id>.pml` in the current directory; any name or path you give is used exactly as-is. Needs scoring (not with `--no-score`) |
| `--details` | include extra per-pair fields (template RMSD, score breakdown) |
| `--no-score` | find + classify only, skip quality scoring |
| `--no-download` | never fetch from RCSB; only use a local file |
| `--cache-dir DIR` | where downloaded structures are cached (default `$RUBRIC_CACHE_DIR` or `~/.cache/rubric`) |

---

## 4. Output

RUBRIC writes one JSON object per structure:

| field | what it holds |
|---|---|
| `pairs`, `n_pairs` | the base pairs (example below) |
| `pairs_score` | mean score of the scored pairs (0 to 100) |
| `tier_summary` | number of pairs per tier: Preferred, Secondary review, Priority review |
| `backbone_suiteness` | mean suiteness of the backbone (0 to 1) |
| `backbone_residues` | the flagged suites, with their tier and the torsions to move |
| `pucker_outliers` | the flagged riboses |
| `clashes`, `n_clashes` | atom pairs that overlap |

One entry of `pairs` (shortened):

```json
{
  "res_id1": "A-A-146", "res_id2": "A-U-157",
  "classification": { "lw_class": "cWW", "is_ambiguous": false },
  "score": 98.2,
  "tier": "Secondary review",
  "issues": ["distance"]
}
```

Complete outputs for 1GID:
[`1GID.json`](https://github.com/YesselmanLab/RUBRIC/blob/main/examples/1GID.json) (`rubric 1GID`) and
[`1GID_details.json`](https://github.com/YesselmanLab/RUBRIC/blob/main/examples/1GID_details.json) (`rubric 1GID --details`).

---

## 5. PyMOL viewer

With PyMOL installed, `rubric-view` opens a structure, scores it and highlights
the pairs that need attention:

```bash
rubric-view 1GID              # a PDB id
rubric-view file.cif          # or a local file
rubric-view 1GID secondary    # which pairs to show: priority (default), secondary or all
rubric-view 1GID backbone     # open in the backbone view
```

### Commands inside PyMOL

A *suite* is the stretch of backbone between two nucleotides.

**Score the structure**

| command | what it does |
|---|---|
| `rubric_score` | Scores the open structure. It turns gray, base pairs scoring below 80 are drawn as red sticks, and riboses with an inconsistent pucker are cyan. |
| `rubric_score 1GID, secondary` | Draws the pairs scoring between 80 and 100 instead, from yellow (near 100) to orange (near 80). |
| `rubric_score 1GID, all` | Draws both: the pairs below 80 in red, and the pairs between 80 and 100 in yellow to orange. |
| `rubric_overview` | Zooms back out to the whole structure. |
| `rubric_clear` | Removes all RUBRIC colors and objects. |

"1GID" stands for the name of your structure in PyMOL. The pairs drawn are the ones `rubric_list` and `rubric_next` go through.

**Go through the flagged base pairs**

| command | what it does |
|---|---|
| `rubric_list` | Prints the flagged pairs, grouped by chain, worst first (25 by default; `rubric_list 50` for more). |
| `rubric_next`, `rubric_prev` (or the → and ← keys) | Zooms in on the next or previous pair in that list and prints its score and the measurements that are off. |
| `rubric_goto 12`, `rubric_goto A-169` | Jumps to pair 12 of the list, or to the pair of any residue, flagged or not. |
| click a residue, then `rubric_info` (or Ctrl-I) | Jumps to the flagged pair that residue belongs to. |
| `rubric_ideal` | Overlays the ideal shape of the current pair in green, to show how it should look. Stays on for the next pairs until `rubric_ideal off`. |

**Riboses and backbone**

| command | what it does |
|---|---|
| `rubric_pucker` | Lists the riboses with an inconsistent pucker (the cyan ones). |
| `rubric_pucker 1`, `rubric_pucker A-199` | Zooms in on ribose 1 of that list, or on the ribose of chain A, residue 199, and draws the measurement behind the flag. `rubric_pucker off` hides the cyan. |
| `rubric_backbone` | Switches to the backbone view: pair colors are hidden and flagged suites are colored (yellow: fixable, orange: review, red: fits no known conformation). Prints the list. `rubric_backbone fixable` (or `review`, `flag_only`) lists one group only. |
| `rubric_backbone 1`, `rubric_backbone A-199` | Zooms in on suite 1 of that list, or on the suite of chain A, residue 199, and labels each torsion with its current and target value. |
| `rubric_backbone off` | Returns to the base-pair view. |

**Save and share**

| command | what it does |
|---|---|
| `rubric_dump` | Saves the results to `<structure>_pairs.json`. |
| `rubric_load file.json` | Opens saved results without running RUBRIC again, for example on another computer (open the structure first). |
| `rubric_keys off` | Turns off the arrow-key shortcuts (`rubric_keys on` turns them back on). |

The colors, what the ribose and backbone flags mean, loading the plug-in into a PyMOL
session that is already open, sharing a scored session and the `--pml` script are
described in [`integrations/pymol/README.md`](https://github.com/YesselmanLab/RUBRIC/blob/main/integrations/pymol/README.md).

---

## 6. Troubleshooting

| Problem | Fix |
|---|---|
| `rubric: command not found` after `pip install` | Python's scripts folder is not on your PATH. Install inside a virtual environment, or run `python -m rubric_rna 1GID` |
| `pip install` compiles RUBRIC instead of downloading it | no ready-made package exists for your system, so pip builds from source; install a C++20 compiler and git (section 2) |
| `download failed` / `HTTP 404` | the PDB ID is not on RCSB, or there is no network; give a local `.cif` / `.pdb` file instead |
| `rubric-view needs bash` on Windows | install [Git for Windows](https://gitforwindows.org), which provides bash |

---

## 7. License

RUBRIC is released under the [MIT License](https://github.com/YesselmanLab/RUBRIC/blob/main/LICENSE).
