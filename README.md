# RUBRIC (Rapid Unified Base-pair, Ribose and backbone Integrity Check)

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
| **Internet** | first build downloads & builds GEMMI (the structure parser); the tool also downloads PDB entries at runtime unless you pass a local file |

All other dependencies are bundled in the repo (`nlohmann/json`, the reference
scoring data) or fetched automatically by CMake (GEMMI). There is nothing to
`apt install` / `brew install` beyond a compiler + CMake.

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

### Where the output goes

| you run | JSON goes to |
|---|---|
| no `--out` | **stdout** (redirect with `> file.json` if you like) |
| `--out` (no name) | **`<PDB_id>.json` in the current directory** - e.g. `rubric 1GID --out` → `./1GID.json` |
| `--out name.json` (bare name) | **`./name.json`** (current directory) |
| `--out path/to/name.json` (has a directory) | exactly that path (created if needed) |

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

A JSON object per structure:

```jsonc
{
  "pdb_id": "1EHZ",
  "candidates_total": 1234,      // pairs considered
  "candidates_valid": 200,       // passed geometric validation
  "n_pairs": 31,                 // final selected base pairs
  "pairs_score": 98.46,          // base-pair geometry quality (0-100)
  "backbone_suiteness": 0.7,     // mean Richardson suiteness (0-1, two decimals as in wwPDB reports)
  "pairs": [
    {
      "res_id1": "A-G-1", "res_id2": "A-C-72",
      "sequence": "GC",
      "pair_category": "base-base",
      "classification": { "lw_class": "cWW", "is_ambiguous": false },
      "geometry":   { "dorg": 2.40, "d_v": 0.11, "plane_angle": 24.08, "dNN": 8.42 },
      "rigid_body": { "shear": 2.22, "stretch": 0.89, "stagger": -0.11, "buckle": -12.78, "propeller": 20.40, "opening": 0.33 },
      "num_hbonds": 2,
      "score": 95.12              // per-pair quality (0-100)
    }
    // ...
  ]
}
```

---

## 5. PyMOL viewer

With PyMOL installed, `rubric-view` opens a structure, scores it and highlights
the pairs that need attention:

```bash
rubric-view 1GID              # a PDB id
rubric-view file.cif          # or a local file
rubric-view 1GID secondary    # optional tier: priority (default) | secondary | all
```

### Commands inside PyMOL

| command | what it does |
|---|---|
| `rubric_score [obj], [tier]` | score the structure + highlight the flagged pairs (colored by pair score) and the pucker-mismatch sugars (cyan); `tier` = `priority` (default) / `secondary` / `all` |
| `rubric_list [n]` | print the ranked worklist of the worst pairs (default 25) |
| `rubric_next` / `rubric_prev` | step through the worklist one pair at a time |
| `rubric_goto <rank \| residue>` | jump to worklist entry `12`, **or** inspect any residue by `A-169` / `A169` / `A-G-169` (even a clean pair not in the worklist) - prints its score and what's out of range |
| `rubric_info [selection]` | inspect a clicked/selected residue → jump the worklist to its pair |
| `rubric_overview` | zoom back out to the whole-structure overview |
| `rubric_ideal [on\|off]` | overlay the **idealized** target geometry of the current pair (green ghost) so you can see how to fix it; sticky toggle, default off |
| `rubric_pucker [list \| on \| off \| n \| residue]` | list the sugar-pucker mismatches, show/hide the cyan layer, or a close-up of one (`rubric_pucker 1`, `rubric_pucker A-199`) with its P-perp measurement |
| `rubric_dump [path]` | write the scored pairs (incl. any coordinate edits) to `<obj>_pairs.json` in the current directory, or a path you give |
| `rubric_load <path>` | rebuild the worklist from a `rubric_dump` JSON - **no binary, no re-scoring** (for sharing a scored session) |
| `rubric_keys [on\|off]` | bind/unbind arrow-key navigation |
| `rubric_clear` | remove every overlay and reset |

Colors, keyboard navigation and the `--pml` script are described in
[`integrations/pymol/README.md`](https://github.com/YesselmanLab/RUBRIC/blob/main/integrations/pymol/README.md).

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
