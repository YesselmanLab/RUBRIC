# parse

A fast, self-contained C++ tool that reads an RNA structure file (mmCIF or PDB),
**finds** base pairs, **classifies** them (Leontis–Westhof), and **scores** their
geometric quality — emitting per-pair and whole-structure quality as JSON.

A full 70S ribosome is annotated and scored in **under a second**.

---

## 1. Requirements

You need a C++ toolchain, CMake, git, and an internet connection for the first
build. No databases or other tools to install.

| Requirement | Notes |
|---|---|
| **C++20 compiler** | Apple Clang ≥ 14, Clang ≥ 13, or GCC ≥ 10 |
| **CMake ≥ 3.28** | check with `cmake --version` |
| **git** | used to fetch a dependency during the first build |
| **curl** | used at runtime to download PDB entries from RCSB (preinstalled on macOS and most Linux) |
| **Internet** | first build downloads & builds GEMMI (the structure parser); the tool also downloads PDB entries at runtime unless you pass a local file |

All other dependencies are bundled in the repo (`nlohmann/json`, the reference
scoring data) or fetched automatically by CMake (GEMMI). There is nothing to
`apt install` / `brew install` beyond a compiler + CMake.

---

## 2. Build

```bash
# 1. clone
git clone https://github.com/YesselmanLab/FastRNA.git
cd FastRNA

# 2. configure (Release = optimized)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# 3. build (uses all cores)
cmake --build build -j
```

The first configure/build takes a few minutes because it fetches and compiles
GEMMI; later builds are fast. When it finishes, the executable is at
`build/parse`, and the reference scoring data is staged next to it
automatically — the tool finds it with no extra setup.

### Install (recommended — run `parse` from anywhere)

So you don't have to type `./build/parse` each time:

```bash
cmake --install build --prefix ~/.local
```

This installs `parse` (and the optional `parse-view` viewer) into
`~/.local/bin`, and the reference data into `~/.local/share/parse`. If
`~/.local/bin` is on your `PATH` (it usually is), you can now run **`parse
1GID`** and **`parse-view 1GID`** from any directory.

If they're not found afterwards, add `~/.local/bin` to your PATH:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc   # or ~/.bashrc
```

(You can also install system-wide with `sudo cmake --install build --prefix
/usr/local`.) The examples below use the installed `parse`; if you skip the
install, just prefix them with `./build/`.

---

## 3. Run

The basic form is `parse <input> [options]`, where `<input>` is either a
4-character **PDB ID** (downloaded from RCSB and cached) or a path to a local
`.cif` / `.pdb` (`.gz` accepted).

```bash
# by PDB ID (downloaded + cached under ~/.cache/parse)
parse 1EHZ

# save to a file — bare --out uses the PDB id → ./1EHZ.json
parse 1EHZ --out

# from a local file (no download) — just the name if it's in the current
# folder, or its path if it's elsewhere
parse mymodel.cif --out              # file in the current directory
parse path/to/mymodel.cif --out      # file in another directory

# a full ribosome — still sub-second
parse 7K00 --out
```

### Where the output goes

| you run | JSON goes to |
|---|---|
| no `--out` | **stdout** (redirect with `> file.json` if you like) |
| `--out` (no name) | **`<PDB_id>.json` in the current directory** — e.g. `parse 1GID --out` → `./1GID.json` |
| `--out name.json` (bare name) | **`./name.json`** (current directory) |
| `--out path/to/name.json` (has a directory) | exactly that path (created if needed) |

The output goes exactly where you point it — no folder is imposed.

Downloaded structures are cached in **`~/.cache/parse/`** — nothing is
written into your working directory unless you ask for it.

### Options

| Option | Effect |
|---|---|
| `--out [FILE]` | write JSON to a file (default: stdout). Bare `--out` → `<PDB_id>.json` in the current directory; any name or path you give is used exactly as-is |
| `--pml [FILE]` | also write a PyMOL script that shows the flagged pairs and pucker outliers (open with `pymol FILE`; no plugin needed). Bare `--pml` → `<PDB_id>.pml` in the current directory; any name or path you give is used exactly as-is. Needs scoring (not with `--no-score`) |
| `--details` | include extra per-pair fields (template RMSD, score breakdown) |
| `--no-score` | find + classify only, skip quality scoring |
| `--no-download` | never fetch from RCSB; only use a local file |
| `--cache-dir DIR` | where downloaded structures are cached (default `$PARSE_CACHE_DIR` or `~/.cache/parse`) |

---

## 4. Output

A JSON object per structure:

```jsonc
{
  "pdb_id": "1EHZ",
  "candidates_total": 1234,      // pairs considered
  "candidates_valid": 200,       // passed geometric validation
  "n_pairs": 31,                 // final selected base pairs
  "pairs_score": 98.46,          // base-pair geometry quality (0–100)
  "backbone_suiteness": 0.696,   // mean Richardson suiteness (0–1, wwPDB convention)
  "pairs": [
    {
      "res_id1": "A-G-1", "res_id2": "A-C-72",
      "sequence": "GC",
      "pair_category": "base-base",
      "classification": { "lw_class": "cWW", "is_ambiguous": false },
      "geometry":   { "dorg": 2.40, "d_v": 0.11, "plane_angle": 24.08, "dNN": 8.42 },
      "rigid_body": { "shear": 2.22, "stretch": 0.89, "stagger": -0.11, "buckle": -12.78, "propeller": 20.40, "opening": 0.33 },
      "num_hbonds": 2,
      "score": 95.12              // per-pair quality (0–100)
    }
    // ...
  ]
}
```

- **`pairs_score` and `backbone_suiteness` are reported separately** — they measure
  different quality axes, so they are not blended into a single "overall" number
  (which would over-claim whole-structure quality; clashscore, bond geometry, etc.
  are not considered here). Note the scales differ by design: `pairs_score` is a
  0–100 geometry score, while `backbone_suiteness` is mean Richardson suiteness as a
  **0–1** fraction (the wwPDB RNA-suiteness measure), reported in its native unit.
- Per-pair `score` and `pairs_score` are 0–100 (higher = better); `backbone_suiteness`
  is 0–1 (higher = better).
- A `classification.lw_class` containing `|` (e.g. `cWH|cWW`) with
  `is_ambiguous: true` means the edge assignment is genuinely undetermined.

### Backbone recommendations (`backbone_residues`)

When scoring is on, the output also carries a **`backbone_residues`** array — one
entry per residue whose backbone is flagged (residues that score fine are
omitted). The suiteness *score* is unchanged; this is a refinement-facing
decomposition of *which* suite torsions are off and *where they should go*:

```jsonc
{
  "res_id": "B-U-259",
  "suiteness": 0.443,
  "tier": "fixable",                 // fixable | review | flag_only
  "target_conformer": "1a",          // Richardson rotamer to aim for
  "deviations": [                    // only the torsions that "fire" (ProSco < 5)
    { "angle": "epsilon_prev", "value": 248.5, "target": 212.3, "gap": 36.3, "prosco": 1.41 }
  ]
}
```

- **`tier`** — `fixable`: low suiteness with ≥1 anomalous torsion → rotate the
  listed `deviations` toward `target_conformer`. `review`: low suiteness but no
  single anomalous torsion (jointly distorted → a full-rebuild candidate;
  `deviations` empty). `flag_only`: an outlier that fits no rotamer (the nearest
  conformer is only a hint).
- **`deviations`** — each suite torsion that is genuinely rare for its conformer
  (`prosco` < 5, lower = more anomalous), most-anomalous first, given as
  `value → target` (the rotation to apply) with the signed `gap` (degrees).
- A flagged torsion is a **statistical anomaly worth reviewing against the map**,
  not a guaranteed error — rare conformations do occur in good structures. The
  `prosco` value is the confidence signal (deep, e.g. < 1, is more likely a real
  error than a borderline 3–5).

---

## 5. Optional: PyMOL viewer

If you have **PyMOL** installed, you can open a structure, score it, and
highlight the pairs that need attention in one command:

```bash
parse-view 1GID              # a PDB id
parse-view file.cif          # or a local file
parse-view 1GID secondary    # optional tier: priority (default) | secondary | all
```

The optional **tier** selects which pairs populate the worklist, by whole-pair
score: `priority` (priority review, score < 80 — the genuinely-poor pairs, the default),
`secondary` (secondary review, 80 ≤ score < 100 — the minor-issue middle, a polishing pass), or
`all` (both); the older words `review` and `acceptable` still work. The structure opens as a whole-view health map — gray, with the
flagged pairs lit up — and the scored JSON is also written to
`1GID_pairs.json` in the current directory.

**Coloring:** pairs are colored by their pair
score: **gray** = preferred (score 100, not flagged), **yellow → orange** =
secondary review (80 ≤ score < 100), **red** = priority review (score < 80).
Separately, a sugar whose pucker is inconsistent (δ vs P-perp) gets a **cyan**
ribose ring, whatever its pair's tier, and `parse_ideal` draws the ideal
template of the current pair in **green**.

### Commands inside PyMOL

Once the viewer is open, drive it with these commands (type them in the PyMOL
command line). Arguments in `[...]` are optional.

| command | what it does |
|---|---|
| `parse_score [obj], [tier]` | score the structure + highlight the flagged pairs (colored by pair score) and the pucker-mismatch sugars (cyan); `tier` = `priority` (default) / `secondary` / `all` |
| `parse_list [n]` | print the ranked worklist of the worst pairs (default 25) |
| `parse_next` / `parse_prev` | step through the worklist one pair at a time |
| `parse_goto <rank \| residue>` | jump to worklist entry `12`, **or** inspect any residue by `A-169` / `A169` / `A-G-169` (even a clean pair not in the worklist) — prints its score and what's out of range |
| `parse_info [selection]` | inspect a clicked/selected residue → jump the worklist to its pair |
| `parse_overview` | zoom back out to the whole-structure overview |
| `parse_ideal [on\|off]` | overlay the **idealized** target geometry of the current pair (green ghost) so you can see how to fix it; sticky toggle, default off |
| `parse_pucker [list \| on \| off \| n \| residue]` | list the sugar-pucker mismatches, show/hide the cyan layer, or a close-up of one (`parse_pucker 1`, `parse_pucker A-199`) with its P-perp measurement |
| `parse_dump [path]` | write the scored pairs (incl. any coordinate edits) to `<obj>_pairs.json` in the current directory, or a path you give |
| `parse_load <path>` | rebuild the worklist from a `parse_dump` JSON — **no binary, no re-scoring** (for sharing a scored session) |
| `parse_keys [on\|off]` | bind/unbind arrow-key navigation |
| `parse_clear` | remove every overlay and reset |
| `parse_set_binary <path>` | point the plugin at a specific `parse` binary |
| `parse_set_ideals <path>` | set the idealized-template dir |

A typical loop: `parse_score` (priority review) → `parse_next` through the flagged
pairs → `parse_ideal` to see the target → fix geometry / re-refine → re-run
`parse_score`; when the priority-review worklist is clear, `parse_score obj, secondary`
to polish the middle. Keyboard navigation and the full color scheme are covered in
[`integrations/pymol/README.md`](integrations/pymol/README.md).

### No plugin? Write a standalone script

`parse 1GID --pml` writes `1GID.pml`, a plain PyMOL script that loads the
structure and shows the same picture (flagged pairs colored by score, cyan
pucker-mismatch riboses) with named selections for each tier and one per
priority-review pair. Open it with `pymol 1GID.pml` on any machine with PyMOL;
neither the plugin nor the binary is needed to view it. Details in
[`integrations/pymol/README.md`](integrations/pymol/README.md#without-the-plugin-a-standalone-script-parse---pml).

PyMOL is only needed for this viewer — the core `parse` tool does not require it.

---

## 6. Notes

- **Internet:** needed once (first build, to fetch GEMMI) and at runtime only
  when you give a bare PDB ID (it downloads via `curl` and caches the file).
  Local-file runs need no network — use `--no-download` to enforce that.
- **Cache:** downloaded structures live in `~/.cache/parse` (override with
  `--cache-dir` or `$PARSE_CACHE_DIR`). Safe to delete anytime.
- **Reference data** (`resources/`) ships with the repo and is staged next to the
  binary at build time; the tool finds it automatically.
- **Performance:** ~4,000 base pairs/second; a 70S ribosome (~2,100 pairs) scores
  in ~0.5 s, the largest RNA structures (~16,000 nt) in ~1.5 s.

---

## 7. Project layout

```
include/pairfinder/   public headers (geometry, core, io, algorithms, scoring)
src/                  implementations
apps/main.cpp         the CLI
resources/            bundled reference data + templates (ship with the binary)
integrations/pymol/   optional PyMOL viewer (parse-view, parse_pymol.py)
CMakeLists.txt        build (fetches GEMMI, vendors nlohmann/json)
```

---

## 8. Troubleshooting

| Symptom | Fix |
|---|---|
| `CMake 3.28 or higher is required` | upgrade CMake (`cmake --version` to check) |
| First build fails fetching GEMMI | check internet / proxy; re-run `cmake -S . -B build` |
| `download failed` / `HTTP 404` at runtime | the PDB ID isn't on RCSB, or no network — pass a local `.cif`/`.pdb`, or check `curl` is installed |
| Want a clean rebuild | `rm -rf build` and repeat the Build steps |
