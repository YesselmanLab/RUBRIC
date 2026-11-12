# RUBRIC → PyMOL

Highlight problematic base pairs in PyMOL during refinement, using the
`rubric` binary. Flagged pairs are shown as sticks colored by their **pair
score**: red = priority review, yellow → orange = secondary review. Sugars whose pucker is inconsistent are drawn with a **cyan
ribose**. Re-run after an edit and a fixed pair (or sugar) disappears. The
original object is never modified.

## Quick start (one command)

From the repo, after building `rubric`:

```bash
./integrations/pymol/rubric-view 1GID                 # a PDB id (fetched from RCSB)
./integrations/pymol/rubric-view /path/to/mymodel.pdb # or a local file
./integrations/pymol/rubric-view 1GID secondary       # optional tier: priority|secondary|all
```

This opens PyMOL, loads the plugin, loads the structure, scores it, and
highlights the issues — all in one go. (Override the PyMOL or rubric
executable with `PYMOL=...` / `RUBRIC_BINARY=...`; `RUBRIC_VIEW_DRYRUN=1`
prints the PyMOL commands without launching.)

To drive it interactively inside an already-open PyMOL session instead, use the
commands below.

## Setup

1. Build `rubric` (see the top-level README).
2. Tell the plugin where the binary is — either:
   - export `RUBRIC_BINARY=/path/to/build/rubric` before launching PyMOL, or
   - run `rubric_set_binary /path/to/build/rubric` inside PyMOL.
3. In PyMOL: `run /path/to/integrations/pymol/rubric_pymol.py`

## Use

```
rubric_score                 # score current object; worklist = priority review (score<80)
rubric_score 6DN1            # score object named 6DN1 (priority review)
rubric_score 6DN1, secondary # the minor-issue middle (80<=score<100), a polish pass
rubric_score 6DN1, all       # every non-Preferred pair
rubric_score 6DN1, all, 0.5  # ...and only issues with color severity >= 0.5 (less noise)
rubric_clear                 # remove the overlay
```

Each `rubric_score` run:
1. saves the object's **current** coordinates to a temp PDB,
2. runs `rubric` on it (sub-second, even on a ribosome),
3. builds a non-destructive overlay object `rubric_overlay` containing only the
   flagged pairs (sticks, colored by pair score), plus a `rubric_problems`
   selection you can `zoom rubric_problems` / iterate over,
4. draws every sugar-pucker mismatch as a gray nucleotide with a cyan ribose
   ring (object `rubric_pucker`; see [Sugar-pucker mismatches](#sugar-pucker-mismatches)),
5. prints a summary to the PyMOL log, ranked by each pair's worst-parameter
   |Z′| (ProSco-gated), followed by the pucker mismatches and a one-line color key.

**Which pairs are flagged** is gated by the whole-pair **score tier** (`tier`
arg): `priority` (score < 80, the default), `secondary` (80 ≤ score < 100), or
`all` (both); the older words `review` and `acceptable` still work. The pair score decides *membership and color* — who a refiner
reviews, and how urgently — while the ProSco/Z′ severity decides the *order within
a tier* and names each pair's worst parameter. So the
default `rubric_score` shows only the genuinely-poor pairs, and `rubric_score obj,
secondary` (or `all`) widens the net for a thorough polishing pass. The printed
header shows the structure's tier counts (e.g. `104 Preferred · 32 secondary review ·
10 priority review`) so you know how big each pass is before you run it.

A typical loop: `rubric_score` (priority review) → fix the flagged pairs → re-run; when the
priority-review worklist is empty, `rubric_score obj, secondary` to polish the middle.

## What gets colored

The whole structure is colored **gray** (unflagged pairs recede), and the flagged
pairs are drawn as colored sticks on top, **colored by the pair score** (the
validated whole-pair quantity; a residue in several flagged pairs takes its lowest).
`rubric_score` prints the key as one line:

| color | meaning | drawn as |
|---|---|---|
| **gray** | *preferred* — pair score 100 | not highlighted |
| **yellow → orange** | *secondary review* — 80 ≤ score < 100 (yellow near 100, orange near 80) | pair sticks |
| **red** | *priority review* — score < 80 | pair sticks |
| **cyan** | *pucker mismatch (δ vs P-perp)* — a flag on one sugar, not a tier | the ribose ring of a gray nucleotide |
| **green** | *ideal template* — the target geometry for the current pair | thin ghost sticks + H-bonds (`rubric_ideal`) |

Gray, yellow → orange and red are one scale, the pair score. Cyan is separate: a
pucker mismatch can sit on a pair of any tier, or on an unpaired nucleotide.

The worst parameter (e.g. `buckle 4.3σ`), its Z′ and its Černý tier still appear in
the printed summary and in the per-pair readout. Until 2026-10 the color came from
the single worst parameter's severity instead, so a pair with one large outlier in a
low-weight parameter could show red while it was in secondary review (e.g. 1GID
A226–U249: propeller 6.3σ, score 85). H-bond-count pairs are no longer drawn blue;
they are colored by their score like any other pair, and still sort to the bottom
of the worklist.

Pass `rubric_score target, priority, 0, 0, 0` (last arg `gray_background=0`) to keep
the object's existing colors instead of graying. Only colors change — coordinates
are never modified.

## Commands

| command | effect |
|---|---|
| `rubric_score [obj], [tier], [min_severity], [zoom]` | score + highlight the flagged pairs; `tier` = `priority` (default) / `secondary` / `all` |
| `rubric_list [n]` | print the ranked worklist (default 25) |
| `rubric_dump [path]` | write the current scored pairs (incl. edits) to JSON — default `<obj>_pairs.json` |
| `rubric_load <path>` | rebuild the worklist from a `rubric_dump` JSON — **no binary, no re-scoring** (for sharing a session) |
| `rubric_next` / `rubric_prev` | step through the worklist one pair at a time |
| `rubric_overview` | zoom back out to the whole-structure overview (undim + re-show all flagged pairs) |
| `rubric_goto <rank \| chain+number>` | jump to worklist entry `12`, **or** inspect any residue by `A-169` / `A169` / `A-G-169` — including a clean pair not in the worklist, printing its score + what's out of range |
| `rubric_info [selection]` | inspect a clicked/selected residue → jump the worklist to its pair |
| `rubric_ideal [on\|off]` | overlay the **idealized** version of the current pair (green ghost = target geometry) so a refiner sees how to fix it; sticky toggle, default off |
| `rubric_pucker [list \| on \| off \| n \| chain+number]` | sugar-pucker mismatches: list them (default), show/hide the cyan layer, or a close-up of one (`rubric_pucker 1`, `rubric_pucker A-199`) with its P-perp measurement |
| `rubric_keys [on\|off]` | bind/unbind the navigation keys |
| `rubric_clear` | remove every overlay + reset |
| `rubric_set_binary <path>` | set the rubric binary path |
| `rubric_set_ideals <path>` | set the idealized-template dir (default `resources/basepair-idealized`) |

## Guided navigation (keys, on-screen status, click-to-inspect)

After `rubric_score`, the integration becomes a guided triage tool — you don't have to
keep typing commands:

- **Keyboard:** **→** = next pair, **←** = previous, **Ctrl-I** = inspect the atom you
  last clicked. (Keys bind automatically; `rubric_keys off` unbinds them.)
- **On-screen status:** a one-line HUD at the bottom of the viewport always shows where
  you are — `RUBRIC 3/45  A226·U249 (chain A)  cWW AU  propeller 6.2σ [Of Concern]` — so
  you never have to scroll the log.
- **Click-to-inspect:** click any flagged stick in 3D, then `rubric_info` (or **Ctrl-I**).
  The worklist jumps to that pair and prints its full breakdown + Phenix selection.
  Click a non-flagged residue and it tells you it's Preferred/unscored.

## Worklist (fix them one at a time)

`rubric_score` lights up *every* flagged pair at once — a health map, good for seeing
where the problem regions cluster, but overwhelming on a large structure. To actually
work through them, use the worklist:

```
rubric_score target      # build the queue + global map
rubric_next              # frame the worst pair: bright sticks, faded 8 A context,
                        #   H-bonds, a label, and its distortion printed to the log
rubric_next              # ... the next one, and so on
rubric_overview          # zoom back out to the whole structure (undim, all pairs shown)
rubric_goto 1A-G-169     # jump straight to a specific pair
```

The queue is **chain-grouped, worst-first within each chain**, so stepping stays in
one region instead of teleporting across the structure. On each stop the rest of the
structure is dimmed (still visible for orientation), and the log prints the failing
parameters (e.g. `shear Z'=+11.3, opening Z'=+10.3` — which tells you *how* it's
distorted) plus a **Phenix-style selection string** for the pair.

The queue is ranked by each pair's **single most-deviant parameter's |Z′|** (the
ProSco-gated Černý severity, *uncapped* so a 6.6σ pair outranks a 5.1σ one) — not
by the 0–100 quality score, which saturates near 100 for most pairs and so makes a
poor rank key. Each entry shows **both**: the pair's score (e.g. `score 65/100`)
*and* its worst parameter (e.g. `propeller 6.2σ`). `rubric_overview` additionally
prints the whole-structure score (overall / pairs / backbone).

**Visual priority:** the worklist is a *visual* triage tool, so `hbond_angles` —
a bond-geometry (chemistry) strain that doesn't make a pair *look* distorted — is
down-weighted (×0.3) in the rank and highlight color. A visibly twisted/splayed
pair therefore outranks one whose only large deviation is H-bond angle, and the
headline names the *visible* defect. The true σ is still reported, and H-bond
*distance* keeps full weight (a long/short bond shows as a visibly splayed pair).
This affects only the PyMOL triage order/color — never the 0–100 score or scorer.

Typical loop: `rubric_next` to a problem → understand it → fix it in Phenix / edit the
coordinates → reload and re-run `rubric_score` → the fixed pair drops off the queue.

## Show the ideal pair (what it should look like)

Knowing a pair is distorted is half the job; a refiner also needs to see the
*target*. `rubric_ideal` overlays the **idealized** version of the current pair —
a perfectly planar, H-bond-optimized template for that Leontis–Westhof class and
sequence — as a thin **green ghost** with its ideal H-bonds:

```
rubric_next
rubric_ideal            # green ghost of the ideal pair appears, fitted on the real one
rubric_ideal on         # ... and stays on for every rubric_next / rubric_goto
rubric_ideal off        # hide it again
```

The ideal is **best-fit** onto the actual pair (least-squares over the shared base
atoms of *both* bases), so each base shows how far it sits from where it should be
— flatten/untwist the real (colored) sticks toward the green ghost. Templates live
in `resources/basepair-idealized/<class>/<seq>.pdb`; point elsewhere with
`rubric_set_ideals <dir>` or the `RUBRIC_IDEALS_DIR` env var. It's a sticky toggle,
**off by default** (so the view stays clean until you ask for the target). If no
template exists for a pair's class+sequence, it says so and shows nothing.

## Sugar-pucker mismatches

The binary reports a sugar under `pucker_outliers` when two readings of its pucker
disagree:

- **δ** (the C5′–C4′–C3′–O3′ torsion) reads the pucker from the ring itself;
- **P-perp**, the distance from the 3′ phosphorus to the line of the glycosidic
  bond extended past C1′, reads it from where the phosphate sits (short P-perp =
  C2′-endo, long = C3′-endo).

Each flagged sugar is drawn whatever the pair tier:
the whole nucleotide as gray sticks, with only its five ribose-ring atoms
(C1′ C2′ C3′ C4′ O4′) in **cyan**. The summary lists them:

```
pucker mismatch (3 sugar(s), cyan; rubric_pucker <n> for a close-up):
    1  A-U-199  delta 90.5° (C3'-endo) vs P-perp 2.31 Å (C2'-endo)
```

```
rubric_pucker              # list the flagged sugars
rubric_pucker 1            # close-up of the first one
rubric_pucker A-199        # ... or by chain + residue number
rubric_pucker off          # hide the cyan layer (sticky; `on` brings it back)
rubric_overview            # back to the whole structure
```

The close-up (as in Figure 7D) shows the nucleotide and its 3′ phosphate, the
glycosidic bond extended past C1′ (light dashes), and the perpendicular from the
3′ P onto it (dark dashes, labelled `P-perp 2.31 Å`). The view faces the plane of
that right angle, tilted 30° so the ring's pucker is visible. The plugin
recomputes P-perp from the coordinates for the drawing and prints it; it matches
the binary's value. When `rubric_next` / `rubric_goto` lands on a pair that contains
a flagged sugar, its ribose is cyan on the focus sticks too and the log notes it.
`rubric_dump` / `rubric_load` carry the flags, so a shared session shows them as well.

The flag is report-only: it does not change any pair score or the worklist.

## Share a scored session (colleague navigates with no binary)

A `.pse` saves the *picture* but not the plugin or the worklist, so a colleague
can't navigate from it alone. To let someone step through a scored structure
**without the `rubric` binary and without re-scoring**, send three files:

1. `1y26.pse`        — the session (`save 1y26.pse` after scoring)
2. `1y26_pairs.json` — the worklist (`rubric_dump`)
3. `rubric_pymol.py`  — the plugin

They then run, in PyMOL:

```
load 1y26.pse           # the structure + highlights
run rubric_pymol.py      # registers the rubric_* commands
rubric_load 1y26_pairs.json   # rebuilds the worklist from the JSON (no binary)
rubric_next              # ... and they can step through exactly as you did
```

`rubric_load` reconstructs the queue from the dumped JSON, so `rubric_next` / `rubric_prev`
/ `rubric_list` / `rubric_overview` / `rubric_dump` all work. Only `rubric_score` itself
(scoring a *new* model) needs the binary. Tip: put those three lines in a `view.pml`
so they only have to type `@view.pml`.

## Without the plugin: a standalone script (`rubric --pml`)

The `rubric` binary can also write the same view as a plain PyMOL script, so
anyone with PyMOL can open it, with no plugin and no binary at view time:

```bash
rubric 1GID --pml                                  # writes ./1GID.pml (JSON still goes to stdout)
rubric model.cif --pml views/model.pml --out model.json
pymol 1GID.pml                                    # or, inside PyMOL: @1GID.pml
```

A bare `--pml` writes `<PDB_id>.pml` in the current directory; any name or path
you give is used exactly as-is (missing directories are created), as for `--out`.

The script loads the structure file that was scored (its absolute path is in the
header) and builds:

| object / selection | what it is |
|---|---|
| `rubric_pairs` | a copy of the flagged pairs as sticks, colored by pair score: yellow → orange = secondary review, red = priority review (the same colors as `rubric_score`) |
| `rubric_priority`, `rubric_secondary` | the flagged residues of each tier |
| `rubric_worklist` → `rubric_pr_001`, `rubric_pr_002`, … | one selection per priority-review pair on the model, worst first (e.g. `zoom rubric_pr_001`) |
| `rubric_pucker`, `rubric_pucker_outliers` | pucker-mismatch nucleotides (δ vs P-perp), gray sticks with a cyan ribose ring |

The header also lists, as comments, the base-pair score, backbone suiteness and
tier counts, then every flagged pair (score, tier, LW class, penalized
parameters) worst first and every pucker outlier, so the file doubles as a
readable worklist.

Notes:

- Like the plugin, the script never modifies coordinates. After editing or
  re-refining the model, re-run `rubric ... --pml` to regenerate it.
- The structure path in the script is absolute. If you send the `.pml` to a
  colleague, send the structure file too and edit the path in the `cmd.load` line.
- Interactive features (`rubric_next`, `rubric_ideal`, the `rubric_pucker` close-up)
  need the plugin; the script gives the static overview and named selections.
- `--pml` needs scoring, so it cannot be combined with `--no-score`.
