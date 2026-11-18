# RUBRIC → PyMOL

Details of the RUBRIC PyMOL plug-in. Installing RUBRIC, starting the viewer
(`rubric-view`) and the list of commands are in the
[main README](../../README.md#5-pymol-viewer).

## In a PyMOL session that is already open

`rubric-view` loads the plug-in for you. To load it into a session you already have
open, print the location of the plug-in file:

```bash
python -c "import rubric_rna, pathlib; print(pathlib.Path(rubric_rna.__file__).parent / 'share/rubric/rubric_pymol.py')"
```

and run that file in PyMOL (`run /path/to/rubric_pymol.py`). In a source checkout
the file is `integrations/pymol/rubric_pymol.py`. If PyMOL cannot find the `rubric`
program, give its location with `rubric_set_binary` (the output of `which rubric`).

## Colors

| color | meaning |
|---|---|
| gray | the rest of the structure, including pairs that score 100 |
| red | a base pair scoring below 80 (priority review) |
| yellow to orange | a base pair scoring between 80 and 100 (secondary review): yellow near 100, orange near 80 |
| cyan | a ribose whose pucker is inconsistent |
| green | the ideal shape of the current pair (`rubric_ideal`) |

In the backbone view, yellow, orange and red mark fixable, review and flag-only
suites instead. Only colors change; coordinates are never modified.

## Ribose flags

RUBRIC reads the pucker of each ribose in two ways: from the ring itself (the δ
torsion), and from where the 3′ phosphate sits (P-perp, the distance from the
phosphorus to the line of the glycosidic bond). When the two disagree, the ribose
is flagged and drawn cyan. The close-up (`rubric_pucker 1`) draws the glycosidic
bond extended past C1′ (light dashes) and P-perp (dark dashes, labelled with its
length). The flag does not change any pair score.

## Backbone flags

A suite is the stretch of backbone between two nucleotides. Each suite is assigned
to its closest known conformation and given a suiteness from 0 (poor) to 1 (good).
Flagged suites fall into three groups:

- **fixable** (yellow): one or more torsions are out of range for the conformation.
  The close-up lists them with the values to move them toward.
- **review** (orange): the suiteness is low but no single torsion is out of range.
  The suite needs rebuilding.
- **flag only** (red): the suite fits no known conformation. The target shown is only
  the nearest conformation.

The backbone flags do not change any pair score.

## Share a scored session

To let someone go through your results without running RUBRIC, send them three files:

1. the PyMOL session (`save 1GID.pse` after scoring),
2. the results (`rubric_dump` writes `1GID_pairs.json`),
3. the plug-in file `rubric_pymol.py`.

They then run in PyMOL:

```
load 1GID.pse
run rubric_pymol.py
rubric_load 1GID_pairs.json
rubric_next
```

All commands except `rubric_score` then work. `rubric_ideal` also needs RUBRIC
installed, because the ideal shapes come with it.

## Without the plug-in: `rubric --pml`

`rubric` can also write a plain PyMOL script that shows the same colors, so anyone
with PyMOL can open it without the plug-in:

```bash
rubric 1GID --pml          # writes 1GID.pml
pymol 1GID.pml             # or, inside PyMOL: @1GID.pml
```

The script colors the flagged pairs and riboses as above and creates named
selections: `rubric_pairs` (all flagged pairs), `rubric_priority` and
`rubric_secondary` (the residues of each group), `rubric_pr_001`, `rubric_pr_002`,
... (one per pair scoring below 80, worst first) and `rubric_pucker` (the flagged
riboses). Its header lists the scores and every flagged pair. It has no interactive
commands and no backbone view. The structure path inside the script is absolute, so
if you send the script to someone, send the structure file too and correct the path
in its `cmd.load` line.
