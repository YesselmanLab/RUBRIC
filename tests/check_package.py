"""Check an installed `rubric` against the reference result for 1GID.

Runs `rubric tests/data/1GID.cif.gz` and compares the base pairs (residues, LW class,
score), the structure scores, the tier counts and the pucker outliers with
tests/data/1GID.expected.json. Scores are compared within a small tolerance, because
compilers on different platforms need not round floating-point results identically.

    python tests/check_package.py            # check
    python tests/check_package.py --update   # rewrite the reference from this build
"""
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
STRUCTURE = DATA / "1GID.cif.gz"
EXPECTED = DATA / "1GID.expected.json"
SCORE_TOL = 0.05
SUITENESS_TOL = 0.002


def summary(report):
    pairs = sorted(
        [p["res_id1"], p["res_id2"], (p.get("classification") or {}).get("lw_class"), p.get("score")]
        for p in report["pairs"]
    )
    return {
        "pdb_id": report["pdb_id"],
        "n_pairs": report["n_pairs"],
        "pairs_score": report["pairs_score"],
        "backbone_suiteness": report["backbone_suiteness"],
        "tier_summary": report["tier_summary"],
        "pucker_outliers": sorted(o["res_id"] for o in report.get("pucker_outliers") or []),
        "pairs": pairs,
    }


def run_rubric():
    exe = shutil.which("rubric")
    if exe is None:
        sys.exit("`rubric` is not on PATH")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "report.json"
        subprocess.run([exe, str(STRUCTURE), "--no-download", "--out", str(out)],
                       check=True, capture_output=True)
        return json.loads(out.read_text())


def close(a, b, tol):
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, abs_tol=tol)


def compare(got, exp):
    problems = []
    for key in ("pdb_id", "n_pairs", "tier_summary", "pucker_outliers"):
        if got[key] != exp[key]:
            problems.append(f"{key}: {got[key]} != expected {exp[key]}")
    if not close(got["pairs_score"], exp["pairs_score"], SCORE_TOL):
        problems.append(f"pairs_score {got['pairs_score']} != {exp['pairs_score']}")
    if not close(got["backbone_suiteness"], exp["backbone_suiteness"], SUITENESS_TOL):
        problems.append(f"backbone_suiteness {got['backbone_suiteness']} != {exp['backbone_suiteness']}")
    got_pairs = {tuple(p[:3]): p[3] for p in got["pairs"]}
    exp_pairs = {tuple(p[:3]): p[3] for p in exp["pairs"]}
    for key in sorted(set(got_pairs) ^ set(exp_pairs)):
        problems.append(f"pair {key} {'unexpected' if key in got_pairs else 'missing'}")
    for key in sorted(set(got_pairs) & set(exp_pairs)):
        if not close(got_pairs[key], exp_pairs[key], SCORE_TOL):
            problems.append(f"pair {key} score {got_pairs[key]} != {exp_pairs[key]}")
    return problems


def main():
    got = summary(run_rubric())
    if "--update" in sys.argv:
        EXPECTED.write_text(json.dumps(got, indent=1) + "\n")
        print(f"wrote {EXPECTED}")
        return
    problems = compare(got, json.loads(EXPECTED.read_text()))
    if problems:
        print("\n".join(problems[:50]))
        sys.exit(f"{len(problems)} difference(s) from the reference result")
    print(f"rubric matches the reference for 1GID ({got['n_pairs']} pairs)")


if __name__ == "__main__":
    main()
