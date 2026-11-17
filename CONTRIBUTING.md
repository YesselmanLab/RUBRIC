# Contributing to RUBRIC

These guidelines follow the cctbx contribution guidelines
(https://github.com/cctbx/cctbx_project/blob/master/CONTRIBUTING.md), adapted to a
C++ program.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

See `tests/README.md` for what the tests cover and how to add one.

## Expectations for a change

- **Tests.** New functionality comes with tests: unit tests for new functions and a
  command-line or regression test for new user-facing behavior. A bug fix comes
  with a test that fails before the fix.
- **All tests pass** on Linux, macOS and Windows, including the sanitizer run, before
  a change is merged. A failing test is fixed first.
- **Results.** If a change alters RUBRIC's results on purpose, update the stored
  panel results (`tests/README.md`) and explain the change in the commit.
- **No new dependencies.** RUBRIC depends only on GEMMI (fetched at build time) and
  nlohmann/json (bundled header). A new dependency needs discussion first. The tests
  add none. If the GEMMI version changes, change it in `conda/recipe.yaml` too.
- **Documentation.** Public functions and classes have a comment saying what they
  do, their inputs and their result. Command-line options are described in
  `--help` and the README.

## Style

- C++20. Follow the style of the file being edited.
- Clear names over short ones; functions that fit on a screen.
- Comments explain why, not what the code already says. No debugging output in
  committed code; diagnostics go to standard error behind an option.
- No compiler warnings (the CI builds with warnings as errors on the current GCC and
  Clang).
- Constants that define results (cutoffs, bands, weights) live in one place and are
  named.

## Releases

1. Set the new version in `CMakeLists.txt` (`project(... VERSION 1.1)`) and in
   `conda/recipe.yaml`.
2. Commit and push to `main`; wait for the CI to pass.
3. On GitHub, publish a release (Releases, Draft a new release) with a new tag such
   as `v1.1`. With the Zenodo integration enabled, the release is also archived on
   Zenodo with its own DOI.
4. The tag starts the release workflow. It builds the packages, installs them on every
   supported Python version, runs the test suite, then waits for approval in the
   `pypi` environment before publishing to PyPI.

A version number can be published to PyPI only once.
