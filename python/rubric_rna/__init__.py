"""RUBRIC: Rapid Unified Base-pair, Ribose and backbone Integrity Check.

The ``rubric`` and ``rubric-view`` commands installed by this package run the
compiled RUBRIC program and its PyMOL viewer, which ship inside the package.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from importlib.metadata import PackageNotFoundError, version
    __version__ = version("rubric-rna")
except PackageNotFoundError:  # running from a source checkout
    __version__ = "unknown"

_PACKAGE = Path(__file__).resolve().parent
_BIN = _PACKAGE / "bin"


def binary_path() -> Path:
    """Path to the bundled ``rubric`` executable."""
    return _BIN / ("rubric.exe" if os.name == "nt" else "rubric")


def _require_binary() -> Path:
    """The bundled executable, or a clear error if the package was not built."""
    exe = binary_path()
    if not exe.is_file():
        raise SystemExit(
            f"The compiled rubric program was not found at {exe}. Install the "
            "package with `pip install rubric-rna` (or `pip install .` from a checkout)."
        )
    return exe


def _run(command, env=None):
    """Replace this process with ``command`` (or wait for it on Windows)."""
    if os.name == "nt":
        raise SystemExit(subprocess.call(command, env=env))
    if env is None:
        os.execv(command[0], command)
    os.execve(command[0], command, env)


def main() -> None:
    """Entry point of the ``rubric`` command."""
    _run([str(_require_binary()), *sys.argv[1:]])


def view_main() -> None:
    """Entry point of the ``rubric-view`` command (needs PyMOL)."""
    env = dict(os.environ)
    env.setdefault("RUBRIC_BINARY", str(_require_binary()))
    bash = shutil.which("bash")
    if bash is None:
        raise SystemExit("rubric-view needs bash, which was not found on PATH")
    _run([bash, str(_BIN / "rubric-view"), *sys.argv[1:]], env)
