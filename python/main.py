"""CLI wrapper that runs the C++ airGM solver via WSL."""

from __future__ import annotations

import shlex
import subprocess
import sys
from pathlib import Path


def windows_path_to_wsl(path: Path) -> str:
    drive = path.drive.rstrip(":").lower()
    tail = path.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{tail}"


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv

    repo_root = Path(__file__).resolve().parent.parent
    repo_wsl = windows_path_to_wsl(repo_root)
    solver_args = " ".join(shlex.quote(token) for token in argv[1:])
    bash_cmd = (
        f"cd {shlex.quote(repo_wsl)} && "
        f"if [ ! -x ./src/airGM2.1 ]; then (cd src && make); fi && "
        f"./src/airGM2.1 {solver_args}"
    )

    try:
        completed = subprocess.run(["wsl.exe", "--", "bash", "-lc", bash_cmd], check=False)
    except FileNotFoundError:
        print("ERROR: wsl.exe not found. Install WSL to run the C++ backend from Python.")
        return 1
    return int(completed.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
