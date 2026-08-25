#!/usr/bin/env python3
"""Stage a compatible CP2K wavefunction as a matched-diagnostic start guess."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("case_list", type=Path)
    parser.add_argument("seed_wfn", type=Path)
    args = parser.parse_args()

    if not args.seed_wfn.is_file():
        raise FileNotFoundError(args.seed_wfn)

    staged = 0
    for line in args.case_list.read_text().splitlines():
        relative = line.strip()
        if not relative or relative.startswith("#"):
            continue
        run_dir = args.root / relative
        input_file = run_dir / "input.inp"
        text = input_file.read_text()
        projects = re.findall(r"^\s*PROJECT_NAME\s+(\S+)\s*$", text, re.MULTILINE)
        if len(projects) != 1:
            raise ValueError(f"Expected one PROJECT_NAME in {input_file}")
        restart_text, count = re.subn(
            r"(?m)^(\s*SCF_GUESS\s+)\S+\s*$", r"\g<1>RESTART", text, count=1
        )
        if count != 1:
            raise ValueError(f"Expected one SCF_GUESS in {input_file}")
        (run_dir / "input.restart.inp").write_text(restart_text)
        shutil.copy2(args.seed_wfn, run_dir / f"{projects[0]}-RESTART.wfn")
        staged += 1

    print(f"Staged {staged} restart guesses from {args.seed_wfn}")


if __name__ == "__main__":
    main()
