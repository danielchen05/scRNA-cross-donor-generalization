from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

from _common import bootstrap_repo


DEFAULT_CONFIGS = [
    Path("configs/preprocessing/kidney.yaml"),
    Path("configs/preprocessing/pancreas.yaml"),
    Path("configs/preprocessing/lung.yaml"),
    Path("configs/preprocessing/pbmc.yaml"),
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run kidney + pancreas + lung + PBMC preprocessing."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite benchmark-ready files",
    )
    args = parser.parse_args()

    root = bootstrap_repo()
    for config in DEFAULT_CONFIGS:
        cmd = [
            sys.executable,
            str(root / "scripts" / "preprocess.py"),
            "--config",
            str(config),
        ]
        if args.force:
            cmd.append("--force")
        print("\n$", " ".join(cmd))
        subprocess.run(cmd, cwd=root, check=True)


if __name__ == "__main__":
    main()
