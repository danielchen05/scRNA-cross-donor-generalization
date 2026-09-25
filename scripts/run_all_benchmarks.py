from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

from _common import bootstrap_repo

DEFAULT_CONFIGS = [
    Path("configs/benchmarks/kidney.yaml"),
    Path("configs/benchmarks/pancreas.yaml"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run kidney + pancreas plot-free benchmarks.")
    parser.add_argument(
        "--clean-output",
        action="store_true",
        help="Delete each dataset's old benchmark result directory first",
    )
    args = parser.parse_args()

    root = bootstrap_repo()
    for config in DEFAULT_CONFIGS:
        cmd = [sys.executable, str(root / "scripts" / "run_benchmark.py"), "--config", str(config)]
        if args.clean_output:
            cmd.append("--clean-output")
        print("\n$", " ".join(cmd))
        subprocess.run(cmd, cwd=root, check=True)


if __name__ == "__main__":
    main()
