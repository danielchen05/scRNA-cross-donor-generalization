from __future__ import annotations

"""Small CLI bootstrap helpers used by repository scripts."""

from pathlib import Path
import sys


def find_repo_root(start: str | Path | None = None) -> Path:
    here = Path(start or __file__).resolve()
    if here.is_file():
        here = here.parent
    for candidate in [here, *here.parents]:
        if (candidate / "src" / "scrna_benchmark").is_dir():
            return candidate
    raise RuntimeError("Could not locate repository root containing src/scrna_benchmark")


def bootstrap_repo() -> Path:
    root = find_repo_root()
    src = root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    return root
