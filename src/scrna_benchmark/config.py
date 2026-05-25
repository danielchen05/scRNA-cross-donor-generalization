# src/scrna_benchmark/config.py

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


def _default_representations() -> dict[str, str]:
    return {
        "hvg": "hvg",
        "pca": "X_pca",
        "harmony": "X_pca_harmony",
        "scvi": "X_scVI",
    }


@dataclass
class DatasetConfig:
    """Configuration for one dataset in the benchmark.

    The notebook should define one DatasetConfig and pass it to
    scrna_benchmark.pipeline.run_dataset_benchmark.
    """

    dataset_name: str
    adata_path: str | Path
    output_dir: str | Path

    celltype_col: str
    donor_col: str
    batch_col: str | None = None
    group_col: str | None = None

    rename_obs: dict[str, str] = field(default_factory=dict)
    representations: dict[str, str] = field(default_factory=_default_representations)

    min_cells: int | None = 200
    min_donors: int | None = 5

    run_random_split: bool = True
    run_donor_cv: bool = True
    run_group_transfer: bool = False
    run_donor_ablation: bool = False

    test_size: float = 0.2
    n_folds: int = 5
    random_state: int = 42

    random_split_n_repeats: int = 5
    random_split_seeds: list[int] | None = None

    group_transfer_pairs: list[tuple[str, str]] = field(default_factory=list)
    donor_ablation_k_values: list[int] = field(default_factory=list)
    donor_ablation_n_repeats: int = 5

    save_filtered_adata: bool = False
    verbose: bool = True

    def __post_init__(self) -> None:
        self.adata_path = Path(self.adata_path)
        self.output_dir = Path(self.output_dir)

        if not self.dataset_name:
            raise ValueError("dataset_name must be non-empty.")
        if self.test_size <= 0 or self.test_size >= 1:
            raise ValueError("test_size must be between 0 and 1.")
        if self.n_folds < 2:
            raise ValueError("n_folds must be at least 2.")
        if self.random_split_n_repeats < 1:
            raise ValueError("random_split_n_repeats must be at least 1.")
        if self.random_split_seeds is not None and len(self.random_split_seeds) == 0:
            raise ValueError("random_split_seeds must be None or a non-empty list.")

    @property
    def dataset_output_dir(self) -> Path:
        return Path(self.output_dir) / self.dataset_name

    @property
    def filtering_output_dir(self) -> Path:
        return self.dataset_output_dir / "filtering"

    @property
    def random_output_dir(self) -> Path:
        return self.dataset_output_dir / "random_split"

    @property
    def donor_cv_output_dir(self) -> Path:
        return self.dataset_output_dir / "donor_cv"

    @property
    def group_transfer_output_dir(self) -> Path:
        return self.dataset_output_dir / "group_transfer"

    @property
    def donor_ablation_output_dir(self) -> Path:
        return self.dataset_output_dir / "donor_ablation"

    def required_obs_columns(self) -> list[str]:
        cols = [self.celltype_col, self.donor_col]
        if self.batch_col is not None:
            cols.append(self.batch_col)
        if self.group_col is not None:
            cols.append(self.group_col)
        return cols

    def to_dict(self) -> dict[str, Any]:
        out = asdict(self)
        out["adata_path"] = str(self.adata_path)
        out["output_dir"] = str(self.output_dir)
        return out
