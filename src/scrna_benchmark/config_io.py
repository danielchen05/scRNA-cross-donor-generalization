from __future__ import annotations

"""YAML configuration loaders for benchmark and preprocessing stages."""

from pathlib import Path

import yaml

from .config import DatasetConfig
from .preprocessing import PreprocessingConfig, load_preprocessing_config


def load_dataset_config(path: str | Path) -> DatasetConfig:
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    if raw.get("group_transfer_pairs") is not None:
        raw["group_transfer_pairs"] = [tuple(x) for x in raw["group_transfer_pairs"]]
    return DatasetConfig(**raw)


__all__ = ["DatasetConfig", "PreprocessingConfig", "load_dataset_config", "load_preprocessing_config"]
