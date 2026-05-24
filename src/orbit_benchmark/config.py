from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class OrbitConfig:
    top_n: int = 150
    top_k: tuple[int, ...] = (10, 25, 50)
    program_k: int = 50
    quality_thresholds: dict[str, float] = field(
        default_factory=lambda: {"moderate": 0.2, "stringent": 0.5}
    )
    matched_contexts: dict[str, list[str]] = field(default_factory=dict)
    core_contexts: list[str] = field(default_factory=list)
    timepoints: tuple[int, ...] = (6, 24)
    stress_lambda: float = 0.35
    stress_high_z: float = 1.5
    bootstrap_iterations: int = 500
    random_seed: int = 2026

    @classmethod
    def from_mapping(cls, values: dict[str, Any]) -> "OrbitConfig":
        data = dict(values)
        if "top_k" in data:
            data["top_k"] = tuple(int(x) for x in data["top_k"])
        if "timepoints" in data:
            data["timepoints"] = tuple(int(x) for x in data["timepoints"])
        return cls(**data)


def load_config(path: str | Path | None) -> OrbitConfig:
    if path is None:
        return OrbitConfig()
    with Path(path).open("r", encoding="utf-8") as handle:
        values = yaml.safe_load(handle) or {}
    return OrbitConfig.from_mapping(values)

