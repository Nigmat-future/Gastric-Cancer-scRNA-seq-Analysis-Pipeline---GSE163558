"""Input/output helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def ensure_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def save_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")

