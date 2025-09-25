"""Shared logging helpers."""

from __future__ import annotations

import logging
from pathlib import Path
from contextlib import contextmanager
import json
import time

_LOGGER_NAME = "gc_scRNA_py"


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a configured logger instance.

    Parameters
    ----------
    name:
        Optional child logger name. Defaults to the root logger for the package.
    """
    logger_name = f"{_LOGGER_NAME}.{name}" if name else _LOGGER_NAME
    logger = logging.getLogger(logger_name)
    if not logging.getLogger(_LOGGER_NAME).handlers:
        _configure_root_logger()
    return logger


def _configure_root_logger() -> None:
    logger = logging.getLogger(_LOGGER_NAME)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def set_log_file(path: Path) -> None:
    """Add a file handler that mirrors the console log output."""
    logger = logging.getLogger(_LOGGER_NAME)
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


@contextmanager
def time_block(task_name: str, *, write_jsonl: Path | None = None):
    """Context manager to time a code block and optionally append JSONL.

    Parameters
    ----------
    task_name:
        A short human-readable name for the task.
    write_jsonl:
        If provided, appends a JSON line with timing info to this file.
    """
    logger = get_logger(__name__)
    start = time.time()
    logger.info("[TIMER] %s | started", task_name)
    try:
        yield
    finally:
        end = time.time()
        elapsed = end - start
        logger.info("[TIMER] %s | completed in %.2f s", task_name, elapsed)
        if write_jsonl is not None:
            try:
                write_jsonl.parent.mkdir(parents=True, exist_ok=True)
                payload = {
                    "task": task_name,
                    "t_start": start,
                    "t_end": end,
                    "elapsed_sec": round(elapsed, 3),
                }
                with write_jsonl.open("a", encoding="utf-8") as f:
                    f.write(json.dumps(payload, ensure_ascii=False) + "\n")
            except Exception as exc:  # best-effort, don't crash main flow
                logger.warning("Failed to write timing JSONL for %s: %s", task_name, exc)

