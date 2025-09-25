"""Entry point for Step 12: cell-cell communication."""

from __future__ import annotations

from ..cell_communication import run_cell_communication
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 12 cell-cell communication analysis")
    run_cell_communication()


if __name__ == "__main__":
    main()
