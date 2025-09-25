"""Entry point for Step 5: composition analysis."""

from __future__ import annotations

from ..proportions import run_step5_proportions
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 5 cell-type proportions")
    run_step5_proportions()


if __name__ == "__main__":
    main()
