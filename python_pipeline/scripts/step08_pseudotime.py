"""Entry point for Step 8: pseudotime analysis."""

from __future__ import annotations

from ..pseudotime import run_pseudotime
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 8 pseudotime")
    run_pseudotime()


if __name__ == "__main__":
    main()
