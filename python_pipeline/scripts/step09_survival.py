"""Entry point for Step 9: TCGA survival analysis."""

from __future__ import annotations

from ..survival import run_survival_analysis
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 9 survival analysis")
    run_survival_analysis()


if __name__ == "__main__":
    main()
