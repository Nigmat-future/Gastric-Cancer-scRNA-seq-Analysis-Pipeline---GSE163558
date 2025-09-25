"""Entry point for Step 7: malignant signatures and enrichment."""

from __future__ import annotations

from ..signature_analysis import run_step7_signatures
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 7 malignant signature analysis")
    run_step7_signatures()


if __name__ == "__main__":
    main()
