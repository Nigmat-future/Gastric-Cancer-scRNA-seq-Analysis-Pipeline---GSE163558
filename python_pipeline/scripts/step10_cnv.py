"""Entry point for Step 10: CNV inference."""

from __future__ import annotations

from ..cnv import run_cnv_inference
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 10 CNV inference")
    run_cnv_inference()


if __name__ == "__main__":
    main()
