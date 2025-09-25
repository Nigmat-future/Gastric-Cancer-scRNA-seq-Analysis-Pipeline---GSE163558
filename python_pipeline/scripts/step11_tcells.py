"""Entry point for Step 11: T-cell annotation."""

from __future__ import annotations

from ..tcell_annotation import run_step11_tcells
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 11 T-cell annotation")
    run_step11_tcells()


if __name__ == "__main__":
    main()
