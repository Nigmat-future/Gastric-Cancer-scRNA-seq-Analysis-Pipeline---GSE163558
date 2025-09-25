"""Entry point for Step 3: cell type annotation."""

from __future__ import annotations

from ..clustering import run_step3_celltype
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 3 cell type annotation")
    run_step3_celltype(resolution=0.5)


if __name__ == "__main__":
    main()
