"""Entry point for Step 6: TCGA-STAD processing and malignant scoring."""

from __future__ import annotations

from ..malignant_scoring import run_step6_malignant
from ..tcga_bulk import run_tcga_pipeline, download_tcga_inputs
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running TCGA-STAD bulk processing")
    # Ensure inputs exist (auto-download if missing)
    download_tcga_inputs()
    run_tcga_pipeline()
    LOGGER.info("Running malignant scoring on epithelial cells")
    run_step6_malignant()


if __name__ == "__main__":
    main()
