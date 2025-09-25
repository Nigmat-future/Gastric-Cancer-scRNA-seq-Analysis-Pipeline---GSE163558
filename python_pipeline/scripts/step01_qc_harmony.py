"""Entry point for Step 1-2: QC and Harmony integration."""

from __future__ import annotations

from pathlib import Path

from .. import config
from ..preprocessing import preprocess_pipeline
from ..utils.io import ensure_dirs
from ..utils.logging import get_logger, set_log_file

LOGGER = get_logger(__name__)


def main() -> None:
    ensure_dirs(config.DEFAULT_OUTPUT_SUBDIRS.values())
    log_path = config.DEFAULT_OUTPUT_SUBDIRS["harmony"] / "step01_qc_harmony.log"
    set_log_file(log_path)
    LOGGER.info("Starting QC + Harmony pipeline")
    adata = preprocess_pipeline()
    LOGGER.info("Completed QC pipeline; final cell count: %d", adata.n_obs)


if __name__ == "__main__":
    main()
