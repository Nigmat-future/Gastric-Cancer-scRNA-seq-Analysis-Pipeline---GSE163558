"""Entry point for Step 4: feature plots and marker heatmaps."""

from __future__ import annotations

from ..marker_analysis import run_step4_feature_heatmap
from ..utils.logging import get_logger

LOGGER = get_logger(__name__)


def main() -> None:
    LOGGER.info("Running Step 4 marker visualisation")
    run_step4_feature_heatmap()


if __name__ == "__main__":
    main()
