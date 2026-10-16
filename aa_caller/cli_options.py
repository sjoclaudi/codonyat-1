"""Shared opt-in CLI flags for research statistics and AAVF."""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from .constants import DEFAULT_MIN_COVERAGE, DEFAULT_MIN_FREQ, DEFAULT_SB_ALPHA


def add_research_output_args(parser: ArgumentParser) -> None:
    """Attach opt-in 1.3.0 flags. Defaults leave TSV/XML unchanged."""
    parser.add_argument(
        "--min-baseq",
        type=float,
        default=None,
        help=(
            "If set, count a codon only when all three bases have Phred-like "
            "quality ≥ this value (e.g. 20). Off by default so TSV/XML stay "
            "byte-identical to 1.2.1."
        ),
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=None,
        help="Optional stats TSV (Wilson 95%% CI, strand-bias p, flags, diversity).",
    )
    parser.add_argument(
        "--aavf",
        type=Path,
        default=None,
        help="Optional AAVF v1.0 path. Research use only; CAT/SRVL are omitted.",
    )
    parser.add_argument(
        "--min-freq",
        type=float,
        default=DEFAULT_MIN_FREQ,
        help="Frequency floor (fraction 0–1) for the lowfreq flag / AAVF FILTER (default: 0.01).",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=DEFAULT_MIN_COVERAGE,
        help="Coverage floor for the lowcov flag / AAVF FILTER (default: 10).",
    )
    parser.add_argument(
        "--sb-alpha",
        type=float,
        default=DEFAULT_SB_ALPHA,
        help="Fisher exact p-value cutoff for the sb (strand-bias) flag (default: 0.05).",
    )
