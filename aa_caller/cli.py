from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .constants import (
    DEFAULT_ENTROPY_THRESHOLD,
    DEFAULT_PROTEIN,
    DEFAULT_RATIO_LOWER,
    DEFAULT_RATIO_UPPER,
)
from .container import SamContainer
from .reference import FullReference
from .validators import (
    parse_amplicons,
    validate_amplicon_file,
    validate_reference_file,
    validate_sam_file,
)

logger = logging.getLogger(__name__)


def main() -> None:
    """Command-line entry point that wires inputs/key outputs together."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(
        description="Codon-aware amino acid variant typer from SAM alignments. Research use only."
    )
    parser.add_argument("sam_file", type=Path)
    parser.add_argument("reference_file", type=Path)
    parser.add_argument("amplicons_file", type=Path)
    parser.add_argument(
        "--annotation",
        type=Path,
        help="Optional TSV annotation file (record, protein, description, intervals, strand)",
    )
    parser.add_argument(
        "--ratio-upper",
        type=float,
        default=DEFAULT_RATIO_UPPER,
        help="Upper bound for strand-ratio balancing (default mirrors Perl).",
    )
    parser.add_argument(
        "--ratio-lower",
        type=float,
        default=DEFAULT_RATIO_LOWER,
        help="Lower bound for strand-ratio balancing (default mirrors Perl).",
    )
    parser.add_argument(
        "--entropy-threshold",
        type=float,
        default=DEFAULT_ENTROPY_THRESHOLD,
        help="Minimum Shannon entropy required to mark a position as diverse.",
    )
    parser.add_argument(
        "--protein",
        type=str,
        default=DEFAULT_PROTEIN,
        help="Target protein name to analyze (default: RT).",
    )
    args = parser.parse_args()

    if not args.sam_file.exists():
        parser.error(f"SAM file {args.sam_file} does not exist")
    if not args.reference_file.exists():
        parser.error(f"Reference file {args.reference_file} does not exist")
    if not args.amplicons_file.exists():
        parser.error(f"Amplicon file {args.amplicons_file} does not exist")
    if args.annotation and not args.annotation.exists():
        parser.error(f"Annotation file {args.annotation} does not exist")

    logger.info("Validating input formats before parsing")
    validate_sam_file(args.sam_file)
    validate_reference_file(args.reference_file, protein_name=args.protein, annotation_path=args.annotation)
    validate_amplicon_file(args.amplicons_file)

    logger.info("Reading amplicon configuration")
    parse_amplicons(args.amplicons_file)

    logger.info("Loading reference")
    reference = FullReference(args.reference_file, annotation_path=args.annotation)

    logger.info("Parsing SAM entries")
    container = SamContainer(
        args.sam_file,
        ratio_upper=args.ratio_upper,
        ratio_lower=args.ratio_lower,
        entropy_threshold=args.entropy_threshold,
    )
    container.load()

    logger.info("Calculating %s variant frequencies", args.protein)
    container.calculate_variant_frequencies(args.protein, reference)

    csv_path = args.sam_file.with_suffix(args.sam_file.suffix + ".tsv")
    xml_path = args.sam_file.with_suffix(args.sam_file.suffix + ".xml")
    logger.info("Writing TSV results")
    container.write_csv(str(args.sam_file), reference, csv_path, protein_name=args.protein)
    logger.info("Writing XML dump")
    container.write_xml(str(args.sam_file), reference, xml_path)
