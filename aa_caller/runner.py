from argparse import Namespace, ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping

from .constants import DEFAULT_ENTROPY_THRESHOLD, DEFAULT_RATIO_LOWER, DEFAULT_RATIO_UPPER
from .container import SamContainer
from .models import Amplicon
from .reference import FullReference
from .validators import parse_amplicons, validate_amplicon_file, validate_reference_file, validate_sam_file


@dataclass
class VariantCallResult:
    """Encapsulates the runner outputs along with the container state."""

    csv_path: Path
    xml_path: Path
    container: SamContainer
    reference: FullReference
    amplicons: Dict[str, Amplicon]


def call_variants(
    sam_path: Path | str,
    reference_path: Path | str,
    amplicons_path: Path | str,
    *,
    annotation_path: Path | str | None = None,
    ratio_upper: float = DEFAULT_RATIO_UPPER,
    ratio_lower: float = DEFAULT_RATIO_LOWER,
    entropy_threshold: float = DEFAULT_ENTROPY_THRESHOLD,
    protein: str = "RT",
    csv_path: Path | str | None = None,
    xml_path: Path | str | None = None,
) -> VariantCallResult:
    """Run the codon-level pipeline and return the generated artifacts."""

    sam_path = Path(sam_path)
    reference_path = Path(reference_path)
    amplicons_path = Path(amplicons_path)
    annotation_path = Path(annotation_path) if annotation_path else None

    validate_sam_file(sam_path)
    validate_reference_file(reference_path, protein_name=protein, annotation_path=annotation_path)
    validate_amplicon_file(amplicons_path)

    amplicons = parse_amplicons(amplicons_path)
    reference = FullReference(reference_path, annotation_path=annotation_path)
    container = SamContainer(
        sam_path,
        ratio_upper=ratio_upper,
        ratio_lower=ratio_lower,
        entropy_threshold=entropy_threshold,
    )
    container.load()
    container.calculate_variant_frequencies(protein, reference)

    csv_path = Path(csv_path) if csv_path else sam_path.with_suffix(sam_path.suffix + ".tsv")
    xml_path = Path(xml_path) if xml_path else sam_path.with_suffix(sam_path.suffix + ".xml")

    container.write_csv(str(sam_path), reference, csv_path, protein_name=protein)
    container.write_xml(str(sam_path), reference, xml_path)

    return VariantCallResult(
        csv_path=csv_path,
        xml_path=xml_path,
        container=container,
        reference=reference,
        amplicons=amplicons,
    )


def call_variants_from_args(
    args: Namespace | Mapping[str, Any], *, ratio_upper: float | None = None, ratio_lower: float | None = None, entropy_threshold: float | None = None,
    protein: str | None = None,
    csv_path: Path | str | None = None,
    xml_path: Path | str | None = None,
) -> VariantCallResult:
    """Run the pipeline using parsed CLI arguments or a plain mapping."""

    def _pick(key: str, default: Any = None) -> Any:
        if isinstance(args, Namespace):
            return getattr(args, key, default)
        if isinstance(args, Mapping):
            return args.get(key, default)
        return default

    sam = _pick("sam_file", _pick("sam_path"))
    reference = _pick("reference_file", _pick("reference_path"))
    amplicons = _pick("amplicons_file", _pick("amplicons_path"))
    annotation = _pick("annotation", _pick("annotation_path"))
    if not (sam and reference and amplicons):
        raise ValueError("args must include sam_file/reference_file/amplicons_file")

    effective_ratio_upper = ratio_upper if ratio_upper is not None else _pick("ratio_upper", DEFAULT_RATIO_UPPER)
    effective_ratio_lower = ratio_lower if ratio_lower is not None else _pick("ratio_lower", DEFAULT_RATIO_LOWER)
    effective_entropy_threshold = entropy_threshold if entropy_threshold is not None else _pick("entropy_threshold", DEFAULT_ENTROPY_THRESHOLD)
    effective_protein = protein if protein is not None else _pick("protein", "RT")

    overriding_csv = csv_path if csv_path is not None else _pick("csv_path")
    overriding_xml = xml_path if xml_path is not None else _pick("xml_path")

    return call_variants(
        sam_path=sam,
        reference_path=reference,
        amplicons_path=amplicons,
        annotation_path=annotation,
        ratio_upper=effective_ratio_upper,
        ratio_lower=effective_ratio_lower,
        entropy_threshold=effective_entropy_threshold,
        protein=effective_protein,
        csv_path=overriding_csv,
        xml_path=overriding_xml,
    )


def runner_cli() -> None:
    """Simple CLI that forwards parsed arguments to `call_variants_from_args`."""

    parser = ArgumentParser(description="Run the codonyat variant caller from Python. Research use only.")
    parser.add_argument("sam_file", type=Path, help="Path to the input SAM file")
    parser.add_argument("reference_file", type=Path, help="FASTA reference with protein annotations")
    parser.add_argument("amplicons_file", type=Path, help="Amplicon TSV/CSV describing primers")
    parser.add_argument("--annotation", type=Path, help="Optional TSV annotation file")
    parser.add_argument("--ratio-upper", type=float, default=DEFAULT_RATIO_UPPER, help="Upper bound for strand balance")
    parser.add_argument("--ratio-lower", type=float, default=DEFAULT_RATIO_LOWER, help="Lower bound for strand balance")
    parser.add_argument("--entropy-threshold", type=float, default=DEFAULT_ENTROPY_THRESHOLD, help="Minimum entropy to mark a position")
    parser.add_argument("--protein", type=str, default="RT", help="Target protein name to analyze (default: RT)")
    parser.add_argument("--csv-path", type=Path, help="Override the TSV output path")
    parser.add_argument("--xml-path", type=Path, help="Override the XML output path")

    args = parser.parse_args()
    result = call_variants_from_args(args, csv_path=args.csv_path, xml_path=args.xml_path)
    print(f"Wrote TSV: {result.csv_path}")
    print(f"Wrote XML: {result.xml_path}")
