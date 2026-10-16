from __future__ import annotations

import re
from pathlib import Path

from Bio import SeqIO

from .models import Amplicon
from .reference import proteins_from_fasta_record


def validate_sam_file(path: Path) -> None:
    """Perform a lightweight sanity check on the SAM file before parsing."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11 or not fields[0] or not fields[3].isdigit():
                raise ValueError(f"SAM file {path} looks malformed: {line.strip()}")
            return
        raise ValueError(f"SAM file {path} contains no alignment entries")


def validate_reference_file(
    path: Path, *, protein_name: str = "RT", annotation_path: Path | None = None
) -> None:
    """Ensure the reference FASTA is valid and has protein annotations.

    Supports multi-record FASTA files. If annotation_path is provided,
    validates that instead of FASTA headers.

    ``protein_name`` may be one name, a comma-separated list, or ``all``
    (which only requires the FASTA to be readable).
    """
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"Reference file {path} is empty or not FASTA")

    if annotation_path:
        return

    annotated_names = {
        protein.name for record in records for protein in proteins_from_fasta_record(record)
    }
    for name in protein_name.split(","):
        name = name.strip()
        if not name or name.lower() == "all":
            continue
        if name not in annotated_names:
            raise ValueError(f"Reference {path} header is missing an {name} protein annotation")


def validate_amplicon_file(path: Path) -> None:
    """Check that the amplicon TSV/CSV has at least a header plus one valid row."""
    with open(path) as fh:
        header = fh.readline()
        if not header or len(header.strip().split("\t")) < 7:
            raise ValueError(f"Amplicon file {path} does not contain the expected header")
        for line in fh:
            if not line.strip():
                continue
            parts = re.split(r"[\t,\s]+", line.strip().replace('"', ""))
            if len(parts) < 7:
                raise ValueError(f"Amplicon line needs 7 columns: {line.strip()}")
            return
        raise ValueError(f"Amplicon file {path} contains no definitions")


def parse_amplicons(path: Path) -> dict[str, Amplicon]:
    """Read the amplicon configuration file and construct Amplicon objects."""
    amplicons: dict[str, Amplicon] = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            if not line.strip():
                continue
            # parse each non-empty line into an Amplicon object
            amplicon = Amplicon.from_string(line)
            amplicons[amplicon.label] = amplicon
    return amplicons
