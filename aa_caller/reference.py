from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

from Bio import SeqIO

from .models import Interval, Protein


class FullReference:
    """Loads the reference sequence(s) and annotated proteins from FASTA.

    Supports multi-record FASTA files and optional TSV annotation sidecar.
    """
    def __init__(self, fasta_path: Path, annotation_path: Path | None = None) -> None:
        """Initialize the reference by parsing FASTA file(s) with protein metadata.

        Args:
            fasta_path: Path to the FASTA file (can contain multiple records)
            annotation_path: Optional path to TSV annotation file
        """
        self.path = fasta_path
        self.annotation_path = annotation_path
        self.sequences: dict[str, str] = {}
        self.record_ids: list[str] = []
        self.proteins: dict[tuple[str, str], Protein] = {}
        self._load_reference()

    @property
    def seq(self) -> str:
        """Return the first sequence for backward compatibility."""
        return self.sequences[self.record_ids[0]] if self.record_ids else ""

    @property
    def id(self) -> str:
        """Return the first record ID for backward compatibility."""
        return self.record_ids[0] if self.record_ids else ""

    def _load_reference(self) -> None:
        """Populate sequences and protein dictionaries from FASTA and annotations."""
        for record in SeqIO.parse(str(self.path), "fasta"):
            self.sequences[record.id] = str(record.seq).upper()
            self.record_ids.append(record.id)

            if self.annotation_path is None:
                for protein in proteins_from_fasta_record(record):
                    self.proteins[(record.id, protein.name)] = protein

        if self.annotation_path:
            self._load_tsv_annotation()

    def _load_tsv_annotation(self) -> None:
        """Load protein annotations from TSV file.

        Format: record, protein, description, intervals, strand
        Intervals: 1-based inclusive, e.g., "266-13468;13468-21555"
        """
        with open(self.annotation_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                record_id = row["record"]
                protein_name = row["protein"]
                description = row.get("description", "")
                intervals_text = row["intervals"]
                strand = row.get("strand", "+")

                if strand not in ["+", "-"]:
                    raise ValueError(f"Invalid strand '{strand}' for protein {protein_name} on {record_id}")

                if strand == "-":
                    raise ValueError(
                        f"Minus-strand genes are not supported: {protein_name} on {record_id} (strand={strand}). "
                        "codonyat requires plus-strand annotations only."
                    )

                intervals = []
                for interval_text in intervals_text.split(";"):
                    interval_text = interval_text.strip()
                    if not interval_text:
                        continue
                    start, end = interval_text.split("-")
                    intervals.append(Interval(start=int(start), end=int(end)))

                protein = Protein(
                    name=protein_name,
                    description=description,
                    record_id=record_id,
                    intervals=intervals,
                    strand=strand
                )
                key = (record_id, protein_name)
                self.proteins[key] = protein

    def get_protein(self, protein_name: str, record_id: str | None = None) -> Protein | None:
        """Get a protein by name, optionally specifying the record.

        For backward compatibility, if record_id is not specified and there's
        only one record, return the protein from that record.
        """
        if record_id:
            return self.proteins.get((record_id, protein_name))

        if len(self.record_ids) == 1:
            return self.proteins.get((self.record_ids[0], protein_name))

        for key, protein in self.proteins.items():
            if key[1] == protein_name:
                return protein

        return None

    def get_seq_at(
        self,
        position: int,
        length: int = 3,
        record_id: str | None = None,
        positions: Sequence[int] | None = None,
    ) -> str:
        """Return a substring of the reference at the provided base coordinate.

        Args:
            position: 1-based position (used when ``positions`` is omitted)
            length: Length of substring to extract (consecutive slice)
            record_id: Optional record ID; defaults to first record
            positions: Actual part coordinates of a codon. When given, those
                bases are joined in order so a junction codon is reconstructed
                even if the parts are not a consecutive genomic slice.
        """
        if record_id is None:
            record_id = self.record_ids[0] if self.record_ids else ""

        seq = self.sequences.get(record_id, "")
        if positions is not None:
            bases: list[str] = []
            for coord in positions:
                idx = int(coord) - 1
                if 0 <= idx < len(seq):
                    bases.append(seq[idx])
            return "".join(bases)
        return seq[position - 1 : position - 1 + length]

    def proteins_in_order(self) -> list[Protein]:
        """Return annotated proteins in FASTA-header / TSV order."""
        return list(self.proteins.values())


def proteins_from_fasta_record(record) -> list[Protein]:
    """Parse protein annotations from a FASTA record header.

    Names are the annotation identifiers (``RT``, ``S``, ``ORF1ab``), not
    substrings of the description line.
    """
    proteins: list[Protein] = []
    header_parts = record.description.split(";")
    for i, part in enumerate(header_parts):
        part = part.strip()
        if not part:
            continue
        if i == 0 and part.startswith(record.id):
            part = part[len(record.id) :].strip()
            if not part:
                continue
        if ":" not in part or "(" not in part:
            continue
        proteins.append(Protein.from_string(part, record_id=record.id))
    return proteins


def parse_protein_spec(spec: str | Sequence[str]) -> list[str]:
    """Split a protein selection (``"RT"``, ``"PR,RT,INT"``, ``"all"`` or a list) into names.

    Order is kept and duplicates are dropped. ``"all"`` is returned as ``["all"]``
    and expanded by :func:`resolve_proteins` once the reference is loaded.
    """
    if isinstance(spec, str):
        items = [part.strip() for part in spec.split(",")]
    else:
        items = [str(part).strip() for part in spec]
    names: list[str] = []
    for item in items:
        if item and item not in names:
            names.append(item)
    if not names:
        raise ValueError("No protein selected")
    if "all" in (name.lower() for name in names) and len(names) > 1:
        raise ValueError("'all' cannot be combined with other protein names")
    return names


def resolve_protein_targets(spec: str | Sequence[str], reference: FullReference) -> list[Protein]:
    """Resolve a protein selection to annotated Protein objects.

    ``"all"`` selects every annotated protein in FASTA-header / TSV order.
    A comma-separated list keeps the requested order. On multi-record
    references, each requested name includes every matching record.
    """
    names = parse_protein_spec(spec)
    if len(names) == 1 and names[0].lower() == "all":
        proteins = reference.proteins_in_order()
        if not proteins:
            raise ValueError(f"Reference {reference.path} has no annotated proteins")
        return proteins
    selected: list[Protein] = []
    missing: list[str] = []
    for name in names:
        matches = [protein for key, protein in reference.proteins.items() if key[1] == name]
        if not matches:
            missing.append(name)
        else:
            selected.extend(matches)
    if missing:
        raise ValueError(f"Protein(s) not found in reference: {', '.join(missing)}")
    return selected


def resolve_proteins(spec: str | Sequence[str], reference: FullReference) -> list[str]:
    """Resolve a protein selection to protein names.

    ``"all"`` selects every annotated protein in FASTA-header / TSV order.
    Names are unique; duplicate names across records are collapsed in
    first-seen order. Use :func:`resolve_protein_targets` when record
    identity must be preserved.
    """
    names = parse_protein_spec(spec)
    if len(names) == 1 and names[0].lower() == "all":
        if not reference.proteins:
            raise ValueError(f"Reference {reference.path} has no annotated proteins")
        seen: list[str] = []
        for _record_id, protein_name in reference.proteins:
            if protein_name not in seen:
                seen.append(protein_name)
        return seen
    missing = [name for name in names if reference.get_protein(name) is None]
    if missing:
        raise ValueError(f"Protein(s) not found in reference: {', '.join(missing)}")
    return names
