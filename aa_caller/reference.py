from __future__ import annotations

import csv
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
                self._parse_fasta_header(record)

        if self.annotation_path:
            self._load_tsv_annotation()

    def _parse_fasta_header(self, record) -> None:
        """Parse protein annotations from FASTA header."""
        header_parts = record.description.split(";")
        for i, part in enumerate(header_parts):
            part = part.strip()
            if not part:
                continue
            if i == 0 and part.startswith(record.id):
                part = part[len(record.id):].strip()
                if not part:
                    continue
            if ":" not in part or "(" not in part:
                continue
            protein = Protein.from_string(part, record_id=record.id)
            key = (record.id, protein.name)
            self.proteins[key] = protein

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

    def get_seq_at(self, position: int, length: int = 3, record_id: str | None = None) -> str:
        """Return a substring of the reference at the provided base coordinate.

        Args:
            position: 1-based position
            length: Length of substring to extract
            record_id: Optional record ID; defaults to first record
        """
        if record_id is None:
            record_id = self.record_ids[0] if self.record_ids else ""

        seq = self.sequences.get(record_id, "")
        return seq[position - 1 : position - 1 + length]
