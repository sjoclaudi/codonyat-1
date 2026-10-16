from __future__ import annotations

__all__ = ["SamEntry"]

import re

from .models import Qual

# Compiled once at module load — reused across every SamEntry instance.
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


class SamEntry:
    """Wraps a SAM record with convenience helpers used during variant aggregation."""

    def __init__(self, line: str) -> None:
        """Parse a single SAM line and cache derived attributes that the pipeline needs."""
        self.raw = line.strip()
        self.fields = self.raw.split("\t")
        self.identifier = self.fields[0]
        self.flag = int(self.fields[1])
        self.reference = self.fields[2]
        self.coordinate = int(self.fields[3])
        self.map_quality = int(self.fields[4])
        self.cigar = self.fields[5]
        self.sequence = self.fields[9]
        self.quality_string = self.fields[10]
        self.options = self.fields[11:]
        self.orientation = self._determine_orientation()
        self.occurences = self._parse_occurences()
        self._reference_covered_length = self._compute_reference_covered_length()
        self._ref_to_read = self._build_ref_to_read_map()
        self._quality = Qual(self.quality_string, "sanger")
        self.amplicon = self._parse_amplicon()

    def _determine_orientation(self) -> str:
        """Return whether the read is mapped to the forward or reverse strand."""
        if self.flag & 16:
            return "R"
        return "F"

    def _parse_occurences(self) -> int:
        """Detect if the identifier encodes a multiplicity suffix and return count."""
        parts = self.identifier.split("_")
        tail = parts[-1]
        if tail.isdigit():
            return max(1, int(tail))
        return 1

    def _compute_reference_covered_length(self) -> int:
        """Sum the CIGAR operations that consume reference bases to get coverage."""
        ops = _CIGAR_RE.findall(self.cigar)
        total = 0
        for length, op in ops:
            if op in "MDN=X":
                total += int(length)
        return total

    def _build_ref_to_read_map(self) -> dict[int, int | None]:
        """Walk the CIGAR string and build a mapping from reference positions to read positions.

        CIGAR operations:
        - M/=/X: consume both reference and query
        - I: consume query only (skip in reference mapping)
        - D/N: consume reference only (map to None = gap)
        - S: consume query only (skip in reference mapping)
        - H/P: consume neither
        """
        ref_to_read: dict[int, Optional[int]] = {}
        ops = _CIGAR_RE.findall(self.cigar)
        ref_pos = self.coordinate
        read_pos = 0
        for length_str, op in ops:
            length = int(length_str)
            if op in "M=X":
                for i in range(length):
                    ref_to_read[ref_pos + i] = read_pos + i
                ref_pos += length
                read_pos += length
            elif op == "I":
                read_pos += length
            elif op in "DN":
                for i in range(length):
                    ref_to_read[ref_pos + i] = None
                ref_pos += length
            elif op == "S":
                read_pos += length
            # H and P consume neither
        return ref_to_read

    def _parse_amplicon(self) -> str:
        """Extract the amplicon tag embedded in the read identifier, if present."""
        match = re.search(r"(Amp_[0-9]+)", self.identifier)
        if match:
            return match.group(1)
        return "Amp_NONE"

    def covers(self, position: int) -> bool:
        """Check whether this read spans the requested reference position."""
        return self.coordinate <= position <= self.coordinate + self._reference_covered_length - 1

    def codon_at(self, position: int) -> str | None:
        """Return the codon at the given reference position using the CIGAR-aware mapping."""
        bases = []
        for offset in range(3):
            ref_pos = position + offset
            read_idx = self._ref_to_read.get(ref_pos)
            if read_idx is None and ref_pos in self._ref_to_read:
                bases.append("-")
            elif read_idx is None:
                return None
            else:
                bases.append(self.sequence[read_idx])
        return "".join(bases)

    def codon_qualities_at(self, position: int) -> list[float] | None:
        """Return Phred-like qualities for the three codon bases, or None if any is missing/gapped."""
        quals: list[float] = []
        for offset in range(3):
            ref_pos = position + offset
            read_idx = self._ref_to_read.get(ref_pos)
            if read_idx is None:
                return None
            quals.append(self._quality.return_quality_by_pos(read_idx))
        return quals

    def is_mapped(self) -> bool:
        """True when the read is not flagged as unmapped."""
        return not (self.flag & 4)
