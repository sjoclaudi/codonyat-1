from __future__ import annotations

import csv
import logging
import math
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from pathlib import Path

from .genetic_code import codon_to_aminoacid
from .models import Protein, Variant
from .reference import FullReference
from .sam import SamEntry

logger = logging.getLogger(__name__)


class SamContainer:
    """Aggregates SAM entries and computes RT variant statistics for an amplicon set."""

    def __init__(
        self,
        sam_path: Path,
        ratio_upper: float,
        ratio_lower: float,
        entropy_threshold: float,
        min_baseq: float | None = None,
        collapsed_read_counts: bool = False,
    ) -> None:
        """Initialize the container with the SAM path, ratio window, and entropy threshold."""
        self.sam_path = sam_path
        self.reads: list[SamEntry] = []
        self.variants: dict[int, dict] = {}
        self.protein_variants: dict[str, dict[int, dict]] = {}
        self._protein_result_list: list[tuple[str, dict[int, dict]]] = []
        self.ratio_upper = ratio_upper
        self.ratio_lower = ratio_lower
        self.entropy_threshold = entropy_threshold
        self.min_baseq = min_baseq
        self.collapsed_read_counts = collapsed_read_counts

    def load(self) -> None:
        """Load every mapped read from the SAM file so they can be aggregated."""
        with open(self.sam_path) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                entry = SamEntry(line, collapsed_read_counts=self.collapsed_read_counts)
                if entry.is_mapped():
                    self.reads.append(entry)

    def calculate_variant_frequencies(
        self, protein_name: str, reference: FullReference, record_id: str | None = None
    ) -> None:
        """Walk every position of the requested protein and build variant histories.

        Equivalent to ``calculate_variant_frequencies_multi([protein_name], reference)``
        when ``record_id`` is omitted. Supports multi-interval proteins
        (frameshift, splicing, circular).
        """
        protein = reference.get_protein(protein_name, record_id)
        if protein is None:
            raise ValueError(f"Protein {protein_name} not found in reference")
        self.calculate_for_proteins([protein], reference)

    def calculate_variant_frequencies_multi(
        self, protein_names: Sequence[str], reference: FullReference
    ) -> None:
        """Count codons for several proteins in a single pass over loaded reads.

        Results are stored per protein in ``self.protein_variants``;
        ``self.variants`` points at the first requested protein (the only one in
        single-protein runs), as in 1.1.0.
        """
        if not protein_names:
            raise ValueError("At least one protein is required")
        proteins: list[Protein] = []
        for name in protein_names:
            protein = reference.get_protein(name)
            if protein is None:
                raise ValueError(f"Protein {name} not found in reference")
            proteins.append(protein)
        self.calculate_for_proteins(proteins, reference)

    def calculate_for_proteins(self, proteins: Sequence[Protein], reference: FullReference) -> None:
        """Count every requested protein in one pass over ``self.reads``.

        Supports multi-record references (reads are matched by RNAME) and
        multi-interval proteins.
        """
        if not proteins:
            raise ValueError("At least one protein is required")

        jobs: list[tuple[Protein, dict[int, dict]]] = []
        for protein in proteins:
            entries: dict[int, dict] = {}
            aa_pos = 0
            for pos, record_id, codon_coords in self._get_codon_positions(protein, reference):
                aa_pos += 1
                entries[pos] = {
                    "variants": {},
                    "fw_depth": 0,
                    "rv_depth": 0,
                    "depth": 0,
                    "amplicons": {},
                    "record_id": record_id,
                    "aa_pos": aa_pos,
                    "codon_coords": codon_coords,
                }
            jobs.append((protein, entries))

        for read in self.reads:
            for protein, entries in jobs:
                if read.reference != protein.record_id:
                    continue
                for pos, entry in entries.items():
                    self._accumulate_read_at_position(read, pos, entry)

        self.protein_variants = {}
        self._protein_result_list = []
        for protein, entries in jobs:
            self.variants = entries
            for pos in self.variants:
                self.calculate_ratios_nucleotide(pos)
                self.calculate_ratios_aminoacid(pos)
                self.calculate_shannon_entropy_single_position(pos)
                self.log_position_stats(pos)
            self.protein_variants[protein.name] = self.variants
            self._protein_result_list.append((protein.name, self.variants))
        self.variants = jobs[0][1]

    def _accumulate_read_at_position(self, read: SamEntry, pos: int, entry: dict) -> None:
        """Add one mapped read's codon to a protein position, if it covers it."""
        codon_coords = entry.get("codon_coords")
        if codon_coords is None:
            codon_coords = (pos, pos + 1, pos + 2)
        if not any(read.covers(coord) for coord in codon_coords):
            return
        codon = read.codon_at(pos, positions=codon_coords)
        if not codon or len(codon) < 3:
            return
        if self.min_baseq is not None:
            quals = read.codon_qualities_at(pos, positions=codon_coords)
            if quals is None or any(qual < self.min_baseq for qual in quals):
                return
        depth_increment = read.occurences
        variants: dict[str, Variant] = entry["variants"]
        if read.orientation == "F":
            entry["fw_depth"] += depth_increment
        else:
            entry["rv_depth"] += depth_increment
        entry["depth"] += depth_increment
        variant = variants.setdefault(codon, Variant(codon=codon))
        variant.count += depth_increment
        if read.orientation == "F":
            variant.fw_reads += depth_increment
        else:
            variant.rv_reads += depth_increment
        amplicon_label = read.amplicon or "Amp_NONE"
        amp_variant = variant.amplicons.setdefault(amplicon_label, {"fw": 0, "rv": 0})
        if read.orientation == "F":
            amp_variant["fw"] += depth_increment
        else:
            amp_variant["rv"] += depth_increment
        amp_position = entry["amplicons"].setdefault(
            amplicon_label, {"fw_depth": 0, "rv_depth": 0, "depth": 0}
        )
        if read.orientation == "F":
            amp_position["fw_depth"] += depth_increment
        else:
            amp_position["rv_depth"] += depth_increment
        amp_position["depth"] += depth_increment

    def _get_codon_positions(self, protein, reference: FullReference) -> list[tuple]:
        """Generate (start, record_id, codon_coords) for every codon in the protein.

        Coding bases are concatenated in interval order so the reading frame is
        carried across parts. A codon that sits on a junction keeps the actual
        part coordinates (for example positions 10, 11, 12 of 1-10;11-30).
        Incomplete trailing bases at the end of the protein are dropped.
        """
        coding_bases: list[int] = []
        for interval in protein.intervals:
            coding_bases.extend(range(interval.start, interval.end + 1))

        positions = []
        for i in range(0, len(coding_bases) // 3 * 3, 3):
            codon_coords = (coding_bases[i], coding_bases[i + 1], coding_bases[i + 2])
            positions.append((codon_coords[0], protein.record_id, codon_coords))
        return positions

    def calculate_ratios_nucleotide(self, position: int) -> None:
        """Determine strand-specific ratios per variant and flag the balanced reads."""
        entry = self.variants.get(position)
        if not entry or entry["fw_depth"] == 0 or entry["rv_depth"] == 0:
            return
        for variant in entry["variants"].values():
            # compute the per-strand frequency for each codon then derive the ratio
            variant.fw_freq = variant.fw_reads / entry["fw_depth"] if entry["fw_depth"] else 0
            variant.rv_freq = variant.rv_reads / entry["rv_depth"] if entry["rv_depth"] else 0
            variant.ratio = (variant.fw_freq / variant.rv_freq) if variant.rv_freq else 0
            # only mark variants within the balanced ratio window
            if self.ratio_lower <= variant.ratio <= self.ratio_upper:
                variant.balanced_fw += variant.fw_reads
                variant.balanced_rv += variant.rv_reads

    def has_balanced_variants(self, position: int) -> bool:
        """Return True when at least one variant falls within the configured strand-ratio window."""
        entry = self.variants.get(position)
        if not entry:
            return False
        return any(self.ratio_lower <= variant.ratio <= self.ratio_upper for variant in entry["variants"].values())

    def calculate_ratios_aminoacid(self, position: int) -> None:
        """Aggregate amino-acid frequencies across codon variants and capture balance metadata."""
        entry = self.variants.get(position)
        if not entry:
            return
        aa_stats: dict[str, dict[str, float]] = {}
        for variant in entry["variants"].values():
            aa = codon_to_aminoacid(variant.codon)
            stats = aa_stats.setdefault(aa, {"count": 0.0, "fw": 0.0, "rv": 0.0})
            stats["count"] += variant.count
            stats["fw"] += variant.fw_reads
            stats["rv"] += variant.rv_reads
        for aa, stats in aa_stats.items():
            fw = stats["fw"]
            rv = stats["rv"]
            ratio = (fw / rv) if rv else 0
            # flag balanced amino-acid counts using the same ratio window as nucleotides
            stats["ratio"] = ratio
            stats["balanced"] = (fw + rv) if self.ratio_lower <= ratio <= self.ratio_upper else 0
        entry["aminoacid_stats"] = aa_stats

    def calculate_shannon_entropy_single_position(self, position: int) -> float:
        """Compute the Shannon entropy for the observed codon mixture at a position."""
        entry = self.variants.get(position)
        if not entry:
            return 0.0
        depth = entry.get("depth", 0)
        if depth == 0:
            entry["shannon_entropy"] = 0.0
            return 0.0
        entropy = 0.0
        for variant in entry["variants"].values():
            if variant.count == 0:
                continue
            p = variant.count / depth
            # standard Shannon formula: sum of -p log(p) for each variant
            entropy -= p * math.log(p)
        entry["shannon_entropy"] = entropy
        entry["entropy_pass"] = entropy >= self.entropy_threshold
        return entropy

    def log_position_stats(self, position: int) -> None:
        """Log strand balance and entropy so downstream debugging is easier."""
        entry = self.variants.get(position)
        if not entry:
            return
        balanced = [v.codon for v in entry["variants"].values() if self.ratio_lower <= v.ratio <= self.ratio_upper]
        entropy = entry.get("shannon_entropy", 0.0)
        logger.debug(
            "Position %s -- depth=%s, fw=%s, rv=%s, entropy=%.3f, balanced=%s",
            position,
            entry["depth"],
            entry["fw_depth"],
            entry["rv_depth"],
            entropy,
            balanced,
        )

    def _entries_for(self, protein_name: str) -> dict[int, dict]:
        """Per-position results for a protein (falls back to ``self.variants``)."""
        return self.protein_variants.get(protein_name, self.variants)

    def _csv_fieldnames(self, reference: FullReference) -> list[str]:
        has_multiple_records = len(reference.record_ids) > 1
        fieldnames = ["FILE", "REFERENCE"]
        if has_multiple_records:
            fieldnames.append("RECORD")
        fieldnames.extend(
            ["PROTEIN", "VARIANT", "POSITION", "FREQ", "FWCOV", "RVCOV", "TOTALCOV", "RATIO"]
        )
        return fieldnames

    def _write_csv_rows(
        self,
        writer: csv.DictWriter,
        sample_name: str,
        reference: FullReference,
        protein_name: str,
        entries: dict[int, dict] | None = None,
    ) -> None:
        has_multiple_records = len(reference.record_ids) > 1
        if entries is None:
            entries = self._entries_for(protein_name)
        for pos, entry in sorted(entries.items()):
            depth = entry["depth"]
            ratio = entry["fw_depth"] / entry["rv_depth"] if entry["rv_depth"] else 0
            record_id = entry.get("record_id", reference.id)
            for variant in entry["variants"].values():
                row = {
                    "FILE": sample_name,
                    "REFERENCE": record_id,
                    "PROTEIN": protein_name,
                    "VARIANT": variant.codon,
                    "POSITION": pos,
                    "FREQ": round((variant.count / depth * 100) if depth else 0, 3),
                    "FWCOV": entry["fw_depth"],
                    "RVCOV": entry["rv_depth"],
                    "TOTALCOV": depth,
                    "RATIO": round(ratio, 3),
                }
                if has_multiple_records:
                    row["RECORD"] = record_id
                writer.writerow(row)

    def write_csv(
        self, sample_name: str, reference: FullReference, output_path: Path, *, protein_name: str = "RT"
    ) -> None:
        """Emit the variant summary as a TSV file that mirrors the Perl outputs.

        For multi-record references, adds a RECORD column after REFERENCE.
        Single-record outputs remain byte-identical to 1.1.0.
        """
        self.write_csv_proteins(sample_name, reference, output_path, [protein_name])

    def write_csv_proteins(
        self,
        sample_name: str,
        reference: FullReference,
        output_path: Path,
        protein_names: Sequence[str],
    ) -> None:
        """Emit one TSV with a single header and the rows of each protein, in the given order."""
        fieldnames = self._csv_fieldnames(reference)
        pairs = self._output_pairs(protein_names)
        with open(output_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for protein_name, entries in pairs:
                self._write_csv_rows(writer, sample_name, reference, protein_name, entries)

    def _output_pairs(self, protein_names: Sequence[str]) -> list[tuple[str, dict[int, dict]]]:
        """Map requested names to counted entries, preserving duplicate-name order."""
        calculated = [name for name, _entries in self._protein_result_list]
        if self._protein_result_list and calculated == list(protein_names):
            return list(self._protein_result_list)
        used: set[int] = set()
        pairs: list[tuple[str, dict[int, dict]]] = []
        for name in protein_names:
            found = False
            for index, (result_name, entries) in enumerate(self._protein_result_list):
                if result_name == name and index not in used:
                    used.add(index)
                    pairs.append((result_name, entries))
                    found = True
                    break
            if not found:
                pairs.append((name, self._entries_for(name)))
        return pairs

    @staticmethod
    def _add_positions(parent: ET.Element, entries: dict[int, dict]) -> None:
        for pos, entry in sorted(entries.items()):
            position_elem = ET.SubElement(parent, "Position", index=str(pos))
            ET.SubElement(position_elem, "Depth").text = str(entry["depth"])
            ET.SubElement(position_elem, "FwCover").text = str(entry["fw_depth"])
            ET.SubElement(position_elem, "RvCover").text = str(entry["rv_depth"])
            variants_elem = ET.SubElement(position_elem, "Variants")
            for variant in entry["variants"].values():
                # describe each variant using XML attributes for diagnostics
                var_elem = ET.SubElement(variants_elem, "Variant", codon=variant.codon)
                var_elem.set("count", str(variant.count))
                var_elem.set("fw_reads", str(variant.fw_reads))
                var_elem.set("rv_reads", str(variant.rv_reads))

    def write_xml(self, sample_name: str, reference: FullReference, output_path: Path) -> None:
        """Dump the per-position state of the (first) protein into XML for diagnostics."""
        root = ET.Element("SamContainer", sample=sample_name, reference=reference.id)
        self._add_positions(root, self.variants)
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding="utf-8", xml_declaration=True)

    def write_xml_proteins(
        self,
        sample_name: str,
        reference: FullReference,
        output_path: Path,
        protein_names: Sequence[str],
    ) -> None:
        """Multi-protein XML: one ``<Protein name=...>`` element per protein wrapping its positions.

        For a single protein use :meth:`write_xml`, whose layout is unchanged.
        """
        root = ET.Element("SamContainer", sample=sample_name, reference=reference.id)
        for protein_name, entries in self._output_pairs(protein_names):
            protein_elem = ET.SubElement(root, "Protein", name=protein_name)
            self._add_positions(protein_elem, entries)
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding="utf-8", xml_declaration=True)

    def write_stats(
        self,
        sample_name: str,
        reference: FullReference,
        output_path: Path,
        protein_names: Sequence[str],
        *,
        min_freq: float,
        min_coverage: int,
        sb_alpha: float,
    ) -> None:
        """Write the opt-in stats TSV (Wilson CI, strand bias, diversity)."""
        from .stats_writer import write_stats_tsv

        write_stats_tsv(
            sample_name,
            reference,
            output_path,
            self._output_pairs(protein_names),
            min_freq=min_freq,
            min_coverage=min_coverage,
            sb_alpha=sb_alpha,
        )

    def write_aavf(
        self,
        sample_name: str,
        reference: FullReference,
        output_path: Path,
        protein_names: Sequence[str],
        *,
        min_freq: float,
        min_coverage: int,
        sb_alpha: float,
        file_date: str | None = None,
    ) -> None:
        """Write the opt-in AAVF v1.0 file. Research use only; no CAT/SRVL."""
        from .aavf import write_aavf as _write_aavf

        _write_aavf(
            sample_name,
            reference,
            output_path,
            self._output_pairs(protein_names),
            min_freq=min_freq,
            min_coverage=min_coverage,
            sb_alpha=sb_alpha,
            file_date=file_date,
        )
