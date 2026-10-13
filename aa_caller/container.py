from __future__ import annotations

import csv
import logging
import math
import xml.etree.ElementTree as ET
from pathlib import Path

from .genetic_code import codon_to_aminoacid
from .models import Variant
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
    ) -> None:
        """Initialize the container with the SAM path, ratio window, and entropy threshold."""
        self.sam_path = sam_path
        self.reads: list[SamEntry] = []
        self.variants: dict[int, dict[str, Variant]] = {}
        self.ratio_upper = ratio_upper
        self.ratio_lower = ratio_lower
        self.entropy_threshold = entropy_threshold

    def load(self) -> None:
        """Load every mapped read from the SAM file so they can be aggregated."""
        with open(self.sam_path) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                entry = SamEntry(line)
                if entry.is_mapped():
                    self.reads.append(entry)

    def calculate_variant_frequencies(self, protein_name: str, reference: FullReference, record_id: str = None) -> None:
        """Walk every position of the requested protein and build variant histories.

        Supports multi-interval proteins (frameshift, splicing, circular).
        """
        protein = reference.get_protein(protein_name, record_id)
        if protein is None:
            raise ValueError(f"Protein {protein_name} not found in reference")

        codon_positions = self._get_codon_positions(protein, reference)

        for pos, record_id in codon_positions:
            variants: dict[str, Variant] = {}
            fw_depth = 0
            rv_depth = 0
            amplicons_info: dict[str, dict[str, int]] = {}
            for read in self.reads:
                if read.reference != record_id:
                    continue
                if not read.covers(pos):
                    continue
                codon = read.codon_at(pos)
                if not codon or len(codon) < 3:
                    continue
                depth_increment = read.occurences
                if read.orientation == "F":
                    fw_depth += depth_increment
                else:
                    rv_depth += depth_increment
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
                amp_position = amplicons_info.setdefault(amplicon_label, {"fw_depth": 0, "rv_depth": 0, "depth": 0})
                if read.orientation == "F":
                    amp_position["fw_depth"] += depth_increment
                else:
                    amp_position["rv_depth"] += depth_increment
                amp_position["depth"] += depth_increment
            entry = {
                "variants": variants,
                "fw_depth": fw_depth,
                "rv_depth": rv_depth,
                "depth": fw_depth + rv_depth,
                "amplicons": amplicons_info,
                "record_id": record_id,
            }
            self.variants[pos] = entry
            self.calculate_ratios_nucleotide(pos)
            self.calculate_ratios_aminoacid(pos)
            self.calculate_shannon_entropy_single_position(pos)
            self.log_position_stats(pos)

    def _get_codon_positions(self, protein, reference: FullReference) -> list[tuple]:
        """Generate list of (position, record_id) tuples for all codons in the protein.

        Handles multi-interval proteins including frameshifts.
        """
        positions = []
        current_offset = 0

        for interval in protein.intervals:
            interval_length = interval.end - interval.start + 1
            pos = interval.start

            while pos <= interval.end:
                remaining_in_interval = interval.end - pos + 1
                if remaining_in_interval >= 3:
                    positions.append((pos, protein.record_id))
                    pos += 3
                else:
                    break

            current_offset += interval_length

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

    def write_csv(self, sample_name: str, reference: FullReference, output_path: Path, *, protein_name: str = "RT") -> None:
        """Emit the variant summary as a TSV file that mirrors the Perl outputs.

        For multi-record references, adds a RECORD column after REFERENCE.
        Single-record outputs remain byte-identical to 1.1.0.
        """
        has_multiple_records = len(reference.record_ids) > 1

        fieldnames = ["FILE", "REFERENCE"]
        if has_multiple_records:
            fieldnames.append("RECORD")
        fieldnames.extend(["PROTEIN", "VARIANT", "POSITION", "FREQ", "FWCOV", "RVCOV", "TOTALCOV", "RATIO"])

        with open(output_path, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for pos, entry in sorted(self.variants.items()):
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

    def write_xml(self, sample_name: str, reference: FullReference, output_path: Path) -> None:
        """Dump the internal state into XML for downstream diagnostics."""
        root = ET.Element("SamContainer", sample=sample_name, reference=reference.id)
        for pos, entry in sorted(self.variants.items()):
            position_elem = ET.SubElement(root, "Position", index=str(pos))
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
        tree = ET.ElementTree(root)
        tree.write(output_path, encoding="utf-8", xml_declaration=True)
