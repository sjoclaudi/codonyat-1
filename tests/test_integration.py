"""Integration tests for the full variant-calling pipeline and runner API."""

from __future__ import annotations

import csv
import math
import xml.etree.ElementTree as ET
from argparse import Namespace
from pathlib import Path

import pytest

from aa_caller.app import (
    DEFAULT_ENTROPY_THRESHOLD,
    DEFAULT_RATIO_LOWER,
    DEFAULT_RATIO_UPPER,
    FullReference,
    SamContainer,
)
from aa_caller.runner import call_variants, call_variants_from_args

# ---------------------------------------------------------------------------
# Helpers to create synthetic test data
# ---------------------------------------------------------------------------

# A short reference: 30 bases = 10 codons, positions 1-30
_REF_SEQ = "ATGCCCTTTAAAGGGCCCAAATTTGGGCCC"

# SAM fields template: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
_SAM_TEMPLATE = "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"


def _make_fasta(path: Path, protein_header: str = "RT(ReverseTranscriptase):1-30") -> Path:
    """Write a minimal FASTA file with protein annotation.

    Uses the standard format: >seqid protein(description):start-end
    where seqid is separated from the protein annotation by a space.
    """
    path.write_text(f">testref {protein_header}\n{_REF_SEQ}\n")
    return path


def _make_amplicon_tsv(path: Path) -> Path:
    """Write a minimal amplicon TSV."""
    lines = [
        "label\tprotein\treference\t5p\t3p\tstart\tend",
        "Amp_1\tRT\ttestref\tATG\tCCC\t1\t30",
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def _qual(length: int) -> str:
    """Generate a quality string of given length (all high-quality)."""
    return "I" * length


def _make_sam(path: Path, reads: list[str] | None = None) -> Path:
    """Write a SAM file.  If *reads* is None, use a default set of 8 reads."""
    if reads is None:
        reads = _default_reads()
    header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:testref\tLN:30\n"
    path.write_text(header + "".join(reads))
    return path


def _default_reads() -> list[str]:
    """8 reads covering positions 1-30 (the full RT protein).

    Mix of forward (flag 0) and reverse (flag 16), all Amp_1.
    All have simple 30M CIGARs matching the full reference.
    Reads 1-4 carry the reference codon ATG at pos 1.
    Reads 5-6 carry a variant TGA at pos 1.
    Reads 7-8 cover the same region but are reverse strand.
    """
    seq_ref = _REF_SEQ  # reference sequence
    seq_var = "TGA" + _REF_SEQ[3:]  # variant at codon 1
    q = _qual(30)
    return [
        # Forward reads with reference codons
        _SAM_TEMPLATE.format(qname="Amp_1_read1", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_ref, qual=q),
        _SAM_TEMPLATE.format(qname="Amp_1_read2", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_ref, qual=q),
        # Reverse reads with reference codons
        _SAM_TEMPLATE.format(qname="Amp_1_read3", flag=16, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_ref, qual=q),
        _SAM_TEMPLATE.format(qname="Amp_1_read4", flag=16, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_ref, qual=q),
        # Forward reads with variant codon at pos 1
        _SAM_TEMPLATE.format(qname="Amp_1_read5", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_var, qual=q),
        # Reverse reads with variant codon at pos 1
        _SAM_TEMPLATE.format(qname="Amp_1_read6", flag=16, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_var, qual=q),
        # Additional forward/reverse with reference
        _SAM_TEMPLATE.format(qname="Amp_1_read7", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_ref, qual=q),
        _SAM_TEMPLATE.format(qname="Amp_1_read8", flag=16, rname="testref", pos=1, mapq=30, cigar="30M", seq=seq_ref, qual=q),
    ]


# ---------------------------------------------------------------------------
# Tests: SamContainer.load()
# ---------------------------------------------------------------------------

class TestSamContainerLoad:
    def test_load_skips_headers_and_unmapped(self, tmp_path: Path) -> None:
        sam_path = _make_sam(tmp_path / "reads.sam", reads=[
            "@HD\tVN:1.6\n",  # extra header line (already have one in _make_sam header)
            _SAM_TEMPLATE.format(qname="mapped", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=_REF_SEQ, qual=_qual(30)),
            _SAM_TEMPLATE.format(qname="unmapped", flag=4, rname="*", pos=0, mapq=0, cigar="*", seq=_REF_SEQ, qual=_qual(30)),
        ])
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        assert len(container.reads) == 1
        assert container.reads[0].identifier == "mapped"

    def test_load_reads_correct_count(self, tmp_path: Path) -> None:
        sam_path = _make_sam(tmp_path / "reads.sam")
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        assert len(container.reads) == 8


# ---------------------------------------------------------------------------
# Tests: calculate_variant_frequencies
# ---------------------------------------------------------------------------

class TestVariantFrequencies:
    def _run(self, tmp_path: Path) -> SamContainer:
        sam_path = _make_sam(tmp_path / "reads.sam")
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        return container

    def test_variant_counts_at_position_1(self, tmp_path: Path) -> None:
        container = self._run(tmp_path)
        entry = container.variants[1]
        variants = entry["variants"]
        # ATG: 3 fw (read1, read2, read7) + 3 rv (read3, read4, read8) = 6
        assert variants["ATG"].count == 6
        # TGA: 1 fw (read5) + 1 rv (read6) = 2
        assert variants["TGA"].count == 2
        assert entry["depth"] == 8

    def test_fw_rv_depth(self, tmp_path: Path) -> None:
        container = self._run(tmp_path)
        entry = container.variants[1]
        assert entry["fw_depth"] == 4  # read1, read2, read5, read7
        assert entry["rv_depth"] == 4  # read3, read4, read6, read8

    def test_all_positions_computed(self, tmp_path: Path) -> None:
        container = self._run(tmp_path)
        # Protein is 1-30, step 3, so positions 1, 4, 7, 10, 13, 16, 19, 22, 25, 28
        assert sorted(container.variants.keys()) == [1, 4, 7, 10, 13, 16, 19, 22, 25, 28]


# ---------------------------------------------------------------------------
# Tests: ratio calculations
# ---------------------------------------------------------------------------

class TestRatios:
    def test_balanced_variant(self, tmp_path: Path) -> None:
        """With equal fw/rv, ratio should be 1.0 — well within the default window."""
        sam_path = _make_sam(tmp_path / "reads.sam")
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        assert container.has_balanced_variants(1)

    def test_unbalanced_variant(self, tmp_path: Path) -> None:
        """All reads forward → ratio undefined (rv=0), should not be balanced."""
        q = _qual(30)
        reads = [
            _SAM_TEMPLATE.format(qname="Amp_1_r1", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=_REF_SEQ, qual=q),
            _SAM_TEMPLATE.format(qname="Amp_1_r2", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=_REF_SEQ, qual=q),
        ]
        sam_path = _make_sam(tmp_path / "reads.sam", reads=reads)
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        assert not container.has_balanced_variants(1)


# ---------------------------------------------------------------------------
# Tests: Shannon entropy
# ---------------------------------------------------------------------------

class TestEntropy:
    def test_entropy_single_variant_is_zero(self, tmp_path: Path) -> None:
        """A single codon at 100% frequency should have entropy 0."""
        q = _qual(30)
        reads = [
            _SAM_TEMPLATE.format(qname="Amp_1_r1", flag=0, rname="testref", pos=1, mapq=30, cigar="30M", seq=_REF_SEQ, qual=q),
            _SAM_TEMPLATE.format(qname="Amp_1_r2", flag=16, rname="testref", pos=1, mapq=30, cigar="30M", seq=_REF_SEQ, qual=q),
        ]
        sam_path = _make_sam(tmp_path / "reads.sam", reads=reads)
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        assert container.variants[1]["shannon_entropy"] == pytest.approx(0.0)

    def test_entropy_two_variants_positive(self, tmp_path: Path) -> None:
        """Two codons should yield positive entropy."""
        sam_path = _make_sam(tmp_path / "reads.sam")  # default reads have ATG + TGA
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        # 6 ATG + 2 TGA = 8 total → entropy = -(6/8)*ln(6/8) - (2/8)*ln(2/8)
        p1, p2 = 6 / 8, 2 / 8
        expected = -(p1 * math.log(p1) + p2 * math.log(p2))
        assert container.variants[1]["shannon_entropy"] == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Tests: CSV output
# ---------------------------------------------------------------------------

class TestCSVOutput:
    def test_csv_columns_and_content(self, tmp_path: Path) -> None:
        sam_path = _make_sam(tmp_path / "reads.sam")
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        csv_path = tmp_path / "output.tsv"
        container.write_csv("sample", ref, csv_path)

        with open(csv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        assert len(rows) > 0
        expected_cols = {"FILE", "REFERENCE", "PROTEIN", "VARIANT", "POSITION", "FREQ", "FWCOV", "RVCOV", "TOTALCOV", "RATIO"}
        assert expected_cols == set(rows[0].keys())

        # Check that position 1 has both ATG and TGA variants
        pos1_rows = [r for r in rows if r["POSITION"] == "1"]
        codons = {r["VARIANT"] for r in pos1_rows}
        assert "ATG" in codons
        assert "TGA" in codons


# ---------------------------------------------------------------------------
# Tests: XML output
# ---------------------------------------------------------------------------

class TestXMLOutput:
    def test_xml_valid_and_has_elements(self, tmp_path: Path) -> None:
        sam_path = _make_sam(tmp_path / "reads.sam")
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        xml_path = tmp_path / "output.xml"
        container.write_xml("sample", ref, xml_path)

        tree = ET.parse(xml_path)
        root = tree.getroot()
        assert root.tag == "SamContainer"
        assert root.attrib["sample"] == "sample"

        positions = root.findall("Position")
        assert len(positions) == 10  # 10 codon positions

        # Check the first position has variants
        first_pos = positions[0]
        variants = first_pos.find("Variants")
        assert variants is not None
        variant_elems = variants.findall("Variant")
        assert len(variant_elems) >= 1


# ---------------------------------------------------------------------------
# Tests: runner.py API
# ---------------------------------------------------------------------------

class TestRunnerAPI:
    def test_call_variants_end_to_end(self, tmp_path: Path) -> None:
        sam_path = _make_sam(tmp_path / "reads.sam")
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        amp_path = _make_amplicon_tsv(tmp_path / "amps.tsv")
        csv_out = tmp_path / "result.tsv"
        xml_out = tmp_path / "result.xml"

        result = call_variants(
            sam_path=sam_path,
            reference_path=fasta_path,
            amplicons_path=amp_path,
            csv_path=csv_out,
            xml_path=xml_out,
        )

        assert result.csv_path == csv_out
        assert result.xml_path == xml_out
        assert csv_out.exists()
        assert xml_out.exists()
        assert len(result.container.reads) == 8
        assert "Amp_1" in result.amplicons

    def test_call_variants_from_args_with_namespace(self, tmp_path: Path) -> None:
        sam_path = _make_sam(tmp_path / "reads.sam")
        fasta_path = _make_fasta(tmp_path / "ref.fasta")
        amp_path = _make_amplicon_tsv(tmp_path / "amps.tsv")
        csv_out = tmp_path / "result.tsv"
        xml_out = tmp_path / "result.xml"

        args = Namespace(
            sam_file=sam_path,
            reference_file=fasta_path,
            amplicons_file=amp_path,
            ratio_upper=DEFAULT_RATIO_UPPER,
            ratio_lower=DEFAULT_RATIO_LOWER,
            entropy_threshold=DEFAULT_ENTROPY_THRESHOLD,
        )
        result = call_variants_from_args(args, csv_path=csv_out, xml_path=xml_out)

        assert result.csv_path == csv_out
        assert result.xml_path == xml_out
        assert csv_out.exists()
        assert xml_out.exists()
        assert len(result.container.variants) == 10

    def test_non_rt_protein(self, tmp_path: Path) -> None:
        """Verify that a non-RT protein can be analyzed via the --protein flag."""
        sam_path = _make_sam(tmp_path / "reads.sam")
        fasta_path = _make_fasta(tmp_path / "ref.fasta", protein_header="POL(Polymerase):1-30")
        amp_path = _make_amplicon_tsv(tmp_path / "amps.tsv")
        csv_out = tmp_path / "result.tsv"
        xml_out = tmp_path / "result.xml"

        result = call_variants(
            sam_path=sam_path,
            reference_path=fasta_path,
            amplicons_path=amp_path,
            protein="POL",
            csv_path=csv_out,
            xml_path=xml_out,
        )

        assert csv_out.exists()
        assert xml_out.exists()
        assert len(result.container.variants) == 10
        # Verify the CSV uses the correct protein name
        with open(csv_out) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                assert row["PROTEIN"] == "POL"
                break
