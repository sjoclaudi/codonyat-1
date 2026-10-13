"""Tests for codonyat 1.2.0 features: multi-record, multi-interval, TSV annotation."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from aa_caller.app import (
    DEFAULT_ENTROPY_THRESHOLD,
    DEFAULT_RATIO_LOWER,
    DEFAULT_RATIO_UPPER,
    FullReference,
    SamContainer,
)
from aa_caller.models import Protein

# ---------------------------------------------------------------------------
# Test multi-interval protein parsing
# ---------------------------------------------------------------------------

class TestMultiIntervalProtein:
    def test_single_interval_backward_compat(self) -> None:
        """Single interval should work as before."""
        protein = Protein.from_string("RT(Reverse Transcriptase):100-200", record_id="ref1")
        assert protein.name == "RT"
        assert protein.description == "Reverse Transcriptase"
        assert protein.record_id == "ref1"
        assert len(protein.intervals) == 1
        assert protein.intervals[0].start == 100
        assert protein.intervals[0].end == 200
        assert protein.start_coordinate == 100
        assert protein.end_coordinate == 200

    def test_multi_interval_frameshift(self) -> None:
        """Parse frameshift with semicolon-separated intervals."""
        protein = Protein.from_string("ORF1ab(Polyprotein):266-13468;13468-21555", record_id="MN908947.3")
        assert protein.name == "ORF1ab"
        assert len(protein.intervals) == 2
        assert protein.intervals[0].start == 266
        assert protein.intervals[0].end == 13468
        assert protein.intervals[1].start == 13468
        assert protein.intervals[1].end == 21555
        assert protein.start_coordinate == 266
        assert protein.end_coordinate == 21555

    def test_multi_interval_circular(self) -> None:
        """Parse circular origin-crossing protein."""
        protein = Protein.from_string("P(Polymerase):2307-3221;1-1623", record_id="NC_003977.2")
        assert protein.name == "P"
        assert len(protein.intervals) == 2
        assert protein.intervals[0].start == 2307
        assert protein.intervals[0].end == 3221
        assert protein.intervals[1].start == 1
        assert protein.intervals[1].end == 1623


# ---------------------------------------------------------------------------
# Test multi-record FASTA
# ---------------------------------------------------------------------------

class TestMultiRecordFASTA:
    def test_single_record_backward_compat(self, tmp_path: Path) -> None:
        """Single-record FASTA should work as before."""
        fasta_path = tmp_path / "ref.fasta"
        fasta_path.write_text(">K03455 RT(Reverse Transcriptase):1-30\nATGCCCTTTAAAGGGCCCAAAT TTGGGCCC\n")
        ref = FullReference(fasta_path)
        assert ref.id == "K03455"
        assert ref.seq == "ATGCCCTTTAAAGGGCCCAAATTTGGGCCC"
        assert len(ref.record_ids) == 1
        protein = ref.get_protein("RT")
        assert protein is not None
        assert protein.name == "RT"

    def test_multi_record_fasta(self, tmp_path: Path) -> None:
        """Multi-record FASTA should load all records."""
        fasta_path = tmp_path / "ref.fasta"
        fasta_path.write_text(
            ">Segment1 GeneA(ProteinA):1-30\nATGCCCTTTAAAGGGCCCAAAT TTGGGCCC\n"
            ">Segment2 GeneB(ProteinB):1-30\nTTTGGGCCCAAAGGGTTTACGTCGATCGAT\n"
        )
        ref = FullReference(fasta_path)
        assert len(ref.record_ids) == 2
        assert "Segment1" in ref.record_ids
        assert "Segment2" in ref.record_ids
        assert len(ref.sequences) == 2

        protein_a = ref.get_protein("GeneA", "Segment1")
        assert protein_a is not None
        assert protein_a.name == "GeneA"
        assert protein_a.record_id == "Segment1"

        protein_b = ref.get_protein("GeneB", "Segment2")
        assert protein_b is not None
        assert protein_b.name == "GeneB"
        assert protein_b.record_id == "Segment2"


# ---------------------------------------------------------------------------
# Test TSV annotation format
# ---------------------------------------------------------------------------

class TestTSVAnnotation:
    def test_tsv_annotation_loading(self, tmp_path: Path) -> None:
        """TSV annotation file should load correctly."""
        fasta_path = tmp_path / "ref.fasta"
        fasta_path.write_text(">MN908947.3\nATGCCCTTTAAAGGGCCCAAAT TTGGGCCCAAAGGGAAACCCGGGTTTTTT\n")

        tsv_path = tmp_path / "annotation.tsv"
        tsv_path.write_text(
            "record\tprotein\tdescription\tintervals\tstrand\n"
            "MN908947.3\tS\tSpike glycoprotein\t1-30\t+\n"
            "MN908947.3\tN\tNucleocapsid\t31-48\t+\n"
        )

        ref = FullReference(fasta_path, annotation_path=tsv_path)
        assert len(ref.record_ids) == 1

        protein_s = ref.get_protein("S", "MN908947.3")
        assert protein_s is not None
        assert protein_s.name == "S"
        assert protein_s.description == "Spike glycoprotein"
        assert len(protein_s.intervals) == 1
        assert protein_s.intervals[0].start == 1
        assert protein_s.intervals[0].end == 30

        protein_n = ref.get_protein("N", "MN908947.3")
        assert protein_n is not None
        assert protein_n.name == "N"

    def test_tsv_annotation_multi_interval(self, tmp_path: Path) -> None:
        """TSV annotation with multi-interval proteins."""
        fasta_path = tmp_path / "ref.fasta"
        fasta_path.write_text(">MN908947.3\n" + "A" * 22000 + "\n")

        tsv_path = tmp_path / "annotation.tsv"
        tsv_path.write_text(
            "record\tprotein\tdescription\tintervals\tstrand\n"
            "MN908947.3\tORF1ab\tPolyprotein\t266-13468;13468-21555\t+\n"
        )

        ref = FullReference(fasta_path, annotation_path=tsv_path)
        protein = ref.get_protein("ORF1ab", "MN908947.3")
        assert protein is not None
        assert len(protein.intervals) == 2
        assert protein.intervals[0].start == 266
        assert protein.intervals[0].end == 13468
        assert protein.intervals[1].start == 13468
        assert protein.intervals[1].end == 21555

    def test_tsv_annotation_rejects_minus_strand(self, tmp_path: Path) -> None:
        """Minus-strand genes should be rejected with clear error."""
        fasta_path = tmp_path / "ref.fasta"
        fasta_path.write_text(">ref\nATGCCCTTTAAAGGGCCCAAAT TTGGGCCC\n")

        tsv_path = tmp_path / "annotation.tsv"
        tsv_path.write_text(
            "record\tprotein\tdescription\tintervals\tstrand\n"
            "ref\tGeneX\tSome gene\t1-30\t-\n"
        )

        with pytest.raises(ValueError, match="Minus-strand genes are not supported"):
            FullReference(fasta_path, annotation_path=tsv_path)


# ---------------------------------------------------------------------------
# Test output format for multi-record references
# ---------------------------------------------------------------------------

class TestMultiRecordOutput:
    def test_single_record_no_record_column(self, tmp_path: Path) -> None:
        """Single-record output should NOT have RECORD column (backward compat)."""
        fasta_path = tmp_path / "ref.fasta"
        fasta_path.write_text(">testref RT(ReverseTranscriptase):1-30\nATGCCCTTTAAAGGGCCCAAAT TTGGGCCC\n")

        sam_path = tmp_path / "reads.sam"
        sam_path.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:testref\tLN:30\n"
            "Amp_1_read1\t0\ttestref\t1\t30\t30M\t*\t0\t0\tATGCCCTTTAAAGGGCCCAAAT TTGGGCCC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        )

        amp_path = tmp_path / "amps.tsv"
        amp_path.write_text("label\tprotein\treference\t5p\t3p\tstart\tend\nAmp_1\tRT\ttestref\tATG\tCCC\t1\t30\n")

        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)

        csv_path = tmp_path / "output.tsv"
        container.write_csv("sample", ref, csv_path, protein_name="RT")

        with open(csv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            header = reader.fieldnames
            assert "RECORD" not in header
            assert header == ["FILE", "REFERENCE", "PROTEIN", "VARIANT", "POSITION", "FREQ", "FWCOV", "RVCOV", "TOTALCOV", "RATIO"]

    def test_multi_record_has_record_column(self, tmp_path: Path) -> None:
        """Multi-record output SHOULD have RECORD column."""
        fasta_path = tmp_path / "ref.fasta"
        fasta_path.write_text(
            ">Seg1 RT(ReverseTranscriptase):1-30\nATGCCCTTTAAAGGGCCCAAAT TTGGGCCC\n"
            ">Seg2 POL(Polymerase):1-30\nTTTGGGCCCAAAGGGTTTACGTCGATCGAT\n"
        )

        sam_path = tmp_path / "reads.sam"
        sam_path.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:Seg1\tLN:30\n"
            "Amp_1_read1\t0\tSeg1\t1\t30\t30M\t*\t0\t0\tATGCCCTTTAAAGGGCCCAAAT TTGGGCCC\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        )

        amp_path = tmp_path / "amps.tsv"
        amp_path.write_text("label\tprotein\treference\t5p\t3p\tstart\tend\nAmp_1\tRT\tSeg1\tATG\tCCC\t1\t30\n")

        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)

        csv_path = tmp_path / "output.tsv"
        container.write_csv("sample", ref, csv_path, protein_name="RT")

        with open(csv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            header = reader.fieldnames
            assert "RECORD" in header
            assert header == ["FILE", "REFERENCE", "RECORD", "PROTEIN", "VARIANT", "POSITION", "FREQ", "FWCOV", "RVCOV", "TOTALCOV", "RATIO"]

            rows = list(reader)
            assert len(rows) > 0
            assert rows[0]["RECORD"] == "Seg1"


# ---------------------------------------------------------------------------
# Backward compatibility regression test
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    def test_single_record_single_interval_byte_identical(self) -> None:
        """Single-record, single-interval output should be byte-identical to v1.1.0.
        
        Uses golden files generated from codonyat==1.1.0 from PyPI on real HIV-1 data.
        """
        import filecmp
        from pathlib import Path
        
        test_data_dir = Path(__file__).parent / "data"
        fasta_path = test_data_dir / "hiv1_ref.fasta"
        sam_path = test_data_dir / "hiv1_reads.sam"
        amp_path = test_data_dir / "hiv1_amps.tsv"
        golden_tsv = test_data_dir / "hiv1_golden_v1.1.0.tsv"
        golden_xml = test_data_dir / "hiv1_golden_v1.1.0.xml"
        
        assert fasta_path.exists(), "Test data missing: hiv1_ref.fasta"
        assert sam_path.exists(), "Test data missing: hiv1_reads.sam"
        assert amp_path.exists(), "Test data missing: hiv1_amps.tsv"
        assert golden_tsv.exists(), "Golden file missing: hiv1_golden_v1.1.0.tsv"
        assert golden_xml.exists(), "Golden file missing: hiv1_golden_v1.1.0.xml"
        
        ref = FullReference(fasta_path)
        container = SamContainer(sam_path, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD)
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            output_tsv = Path(tmpdir) / "output.tsv"
            output_xml = Path(tmpdir) / "output.xml"
            
            container.write_csv("hiv1_reads.sam", ref, output_tsv, protein_name="RT")
            container.write_xml("hiv1_reads.sam", ref, output_xml)
            
            tsv_match = filecmp.cmp(output_tsv, golden_tsv, shallow=False)
            xml_match = filecmp.cmp(output_xml, golden_xml, shallow=False)
            
            if not tsv_match:
                print(f"\n=== Expected TSV (v1.1.0) ===\n{golden_tsv.read_text()}")
                print(f"\n=== Actual TSV (v1.2.0) ===\n{output_tsv.read_text()}")
            
            if not xml_match:
                print(f"\n=== Expected XML (v1.1.0) ===\n{golden_xml.read_text()[:500]}")
                print(f"\n=== Actual XML (v1.2.0) ===\n{output_xml.read_text()[:500]}")
            
            assert tsv_match, "TSV output is not byte-identical to v1.1.0"
            assert xml_match, "XML output is not byte-identical to v1.1.0"
