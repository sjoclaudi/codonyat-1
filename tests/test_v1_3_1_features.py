"""1.3.1: collapsed-read-counts opt-in, junction codons, exact protein names."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

from aa_caller.app import (
    DEFAULT_ENTROPY_THRESHOLD,
    DEFAULT_RATIO_LOWER,
    DEFAULT_RATIO_UPPER,
    FullReference,
    SamContainer,
    SamEntry,
    validate_reference_file,
)
from aa_caller.models import Protein
from aa_caller.runner import call_variants


def _sam_line(qname: str, seq: str = "ATGCCCTTTAAAGGGCCCAAATTTGGGCCC") -> str:
    return (
        f"{qname}\t0\tref\t1\t30\t{len(seq)}M\t*\t0\t0\t{seq}\t{'I' * len(seq)}\n"
    )


class TestCollapsedReadCounts:
    def test_sample_123_counts_once_by_default(self) -> None:
        entry = SamEntry(_sam_line("sample_123"))
        assert entry.occurences == 1

    def test_sample_123_counts_as_123_with_option(self) -> None:
        entry = SamEntry(_sam_line("sample_123"), collapsed_read_counts=True)
        assert entry.occurences == 123

    def test_pipeline_depth_respects_collapsed_read_counts(self, tmp_path: Path) -> None:
        seq = "ATGCCCTTTAAAGGGCCCAAATTTGGGCCC"
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(f">ref RT(ReverseTranscriptase):1-30\n{seq}\n")
        sam = tmp_path / "reads.sam"
        sam.write_text("@HD\tVN:1.6\n@SQ\tSN:ref\tLN:30\n" + _sam_line("sample_123", seq))
        amps = tmp_path / "amps.tsv"
        amps.write_text("label\tprotein\treference\t5p\t3p\tstart\tend\nAmp_1\tRT\tref\tATG\tCCC\t1\t30\n")

        default = call_variants(
            sam,
            fasta,
            amps,
            csv_path=tmp_path / "default.tsv",
            xml_path=tmp_path / "default.xml",
        )
        collapsed = call_variants(
            sam,
            fasta,
            amps,
            collapsed_read_counts=True,
            csv_path=tmp_path / "collapsed.tsv",
            xml_path=tmp_path / "collapsed.xml",
        )
        assert default.container.variants[1]["depth"] == 1
        assert collapsed.container.variants[1]["depth"] == 123
        default_row = next(
            row
            for row in csv.DictReader((tmp_path / "default.tsv").open(), delimiter="\t")
            if row["POSITION"] == "1"
        )
        collapsed_row = next(
            row
            for row in csv.DictReader((tmp_path / "collapsed.tsv").open(), delimiter="\t")
            if row["POSITION"] == "1"
        )
        assert default_row["TOTALCOV"] == "1"
        assert collapsed_row["TOTALCOV"] == "123"


class TestJunctionCodons:
    def test_codon_crossing_parts_1_10_and_11_30(self, tmp_path: Path) -> None:
        # 1-9 ATGCCCTTT, 10-12 CAT (junction), 13-30 AAAGGGCCCAAATTTGGA
        seq = "ATGCCCTTTCATAAAGGGCCCAAATTTGGA"
        assert len(seq) == 30
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(f">ref PX(Parted):1-10,11-30\n{seq}\n")
        sam = tmp_path / "reads.sam"
        sam.write_text("@HD\tVN:1.6\n@SQ\tSN:ref\tLN:30\n" + _sam_line("read1", seq))
        amps = tmp_path / "amps.tsv"
        amps.write_text("label\tprotein\treference\t5p\t3p\tstart\tend\nAmp_1\tPX\tref\tATG\tGGA\t1\t30\n")

        result = call_variants(
            sam,
            fasta,
            amps,
            protein="PX",
            csv_path=tmp_path / "out.tsv",
            xml_path=tmp_path / "out.xml",
        )
        entries = result.container.variants
        assert 10 in entries
        assert 11 not in entries
        assert sorted(entries) == [1, 4, 7, 10, 13, 16, 19, 22, 25, 28]
        assert entries[10]["codon_coords"] == (10, 11, 12)
        assert "CAT" in entries[10]["variants"]

        ref = FullReference(fasta)
        assert ref.get_seq_at(10, 3, positions=(10, 11, 12)) == "CAT"
        read = SamEntry(_sam_line("read1", seq))
        assert read.codon_at(10, positions=(10, 11, 12)) == "CAT"

    def test_orf1ab_still_has_7097_codons(self, tmp_path: Path) -> None:
        protein = Protein.from_string(
            "ORF1ab(Polyprotein):266-13468,13468-21555", record_id="MN908947.3"
        )
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">MN908947.3 ORF1ab(Polyprotein):266-13468,13468-21555\nA\n")
        dummy = tmp_path / "reads.sam"
        dummy.write_text("@HD\tVN:1.6\nread1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")
        container = SamContainer(
            dummy, DEFAULT_RATIO_UPPER, DEFAULT_RATIO_LOWER, DEFAULT_ENTROPY_THRESHOLD
        )
        ref = FullReference(fasta)
        coords = container._get_codon_positions(protein, ref)
        assert len(coords) == 7097
        # Frameshift: last codon of the first part ends at 13468; next codon starts there.
        starts = [item[0] for item in coords]
        assert 13466 in starts
        assert 13468 in starts
        first_shift = next(item for item in coords if item[0] == 13468)
        assert first_shift[2] == (13468, 13469, 13470)


class TestExactProteinNameMatch:
    def test_protein_s_fails_on_header_without_s(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">ref NS(Nonstructural):1-99\n" + "A" * 99 + "\n")
        with pytest.raises(ValueError, match="missing an S protein annotation"):
            validate_reference_file(fasta, protein_name="S")

    def test_protein_s_still_accepted_when_annotated(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">ref S(Spike):1-99\n" + "A" * 99 + "\n")
        validate_reference_file(fasta, protein_name="S")

    def test_cli_protein_s_fails_on_header_without_s(self, tmp_path: Path) -> None:
        seq = "ATGCCCTTTAAAGGGCCCAAATTTGGGCCC"
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(f">ref NS(Nonstructural):1-30\n{seq}\n")
        sam = tmp_path / "reads.sam"
        sam.write_text("@HD\tVN:1.6\n@SQ\tSN:ref\tLN:30\n" + _sam_line("read1", seq))
        amps = tmp_path / "amps.tsv"
        amps.write_text(
            "label\tprotein\treference\t5p\t3p\tstart\tend\nAmp_1\tNS\tref\tATG\tCCC\t1\t30\n"
        )
        proc = subprocess.run(
            [sys.executable, "-m", "aa_caller", str(sam), str(fasta), str(amps), "--protein", "S"],
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode != 0
        assert "missing an S protein annotation" in proc.stderr
