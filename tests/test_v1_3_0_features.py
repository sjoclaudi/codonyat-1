"""Opt-in 1.3.0 stats / AAVF / min-baseq tests. Default TSV/XML stay 1.2.1."""

from __future__ import annotations

import csv
import filecmp
from pathlib import Path

import pytest

from aa_caller.app import (
    DEFAULT_ENTROPY_THRESHOLD,
    DEFAULT_RATIO_LOWER,
    DEFAULT_RATIO_UPPER,
    FullReference,
    SamContainer,
)
from aa_caller.constants import PACKAGE_VERSION
from aa_caller.runner import call_variants

DATA = Path(__file__).parent / "data"


class TestDefaultOutputUnchanged:
    def test_hiv1_rt_tsv_xml_byte_identical_to_1_2_1(self, tmp_path: Path) -> None:
        result = call_variants(
            sam_path=DATA / "hiv1_reads.sam",
            reference_path=DATA / "hiv1_ref.fasta",
            amplicons_path=DATA / "hiv1_amps.tsv",
            protein="RT",
            csv_path=tmp_path / "out.tsv",
            xml_path=tmp_path / "out.xml",
        )
        # FILE column uses the path string; rewrite with the golden sample name.
        ref = FullReference(DATA / "hiv1_ref.fasta")
        container = SamContainer(
            DATA / "hiv1_reads.sam",
            DEFAULT_RATIO_UPPER,
            DEFAULT_RATIO_LOWER,
            DEFAULT_ENTROPY_THRESHOLD,
        )
        container.load()
        container.calculate_variant_frequencies("RT", ref)
        container.write_csv("hiv1_reads.sam", ref, tmp_path / "named.tsv", protein_name="RT")
        container.write_xml("hiv1_reads.sam", ref, tmp_path / "named.xml")
        assert filecmp.cmp(tmp_path / "named.tsv", DATA / "hiv1_golden_v1.1.0.tsv", shallow=False)
        assert filecmp.cmp(tmp_path / "named.xml", DATA / "hiv1_golden_v1.1.0.xml", shallow=False)
        assert result.stats_path is None
        assert result.aavf_path is None

    def test_writing_stats_does_not_change_tsv_xml(self, tmp_path: Path) -> None:
        plain = call_variants(
            sam_path=DATA / "hiv1_reads.sam",
            reference_path=DATA / "hiv1_ref.fasta",
            amplicons_path=DATA / "hiv1_amps.tsv",
            csv_path=tmp_path / "plain.tsv",
            xml_path=tmp_path / "plain.xml",
        )
        extra = call_variants(
            sam_path=DATA / "hiv1_reads.sam",
            reference_path=DATA / "hiv1_ref.fasta",
            amplicons_path=DATA / "hiv1_amps.tsv",
            csv_path=tmp_path / "extra.tsv",
            xml_path=tmp_path / "extra.xml",
            stats_path=tmp_path / "extra.stats.tsv",
            aavf_path=tmp_path / "extra.aavf",
            aavf_file_date="20261016",
        )
        assert plain.csv_path.read_bytes() == extra.csv_path.read_bytes()
        assert plain.xml_path.read_bytes() == extra.xml_path.read_bytes()
        assert extra.stats_path.exists()
        assert extra.aavf_path.exists()

    def test_hiv1_all_proteins_byte_identical(self, tmp_path: Path, monkeypatch) -> None:
        for name in ("hiv1_reads.sam", "hiv1_pol_ref.fasta", "hiv1_amps.tsv"):
            (tmp_path / name).write_bytes((DATA / name).read_bytes())
        monkeypatch.chdir(tmp_path)
        result = call_variants(
            sam_path="hiv1_reads.sam",
            reference_path="hiv1_pol_ref.fasta",
            amplicons_path="hiv1_amps.tsv",
            protein="all",
        )
        assert filecmp.cmp(result.csv_path, DATA / "hiv1_pol_all_v1.1.0.tsv", shallow=False)
        assert filecmp.cmp(result.xml_path, DATA / "hiv1_pol_all_v1.1.0.xml", shallow=False)


class TestMinBaseq:
    def test_low_quality_codon_dropped(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">testref RT(ReverseTranscriptase):1-30\nATGCCCTTTAAAGGGCCCAAATTTGGGCCC\n")
        high = "I" * 30
        low = "!" * 30
        sam = tmp_path / "reads.sam"
        sam.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:testref\tLN:30\n"
            f"Amp_1_good\t0\ttestref\t1\t30\t30M\t*\t0\t0\tATGCCCTTTAAAGGGCCCAAATTTGGGCCC\t{high}\n"
            f"Amp_1_bad\t0\ttestref\t1\t30\t30M\t*\t0\t0\tTTGCCCTTTAAAGGGCCCAAATTTGGGCCC\t{low}\n"
        )
        amps = tmp_path / "amps.tsv"
        amps.write_text("label\tprotein\treference\t5p\t3p\tstart\tend\nAmp_1\tRT\ttestref\tATG\tCCC\t1\t30\n")

        unfiltered = call_variants(
            sam_path=sam,
            reference_path=fasta,
            amplicons_path=amps,
            csv_path=tmp_path / "plain.tsv",
            xml_path=tmp_path / "plain.xml",
        )
        filtered = call_variants(
            sam_path=sam,
            reference_path=fasta,
            amplicons_path=amps,
            min_baseq=20,
            csv_path=tmp_path / "q20.tsv",
            xml_path=tmp_path / "q20.xml",
        )
        plain_codons = {row["VARIANT"] for row in csv.DictReader(unfiltered.csv_path.open(), delimiter="\t") if row["POSITION"] == "1"}
        q20_codons = {row["VARIANT"] for row in csv.DictReader(filtered.csv_path.open(), delimiter="\t") if row["POSITION"] == "1"}
        assert "ATG" in plain_codons and "TTG" in plain_codons
        assert q20_codons == {"ATG"}


class TestStatsAndAavf:
    def test_stats_has_variant_codon_gene_and_reference_values(self, tmp_path: Path) -> None:
        result = call_variants(
            sam_path=DATA / "hiv1_reads.sam",
            reference_path=DATA / "hiv1_ref.fasta",
            amplicons_path=DATA / "hiv1_amps.tsv",
            protein="RT",
            csv_path=tmp_path / "out.tsv",
            xml_path=tmp_path / "out.xml",
            stats_path=tmp_path / "out.stats.tsv",
        )
        rows = list(csv.DictReader(result.stats_path.open(), delimiter="\t"))
        levels = {row["LEVEL"] for row in rows}
        assert levels == {"variant", "codon", "gene"}
        variants = [row for row in rows if row["LEVEL"] == "variant"]
        assert variants
        # Two balanced reads, one codon: freq=1, Wilson of 2/2 is not (1,1)
        first = variants[0]
        assert first["FREQ"] == "1"
        assert float(first["FREQ_CI_LOW"]) < 1.0
        assert float(first["FREQ_CI_HIGH"]) == 1.0
        assert first["FLAGS"] == "lowcov"  # coverage 2 < default 10
        assert first["AA_POS"] == "1"
        gene = next(row for row in rows if row["LEVEL"] == "gene")
        assert gene["PROTEIN"] == "RT"
        assert gene["N_CODONS"] == "10"
        assert gene["SHANNON"] == "0"
        assert gene["PI_NT"] == "0"

    def test_aavf_header_and_no_resistance_fields(self, tmp_path: Path) -> None:
        result = call_variants(
            sam_path=DATA / "hiv1_reads.sam",
            reference_path=DATA / "hiv1_ref.fasta",
            amplicons_path=DATA / "hiv1_amps.tsv",
            protein="RT",
            csv_path=tmp_path / "out.tsv",
            xml_path=tmp_path / "out.xml",
            aavf_path=tmp_path / "out.aavf",
            aavf_file_date="20261016",
        )
        text = result.aavf_path.read_text()
        assert text.startswith("##fileformat=AAVFv1.0\n")
        assert f"##source=codonyat-{PACKAGE_VERSION} (research use only)" in text
        assert "CAT" not in text and "SRVL" not in text
        assert "resistance" not in text.lower()
        lines = text.splitlines()
        header = next(line for line in lines if line.startswith("#CHROM"))
        assert header.split("\t") == [
            "#CHROM",
            "GENE",
            "POS",
            "REF",
            "ALT",
            "FILTER",
            "ALT_FREQ",
            "COVERAGE",
            "INFO",
        ]
        data = [line for line in lines if not line.startswith("#")]
        assert data
        chrom, gene, pos, ref, alt, filt, freq, cov, info = data[0].split("\t")
        assert chrom == "K03455"
        assert gene == "RT"
        assert pos == "1"
        assert ref == "M" and alt == "M"
        assert "lowcov" in filt
        assert "RC=atg" in info
        assert "FWC=" in info and "RVC=" in info and "SR=" in info

    def test_aavf_groups_codons_by_amino_acid(self, tmp_path: Path) -> None:
        fasta = tmp_path / "ref.fasta"
        fasta.write_text(">testref RT(ReverseTranscriptase):1-30\nATGCCCTTTAAAGGGCCCAAATTTGGGCCC\n")
        seq_ref = "ATGCCCTTTAAAGGGCCCAAATTTGGGCCC"
        seq_var = "TTG" + seq_ref[3:]  # M→L
        q = "I" * 30
        sam = tmp_path / "reads.sam"
        sam.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:testref\tLN:30\n"
            f"Amp_1_r1\t0\ttestref\t1\t30\t30M\t*\t0\t0\t{seq_ref}\t{q}\n"
            f"Amp_1_r2\t16\ttestref\t1\t30\t30M\t*\t0\t0\t{seq_ref}\t{q}\n"
            f"Amp_1_r3\t0\ttestref\t1\t30\t30M\t*\t0\t0\t{seq_var}\t{q}\n"
            f"Amp_1_r4\t16\ttestref\t1\t30\t30M\t*\t0\t0\t{seq_var}\t{q}\n"
        )
        amps = tmp_path / "amps.tsv"
        amps.write_text("label\tprotein\treference\t5p\t3p\tstart\tend\nAmp_1\tRT\ttestref\tATG\tCCC\t1\t30\n")
        result = call_variants(
            sam_path=sam,
            reference_path=fasta,
            amplicons_path=amps,
            csv_path=tmp_path / "out.tsv",
            xml_path=tmp_path / "out.xml",
            aavf_path=tmp_path / "out.aavf",
            stats_path=tmp_path / "out.stats.tsv",
            aavf_file_date="20261016",
            min_coverage=1,
            min_freq=0.01,
        )
        rows = [line.split("\t") for line in result.aavf_path.read_text().splitlines() if not line.startswith("#")]
        pos1 = [row for row in rows if row[2] == "1"]
        alts = {row[4]: row for row in pos1}
        assert set(alts) == {"M", "L"}
        assert alts["M"][3] == "M"
        assert alts["L"][3] == "M"
        assert float(alts["M"][6]) == pytest.approx(0.5)
        assert float(alts["L"][6]) == pytest.approx(0.5)
        assert "AC=Ttg" in alts["L"][8] or "AC=ttg" in alts["L"][8].lower()

        stats = list(csv.DictReader(result.stats_path.open(), delimiter="\t"))
        codon1 = next(row for row in stats if row["LEVEL"] == "codon" and row["POSITION"] == "1")
        # 2 ATG + 2 TTG → Shannon of two equal types = ln(2)
        import math

        assert float(codon1["SHANNON"]) == pytest.approx(math.log(2), rel=1e-6)
        assert float(codon1["SHANNON_NORM"]) == pytest.approx(1.0, rel=1e-6)
        assert float(codon1["GINI_SIMPSON"]) == pytest.approx(0.5, rel=1e-6)
        assert float(codon1["PI_NT"]) == pytest.approx((4 / 3 * 0.5) / 3.0, rel=1e-6)
        assert float(codon1["PI_AA"]) == pytest.approx(4 / 3 * 0.5, rel=1e-6)
