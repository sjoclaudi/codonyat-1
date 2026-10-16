"""Non-HIV fixtures: SARS-CoV-2 (single record, frameshifted ORF1ab) and
influenza A/PR/8/34 (multi-record HA + NA). Reads are synthetic with known
codon changes, so frequencies are exact. See tests/data/non_hiv/README.md."""

from __future__ import annotations

import csv
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from aa_caller.reference import FullReference
from aa_caller.runner import call_variants

DATA = Path(__file__).parent / "data" / "non_hiv"


def _run(tmp_path: Path, virus: str, protein: str = "all"):
    d = DATA / virus
    res = call_variants(
        d / "reads.sam",
        d / "reference.fasta",
        d / "amplicons.tsv",
        annotation_path=d / "annotation.tsv",
        protein=protein,
        csv_path=tmp_path / f"{virus}.tsv",
        xml_path=tmp_path / f"{virus}.xml",
    )
    with open(res.csv_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return res, rows


def _row(rows, protein, pos, codon):
    hits = [r for r in rows if r["PROTEIN"] == protein and r["POSITION"] == str(pos) and r["VARIANT"] == codon]
    assert len(hits) == 1, (protein, pos, codon, hits)
    return hits[0]


# ---------------------------------------------------------------- SARS-CoV-2


def test_sarscov2_reference_and_annotation():
    ref = FullReference(DATA / "sarscov2" / "reference.fasta", DATA / "sarscov2" / "annotation.tsv")
    assert ref.record_ids == ["MN908947.3"]
    assert len(ref.seq) == 29903
    assert [p.name for p in ref.proteins_in_order()] == ["ORF1ab", "S", "E", "N"]
    orf1ab = ref.get_protein("ORF1ab")
    assert [(i.start, i.end) for i in orf1ab.intervals] == [(266, 13468), (13468, 21555)]
    # Reference codons at the sites mutated in the fixture.
    assert ref.get_seq_at(23402) == "GAT"  # S D614
    assert ref.get_seq_at(28880) + ref.get_seq_at(28883) == "AGGGGA"  # N R203, G204


def test_sarscov2_all_proteins(tmp_path):
    res, rows = _run(tmp_path, "sarscov2")
    assert res.proteins == ["ORF1ab", "S", "E", "N"]
    # Only covered codons reach the TSV: 29 S (read starts mid-codon) + 30 N,
    # plus one extra row per mutated codon (1 in S, 2 in N).
    assert len(rows) == 62
    assert {r["PROTEIN"] for r in rows} == {"S", "N"}
    assert "RECORD" not in rows[0]  # single-record layout unchanged
    d614g = _row(rows, "S", 23402, "GGT")
    assert float(d614g["FREQ"]) == pytest.approx(70.0)
    assert (d614g["FWCOV"], d614g["RVCOV"], d614g["TOTALCOV"]) == ("10", "10", "20")
    assert float(_row(rows, "S", 23402, "GAT")["FREQ"]) == pytest.approx(30.0)
    assert float(_row(rows, "N", 28880, "AAA")["FREQ"]) == pytest.approx(40.0)  # R203K
    assert float(_row(rows, "N", 28883, "CGA")["FREQ"]) == pytest.approx(40.0)  # G204R


def test_sarscov2_xml_codon_counts(tmp_path):
    res, _ = _run(tmp_path, "sarscov2")
    root = ET.parse(res.xml_path).getroot()
    counts = {p.get("name"): len(p.findall("Position")) for p in root.findall("Protein")}
    # ORF1ab across the -1 frameshift: 4401 + 2696 codons (stop included).
    assert counts == {"ORF1ab": 7097, "S": 1274, "E": 76, "N": 420}


def test_sarscov2_single_protein(tmp_path):
    _, rows = _run(tmp_path, "sarscov2", protein="N")
    assert {r["PROTEIN"] for r in rows} == {"N"}
    assert len(rows) == 32


# ------------------------------------------------------------- influenza A


def test_influenza_multi_record_reference():
    ref = FullReference(DATA / "influenza" / "reference.fasta", DATA / "influenza" / "annotation.tsv")
    assert ref.record_ids == ["NC_002017.1", "NC_002018.1"]
    assert [len(ref.sequences[r]) for r in ref.record_ids] == [1778, 1413]
    assert ref.get_protein("HA").record_id == "NC_002017.1"
    assert ref.get_protein("NA").record_id == "NC_002018.1"
    assert ref.get_seq_at(45, record_id="NC_002017.1") == "CTA"
    assert ref.get_seq_at(27, record_id="NC_002018.1") == "CCA"


def test_influenza_reads_stay_on_their_segment(tmp_path):
    res, rows = _run(tmp_path, "influenza")
    assert res.proteins == ["HA", "NA"]
    assert len(rows) == 62  # 30 HA + 30 NA codons, +1 per mutated codon
    assert all(r["RECORD"] == r["REFERENCE"] for r in rows)
    ha = [r for r in rows if r["PROTEIN"] == "HA"]
    na = [r for r in rows if r["PROTEIN"] == "NA"]
    assert {r["RECORD"] for r in ha} == {"NC_002017.1"}
    assert {r["RECORD"] for r in na} == {"NC_002018.1"}
    # Depths prove NA reads are not counted on HA and vice versa.
    assert {r["TOTALCOV"] for r in ha} == {"16"}
    assert {r["TOTALCOV"] for r in na} == {"10"}
    assert float(_row(rows, "HA", 45, "GTA")["FREQ"]) == pytest.approx(25.0)  # L5V
    assert float(_row(rows, "NA", 27, "GCA")["FREQ"]) == pytest.approx(60.0)  # P3A


def test_influenza_cli(tmp_path):
    d = DATA / "influenza"
    sam = tmp_path / "reads.sam"
    sam.write_text((d / "reads.sam").read_text())
    subprocess.run(
        [sys.executable, "-m", "aa_caller", str(sam), str(d / "reference.fasta"), str(d / "amplicons.tsv"),
         "--annotation", str(d / "annotation.tsv"), "--protein", "NA"],
        check=True, capture_output=True,
    )
    with open(tmp_path / "reads.sam.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 31 and {r["PROTEIN"] for r in rows} == {"NA"}
