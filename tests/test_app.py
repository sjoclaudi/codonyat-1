from pathlib import Path

import pytest

from aa_caller import app


def _write_fasta(path: Path, header: str) -> None:
    path.write_text(f">{header}\nATGCCC\n")


def _write_amplicon_csv(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines))


def test_codon_to_aminoacid_translates_known_codon() -> None:
    assert app.codon_to_aminoacid("ATG") == "M"
    assert app.codon_to_aminoacid("tca") == "S"


def test_codon_to_aminoacid_unknown_returns_x() -> None:
    assert app.codon_to_aminoacid("ZZZ") == "X"


def test_protein_from_string_parses_interval() -> None:
    text = "RT(Reverse Transcriptase):1-99"
    protein = app.Protein.from_string(text)
    assert protein.name == "RT"
    assert protein.description == "Reverse Transcriptase"
    assert protein.start_coordinate == 1
    assert protein.end_coordinate == 99


def test_amplicon_from_string_and_parse(tmp_path: Path) -> None:
    line = "Amp_1\tRT\tref\tACT\tGCA\t10\t20"
    amplicon = app.Amplicon.from_string(line)
    assert amplicon.label == "Amp_1"
    assert amplicon.protein == "RT"
    assert amplicon.start_coordinate == 10
    assert amplicon.end_coordinate == 20

    csv_file = tmp_path / "amps.tsv"
    _write_amplicon_csv(csv_file, ["label\tprotein\treference\t5p\t3p\tstart\tend", line])
    parsed = app.parse_amplicons(csv_file)
    assert "Amp_1" in parsed


def test_amplicon_from_string_requires_columns() -> None:
    with pytest.raises(ValueError, match="requires at least 7 columns"):
        app.Amplicon.from_string("Amp_1\tRT\tref")


def test_validate_reference_file_accepts_rt(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    _write_fasta(fasta, "ref RT(Reverse Transcriptase):1-99")
    app.validate_reference_file(fasta)


def test_validate_reference_file_rejects_missing_rt(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    _write_fasta(fasta, "ref POL(Polymerase):1-99")
    with pytest.raises(ValueError, match="missing a RT protein annotation"):
        app.validate_reference_file(fasta)


def test_validate_amplicon_file_ok(tmp_path: Path) -> None:
    amplicon_file = tmp_path / "amps.tsv"
    header = "label\tprotein\treference\t5p\t3p\tstart\tend"
    line = "Amp_1\tRT\tref\tACT\tGCA\t100\t200"
    _write_amplicon_csv(amplicon_file, [header, line])
    app.validate_amplicon_file(amplicon_file)


def test_validate_amplicon_file_header_too_short(tmp_path: Path) -> None:
    amplicon_file = tmp_path / "amps.tsv"
    _write_amplicon_csv(amplicon_file, ["label\tprotein	reference"])
    with pytest.raises(ValueError, match="does not contain the expected header"):
        app.validate_amplicon_file(amplicon_file)


def test_validate_amplicon_file_no_rows(tmp_path: Path) -> None:
    amplicon_file = tmp_path / "amps.tsv"
    header = "label\tprotein\treference\t5p\t3p\tstart\tend"
    _write_amplicon_csv(amplicon_file, [header])
    with pytest.raises(ValueError, match="contains no definitions"):
        app.validate_amplicon_file(amplicon_file)


def test_validate_sam_file_accepts_minimal(tmp_path: Path) -> None:
    sam_file = tmp_path / "reads.sam"
    sam_file.write_text("@HD\nread1\t0\tref\t5\t30\t3M\t*\t0\t0\tACT\t!!!\n")
    app.validate_sam_file(sam_file)


def test_validate_sam_file_rejects_bad_row(tmp_path: Path) -> None:
    sam_file = tmp_path / "reads.sam"
    sam_file.write_text("\t\t\t\n")
    with pytest.raises(ValueError, match="looks malformed"):
        app.validate_sam_file(sam_file)
