"""Regenerate the tiny non-HIV fixtures (deterministic, no randomness).

Inputs: the NCBI FASTA files fetched with
  efetch -db nuccore -id MN908947.3,NC_002017.1,NC_002018.1 -format fasta
saved under --src. Reads are synthetic (simulated from the reference with
known codon changes), so the expected frequencies are exact.
Run: python make_fixtures.py --src /path/to/fastas
"""
from __future__ import annotations

import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
READ_LEN = 90


def read_fasta(path: Path) -> tuple[str, str]:
    lines = path.read_text().splitlines()
    return lines[0][1:].split()[0], "".join(lines[1:]).upper()


def mutate(seq: str, changes: dict[int, str]) -> str:
    s = list(seq)
    for pos, base in changes.items():
        s[pos - 1] = base
    return "".join(s)


def reads_for(rname, ref, amp, start, n_per_strand, n_mut_per_strand, changes):
    """n reads per strand covering start..start+READ_LEN-1; first n_mut carry changes."""
    mut = mutate(ref, changes)
    out = []
    for strand, flag in (("F", 0), ("R", 16)):
        for i in range(n_per_strand):
            src = mut if i < n_mut_per_strand else ref
            seq = src[start - 1 : start - 1 + READ_LEN]
            out.append(
                f"{amp}_{strand}{i}\t{flag}\t{rname}\t{start}\t60\t{READ_LEN}M\t*\t0\t0\t{seq}\t{'I' * READ_LEN}"
            )
    return out


def write_sam(path, records, reads):
    hdr = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{r}\tLN:{len(s)}" for r, s in records]
    path.write_text("\n".join(hdr + reads) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True)
    src = ap.parse_args().src

    # --- SARS-CoV-2 (MN908947.3), single record, TSV annotation -------------
    rid, seq = read_fasta(src / "MN908947.3.fa")
    d = HERE / "sarscov2"
    d.mkdir(exist_ok=True)
    (d / "reference.fasta").write_text(f">{rid}\n" + "\n".join(seq[i:i+70] for i in range(0, len(seq), 70)) + "\n")
    (d / "annotation.tsv").write_text(
        "record\tprotein\tdescription\tintervals\tstrand\n"
        f"{rid}\tORF1ab\tORF1ab polyprotein (-1 frameshift)\t266-13468;13468-21555\t+\n"
        f"{rid}\tS\tsurface glycoprotein\t21563-25384\t+\n"
        f"{rid}\tE\tenvelope protein\t26245-26472\t+\n"
        f"{rid}\tN\tnucleocapsid phosphoprotein\t28274-29533\t+\n"
    )
    (d / "amplicons.tsv").write_text(
        "label\tprotein\treference\t5p\t3p\tstart\tend\n"
        f"Amp_1\tS\t{rid}\t{seq[23360:23380]}\t{seq[23430:23450]}\t23361\t23450\n"
        f"Amp_2\tN\t{rid}\t{seq[28840:28860]}\t{seq[28910:28930]}\t28841\t28930\n"
    )
    # S D614G: codon 23402-23404 GAT->GGT (A23403G), 7/10 reads per strand.
    # N R203K/G204R: GGG->AAC at 28881-28883 (codons 28880 AGG->AAA, 28883 GGA->CGA), 4/10.
    reads = reads_for(rid, seq, "Amp_1", 23361, 10, 7, {23403: "G"})
    reads += reads_for(rid, seq, "Amp_2", 28841, 10, 4, {28881: "A", 28882: "A", 28883: "C"})
    write_sam(d / "reads.sam", [(rid, seq)], reads)

    # --- Influenza A/PR/8/34, multi-record (HA seg 4 + NA seg 6) -------------
    ha_id, ha = read_fasta(src / "NC_002017.1.fa")
    na_id, na = read_fasta(src / "NC_002018.1.fa")
    d = HERE / "influenza"
    d.mkdir(exist_ok=True)
    (d / "reference.fasta").write_text(f">{ha_id}\n{ha}\n>{na_id}\n{na}\n")
    (d / "annotation.tsv").write_text(
        "record\tprotein\tdescription\tintervals\tstrand\n"
        f"{ha_id}\tHA\thaemagglutinin\t33-1733\t+\n"
        f"{na_id}\tNA\tneuraminidase\t21-1385\t+\n"
    )
    (d / "amplicons.tsv").write_text(
        "label\tprotein\treference\t5p\t3p\tstart\tend\n"
        f"Amp_1\tHA\t{ha_id}\t{ha[32:52]}\t{ha[102:122]}\t33\t122\n"
        f"Amp_2\tNA\t{na_id}\t{na[20:40]}\t{na[90:110]}\t21\t110\n"
    )
    # HA codon 5 (pos 45-47) third-position change; NA codon 3 (pos 27-29) first base.
    ha_codon5 = ha[44:47]
    na_codon3 = na[26:29]
    ha_change = {45: "C" if ha_codon5[0] != "C" else "G"}
    na_change = {27: "C" if na_codon3[0] != "C" else "G"}
    reads = reads_for(ha_id, ha, "Amp_1", 33, 8, 2, ha_change)
    reads += reads_for(na_id, na, "Amp_2", 21, 5, 3, na_change)
    write_sam(d / "reads.sam", [(ha_id, ha), (na_id, na)], reads)


if __name__ == "__main__":
    main()
