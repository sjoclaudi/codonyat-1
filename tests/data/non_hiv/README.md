# Non-HIV test fixtures

Small fixtures that exercise codonyat outside HIV-1: TSV annotation, a
multi-interval (frameshifted) protein and a multi-record reference.
Research use only.

| Dir | Reference | Records | Proteins (annotation.tsv) |
|---|---|---|---|
| `sarscov2/` | MN908947.3 (SARS-CoV-2 Wuhan-Hu-1), 29,903 nt | 1 | ORF1ab `266-13468;13468-21555`, S `21563-25384`, E `26245-26472`, N `28274-29533` |
| `influenza/` | NC_002017.1 (A/PR/8/34 segment 4, HA, 1,778 nt) + NC_002018.1 (segment 6, NA, 1,413 nt) | 2 | HA `33-1733` on NC_002017.1, NA `21-1385` on NC_002018.1 |

Sequences: NCBI GenBank/RefSeq (INSDC). NCBI places no restrictions on use or
distribution of GenBank data; cite the accessions. Coordinates come from the
NCBI feature tables of the same accessions.

Reads (`reads.sam`) are synthetic, written by `make_fixtures.py` (no
randomness): 90 nt perfect-match reads, equal numbers on each strand, with a
known fraction carrying a codon change.

| Fixture | Amplicon / region | Reads (F+R) | Change | Expected FREQ |
|---|---|---|---|---|
| sarscov2 | Amp_1, S 23361-23450 | 10+10 | S D614G, A23403G (codon 23402 GAT>GGT) | 70% |
| sarscov2 | Amp_2, N 28841-28930 | 10+10 | N R203K/G204R, GGG>AAC 28881-28883 | 40% (codons 28880, 28883) |
| influenza | Amp_1, HA 33-122 | 8+8 | HA codon 5 (pos 45) CTA>GTA (L5V) | 25% |
| influenza | Amp_2, NA 21-110 | 5+5 | NA codon 3 (pos 27) CCA>GCA (P3A) | 60% |

Regenerate: fetch the three FASTA files from NCBI efetch
(`db=nuccore&rettype=fasta`) into a folder and run
`python make_fixtures.py --src <folder>`.

The SAM files are plain SAM because codonyat reads SAM only. To make a BAM for
pipeline tests: `samtools sort -o reads.bam reads.sam && samtools index reads.bam`.
