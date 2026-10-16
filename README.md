# codonyat

codon-yat — A codon-aware amino acid variant typer from SAM alignments of viral NGS data.

A Python package for codon-level amino acid variant calling from SAM alignments of viral amplicon sequencing, including:

- SAM parsing that honors amplicon labels and strand orientation.
- CIGAR-aware codon extraction that correctly handles insertions, deletions, and soft clips.
- Ratio balancing, entropy tracking, and TSV/XML exporters mirroring the legacy Perl outputs.
- Configurable CLI flags for strand ratio bounds, entropy sensitivity, and target protein.
- Opt-in research statistics (Wilson CI, strand-bias test, diversity indices) and AAVF v1.0. Research use only.

## Highlights

- Produces consistent TSV/XML diagnostics while exposing Python objects (`SamContainer`, `FullReference`, etc.) that downstream tooling can import.
- Validates every input file before parsing to surface malformed SAM/FASTA/amplicon data early.
- Logs per-position entropy and strand balance for easier debugging in CI or local runs.
- Supports any protein target via `--protein` (one name, a comma-separated list, or `all`; defaults to RT).

## Installation

```bash
pip install codonyat
```

Requires Python 3.10 or newer. The only runtime dependency is Biopython.

## Running it in a workflow

This repository contains only the Python package. To run codonyat on many samples
(from FASTQ or SAM, with read QC, trimming, alignment and a merged report), use the
separate Nextflow pipeline **[nf-codonyat](https://github.com/sjoclaudi/nf-codonyat)**,
which installs this package from PyPI.

## CLI usage

Once installed, the `codonyat` entry point is available:

```bash
codonyat sample.sam reference.fasta amplicons.tsv
```

For multi-record references or complex annotations:

```bash
codonyat sample.sam reference.fasta amplicons.tsv --annotation annotations.tsv --protein S
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--annotation` | None | Optional TSV annotation file (record, protein, description, intervals, strand) |
| `--protein` | `RT` | Protein(s) to analyse: one name, a comma-separated list (e.g. `PR,RT,INT`, counted in a single pass) or `all` (every protein annotated in the reference header or TSV). Names must match the FASTA header annotations |
| `--ratio-upper` | `3.162` | Upper bound for strand-ratio balancing |
| `--ratio-lower` | `0.316` | Lower bound for strand-ratio balancing |
| `--entropy-threshold` | `0.0` | Minimum Shannon entropy to mark a position as diverse |
| `--min-baseq` | off | Count a codon only if all 3 bases have Phred-like quality ≥ this value |
| `--stats` | off | Optional stats TSV (CI, strand-bias, flags, diversity) |
| `--aavf` | off | Optional AAVF v1.0 path (no CAT/SRVL) |
| `--min-freq` | `0.01` | Frequency floor (fraction) for `lowfreq` / AAVF FILTER |
| `--min-coverage` | `10` | Coverage floor for `lowcov` / AAVF FILTER |
| `--sb-alpha` | `0.05` | Fisher p-value cutoff for the `sb` flag |

### Output

The CLI writes two files alongside the input SAM:

- `[sam-file].tsv` — TSV with columns: FILE, REFERENCE, [RECORD], PROTEIN, VARIANT, POSITION, FREQ, FWCOV, RVCOV, TOTALCOV, RATIO
  - The RECORD column is added only for multi-record references (v1.2.0+)
  - Single-record outputs remain byte-identical to v1.1.0
  - Rows are ordered by genomic position (codon index) along the protein's intervals
  - For multi-interval proteins, this follows the biological reading order through all segments
- `[sam-file].xml` — per-position `<Depth>`, `<FwCover>`, `<RvCover>`, `<Variants>`
  - With several proteins, each protein's positions are wrapped in `<Protein name="…">`; single-protein XML is unchanged

Default TSV/XML are byte-identical to 1.2.1 when the opt-in flags are omitted.

### Research statistics (opt-in)

`--stats PATH` writes a sidecar TSV with `LEVEL` = `variant`, `codon`, or `gene`:

- **variant**: frequency as a fraction, Wilson 95% CI, Fisher strand-bias p, flags (`lowfreq`, `sb`, `lowcov`)
- **codon / gene**: Shannon entropy (raw and normalised), Gini–Simpson, nucleotide π, amino-acid π

`--aavf PATH` writes AAVF v1.0. FILTER is taken from the flags above. The file
does not contain resistance or surveillance fields. Formulas and citations:
[docs/methods.md](docs/methods.md).

```bash
codonyat-runner reads.sam reference.fasta amps.tsv \
  --protein all --stats sample.stats.tsv --aavf sample.aavf
```

### Reference FASTA format

The reference FASTA supports both single-record and multi-record files. Protein annotations can be provided either in the FASTA header or via a separate TSV file.

#### FASTA header annotation (backward compatible)

Protein annotations in the format `Name(Description):start-end`, separated by semicolons:

```
>K03455|HIVHXB2CG PR(Protease):2253-2549;RT(Reverse Transcriptase):2550-3869;INT(Integrase):4230-5093
TGGAAGGGCTAATTCACTCCCAACGAAGAC...
```

#### Multi-interval proteins

For proteins spanning multiple intervals (frameshift, splicing, or circular origin-crossing), use comma-separated intervals. All coordinates are **1-based inclusive**:

```
>MN908947.3 ORF1ab(Polyprotein):266-13468,13468-21555
ATTAAAGGTTTATACCTTCCCAGGTAACAA...
```

The interval `266-13468,13468-21555` represents a -1 programmed ribosomal frameshift where base 13468 is read twice. Semicolons separate different proteins in the header.

#### TSV annotation format (v1.2.0+)

For complex annotations, use a TSV file with the `--annotation` flag:

```
record	protein	description	intervals	strand
MN908947.3	S	Spike glycoprotein	21563-25384	+
MN908947.3	ORF1ab	Polyprotein	266-13468;13468-21555	+
NC_003977.2	P	DNA polymerase	2307-3221;1-1623	+
```

**Important**: Minus-strand genes are not supported and will be rejected with a clear error message.

## Package API

Import `aa_caller` to reuse the core objects:

```python
from aa_caller import SamContainer, FullReference, parse_amplicons
```

### Python wrapper

Use `call_variants` to run variant calling from Python without touching the CLI:

```python
from pathlib import Path
from aa_caller import call_variants

result = call_variants(
    sam_path=Path("reads.sam"),
    reference_path=Path("reference.fasta"),
    amplicons_path=Path("amps.tsv"),
    protein="PR,RT,INT",  # or "all", "RT", or ["PR", "RT"]
    ratio_upper=3.2,
)

print(result.csv_path, result.xml_path)
print(result.container.variants.keys())
```

If you already have an `argparse.Namespace` or mapping of the CLI arguments, `call_variants_from_args` adapts them directly:

```python
from argparse import Namespace

args = Namespace(
    sam_file="reads.sam",
    reference_file="reference.fasta",
    amplicons_file="amps.tsv",
    protein="RT",
    ratio_upper=3.2,
)

call_variants_from_args(args)
```

### CLI wrapper

A lightweight runner at `codonyat-runner` exposes the same arguments as `call_variants_from_args`:

```bash
codonyat-runner reads.sam reference.fasta amps.tsv --protein PR,RT,INT --ratio-upper 3.1 --csv-path results.tsv
```

## Development

Install the repository with the optional dev tooling so your local environment matches CI:

```bash
pip install --upgrade pip
pip install -e .[dev]
```

Now you can run the same checks as [.github/workflows/python-tests.yml](.github/workflows/python-tests.yml):

```bash
ruff check .
python -m pytest
python -m build && twine check dist/*
```

## Project structure

```
codonyat/
├── aa_caller/            # Python package
│   ├── cli.py            # CLI entry point (codonyat)
│   ├── runner.py         # Python API (call_variants, call_variants_from_args) and codonyat-runner
│   ├── container.py      # SamContainer — variant aggregation engine
│   ├── sam.py            # SamEntry — CIGAR-aware SAM record parser
│   ├── reference.py      # FullReference — FASTA + protein annotation loader
│   ├── models.py         # Protein, Amplicon, Variant, Qual dataclasses
│   ├── validators.py     # Input file validation
│   ├── genetic_code.py   # Codon translation table
│   ├── constants.py      # Default thresholds
│   ├── metrics.py        # Wilson CI, Fisher, diversity indices
│   ├── stats_writer.py   # --stats TSV
│   └── aavf.py           # AAVF v1.0 writer
├── tests/                # pytest suite (unit + integration)
├── docs/index.html       # Project page (GitHub Pages)
├── docs/methods.md       # Formulas and citations for 1.3.0 outputs
└── pyproject.toml        # Package metadata
```

## Releasing

Publishing to PyPI is done by `.github/workflows/publish-pypi.yml` when a GitHub
release is published:

1. Set `version` in `pyproject.toml` (e.g. `1.0.3`) and merge it to `main`.
2. Create a tag `v<version>` and a GitHub release from it.
3. The workflow builds the sdist and wheel, runs `twine check`, checks that the tag
   matches the package version, and publishes via PyPI trusted publishing
   (GitHub environment `pypi`).

## License

MIT
