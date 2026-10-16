# Methods (codonyat 1.3.0)

Research use only. These formulas apply to the opt-in `--stats` and `--aavf` outputs.
Default TSV/XML writing is unchanged from 1.2.1.

## Error model and confidence

### Minimum base quality

When `--min-baseq Q` is set, a codon is counted only if all three aligned bases
have a Phred-like quality ≥ Q (Sanger table already used by `Qual`). Gapped or
soft-clipped codon positions fail the check. The flag is off by default so
existing TSV/XML counts stay byte-identical to 1.2.1.

Quality-aware counting of minority variants follows the same motivation as
LoFreq (Wilm et al. 2012, *Nucleic Acids Research*,
https://doi.org/10.1093/nar/gks918) and iVar (Grubaugh et al. 2019,
*Genome Biology*, https://doi.org/10.1186/s13059-018-1618-7). See also
McCrone & Lauring 2016, *Journal of Virology*,
https://doi.org/10.1128/JVI.00667-16.

### Wilson 95% confidence interval

Each variant frequency `p = k / n` (k = codon count, n = codon coverage) is
reported with a Wilson score interval at 95% (`z = Φ⁻¹(0.975)`):

Wilson, E. B. (1927). Probable inference, the law of succession, and statistical
inference. *Journal of the American Statistical Association* 22:209–212.
https://doi.org/10.1080/01621459.1927.10502953

Newcombe, R. G. (1998). Two-sided confidence intervals for the single
proportion. *Statistics in Medicine* 17:857–872.

### Strand-bias test

A two-sided Fisher's exact test is applied to the 2×2 table of
variant vs other × forward vs reverse counts (Fisher 1935). The `sb` flag is
set when `p < --sb-alpha` (default 0.05).

### Flags

| Code | Meaning |
|------|---------|
| `lowfreq` | variant (or AAVF ALT) frequency `< --min-freq` (default 0.01) |
| `lowcov` | codon coverage `< --min-coverage` (default 10) |
| `sb` | Fisher p-value `< --sb-alpha` (default 0.05) |

Flags are written only to `--stats` and to AAVF `FILTER` (semicolon-separated,
or `PASS`). They never alter the default TSV/XML columns.

## Diversity indices

Computed on unambiguous ACGT codons only. Codons containing `N` or `-` are
excluded (documented rule). Per-gene values are unweighted means of per-codon
values over codon sites that retain at least one valid observation.

| Index | Formula | Notes |
|-------|---------|--------|
| Shannon raw | `H = −Σ pᵢ ln pᵢ` | nats; same convention as the internal per-position entropy |
| Shannon normalised | `H / ln(k)` | k = number of observed types; 0 if k ≤ 1. Independent of log base (Gregori et al. 2016) |
| Gini–Simpson | `1 − Σ pᵢ²` | |
| Nucleotide π | mean Nei–Li π over the three codon sites | `n/(n−1) · (1 − Σ f²)` per site |
| Amino-acid π | Nei–Li π of translated ACGT codons | one site = one codon |

Gregori, J. et al. (2016). Viral quasispecies complexity measures. *Virology*
493:227–237. https://doi.org/10.1016/j.virol.2016.03.017

Nei, M. & Li, W. H. (1979). Mathematical model for studying genetic variation
in terms of restriction endonucleases. *PNAS* 76:5269–5273.
https://doi.org/10.1073/pnas.76.10.5269

Zanini, F. et al. (2015). Population genomics of intrapatient HIV-1
evolution. *eLife* 4:e11282. https://doi.org/10.7554/eLife.11282

## AAVF v1.0

Opt-in `--aavf PATH` writes Amino Acid Variant Format v1.0
(https://github.com/winhiv/aavf-spec, Apache-2.0):

- First line `##fileformat=AAVFv1.0`
- `##source=codonyat-1.3.0 (research use only)`
- `##reference` is the reference FASTA filename; one `##contig=<ID=...,length=...>` line per FASTA record
- Columns `#CHROM GENE POS REF ALT FILTER ALT_FREQ COVERAGE INFO`
- `POS` is the 1-based amino-acid index along the protein's intervals
- `REF`/`ALT` are IUPAC amino acids; stop is `*`; partial or gapped codons map to `X`
- `FILTER` uses the flags above (not `CAT` / `SRVL`)
- Reserved INFO: `RC`, `AC`, `ACC`, `ACF`; extra `FWC`, `RVC`, `SR`
- Synonymous observations are emitted as `REF=ALT` with every observed codon in `AC`
- True amino-acid indels are not called

This AAVF file carries **no** resistance or surveillance categories.
