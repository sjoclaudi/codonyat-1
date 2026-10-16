"""AAVF v1.0 writer. Research use only. Never emits CAT or SRVL."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from datetime import date
from pathlib import Path

from .constants import PACKAGE_VERSION
from .genetic_code import codon_to_aminoacid
from .metrics import classify_flags, is_acgt_codon
from .models import Variant
from .reference import FullReference

_AAVF_AA = frozenset("ACDEFGHIKLMNPQRSTVWYX*")


def aavf_amino_acid(codon: str) -> str:
    """Map a codon to an AAVF amino-acid letter (``*`` = stop, ``X`` = other)."""
    aa = codon_to_aminoacid(codon)
    if aa == "_":
        return "*"
    if aa in _AAVF_AA:
        return aa
    return "X"


def mark_changed_bases(ref_codon: str, alt_codon: str) -> str:
    """Lowercase matches to the reference codon; uppercase substitutions.

    Matches the AAVF v1.0 examples (``RC=aaa;AC=aaC``).
    """
    ref = ref_codon.lower()
    alt = alt_codon.lower()
    marked = []
    for i, base in enumerate(alt):
        if i < len(ref) and base == ref[i]:
            marked.append(base)
        else:
            marked.append(base.upper())
    return "".join(marked)


def _fmt_float(value: float) -> str:
    return f"{value:.6g}"


def write_aavf(
    sample_name: str,
    reference: FullReference,
    output_path: Path,
    protein_pairs: Sequence[tuple[str, dict[int, dict]]],
    *,
    min_freq: float,
    min_coverage: int,
    sb_alpha: float,
    file_date: str | None = None,
) -> None:
    """Write an AAVF v1.0 file for the counted proteins.

    FILTER codes are the item-1 flags ``lowfreq``, ``sb`` and ``lowcov``.
    CAT and SRVL are never written. Partial / gapped codons map to ``X``.
    True amino-acid indels are not called.
    """
    dated = file_date or date.today().strftime("%Y%m%d")
    fasta_name = Path(reference.path).name
    contig_lines = [
        f"##contig=<ID={record_id.replace(':', '|')},length={len(reference.sequences[record_id])}>"
        for record_id in reference.record_ids
    ]
    lines = [
        "##fileformat=AAVFv1.0",
        f"##fileDate={dated}",
        f"##source=codonyat-{PACKAGE_VERSION} (research use only)",
        f"##reference={fasta_name}",
        *contig_lines,
        '##INFO=<ID=RC,Number=1,Type=String,Description="Reference codon">',
        '##INFO=<ID=AC,Number=.,Type=String,Description="Alternate codon(s) encoding ALT">',
        '##INFO=<ID=ACC,Number=.,Type=Integer,Description="Alternate codon counts, same order as AC">',
        '##INFO=<ID=ACF,Number=.,Type=Float,Description="Alternate codon frequencies, same order as AC">',
        '##INFO=<ID=FWC,Number=1,Type=Integer,Description="Forward-strand coverage at the codon">',
        '##INFO=<ID=RVC,Number=1,Type=Integer,Description="Reverse-strand coverage at the codon">',
        '##INFO=<ID=SR,Number=1,Type=Float,Description="Forward/reverse coverage ratio at the codon">',
        '##FILTER=<ID=lowfreq,Description="ALT frequency below --min-freq">',
        '##FILTER=<ID=sb,Description="Strand-bias Fisher p-value below --sb-alpha">',
        '##FILTER=<ID=lowcov,Description="Codon coverage below --min-coverage">',
        "#CHROM\tGENE\tPOS\tREF\tALT\tFILTER\tALT_FREQ\tCOVERAGE\tINFO",
    ]

    for protein_name, entries in protein_pairs:
        # aa_pos then genomic position; keeps GENE blocks contiguous.
        ordered = sorted(
            entries.items(),
            key=lambda item: (item[1].get("aa_pos", item[0]), item[0]),
        )
        for pos, entry in ordered:
            depth = int(entry["depth"])
            if depth <= 0:
                continue
            record_id = str(entry.get("record_id", reference.id)).replace(":", "|")
            gene = protein_name.replace(":", "|")
            aa_pos = int(entry.get("aa_pos", 0)) or "."
            ref_codon = reference.get_seq_at(pos, 3, entry.get("record_id"))
            ref_aa = aavf_amino_acid(ref_codon) if ref_codon else "X"
            fw_depth = int(entry["fw_depth"])
            rv_depth = int(entry["rv_depth"])
            sr = (fw_depth / rv_depth) if rv_depth else 0.0

            by_aa: dict[str, list[Variant]] = defaultdict(list)
            for variant in entry["variants"].values():
                by_aa[aavf_amino_acid(variant.codon)].append(variant)

            for alt_aa, variants in by_aa.items():
                alt_count = sum(v.count for v in variants)
                alt_fw = sum(v.fw_reads for v in variants)
                alt_rv = sum(v.rv_reads for v in variants)
                alt_freq = (alt_count / depth) if depth else 0.0
                flags, _p = classify_flags(
                    freq=alt_freq,
                    coverage=depth,
                    fw=alt_fw,
                    rv=alt_rv,
                    other_fw=fw_depth - alt_fw,
                    other_rv=rv_depth - alt_rv,
                    min_freq=min_freq,
                    min_coverage=min_coverage,
                    sb_alpha=sb_alpha,
                )
                filter_field = ";".join(flags) if flags else "PASS"
                ac_parts = []
                acc_parts = []
                acf_parts = []
                for variant in variants:
                    codon = variant.codon
                    labelled = mark_changed_bases(ref_codon, codon) if is_acgt_codon(codon) else codon.lower()
                    ac_parts.append(labelled)
                    acc_parts.append(str(int(variant.count)))
                    acf_parts.append(_fmt_float(variant.count / depth if depth else 0.0))
                info = (
                    f"RC={ref_codon.lower()};"
                    f"AC={','.join(ac_parts)};"
                    f"ACC={','.join(acc_parts)};"
                    f"ACF={','.join(acf_parts)};"
                    f"FWC={fw_depth};RVC={rv_depth};SR={_fmt_float(sr)}"
                )
                lines.append(
                    "\t".join(
                        [
                            record_id,
                            gene,
                            str(aa_pos),
                            ref_aa,
                            alt_aa,
                            filter_field,
                            _fmt_float(alt_freq),
                            str(depth),
                            info,
                        ]
                    )
                )

    output_path.write_text("\n".join(lines) + "\n")
