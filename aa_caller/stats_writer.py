"""Write the opt-in ``--stats`` TSV (error model + diversity)."""

from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

from .genetic_code import codon_to_aminoacid
from .metrics import (
    classify_flags,
    diversity_from_codon_counts,
    mean_or_zero,
    wilson_ci,
)
from .reference import FullReference


def _blank() -> str:
    return ""


def _fmt(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


STATS_FIELDS = [
    "LEVEL",
    "FILE",
    "REFERENCE",
    "RECORD",
    "PROTEIN",
    "POSITION",
    "AA_POS",
    "VARIANT",
    "AA",
    "COUNT",
    "FREQ",
    "FREQ_CI_LOW",
    "FREQ_CI_HIGH",
    "FW",
    "RV",
    "SB_P",
    "FLAGS",
    "COVERAGE",
    "SHANNON",
    "SHANNON_NORM",
    "GINI_SIMPSON",
    "PI_NT",
    "PI_AA",
    "N_TYPES",
    "N_CODONS",
]


def write_stats_tsv(
    sample_name: str,
    reference: FullReference,
    output_path: Path,
    protein_pairs: Sequence[tuple[str, dict[int, dict]]],
    *,
    min_freq: float,
    min_coverage: int,
    sb_alpha: float,
) -> None:
    """Write variant / codon / gene rows for every requested protein."""
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATS_FIELDS, delimiter="\t")
        writer.writeheader()
        for protein_name, entries in protein_pairs:
            codon_divs: list[dict[str, float]] = []
            ordered = sorted(
                entries.items(),
                key=lambda item: (item[1].get("aa_pos", item[0]), item[0]),
            )
            for pos, entry in ordered:
                record_id = entry.get("record_id", reference.id)
                depth = int(entry["depth"])
                fw_depth = int(entry["fw_depth"])
                rv_depth = int(entry["rv_depth"])
                aa_pos = entry.get("aa_pos", "")
                codon_counts = {variant.codon: float(variant.count) for variant in entry["variants"].values()}
                diversity = diversity_from_codon_counts(codon_counts)
                if diversity["n_valid"] > 0:
                    codon_divs.append(diversity)

                for variant in entry["variants"].values():
                    freq = (variant.count / depth) if depth else 0.0
                    lo, hi = wilson_ci(variant.count, depth)
                    flags, p_value = classify_flags(
                        freq=freq,
                        coverage=depth,
                        fw=variant.fw_reads,
                        rv=variant.rv_reads,
                        other_fw=fw_depth - variant.fw_reads,
                        other_rv=rv_depth - variant.rv_reads,
                        min_freq=min_freq,
                        min_coverage=min_coverage,
                        sb_alpha=sb_alpha,
                    )
                    writer.writerow(
                        {
                            "LEVEL": "variant",
                            "FILE": sample_name,
                            "REFERENCE": record_id,
                            "RECORD": record_id,
                            "PROTEIN": protein_name,
                            "POSITION": pos,
                            "AA_POS": aa_pos,
                            "VARIANT": variant.codon,
                            "AA": codon_to_aminoacid(variant.codon),
                            "COUNT": variant.count,
                            "FREQ": _fmt(freq),
                            "FREQ_CI_LOW": _fmt(lo),
                            "FREQ_CI_HIGH": _fmt(hi),
                            "FW": variant.fw_reads,
                            "RV": variant.rv_reads,
                            "SB_P": _fmt(p_value),
                            "FLAGS": ";".join(flags),
                            "COVERAGE": depth,
                            "SHANNON": _blank(),
                            "SHANNON_NORM": _blank(),
                            "GINI_SIMPSON": _blank(),
                            "PI_NT": _blank(),
                            "PI_AA": _blank(),
                            "N_TYPES": _blank(),
                            "N_CODONS": _blank(),
                        }
                    )

                writer.writerow(
                    {
                        "LEVEL": "codon",
                        "FILE": sample_name,
                        "REFERENCE": record_id,
                        "RECORD": record_id,
                        "PROTEIN": protein_name,
                        "POSITION": pos,
                        "AA_POS": aa_pos,
                        "VARIANT": _blank(),
                        "AA": _blank(),
                        "COUNT": _blank(),
                        "FREQ": _blank(),
                        "FREQ_CI_LOW": _blank(),
                        "FREQ_CI_HIGH": _blank(),
                        "FW": fw_depth,
                        "RV": rv_depth,
                        "SB_P": _blank(),
                        "FLAGS": _blank(),
                        "COVERAGE": depth,
                        "SHANNON": _fmt(diversity["shannon"]),
                        "SHANNON_NORM": _fmt(diversity["shannon_norm"]),
                        "GINI_SIMPSON": _fmt(diversity["gini_simpson"]),
                        "PI_NT": _fmt(diversity["pi_nt"]),
                        "PI_AA": _fmt(diversity["pi_aa"]),
                        "N_TYPES": int(diversity["n_types"]),
                        "N_CODONS": _blank(),
                    }
                )

            first_record = next(iter(entries.values()), {}).get("record_id", reference.id)
            writer.writerow(
                {
                    "LEVEL": "gene",
                    "FILE": sample_name,
                    "REFERENCE": first_record,
                    "RECORD": first_record,
                    "PROTEIN": protein_name,
                    "POSITION": _blank(),
                    "AA_POS": _blank(),
                    "VARIANT": _blank(),
                    "AA": _blank(),
                    "COUNT": _blank(),
                    "FREQ": _blank(),
                    "FREQ_CI_LOW": _blank(),
                    "FREQ_CI_HIGH": _blank(),
                    "FW": _blank(),
                    "RV": _blank(),
                    "SB_P": _blank(),
                    "FLAGS": _blank(),
                    "COVERAGE": _blank(),
                    "SHANNON": _fmt(mean_or_zero(d["shannon"] for d in codon_divs)),
                    "SHANNON_NORM": _fmt(mean_or_zero(d["shannon_norm"] for d in codon_divs)),
                    "GINI_SIMPSON": _fmt(mean_or_zero(d["gini_simpson"] for d in codon_divs)),
                    "PI_NT": _fmt(mean_or_zero(d["pi_nt"] for d in codon_divs)),
                    "PI_AA": _fmt(mean_or_zero(d["pi_aa"] for d in codon_divs)),
                    "N_TYPES": _blank(),
                    "N_CODONS": len(codon_divs),
                }
            )
