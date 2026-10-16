"""Wilson CI, Fisher's exact strand-bias test, and diversity indices.

All formulas are implemented with the standard library. Citations are
collected in docs/methods.md. Research use only.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence

from .constants import WILSON_Z_95
from .genetic_code import codon_to_aminoacid

_ACGT = frozenset("ACGT")


def is_acgt_codon(codon: str) -> bool:
    """True when *codon* is three unambiguous DNA bases (no N or gap)."""
    return len(codon) == 3 and all(base in _ACGT for base in codon.upper())


def wilson_ci(
    successes: int | float,
    n: int | float,
    z: float = WILSON_Z_95,
) -> tuple[float, float]:
    """95% Wilson score interval for a binomial proportion (Wilson 1927).

    Returns ``(low, high)`` clipped to ``[0, 1]``. An empty sample yields
    ``(0.0, 1.0)``.
    """
    n = float(n)
    if n <= 0:
        return (0.0, 1.0)
    successes = float(successes)
    if successes < 0 or successes > n:
        raise ValueError("successes must lie between 0 and n")
    p = successes / n
    z2 = z * z
    denom = 1.0 + z2 / n
    centre = (p + z2 / (2.0 * n)) / denom
    adj = z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n) / denom
    return (max(0.0, centre - adj), min(1.0, centre + adj))


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_pvalue(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher's exact p-value for a 2×2 table (Fisher 1935).

    The table is::

        variant     a   b
        other       c   d

    Tables with a zero margin are defined as p = 1. The two-sided p-value
    sums hypergeometric probabilities of every table (same margins) whose
    probability is less than or equal to the observed table.
    """
    a, b, c, d = (int(a), int(b), int(c), int(d))
    if min(a, b, c, d) < 0:
        raise ValueError("Fisher table counts must be non-negative")
    if (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
        return 1.0

    n_total = a + b + c + d
    k_col = a + c
    n_row = a + b
    lo = max(0, n_row + k_col - n_total)
    hi = min(n_row, k_col)
    log_n = _log_comb(n_total, n_row)

    def logp(x: int) -> float:
        return _log_comb(k_col, x) + _log_comb(n_total - k_col, n_row - x) - log_n

    observed = logp(a)
    total = 0.0
    for x in range(lo, hi + 1):
        lp = logp(x)
        if lp <= observed + 1e-12:
            total += math.exp(lp)
    return min(1.0, total)


def shannon_entropy(counts: Sequence[float]) -> float:
    """Shannon entropy in nats: ``H = -Σ p_i ln p_i`` (zero-count types skipped)."""
    total = float(sum(counts))
    if total <= 0:
        return 0.0
    entropy = 0.0
    for count in counts:
        if count <= 0:
            continue
        p = count / total
        entropy -= p * math.log(p)
    return entropy


def shannon_normalised(counts: Sequence[float]) -> float:
    """``H / ln(k)`` for k observed (positive-count) types; 0 when k ≤ 1.

    Equal to Gregori et al. (2016) normalised Shannon regardless of log base.
    """
    positive = [count for count in counts if count > 0]
    if len(positive) <= 1:
        return 0.0
    return shannon_entropy(positive) / math.log(len(positive))


def gini_simpson(counts: Sequence[float]) -> float:
    """Gini–Simpson index ``1 - Σ p_i²``."""
    total = float(sum(counts))
    if total <= 0:
        return 0.0
    return 1.0 - sum((count / total) ** 2 for count in counts if count > 0)


def pairwise_pi(counts: Sequence[float]) -> float:
    """Nei–Li pairwise diversity ``n/(n-1) · (1 - Σ p_i²)``; 0 when n < 2."""
    total = float(sum(counts))
    if total < 2:
        return 0.0
    return (total / (total - 1.0)) * (1.0 - sum((count / total) ** 2 for count in counts if count > 0))


def nucleotide_pi(codon_counts: Mapping[str, float]) -> float:
    """Mean Nei–Li π over the three codon sites (ACGT codons only)."""
    valid = {codon.upper(): float(count) for codon, count in codon_counts.items() if is_acgt_codon(codon) and count > 0}
    if not valid:
        return 0.0
    site_pis: list[float] = []
    for site in range(3):
        bases: dict[str, float] = {}
        for codon, count in valid.items():
            bases[codon[site]] = bases.get(codon[site], 0.0) + count
        site_pis.append(pairwise_pi(list(bases.values())))
    return sum(site_pis) / 3.0


def aminoacid_pi(codon_counts: Mapping[str, float]) -> float:
    """Nei–Li pairwise amino-acid diversity at one codon site (ACGT codons only)."""
    aa_counts: dict[str, float] = {}
    for codon, count in codon_counts.items():
        if not is_acgt_codon(codon) or count <= 0:
            continue
        aa = codon_to_aminoacid(codon)
        aa_counts[aa] = aa_counts.get(aa, 0.0) + float(count)
    return pairwise_pi(list(aa_counts.values()))


def diversity_from_codon_counts(codon_counts: Mapping[str, float]) -> dict[str, float]:
    """Shannon (raw + normalised), Gini–Simpson, π_nt and π_aa for one codon site."""
    valid_counts = [float(count) for codon, count in codon_counts.items() if is_acgt_codon(codon) and count > 0]
    return {
        "shannon": shannon_entropy(valid_counts),
        "shannon_norm": shannon_normalised(valid_counts),
        "gini_simpson": gini_simpson(valid_counts),
        "pi_nt": nucleotide_pi(codon_counts),
        "pi_aa": aminoacid_pi(codon_counts),
        "n_types": float(len(valid_counts)),
        "n_valid": float(sum(valid_counts)),
    }


def mean_or_zero(values: Iterable[float]) -> float:
    items = list(values)
    if not items:
        return 0.0
    return sum(items) / len(items)


def classify_flags(
    *,
    freq: float,
    coverage: int,
    fw: int,
    rv: int,
    other_fw: int,
    other_rv: int,
    min_freq: float,
    min_coverage: int,
    sb_alpha: float,
) -> tuple[list[str], float]:
    """Return ``(flags, fisher_p)`` using the opt-in error-model thresholds."""
    flags: list[str] = []
    if coverage < min_coverage:
        flags.append("lowcov")
    if freq < min_freq:
        flags.append("lowfreq")
    p_value = fisher_exact_pvalue(fw, rv, other_fw, other_rv)
    if p_value < sb_alpha:
        flags.append("sb")
    return flags, p_value
