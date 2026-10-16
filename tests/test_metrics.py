"""Reference-value tests for Wilson CI, Fisher, and diversity indices."""

from __future__ import annotations

import math

import pytest

from aa_caller.aavf import aavf_amino_acid, mark_changed_bases
from aa_caller.metrics import (
    aminoacid_pi,
    classify_flags,
    diversity_from_codon_counts,
    fisher_exact_pvalue,
    gini_simpson,
    is_acgt_codon,
    nucleotide_pi,
    pairwise_pi,
    shannon_entropy,
    shannon_normalised,
    wilson_ci,
)


class TestWilsonCI:
    def test_zero_of_twenty_matches_published_bound(self) -> None:
        # Brown, Cai & DasGupta 2001: n=20, x=0, 95% Wilson upper ≈ 0.1611
        low, high = wilson_ci(0, 20)
        assert low == 0.0
        assert high == pytest.approx(0.161132, rel=1e-4)

    def test_six_of_eight(self) -> None:
        # Independent evaluation of the Wilson formula (z = Φ^{-1}(0.975))
        low, high = wilson_ci(6, 8)
        z = 1.959963984540054
        p, n = 0.75, 8.0
        z2 = z * z
        denom = 1.0 + z2 / n
        centre = (p + z2 / (2.0 * n)) / denom
        adj = z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n) / denom
        assert low == pytest.approx(centre - adj, rel=1e-12)
        assert high == pytest.approx(centre + adj, rel=1e-12)
        assert 0.0 < low < 0.75 < high < 1.0

    def test_empty_sample(self) -> None:
        assert wilson_ci(0, 0) == (0.0, 1.0)


class TestFisherExact:
    def test_balanced_3_1_1_3(self) -> None:
        # N=8, K=4, n=4; P(X=3)=16/70. Two-sided sum of tables with P≤16/70 is 34/70.
        p = fisher_exact_pvalue(3, 1, 1, 3)
        assert p == pytest.approx(34 / 70, rel=1e-9)

    def test_tea_tasting_style_table(self) -> None:
        # Fisher's lady-tasting-tea margins analogue [[3,1],[1,3]] already
        # covered; this is the classic [[1,9],[11,3]] table (p ≈ 0.002759)
        p = fisher_exact_pvalue(1, 9, 11, 3)
        assert p == pytest.approx(0.002759456, rel=1e-5)

    def test_zero_margin_is_one(self) -> None:
        assert fisher_exact_pvalue(5, 0, 10, 0) == 1.0
        assert fisher_exact_pvalue(0, 0, 0, 0) == 1.0

    def test_identical_strands_not_significant(self) -> None:
        p = fisher_exact_pvalue(3, 3, 1, 1)
        assert p == pytest.approx(1.0)


class TestDiversity:
    def test_shannon_two_types_matches_hand_value(self) -> None:
        # 6 + 2 observations: p = 0.75, 0.25
        expected = -(0.75 * math.log(0.75) + 0.25 * math.log(0.25))
        assert shannon_entropy([6, 2]) == pytest.approx(expected, rel=1e-12)
        assert shannon_normalised([6, 2]) == pytest.approx(expected / math.log(2), rel=1e-12)

    def test_shannon_single_type_is_zero(self) -> None:
        assert shannon_entropy([8]) == 0.0
        assert shannon_normalised([8]) == 0.0

    def test_gini_simpson_hand_value(self) -> None:
        # 1 - (0.75^2 + 0.25^2) = 1 - (0.5625 + 0.0625) = 0.375
        assert gini_simpson([6, 2]) == pytest.approx(0.375, rel=1e-12)

    def test_pairwise_pi_finite_sample(self) -> None:
        # n/(n-1) * 0.375 = 8/7 * 0.375
        assert pairwise_pi([6, 2]) == pytest.approx(8 / 7 * 0.375, rel=1e-12)
        assert pairwise_pi([1]) == 0.0

    def test_nucleotide_pi_one_substitution(self) -> None:
        # ATG×6, TTG×2: only site 1 varies (A vs T)
        site1 = 8 / 7 * 0.375
        assert nucleotide_pi({"ATG": 6, "TTG": 2}) == pytest.approx(site1 / 3.0, rel=1e-12)

    def test_aminoacid_pi_matches_pairwise(self) -> None:
        # ATG→M, TTG→L
        assert aminoacid_pi({"ATG": 6, "TTG": 2}) == pytest.approx(8 / 7 * 0.375, rel=1e-12)

    def test_excludes_gaps_and_n(self) -> None:
        div = diversity_from_codon_counts({"ATG": 6, "---": 2, "NNA": 1, "TTG": 2})
        assert div["n_types"] == 2
        assert div["n_valid"] == 8
        assert div["shannon"] == pytest.approx(shannon_entropy([6, 2]))

    def test_is_acgt_codon(self) -> None:
        assert is_acgt_codon("ATG")
        assert not is_acgt_codon("---")
        assert not is_acgt_codon("ATN")
        assert not is_acgt_codon("AT")


class TestFlags:
    def test_lowfreq_and_lowcov(self) -> None:
        flags, _p = classify_flags(
            freq=0.005,
            coverage=4,
            fw=1,
            rv=1,
            other_fw=1,
            other_rv=1,
            min_freq=0.01,
            min_coverage=10,
            sb_alpha=0.05,
        )
        assert "lowfreq" in flags
        assert "lowcov" in flags

    def test_strand_bias_flag(self) -> None:
        flags, p_value = classify_flags(
            freq=0.5,
            coverage=20,
            fw=10,
            rv=0,
            other_fw=0,
            other_rv=10,
            min_freq=0.01,
            min_coverage=10,
            sb_alpha=0.05,
        )
        assert p_value < 0.05
        assert "sb" in flags
        assert "lowfreq" not in flags
        assert "lowcov" not in flags


class TestAavfHelpers:
    def test_stop_and_partial(self) -> None:
        assert aavf_amino_acid("TAA") == "*"
        assert aavf_amino_acid("---") == "X"
        assert aavf_amino_acid("ATG") == "M"

    def test_mark_changed_bases(self) -> None:
        assert mark_changed_bases("AAA", "AAC") == "aaC"
        assert mark_changed_bases("TCA", "TAA") == "tAa"
