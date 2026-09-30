"""Recall estimator module.

Expected genotype counts under Hardy-Weinberg equilibrium (HWE) for a given
alternative-allele frequency and number of genotyped participants.

These are expected counts before consent, eligibility, contactability and
response, so they are an upper bound on recall yield, not a forecast.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class RecallEstimate:
    """Expected genotype counts for one SNP."""

    snp_id: str
    allele_frequency: float = 0.0
    cohort_size: int = 0
    expected_carriers: int = 0  # at least one alternative allele: N(1 - p^2)
    expected_homozygotes: int = 0  # two alternative alleles: N q^2
    arrays_available: list[str] = field(default_factory=list)
    expected_heterozygotes: int = 0  # exactly one alternative allele: N 2pq


def _check_inputs(freq: float, n: int) -> None:
    if not 0.0 <= freq <= 1.0:
        raise ValueError(f"allele frequency must be between 0 and 1, got {freq}")
    if n < 0:
        raise ValueError(f"cohort size must be non-negative, got {n}")


class RecallEstimator:
    """Estimate expected carrier and homozygote counts under HWE.

    ``allele_frequency`` is the frequency q of the allele of interest
    (usually the alternative allele); p = 1 - q. Counts are rounded to the
    nearest whole participant.

    Parameters
    ----------
    default_cohort_size : int
        Number of genotyped participants used when no size is given per SNP.
    """

    def __init__(self, default_cohort_size: int = 50000) -> None:
        self.default_cohort_size = default_cohort_size

    @staticmethod
    def _hwe_carriers(freq: float, n: int) -> int:
        """Expected carriers of at least one copy: N(2pq + q^2) = N(1 - p^2)."""
        _check_inputs(freq, n)
        p = 1.0 - freq
        return round((1.0 - p * p) * n)

    @staticmethod
    def _hwe_homozygotes(freq: float, n: int) -> int:
        """Expected homozygotes for the allele of interest: N q^2."""
        _check_inputs(freq, n)
        return round(freq * freq * n)

    @staticmethod
    def _hwe_heterozygotes(freq: float, n: int) -> int:
        """Expected heterozygotes: N 2pq."""
        _check_inputs(freq, n)
        return round(2.0 * (1.0 - freq) * freq * n)

    def estimate(
        self,
        snp_id: str,
        allele_frequency: float,
        cohort_size: int | None = None,
        arrays_available: list[str] | None = None,
    ) -> RecallEstimate:
        """Estimate expected counts for a single SNP.

        ``cohort_size`` should be the number of participants genotyped for
        this SNP (see ``SNPCoverage.genotyped_samples``), not the whole
        cohort, when the SNP is missing from some arrays.
        """
        n = self.default_cohort_size if cohort_size is None else cohort_size
        return RecallEstimate(
            snp_id=snp_id,
            allele_frequency=allele_frequency,
            cohort_size=n,
            expected_carriers=self._hwe_carriers(allele_frequency, n),
            expected_homozygotes=self._hwe_homozygotes(allele_frequency, n),
            arrays_available=arrays_available or [],
            expected_heterozygotes=self._hwe_heterozygotes(allele_frequency, n),
        )

    def estimate_batch(
        self,
        snp_frequencies: dict[str, float],
        cohort_size: int | None = None,
    ) -> list[RecallEstimate]:
        """Estimate expected counts for multiple SNPs."""
        return [self.estimate(snp_id, freq, cohort_size) for snp_id, freq in snp_frequencies.items()]
