"""Recall estimator module.

Estimates achievable recall sample sizes based on SNP availability,
allele frequencies, and cohort sizes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List


@dataclass
class RecallEstimate:
    """Estimated recall study yield."""

    snp_id: str
    allele_frequency: float = 0.0
    cohort_size: int = 0
    expected_carriers: int = 0
    expected_homozygotes: int = 0
    arrays_available: List[str] = field(default_factory=list)


class RecallEstimator:
    """Estimate recall study yield from allele frequencies.

    Applies Hardy-Weinberg equilibrium to estimate carrier and
    homozygote counts in a given cohort.

    Parameters
    ----------
    default_cohort_size : int
        Default cohort size when not specified per-SNP.
    """

    def __init__(self, default_cohort_size: int = 50000) -> None:
        self.default_cohort_size = default_cohort_size

    @staticmethod
    def _hwe_carriers(freq: float, n: int) -> int:
        """Estimated heterozygous + homozygous alt carriers (2pq + q²) × N."""
        q = freq
        p = 1.0 - q
        carrier_freq = 2 * p * q + q * q
        return int(carrier_freq * n)

    @staticmethod
    def _hwe_homozygotes(freq: float, n: int) -> int:
        """Estimated homozygous alt count (q²) × N."""
        return int(freq * freq * n)

    def estimate(
        self,
        snp_id: str,
        allele_frequency: float,
        cohort_size: int | None = None,
        arrays_available: List[str] | None = None,
    ) -> RecallEstimate:
        """Estimate yield for a single SNP."""
        n = cohort_size or self.default_cohort_size
        return RecallEstimate(
            snp_id=snp_id,
            allele_frequency=allele_frequency,
            cohort_size=n,
            expected_carriers=self._hwe_carriers(allele_frequency, n),
            expected_homozygotes=self._hwe_homozygotes(allele_frequency, n),
            arrays_available=arrays_available or [],
        )

    def estimate_batch(
        self,
        snp_frequencies: Dict[str, float],
        cohort_size: int | None = None,
    ) -> List[RecallEstimate]:
        """Estimate yield for multiple SNPs."""
        return [
            self.estimate(snp_id, freq, cohort_size)
            for snp_id, freq in snp_frequencies.items()
        ]
