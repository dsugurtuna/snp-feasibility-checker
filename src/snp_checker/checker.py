"""Feasibility checker module.

Checks whether target SNPs are on each registered genotyping array and, when
per-array sample counts are known, how many participants have them typed.
Each participant is assumed to be genotyped on exactly one array; if some are
typed on several, ``genotyped_samples`` over-counts them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .catalogue import ArrayCatalogue


@dataclass
class SNPCoverage:
    """Coverage summary for a single SNP."""

    snp_id: str
    present_on: list[str] = field(default_factory=list)
    missing_from: list[str] = field(default_factory=list)
    genotyped_samples: int = 0  # participants on arrays that carry the SNP
    total_samples: int = 0  # participants across all registered arrays

    @property
    def is_available(self) -> bool:
        return len(self.present_on) > 0

    @property
    def coverage_fraction(self) -> float:
        """Fraction of registered arrays that carry the SNP."""
        total = len(self.present_on) + len(self.missing_from)
        if total == 0:
            return 0.0
        return len(self.present_on) / total

    @property
    def sample_fraction(self) -> float:
        """Fraction of participants genotyped for the SNP (needs sample counts)."""
        if self.total_samples == 0:
            return 0.0
        return self.genotyped_samples / self.total_samples


@dataclass
class FeasibilityReport:
    """Feasibility assessment for a set of target SNPs."""

    target_snps: int = 0
    available_count: int = 0
    unavailable_snps: list[str] = field(default_factory=list)
    coverage_details: list[SNPCoverage] = field(default_factory=list)
    array_summary: dict[str, int] = field(default_factory=dict)

    @property
    def feasibility_rate(self) -> float:
        if self.target_snps == 0:
            return 0.0
        return self.available_count / self.target_snps


class FeasibilityChecker:
    """Check SNP availability across genotyping arrays.

    Parameters
    ----------
    catalogue : ArrayCatalogue
        Catalogue of registered arrays.
    """

    def __init__(self, catalogue: ArrayCatalogue) -> None:
        self.catalogue = catalogue

    def check(self, target_snps: list[str]) -> FeasibilityReport:
        """Assess feasibility for a list of target SNPs."""
        report = FeasibilityReport(target_snps=len(target_snps))
        array_names = self.catalogue.array_names
        array_hit_count: dict[str, int] = {a: 0 for a in array_names}
        samples = {a: rec.sample_count for a in array_names if (rec := self.catalogue.get_array(a))}
        total_samples = sum(samples.values())

        for snp_id in target_snps:
            arrays_with = self.catalogue.find_arrays_containing(snp_id)
            arrays_without = [a for a in array_names if a not in arrays_with]

            cov = SNPCoverage(
                snp_id=snp_id,
                present_on=arrays_with,
                missing_from=arrays_without,
                genotyped_samples=sum(samples[a] for a in arrays_with),
                total_samples=total_samples,
            )
            report.coverage_details.append(cov)

            if cov.is_available:
                report.available_count += 1
                for a in arrays_with:
                    array_hit_count[a] += 1
            else:
                report.unavailable_snps.append(snp_id)

        report.array_summary = array_hit_count
        return report

    def check_overlap(
        self,
        target_snps: list[str],
        array_name: str,
    ) -> set[str]:
        """Return the intersection of target SNPs with a specific array."""
        arr = self.catalogue.get_array(array_name)
        if arr is None:
            return set()
        return set(target_snps) & arr.snp_set
