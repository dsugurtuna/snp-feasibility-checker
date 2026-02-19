"""Feasibility checker module.

Checks whether target SNPs are available on genotyping arrays
and reports coverage across study batches.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set

from .catalogue import ArrayCatalogue


@dataclass
class SNPCoverage:
    """Coverage summary for a single SNP."""

    snp_id: str
    present_on: List[str] = field(default_factory=list)
    missing_from: List[str] = field(default_factory=list)

    @property
    def is_available(self) -> bool:
        return len(self.present_on) > 0

    @property
    def coverage_fraction(self) -> float:
        total = len(self.present_on) + len(self.missing_from)
        if total == 0:
            return 0.0
        return len(self.present_on) / total


@dataclass
class FeasibilityReport:
    """Feasibility assessment for a set of target SNPs."""

    target_snps: int = 0
    available_count: int = 0
    unavailable_snps: List[str] = field(default_factory=list)
    coverage_details: List[SNPCoverage] = field(default_factory=list)
    array_summary: Dict[str, int] = field(default_factory=dict)

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

    def check(self, target_snps: List[str]) -> FeasibilityReport:
        """Assess feasibility for a list of target SNPs."""
        report = FeasibilityReport(target_snps=len(target_snps))
        array_names = self.catalogue.array_names
        array_hit_count: Dict[str, int] = {a: 0 for a in array_names}

        for snp_id in target_snps:
            arrays_with = self.catalogue.find_arrays_containing(snp_id)
            arrays_without = [a for a in array_names if a not in arrays_with]

            cov = SNPCoverage(
                snp_id=snp_id,
                present_on=arrays_with,
                missing_from=arrays_without,
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
        target_snps: List[str],
        array_name: str,
    ) -> Set[str]:
        """Return the intersection of target SNPs with a specific array."""
        arr = self.catalogue.get_array(array_name)
        if arr is None:
            return set()
        return set(target_snps) & arr.snp_set
