"""Offline demo on synthetic data.

The arrays, SNP IDs, sample counts and allele frequencies are all made up.
Run from the repository root:  python examples/demo.py
"""

from __future__ import annotations

from pathlib import Path

from snp_checker import ArrayCatalogue, FeasibilityChecker, RecallEstimator

DATA = Path(__file__).parent / "data"

# Synthetic: participants genotyped on each array, and allele frequencies.
SAMPLES = {"ARRAY_A": 8000, "ARRAY_B": 2000}
FREQUENCIES = {"rs9100001": 0.15, "rs9100002": 0.02, "rs9100004": 0.30, "rs9100009": 0.10}


def main() -> None:
    cat = ArrayCatalogue()
    cat.load_manifest_csv("ARRAY_A", DATA / "array_a_manifest.csv", sample_count=SAMPLES["ARRAY_A"])
    cat.load_manifest_csv("ARRAY_B", DATA / "array_b_manifest.csv", sample_count=SAMPLES["ARRAY_B"])

    report = FeasibilityChecker(cat).check(list(FREQUENCIES))
    print(f"SNPs on at least one array: {report.available_count} of {report.target_snps}")
    print()
    print(f"{'snp':<10} {'arrays':<16} {'typed':>6} {'hom':>5} {'het':>5} {'carriers':>8}")
    estimator = RecallEstimator()
    for cov in report.coverage_details:
        est = estimator.estimate(
            cov.snp_id,
            FREQUENCIES[cov.snp_id],
            cohort_size=cov.genotyped_samples,
            arrays_available=cov.present_on,
        )
        arrays = ",".join(cov.present_on) or "-"
        print(
            f"{cov.snp_id:<10} {arrays:<16} {cov.genotyped_samples:>6} "
            f"{est.expected_homozygotes:>5} {est.expected_heterozygotes:>5} "
            f"{est.expected_carriers:>8}"
        )
    print()
    print("Expected counts under HWE, before consent, eligibility and response.")


if __name__ == "__main__":
    main()
