"""Tests for RecallEstimator."""

from snp_checker.estimator import RecallEstimator


class TestRecallEstimator:
    def test_hwe_carriers(self):
        # freq=0.1, N=10000 => carriers = (2*0.9*0.1 + 0.01)*10000 = 1900
        count = RecallEstimator._hwe_carriers(0.1, 10000)
        assert count == 1900

    def test_hwe_homozygotes(self):
        # freq=0.1, N=10000 => homozygotes = 0.01*10000 = 100
        count = RecallEstimator._hwe_homozygotes(0.1, 10000)
        assert count == 100

    def test_estimate_single(self):
        est = RecallEstimator(default_cohort_size=50000)
        result = est.estimate("rs429358", 0.15)
        assert result.snp_id == "rs429358"
        assert result.cohort_size == 50000
        assert result.expected_carriers > 0
        assert result.expected_homozygotes > 0

    def test_estimate_custom_cohort(self):
        est = RecallEstimator()
        result = est.estimate("rs123", 0.05, cohort_size=1000)
        assert result.cohort_size == 1000

    def test_estimate_batch(self):
        est = RecallEstimator(default_cohort_size=10000)
        results = est.estimate_batch({"rs1": 0.1, "rs2": 0.2})
        assert len(results) == 2

    def test_rare_variant(self):
        est = RecallEstimator(default_cohort_size=100000)
        result = est.estimate("rs_rare", 0.001)
        assert result.expected_homozygotes < result.expected_carriers


class TestEstimatorCorrectness:
    def test_counts_are_rounded_not_truncated(self):
        # Exact value is 396 (N * (1 - 0.98**2)); int() truncation gave 395.
        assert RecallEstimator._hwe_carriers(0.02, 10000) == 396

    def test_genotype_classes_sum_to_cohort(self):
        est = RecallEstimator().estimate("rsX", 0.3, cohort_size=1000)
        hom_ref = round(0.7 * 0.7 * 1000)
        assert est.expected_heterozygotes == 420
        assert est.expected_homozygotes == 90
        assert est.expected_carriers == 510
        assert hom_ref + est.expected_heterozygotes + est.expected_homozygotes == 1000

    def test_frequency_out_of_range(self):
        import pytest

        with pytest.raises(ValueError):
            RecallEstimator().estimate("rsX", 1.2)

    def test_zero_cohort_is_not_replaced_by_default(self):
        est = RecallEstimator(default_cohort_size=50000).estimate("rsX", 0.1, cohort_size=0)
        assert est.cohort_size == 0
        assert est.expected_carriers == 0
