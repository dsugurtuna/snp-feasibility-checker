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
