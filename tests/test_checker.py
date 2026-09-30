"""Tests for ArrayCatalogue and FeasibilityChecker."""

from snp_checker.catalogue import ArrayCatalogue, ArrayRecord
from snp_checker.checker import FeasibilityChecker, SNPCoverage


def _make_catalogue():
    cat = ArrayCatalogue()
    cat.register(ArrayRecord("GSA_v3", 650000, frozenset(["rs429358", "rs7412", "rs1234"])))
    cat.register(ArrayRecord("Axiom_UKB", 800000, frozenset(["rs429358", "rs5678"])))
    return cat


class TestArrayCatalogue:
    def test_register_and_lookup(self):
        cat = _make_catalogue()
        assert "GSA_v3" in cat.array_names
        assert cat.get_array("GSA_v3").snp_count == 650000

    def test_find_arrays_containing(self):
        cat = _make_catalogue()
        arrays = cat.find_arrays_containing("rs429358")
        assert "GSA_v3" in arrays
        assert "Axiom_UKB" in arrays

    def test_find_arrays_missing(self):
        cat = _make_catalogue()
        arrays = cat.find_arrays_containing("rs999999")
        assert len(arrays) == 0

    def test_total_unique_snps(self):
        cat = _make_catalogue()
        # rs429358, rs7412, rs1234, rs5678 = 4 unique
        assert cat.total_unique_snps == 4

    def test_load_manifest_csv(self, tmp_path):
        csv_file = tmp_path / "manifest.csv"
        csv_file.write_text("snp_id,chr,pos\nrs111,1,100\nrs222,2,200\n")
        cat = ArrayCatalogue()
        rec = cat.load_manifest_csv("TestArray", csv_file)
        assert rec.snp_count == 2
        assert "rs111" in rec.snp_ids


class TestFeasibilityChecker:
    def test_all_available(self):
        cat = _make_catalogue()
        checker = FeasibilityChecker(cat)
        report = checker.check(["rs429358"])
        assert report.available_count == 1
        assert report.feasibility_rate == 1.0

    def test_partial_availability(self):
        cat = _make_catalogue()
        checker = FeasibilityChecker(cat)
        report = checker.check(["rs429358", "rs999999"])
        assert report.available_count == 1
        assert "rs999999" in report.unavailable_snps

    def test_overlap(self):
        cat = _make_catalogue()
        checker = FeasibilityChecker(cat)
        overlap = checker.check_overlap(["rs429358", "rs7412", "rs5678"], "GSA_v3")
        assert overlap == {"rs429358", "rs7412"}

    def test_overlap_missing_array(self):
        cat = _make_catalogue()
        checker = FeasibilityChecker(cat)
        overlap = checker.check_overlap(["rs429358"], "NonExistent")
        assert len(overlap) == 0

    def test_coverage_fraction(self):
        cov = SNPCoverage(snp_id="rs1", present_on=["A", "B"], missing_from=["C"])
        assert abs(cov.coverage_fraction - 2 / 3) < 0.01


class TestSampleCoverage:
    def test_genotyped_samples_follow_arrays(self):
        cat = ArrayCatalogue()
        cat.register(ArrayRecord("A", 3, frozenset(["rs1", "rs2"]), sample_count=1000))
        cat.register(ArrayRecord("B", 2, frozenset(["rs1"]), sample_count=250))
        report = FeasibilityChecker(cat).check(["rs1", "rs2", "rs3"])
        by_id = {c.snp_id: c for c in report.coverage_details}
        assert by_id["rs1"].genotyped_samples == 1250
        assert by_id["rs2"].genotyped_samples == 1000
        assert by_id["rs2"].sample_fraction == 0.8
        assert by_id["rs3"].genotyped_samples == 0

    def test_manifest_wrong_column_raises(self, tmp_path):
        import pytest

        csv_file = tmp_path / "manifest.csv"
        csv_file.write_text("Name,Chr\nrs111,1\n")
        with pytest.raises(ValueError, match="not found"):
            ArrayCatalogue().load_manifest_csv("X", csv_file)

    def test_manifest_with_bom_and_custom_column(self, tmp_path):
        csv_file = tmp_path / "manifest.csv"
        csv_file.write_bytes("﻿Name,Chr\nrs111,1\n".encode())
        rec = ArrayCatalogue().load_manifest_csv("X", csv_file, snp_column="Name")
        assert rec.snp_ids == frozenset({"rs111"})
