# SNP Feasibility Checker

[![CI](https://github.com/dsugurtuna/snp-feasibility-checker/actions/workflows/ci.yml/badge.svg)](https://github.com/dsugurtuna/snp-feasibility-checker/actions/workflows/ci.yml)

**Assess SNP availability across genotyping arrays and estimate recall study yield using Hardy-Weinberg equilibrium.**

Manages genotyping array manifests, checks target SNP coverage across arrays, and estimates achievable carrier/homozygote counts for recall-by-genotype study design.

> **Portfolio project.** Demonstrates generalised SNP feasibility workflows. No real array manifests or participant data are included.

---

## Architecture

```
src/snp_checker/
    __init__.py      # Public API exports
    catalogue.py     # Array manifest catalogue (ArrayCatalogue)
    checker.py       # SNP feasibility assessment (FeasibilityChecker)
    estimator.py     # HWE-based recall yield estimation (RecallEstimator)
tests/
    test_checker.py  # Catalogue and feasibility tests
    test_estimator.py # HWE estimation tests
```

---

## Quick start

```bash
pip install -e ".[dev]"
pytest -v
```

### Python API

```python
from snp_checker import ArrayCatalogue, FeasibilityChecker, RecallEstimator
from snp_checker.catalogue import ArrayRecord

# Register arrays
cat = ArrayCatalogue()
cat.register(ArrayRecord("GSA_v3", 650000, frozenset(["rs429358", "rs7412"])))
cat.load_manifest_csv("Axiom_UKB", "axiom_manifest.csv")

# Check feasibility
checker = FeasibilityChecker(cat)
report = checker.check(["rs429358", "rs7412", "rs999999"])
print(f"Feasibility rate: {report.feasibility_rate:.0%}")

# Estimate recall yield
estimator = RecallEstimator(default_cohort_size=50000)
estimate = estimator.estimate("rs429358", allele_frequency=0.15)
print(f"Expected carriers: {estimate.expected_carriers}")
```

---

## Key features

| Feature | Detail |
| :--- | :--- |
| **Array catalogue** | Register and query genotyping array manifests |
| **CSV manifest loading** | Bulk-load SNP content from manifest files |
| **Feasibility checking** | Per-SNP coverage across all registered arrays |
| **Array overlap** | Intersection of target SNPs with specific arrays |
| **HWE estimation** | Carrier and homozygote count estimation |
| **Batch estimation** | Yield estimates for multiple SNPs at once |

## Development

```bash
make dev        # install with dev dependencies
make test       # run pytest
make lint       # run ruff
make clean      # remove build artefacts
```

## Jira provenance

| Ticket | Description |
| :--- | :--- |
| BIOIN-298 | SNP feasibility assessment for recall study planning |

---

*Created by [dsugurtuna](https://github.com/dsugurtuna)*
