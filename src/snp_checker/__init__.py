"""SNP Feasibility Checker — assess SNP availability across genotyping arrays."""

__version__ = "1.0.0"

from .catalogue import ArrayCatalogue, ArrayRecord
from .checker import FeasibilityChecker, FeasibilityReport
from .estimator import RecallEstimate, RecallEstimator

__all__ = [
    "ArrayCatalogue",
    "ArrayRecord",
    "FeasibilityChecker",
    "FeasibilityReport",
    "RecallEstimate",
    "RecallEstimator",
]
