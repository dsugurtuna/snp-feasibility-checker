"""Array catalogue module.

Manages genotyping array manifests and SNP content lookups.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, List, Set


@dataclass
class ArrayRecord:
    """Manifest record for a single genotyping array."""

    array_name: str
    snp_count: int = 0
    snp_ids: FrozenSet[str] = field(default_factory=frozenset)

    @property
    def snp_set(self) -> Set[str]:
        return set(self.snp_ids)


class ArrayCatalogue:
    """Catalogue of genotyping array manifests.

    Loads array manifests from CSV files and provides SNP lookup
    across all registered arrays.

    Parameters
    ----------
    arrays : list of ArrayRecord, optional
        Pre-loaded array records.
    """

    def __init__(self, arrays: List[ArrayRecord] | None = None) -> None:
        self._arrays: Dict[str, ArrayRecord] = {}
        if arrays:
            for arr in arrays:
                self._arrays[arr.array_name] = arr

    def register(self, record: ArrayRecord) -> None:
        """Register a genotyping array."""
        self._arrays[record.array_name] = record

    def load_manifest_csv(
        self,
        array_name: str,
        csv_path: str | Path,
        snp_column: str = "snp_id",
    ) -> ArrayRecord:
        """Load array manifest from CSV.

        Parameters
        ----------
        array_name : str
            Name to register the array under.
        csv_path : str or Path
            Path to the manifest CSV.
        snp_column : str
            Column header containing SNP identifiers.
        """
        snps: Set[str] = set()
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                sid = row.get(snp_column, "").strip()
                if sid:
                    snps.add(sid)
        record = ArrayRecord(
            array_name=array_name,
            snp_count=len(snps),
            snp_ids=frozenset(snps),
        )
        self._arrays[array_name] = record
        return record

    def get_array(self, name: str) -> ArrayRecord | None:
        return self._arrays.get(name)

    @property
    def array_names(self) -> List[str]:
        return sorted(self._arrays.keys())

    def find_arrays_containing(self, snp_id: str) -> List[str]:
        """Return names of all arrays that contain a given SNP."""
        return [
            name for name, rec in self._arrays.items()
            if snp_id in rec.snp_ids
        ]

    @property
    def total_unique_snps(self) -> int:
        """Total unique SNPs across all arrays."""
        union: Set[str] = set()
        for rec in self._arrays.values():
            union |= rec.snp_ids
        return len(union)
