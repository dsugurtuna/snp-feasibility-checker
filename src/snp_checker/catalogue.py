"""Array catalogue module.

Manages genotyping array manifests and SNP content lookups.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ArrayRecord:
    """Manifest record for a single genotyping array."""

    array_name: str
    snp_count: int = 0
    snp_ids: frozenset[str] = field(default_factory=frozenset)
    sample_count: int = 0  # participants genotyped on this array, if known

    @property
    def snp_set(self) -> set[str]:
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

    def __init__(self, arrays: list[ArrayRecord] | None = None) -> None:
        self._arrays: dict[str, ArrayRecord] = {}
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
        sample_count: int = 0,
    ) -> ArrayRecord:
        """Load an array manifest from CSV.

        Parameters
        ----------
        array_name : str
            Name to register the array under.
        csv_path : str or Path
            Path to the manifest CSV.
        snp_column : str
            Column header containing SNP identifiers (rsIDs).
        sample_count : int
            Participants genotyped on this array, if known.

        Raises
        ------
        ValueError
            If ``snp_column`` is not in the header. Without this check a wrong
            column name loads an empty array and every SNP looks unavailable.
        """
        snps: set[str] = set()
        with open(csv_path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            if snp_column not in (reader.fieldnames or []):
                raise ValueError(f"{csv_path}: column {snp_column!r} not found; header is {reader.fieldnames}")
            for row in reader:
                sid = (row.get(snp_column) or "").strip()
                if sid:
                    snps.add(sid)
        record = ArrayRecord(
            array_name=array_name,
            snp_count=len(snps),
            snp_ids=frozenset(snps),
            sample_count=sample_count,
        )
        self._arrays[array_name] = record
        return record

    def get_array(self, name: str) -> ArrayRecord | None:
        return self._arrays.get(name)

    @property
    def array_names(self) -> list[str]:
        return sorted(self._arrays.keys())

    def find_arrays_containing(self, snp_id: str) -> list[str]:
        """Return names of all arrays that contain a given SNP."""
        return sorted(name for name, rec in self._arrays.items() if snp_id in rec.snp_ids)

    @property
    def total_unique_snps(self) -> int:
        """Total unique SNPs across all arrays."""
        union: set[str] = set()
        for rec in self._arrays.values():
            union |= rec.snp_ids
        return len(union)
