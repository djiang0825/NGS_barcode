"""
Barcadia: High-performance DNA barcode generation and validation for NGS applications.

This package provides efficient algorithms for generating and validating DNA barcodes
with configurable quality filters including GC content, homopolymer repeats, and
minimum edit distance constraints.

Public API:
    generate_barcodes_core: Generate DNA barcodes with iterative growth algorithm
    validate_barcodes_core: Validate DNA barcodes against quality filters
"""

# Public API - only expose the core functions
from .filter_utils import Filter
from .generate_barcodes import generate_barcodes_core
from .validate_barcodes import validate_barcodes_core

__all__ = [
    "Filter",
    "generate_barcodes_core",
    "validate_barcodes_core",
]
