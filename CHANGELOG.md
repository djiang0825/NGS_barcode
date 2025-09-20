# Changelog

## Version 3.2.2 – Minor fix in documentation
- Updated installation instructions

---

## Version 3.2.1 – Minor fix in documentation
- Updated URL

---

## Version 3.2.0 – PyPI initial release
- Prepared and published the package to PyPI
- Added citation information
- Updated installation instructions in README
- Separated changelog into dedicated CHANGELOG.md file

---

## Version 3.1 – Minor improvements
- Added `--version` flag to CLI
- Enhanced help system with better descriptions and examples

---

## Version 3.0 – First packaged release
- Major refactoring: Restructured as proper Python package with unified CLI
- New CLI: `barcadia generate` and `barcadia validate` commands
- Installation: Now installable via `pip install -e .`
- Enhanced documentation and code organization

---

## Version 2.4
- Enhanced documentation and code organization (added ExistingSequenceSet class to facilitate file loading for generation and validation)

---

## Version 2.3
- Unified parallel processing logic between generation and validation
- Enhanced documentation and code organization

---

## Version 2.2
- Improved adaptive distance method selection criteria for generation and validation
- Optimized argument validation logic for generation (added GV and Hamming bounds calculation when no seeds or equal-length) and validation
- Updated README with new benchmarking results for generation and validation
- Enhanced documentation and code organization

---

## Version 2.1
- Implemented adaptive generation algorithm selection (neighbor enumeration vs pairwise) with significantly improved performance
- Updated README with new benchmarking results for generation
- Enhanced documentation and code organization

---

## Version 2.0
- Enhanced paired barcode generation with seed files (added `--paired-seed1`, `--paired-seed2`)
- Implemented adaptive validation algorithm selection (neighbor enumeration vs pairwise) with significantly improved performance
- Updated README with new benchmarking results for validation
- Optimized progress logging and violation reporting for validation
- Implemented early-stopping optimization for Levenshtein distance calculation
- Enhanced documentation and code organization

---

## Version 1.1
- Initial release with updated readme file
