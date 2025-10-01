# Barcadia (v3.3.0)  
*Best-in-class toolkit for large-scale NGS barcode generation and validation* 

![version](https://img.shields.io/badge/version-3.3.0-blue)  
![license](https://img.shields.io/badge/license-Apache%202.0-brightgreen)  
![platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey) 

---

Barcadia is a **fast, memory-efficient, open-source toolkit** for generating and validating massive DNA barcode libraries for next-generation sequencing (NGS) applications. Designed for speed and scalability, it **outperforms existing tools** for both small- and large-scale operations.

- **High performance & scalability** – Generates 100K barcodes in under 3 minutes and 1M in 40 minutes, largely outperforming existing tools limited to ~100K barcodes and requiring hours of compute. 
- **Memory & compute efficient** – Runs on standard laptops with minimal resources (under 1 GB RAM used for generating 1M barcodes); no multi-core processing required.  
- **Extended functionality** – Supports paired barcode generation for dual-indexing, extension from seed lists, validation of existing barcode sets, and estimation of library size limits — features not found in other tools.  

Barcadia makes it easy to design small or large NGS barcode sets that are optimized for **robust performance in high-throughput sequencing workflows**.   

## Table of Contents
- [Background](#background)
  - [Problem Statement](#problem-statement)
  - [Existing Methods and their Limitations](#existing-methods-and-their-limitations)
  - [This Toolkit and its Advantages](#this-toolkit-and-its-advantages)
- [Default Filter Parameters](#default-filter-parameters)
- [Theoretical Bounds](#theoretical-bounds)
- [Benchmarking Highlights](#benchmarking-highlights)
- [Installation](#installation)
  - [Requirements](#requirements)
  - [Setup](#setup)
- [Quick Start](#quick-start)
- [Package Overview](#package-overview)
  - [Project Structure](#project-structure)
  - [Main Commands](#main-commands)
  - [Shared Modules](#shared-modules)
  - [Utility Scripts](#utility-scripts)
- [Citation](#citation)
- [Changelog](#changelog)

## Background

### Problem Statement

**Next-generation sequencing (NGS) is a high-throughput method that enables millions of DNA fragments to be sequenced in parallel, serving as the core technology for decoding genomes across all living organisms**. In this process, researchers use DNA barcodes to uniquely label and track individual biomolecules. Common examples include multiplex sample indexes and unique molecular identifiers (UMIs). 

The practical utility of a DNA barcode library depends on controlling key features: **GC content** (the percentage of G and C nucleotides), **homopolymer repeats** (the length of the longest stretch of identical nucleotides), and **edit distance** (the measure of dissimilarity among sequences). GC content and homopolymer repeats are **within-sequence** criteria that ensure molecular stability during sequencing and are computationally inexpensive to evaluate. **In contrast, edit distance is a between-sequence constraint that provides error tolerance during analysis, but is computationally demanding to assess for large datasets**.

### Existing Methods and their Limitations

For a set of n sequences, the number of pairwise distance comparisons grows quadratically:

$$\binom{n}{2} = \frac{n(n-1)}{2}$$

**This makes brute-force approaches computationally prohibitive for large datasets**. Apart from the approach described by [Lyons et al. (2017)](https://www.nature.com/articles/s41598-017-12825-2), existing tools are not built to accommodate large-scale barcode library design (≥100K sequences). Lyons et al. circumvented computational limitations using a probabilistic Markov chain model; however, their resulting barcode sets are no longer accessible (invalid URLs in the paper), and the underlying code was not released (likely due to proprietary restrictions). 

### This Toolkit and its Advantages

**Here, I introduce Barcadia, a toolkit for efficient large-scale NGS barcode generation that integrates modern computational optimization with novel distance-constrained algorithms I developed to deliver best-in-class scalability and speed**. The software is openly available to promote reproducibility and is designed to run efficiently on minimal computing resources (e.g., standard laptops), ensuring broad accessibility.

In comparison with [TagGD](https://doi.org/10.1371/journal.pone.0057521)—the only other open-source software reported to support barcode generation at the scale of up to 100,000 sequences—Barcadia generated 20,000 18-bp barcodes in just **1 minute** using only 100 MB of RAM on a comparable 8-core laptop, as opposed to the **5 minutes** highlighted in the TagGD abstract. For 100,000 18-bp barcodes, Barcadia completed the process in **under 5 minutes** with 180 MB of RAM, which is significantly faster than the **1.5 hours** reported in Table 1 of the TagGD paper. 

**Notably, Barcadia can handle the generation of million-scale barcodes within reasonable time (i.e., 40 minutes) on minimal compute setups (e.g., standard laptops; <1 GB RAM and no multi-core processing required), far exceeding the capacity of TagGD and other existing tools**. More detailed benchmarking results are presented in a later section of this document. 

**Additionally, it offers unique features not found in other tools**: paired barcode generation for dual-indexing applications, extension from user-provided seed sequences, and a comprehensive validation script for assessing the quality of existing barcode sets (which together enable the estimation of realistic bounds for achievable library sizes under all specified constraints).

## Default Filter Parameters

Barcadia uses carefully chosen default parameters that **optimize synthetic stability and sequencing reliability** in NGS workflows. Users can configure these based on their specific needs.

### GC Content: 40-60%
- **Rationale**: Sequences with extreme GC content (very low or very high) can form secondary structures and exhibit poor amplification efficiency during PCR
- **Impact**: Moderate GC content ensures reliable sequencing performance across different platforms

### Maximum Homopolymer Length: 2
- **Rationale**: Long homopolymer runs (e.g., AAAA, TTTT) cause sequencing errors on most NGS platforms, particularly Illumina and Ion Torrent
- **Impact**: Short homopolymers prevent read quality degradation and reduce sequencing artifacts

### Minimum Hamming Distance: 3
- **Rationale**: Distance ≥3 allows correction of single-base sequencing errors while maintaining sequence distinguishability
- **Impact**: Provides error tolerance for typical NGS error rates (~0.1-1% per base)

## Theoretical Bounds

Leveraging established results from coding theory, **one can calculate lower and upper bounds on the number of valid barcode sequences under specified edit distance constraints**, assuming equal-length barcodes and applying the Hamming distance metric, which counts the base-pair mismatches.

### Gilbert-Varshamov Bound (Lower Bound)

The Gilbert-Varshamov bound provides a lower bound guarantee that a code of at least this size exists. For DNA sequences (i.e., alphabet of 4 characters: A, T, G, C), it states:

$$M \geq \frac{4^n}{V(n,d-1)}$$

where `V(n,d-1)` is the volume of a Hamming sphere of radius `d-1`:

$$V(n,d-1) = \sum_{i=0}^{d-1} \binom{n}{i} \cdot 3^i$$

### Hamming Bound (Upper Bound)

The Hamming bound provides the theoretical maximum number of codewords that can exist for a given sequence length and minimum distance constraint. For DNA sequences of length `n` with minimum Hamming distance `d`, the bound is:

$$M \leq \frac{4^n}{V(n,t)}$$

where `t = ⌊(d-1)/2⌋` and `V(n,t)` is the volume of a Hamming sphere of radius `t`:

$$V(n,t) = \sum_{i=0}^{t} \binom{n}{i} \cdot 3^i$$

**As the sequence space is exhausted in search of valid barcodes, approaching the theoretical upper bound causes the search to slow down progressively**. In particular, a significant—often exponential—slowdown can be expected once the target size surpasses the Gilbert-Varshamov (GV) bound.

For typical barcode applications (6-16 bp, as longer sequences may introduce molecular complexity) using the default minimum distance of 3, the theoretical bounds are summarized below:

<div align="center">

| Length (bp) | GV Bound (Lower) | Hamming Bound (Upper) |
|:-----------:|:----------------:|:---------------------:|
| 6           | 26               | 215                   |
| 7           | 77               | 744                   |
| 8           | 236              | 2.6K                  |
| 9           | 744              | 9.4K                  |
| 10          | 2.4K             | 34K                   |
| 11          | 7.9K             | 123K                  |
| 12          | 27K              | 453K                  |
| 13          | 90K              | 1.7M                  |
| 14          | 311K             | 6.2M                  |
| 15          | 1.1M             | 23M                   |
| 16          | 3.8M             | 88M                   |

</div>

When no seeds are provided—or when all provided seeds are the same length as the new barcodes to generate—**Barcadia will automatically calculate and display the theoretical bounds**. For the specified barcode length and minimum distance, it issues a warning if the requested count exceeds the GV lower bound (performance slowdown expected), or raises an error if the request is beyond the Hamming upper bound (not achievable).

However, these theoretical bounds only capture the minimum distance constraints. **In practice, within-sequence biological filters (GC content and homopolymer restrictions) will further reduce the achievable library size**. Notably, Barcadia can be used to estimate more realistic bounds that account for all three constraints, following the procedures outlined below:

1. **Generate sequences with only distance constraints** (no biological filters):
   ```bash
   barcadia generate --count <large_number> --length <your_length> --gc-min 0 --gc-max 1 --homopolymer-max <your_length> --output-dir <your_output_dir>
   ```
   *Note: Use a count near the GV bound for your parameters, and set homopolymer-max to your barcode length. This overrides the default biological filters (gc-min=0.4, gc-max=0.6, homopolymer-max=2).*

2. **Check how many pass biological filters** (skipping distance validation):
   ```bash
   barcadia validate --input <your_output_dir>/barcodes.txt --skip-distance --output-dir <your_output_dir>
   ```
   *Note: By default, this validation step uses the biological filters (gc-min=0.4, gc-max=0.6, homopolymer-max=2) to check how many sequences pass. You can also set custom filter values if desired.*

3. **Calculate empirical bounds**: Multiply the theoretical bounds by the observed pass rate from step 2 to get realistic achievable library sizes under all constraints.

## Benchmarking Highlights

Below are performance benchmarks for **barcode generation** on a MacBook Pro 2019 (8-core, 16GB RAM) using default parameters with target sizes near the Gilbert-Varshamov bound (generating equal-length libraries from scratch without seed sequences provided):

<div align="center">

| Length (bp) | Library Size | Time (Peak Memory)    |
|:-----------:|:------------:|:---------------------:|
| 6           | 10           | 0.2 seconds (82 MB)   |
| 8           | 100          | 0.2 seconds (82 MB)   |
| 10          | 1,000        | 0.6 seconds (82 MB)   |
| 12          | 10,000       | 12.3 seconds (88 MB)  |
| 14          | 100,000      | 2.8 minutes (160 MB)  |
| 16          | 1,000,000    | 40 minutes (950 MB)   |

</div>

Below are performance benchmarks for **barcode validation** on a MacBook Pro 2019 (8-core, 16GB RAM) using default parameters with pool sizes near the Gilbert-Varshamov bound (validating barcodes that pass all the within-sequence and between-sequence filters):

<div align="center">

| Length (bp) | Library Size | Time (Peak Memory)    |
|:-----------:|:------------:|:---------------------:|
| 6           | 10           | 0.2 seconds (81 MB)   |
| 8           | 100          | 0.2 seconds (81 MB)   |
| 10          | 1,000        | 0.5 seconds (81 MB)   |
| 12          | 10,000       | 6.1 seconds (88 MB)   |
| 14          | 100,000      | 1.5 minutes (160 MB)  |
| 16          | 1,000,000    | 20 minutes (943 MB)   |

</div>

## Installation

### Requirements
- Python 3.12+
- Dependencies:
  - numpy==2.2.6
  - numba==0.61.2
  - llvmlite==0.44.0
  - psutil==7.0.0

*(For development, Ruff is used for linting and formatting to ensure code quality.)*

### Setup

#### Optional: Virtual Environment Setup

```bash
# Create and activate virtual environment
python -m venv barcode-env
source barcode-env/bin/activate
```

#### Method 1: Install from PyPI (Recommended)
```bash
pip install barcadia
```

#### Method 2: Install from Source (Development)
```bash
git clone https://github.com/djiang0825/NGS_barcode.git
cd NGS_barcode
pip install -e .
```

This installs Barcadia as a Python package with the `barcadia` command-line tool. Dependencies are automatically installed via pyproject.toml for both installation methods. The complete installation (including all dependencies) requires approximately 260MB of disk space.

## Quick Start

```bash
# Generate 1000 barcodes of length 12
barcadia generate --count 1000 --length 12

# Validate existing barcodes
barcadia validate --input test/barcodes.txt

# Check version
barcadia --version

# Get help for any command
barcadia --help
barcadia generate --help
barcadia validate --help
```

## Package Overview

### Project Structure

```
barcadia/                                       # Main package
├── __init__.py                                 # Public API
├── cli.py                                      # Command-line interface
├── generate_barcodes.py                        # Barcode generation
├── validate_barcodes.py                        # Barcode validation
├── config_utils.py                             # Configuration utilities
├── filter_utils.py                             # Core filtering utilities
└── tools/                                      # Utility scripts
    ├── generate_random_sequences.py            # Random sequence generator
    └── memory_benchmark.py                     # Performance monitoring
```

### Main Commands

#### 1. `barcadia generate` - Barcode Generation

**Purpose**: Generate high-performance NGS barcodes efficiently using a novel iterative growth algorithm (extension from seeds and paired mode supported).

**Algorithm Overview**:
1. Load seed sequence files as existing pool and report length distribution (will generate from scratch if no seeds are provided)
2. Generate a batch of unique random sequence candidates passing biological filters (GC content, homopolymer repeats)
3. **Two-step distance filtering with adaptive method selection (neighbor enumeration or pairwise)**:
   - Filter candidates against existing pool (if pairwise, parallelized for large sets when multiple CPUs available)
   - Filter remaining candidates against each other within the current batch (always sequential)
4. Add verified batch of passing candidates to existing pool and repeat until target count reached

**Key Features**:
- Uses Hamming distance for generation when no seeds are provided (always equal-length)
- Supports **extension from seed sequence** file(s) (when `--seeds` is specified) and can accommodate: 
  - Multiple seed files of equal or variable lengths (concatenated automatically as seed pool)
  - Differences in lengths between seed pool and newly-generated sequences (NOTE: new sequences are always the same length as specified by `--length`)
  - Uses Hamming distance for equal-length comparisons, Levenshtein distance for mixed lengths (compared to seed pool)
  - Incompatible with `--paired` mode
- Supports **paired barcode generation** (when `--paired` flag is on) by generating 2x the target count and splitting output into 2 files with suffixes `_paired1.txt` and `_paired2.txt`
- Supports paired barcode generation with **paired seed** files (when`--paired-seed1` and `--paired-seed2` are specified), but is more restrictive than `--seeds` in unpaired mode:
  - Only accepts one file each and both must be specified
  - Both files must have the same number of sequences with all sequences being the same length within and across both files
  - Similar to `--seeds`, can accommodate differences in lengths between seed pool and newly-generated sequences (Hamming for equal/Levenshtein for mixed)
  - Incompatible with `--seeds`
- **Adaptive method selection for distance filtering**:
  - Small barcode sets (<10K sequences including seeds): Pairwise distance checking (fast for small sets)
  - Large mixed-length (within seeds and/or between seeds and new barcodes): Pairwise distance checking (neighbor enumeration requires complex Levenshtein handling)
  - Large equal-length (no seeds or everything equal-length): Choose between pairwise and neighbor enumeration based on min_distance
    * Pairwise distance checking: when min_distance > 4 (large number of neighbors to check)
    * Neighbor enumeration: when min_distance <= 4 (limited number of neighbors to check)

**Basic Usage**:
```bash
# Generate 1000 barcodes of length 12 from scratch (no seeds)
barcadia generate --count 1000 --length 12

# Build from a seed sequence file
barcadia generate --count 1000 --length 12 --seeds seed.txt

# Build from multiple seed sequence files (concatenated automatically)
barcadia generate --count 1000 --length 12 --seeds seed_file1.txt seed_file2.txt

# Generate paired barcodes from scratch (no seeds)
barcadia generate --count 1000 --length 12 --paired

# Generate paired barcodes with paired seed files
barcadia generate --count 1000 --length 12 --paired --paired-seed1 seed_paired1.txt --paired-seed2 seed_paired2.txt
```

**Required Arguments**:
- `--count`: Number of barcodes or barcode pairs to generate
- `--length`: Length of barcodes or barcode pairs to generate

**Optional Arguments**:
- `--gc-min`: Minimum GC content (default: 0.4)
- `--gc-max`: Maximum GC content (default: 0.6)
- `--homopolymer-max`: Maximum homopolymer repeat length (default: 2)
- `--min-distance`: Minimum edit distance between sequences (default: 3)
- `--cpus`: Number of CPU cores to use during the parallel filtering step (default: all available)
- `--seeds`: Seed sequence files (any number of files, one sequence per line as .txt; if not provided, will generate from scratch; incompatible with --paired mode; default: None)
- `--paired`: Generate paired barcodes (doubles target count, splits output randomly into two equal parts; incompatible with --seeds; default: off)
- `--paired-seed1`: Paired seed sequence file 1 (used only with --paired and --paired-seed2, only one file is accepted, all sequences must be same length and match count/length of --paired-seed2; default: None)
- `--paired-seed2`: Paired seed sequence file 2 (used only with --paired and --paired-seed1, only one file is accepted, all sequences must be same length and match count/length of --paired-seed1; default: None)

- `--output-dir`: Output directory (default: test)
- `--output-prefix`: Output filename prefix (default: barcodes)

**Output Files**:
- `{prefix}.txt` or `{prefix}_paired1.txt` & `{prefix}_paired2.txt`: Generated barcodes
- `generate_barcodes_{timestamp}.log`: Detailed generation log

**Important Notes**:
- **Seed sequences are not validated**: If using seed files (paired or unpaired), run `barcadia validate` first to ensure they pass all filters
- In paired mode with seeds, both `--paired-seed1` and `--paired-seed2` must be provided and have the same count/length
- Seeds are preserved in the output files (paired or unpaired)

---

#### 2. `barcadia validate` - Barcode Validation

**Purpose**: Validate existing barcode lists against quality filters with support for variable-length sequences.

**Algorithm Overview**:
1. Load and parse input file(s) and report lengths distribution
2. Check if sequences fail both biological filters (GC content, homopolymer repeats)
3. **Computational complexity-optimized distance validation**:
   - **Method selection**: Automatically chooses optimal validation approach based on barcode set and distance filter characteristics
   - **Early stopping**: Terminates on first violation found (prevents unnecessary computation)
4. Generate detailed validation report with comprehensive violation details

**Key Features**:
- Supports multiple input files (automatically concatenated) with variable lengths
- Uses Hamming distance for equal-length sequences, Levenshtein for mixed lengths (for sequences passing biological filters)
- **Adaptive method selection for distance validation**:
  - Small barcode sets (<10K sequences): Sequential pairwise validation
  - Large barcode sets (≥10K sequences) with mixed lengths and/or min_distance > 4: Parallel pairwise validation (when multiple CPUs available)
  - Large equal-length barcode sets (≥10K sequences) with min_distance ≤ 4: Neighbor enumeration
  - Early stopping on first violation
  - Can be skipped when `--skip-distance` flag is on
- Comprehensive reporting with violation cases (for both biological and distance constraints)

**Basic Usage**:
```bash
# Validate a single file
barcadia validate --input test/barcodes.txt

# Validate multiple files (automatically concatenated)
barcadia validate --input file1.txt file2.txt file3.txt

# Skip distance validation entirely (biological filters only)
barcadia validate --input test/barcodes.txt --skip-distance
```

**Required Arguments**:
- `--input`: Input file(s) containing DNA barcodes (one per line)

**Optional Arguments**:
- `--gc-min`: Minimum GC content (default: 0.4)
- `--gc-max`: Maximum GC content (default: 0.6)
- `--homopolymer-max`: Maximum homopolymer repeat length (default: 2)
- `--min-distance`: Minimum edit distance between sequences (default: 3)
- `--skip-distance`: Skip distance validation entirely (default: off)
- `--cpus`: Number of CPUs for pairwise parallel validation (default: all available)
- `--output-dir`: Output directory for logs and reports (default: test)

**Output Files**:
- `validation_report_{timestamp}.txt`: Detailed validation report
- `validate_barcodes_{timestamp}.log`: Validation process log

---

### Shared Modules

#### `config_utils.py`
Configuration utilities and DNA encoding/decoding:
- DNA bases encoded as 0-3 integers (A=0, T=1, G=2, C=3) for enhanced efficiency
- Convert between DNA strings and integer arrays for optimized processing
- Logging setup and file reading utilities (ExistingSequenceSet class)

#### `filter_utils.py`
High-performance filtering algorithms with Numba JIT compilation:
- Simple validation on filter arguments (GC content, homopolymer, min-distance)
- GC content and homopolymer repeat filtering with early stopping
- Hamming distance for equal-length sequences with early stopping 
- Levenshtein distance for variable-length sequences with early stopping
- Adaptive method selection between pairwise and neighbor enumeration approaches
- Neighbor enumeration for efficient distance constraint checking

### Utility Scripts

#### 1. `generate_random_sequences.py` - Test Data Generation

**Purpose**: Generate random DNA sequences for testing validation scripts.

**Usage**:
```bash
python -m barcadia.tools.generate_random_sequences --count <num> --lengths <length1> [length2...] [--output <file>]
```

**Output**: Auto-generated filename in `test/` directory based on `count` and `length` if `--output` not specified.

---

#### 2. `memory_benchmark.py` - Performance Monitoring

**Purpose**: Monitor memory usage and performance of the main scripts.

**Usage**:
```bash
# General usage
python -m barcadia.tools.memory_benchmark [--mem-output-dir <dir>] <command> [args...]

# Examples
python -m barcadia.tools.memory_benchmark barcadia generate --args
python -m barcadia.tools.memory_benchmark barcadia validate --args
```

**Output**: Memory usage report with peak memory consumption and execution time. Log saved to specified directory (default: `test/`).

## Core APIs

For programmatic usage, the package provides clean APIs for barcode generation and validation:

```python
from barcadia import Filter, generate_barcodes_core, validate_barcodes_core

# Configure filtering criteria
filter_params = Filter(gc_min=0.4, gc_max=0.6, homopolymer_max=2, min_distance=3)

# Generate 1000 barcodes of length 12
barcodes = generate_barcodes_core(target_count=1000, length=12, filter_params=filter_params)

# Validate the generated barcodes
result = validate_barcodes_core(sequences=barcodes, filter_params=filter_params)
```

## Citation

If you use Barcadia in your research, please cite:

```bibtex
@software{barcadia2025,
  title={Barcadia: a high-performance, memory-efficient toolkit for fast generation and validation of large-scale NGS barcodes},
  author={Jiang, Danting},
  year={2025},
  date={2025-09-12},
  url={https://pypi.org/project/barcadia/},
  note={Code repository: https://github.com/djiang0825/NGS_barcode},
  version={3.3.0}
}
```

A preprint will be posted soon. Citation information will be updated with a DOI once available.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for detailed version history.