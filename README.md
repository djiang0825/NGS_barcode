# NGS Barcode Toolkit (v1)

A high-performance toolkit for generating and validating large-scale DNA barcodes for Next-Generation Sequencing (NGS) applications. 

## Features

- **Efficient Generation**: Iterative growth algorithm with early stopping and parallel processing
- **Biological Filters**: GC content and homopolymer length constraints
- **Distance Validation**: Hamming/Levenshtein distance checking with configurable thresholds
- **Scalable Design**: Handles both small and large datasets with memory optimization
- **Performance Monitoring**: Built-in memory benchmarking and detailed logging
- **Flexible Input**: Support for seed sequences and variable-length barcodes

## Background

### Problem Statement

Next-generation sequencing (NGS) is a high-throughput method that enables millions of DNA fragments to be sequenced in parallel, serving as the core technology for decoding genomes across all living organisms. In this process, researchers use DNA barcodes to uniquely label and track individual biomolecules. Common examples include multiplex sample indexes and unique molecular identifiers (UMIs). 

The practical utility of a DNA barcode library depends on controlling key features: **GC content** (the percentage of G and C nucleotides), **homopolymer repeats** (the length of the longest stretch of identical nucleotides), and **edit distance** (the measure of dissimilarity among sequences). GC content and homopolymer repeats are within-sequence criteria that ensure molecular stability during sequencing and are computationally inexpensive to evaluate. In contrast, edit distance is a between-sequence constraint that provides error tolerance during analysis, but is computationally demanding to assess for large datasets.

### Existing Methods and their Limitations

For a set of n sequences, the number of pairwise distance comparisons grows quadratically:

$$\binom{n}{2} = \frac{n(n-1)}{2}$$

This makes brute-force approaches computationally prohibitive for large datasets. Apart from the approach described by [Lyons et al. (2017)](https://www.nature.com/articles/s41598-017-12825-2), existing tools are not built to accommodate large-scale barcode library design (≥100K sequences). Lyons et al. circumvented computational limitations using a probabilistic Markov chain model; however, their resulting barcode sets are no longer accessible (invalid URLs in the paper), and the underlying code was not released (likely due to proprietary restrictions). 

### This Toolkit and its Advantages

Here, I introduce a toolkit for large-scale NGS barcode generation that utilizes modern computational optimization methods and incorporates novel algorithms I have developed. The toolkit is openly available to promote reproducibility and is specifically designed to run efficiently on minimal computing resources (e.g., standard laptops) for broad accessibility.

In comparison with [TagGD](https://doi.org/10.1371/journal.pone.0057521)—the only other method reported to support barcode generation at the scale of up to 100,000 sequences—this toolkit generated 20,000 18-bp barcodes in just 1 minute using only 100 MB of RAM on a comparable 8-core laptop, as opposed to the 5 minutes reported in the TagGD abstract. Notably, for 100,000 18-bp barcodes, this toolkit completed the process in 15 minutes with 175 MB of RAM, which is significantly faster than the 1.5 hours reported in the TagGD paper (Table 1). Additionally, this toolkit can handle the generation of million-scale barcodes within reasonable time (i.e., hours) on standard compute setups (e.g., laptops), far exceeding the capacity of TagGD and other existing tools. More detailed benchmarking results are presented in a later section of this document. 

Additionally, this toolkit offers unique features not available in other tools: paired barcode generation for dual-indexing applications, extension from user-provided seed sequences, and a comprehensive validation script for quality assessment of existing barcode sets.

Together, these algorithmic advances, performance improvements, and extensive features make the toolkit ideally suited for researchers who require large-scale, high-quality barcode libraries for state-of-the-art NGS applications.

## Default Filter Parameters

This toolkit uses carefully chosen default parameters that balance biological stability with computational efficiency (users can configure these based on their specific needs):

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

Leveraging established results from coding theory, one can calculate lower and upper bounds on the number of valid barcode sequences under specified edit distance constraints, assuming equal-length barcodes and applying the Hamming distance metric.

### Gilbert-Varshamov Bound (Lower Bound)

The Gilbert-Varshamov bound provides a lower bound guarantee that a code of at least this size exists. For DNA sequences, it states:

$$M \geq \frac{4^n}{V(n,d-1)}$$

where `V(n,d-1)` is the volume of a Hamming sphere of radius `d-1`:

$$V(n,d-1) = \sum_{i=0}^{d-1} \binom{n}{i} \cdot 3^i$$

### Hamming Bound (Upper Bound)

The Hamming bound provides the theoretical maximum number of codewords that can exist for a given sequence length and minimum distance constraint. For DNA sequences of length `n` with minimum Hamming distance `d`, the bound is:

$$M \leq \frac{4^n}{V(n,t)}$$

where `t = ⌊(d-1)/2⌋` and `V(n,t)` is the volume of a Hamming sphere of radius `t`:

$$V(n,t) = \sum_{i=0}^{t} \binom{n}{i} \cdot 3^i$$

As one exhausts the sequence space in search of valid barcodes that satisfy the distance constraint for all members of the set, the closer the target size approaches the theoretical upper bound, the longer the search will take. In particular, a significant—often exponential—slowdown can be expected once the target size surpasses the Gilbert-Varshamov (GV) bound.

For typical barcode applications (6-16 bp, as longer sequences may introduce molecular complexity) using the default minimum distance of 3, the theoretical bounds are summarized below:

| Length (bp) | GV Bound (Lower) | Hamming Bound (Upper) |
|-------------|------------------|-----------------------|
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

The bounds shown above consider only distance constraints. In practice, additional biological filters (GC content and homopolymer restrictions) further reduce the achievable library sizes. It is possible to estimate realistic bounds that satisfy all three constraints using simulation and sampling methods. A future version of this toolkit will include functionality to calculate these practical bounds for any given parameter set (stay tuned!).

## Benchmarking Highlights

Performance benchmarks on a MacBook Pro 2019 (16-core, 16GB RAM) using default parameters with target sizes near the Gilbert-Varshamov bound:

| Length (bp) | Target Size | Time (Peak Memory) |
|-------------|-------------|-------------------|
| 6           | 10          | 0.3 seconds       |
| 8           | 100         | 0.4 seconds       |
| 10          | 1,000       | 0.7 seconds       |
| 12          | 10,000      | 28 seconds (90 MB) |
| 14          | 100,000     | 16 minutes (170 MB) |
| 16          | 1,000,000   | 10 hours (5 GB)    |

These results demonstrate excellent scalability for small to medium libraries, with generation times remaining practical even for million-scale datasets. Memory usage scales efficiently, remaining under 200 MB for libraries up to 100,000 sequences.

## Installation

### Requirements

- Python 3.12+ (tested with Python 3.12)
- Required packages (install via `pip install -r requirement.txt`):
  - numpy==2.2.6
  - numba==0.61.2 (for JIT compilation acceleration)
  - llvmlite==0.44.0
  - psutil==7.0.0

### Setup

```bash
git clone https://github.com/djiang0825/NGS_barcode.git
cd NGS_barcode
pip install -r requirement.txt
```

## Project Structure

```
Barcode/
├── src/
│   ├── generate_barcodes.py                     # Main barcode generation script
│   ├── validate_barcodes.py                     # Barcode validation script
│   └── utils/
│       ├── dna_utils.py                         # DNA encoding/decoding utilities
│       ├── filter_utils.py                      # Core filtering algorithms (Numba JIT)
│       ├── generate_random_sequences.py         # Random sequence generator for validation testing
│       └── memory_benchmark.py                  # Performance monitoring utility
└── requirement.txt                              # Python dependencies
```

## Scripts Overview

### Main Scripts (User-Facing)

#### 1. `generate_barcodes.py` - Barcode Generation

**Purpose**: Generate diverse NGS barcodes from scratch or extend existing seed sequences using an iterative growth algorithm (paired mode supported).

**Algorithm Overview**:
1. Generate random sequences passing biological filters (GC content, homopolymer length)
2. Filter candidates for minimum distance constraints vs existing pool (parallel)
3. Filter candidates for distance constraints within batch (sequential)
4. Add verified batch to pool and repeat until target count reached

**Key Features**:
- Guarantees no duplicate sequences
- Satisfies all constraints in final pool
- Supports paired barcode generation (when --paired flag is on, will generate 2x the amount of target counts and split output into 2 files)
- Can build from seed sequence files (when --seeds flag is on) and accomodate difference in lengths between seed and new sequences (new sequences are all the same length)
- Parallel processing with configurable CPU usage and intelligent chunk sizing
- Memory-efficient integer encoding

**Basic Usage**:
```bash
# Generate 1000 barcodes of length 12
python src/generate_barcodes.py --count 1000 --length 12

# Generate paired barcodes
python src/generate_barcodes.py --count 1000 --length 12 --paired

# Build from seed sequences
python src/generate_barcodes.py --count 1000 --length 12 \
    --seeds seed_file1.txt seed_file2.txt
```

**Required Arguments**:
- `--count`: Number of barcodes to generate
- `--length`: Length of each barcode sequence

**Optional Arguments**:
- `--gc-min`: Minimum GC content (default: 0.4)
- `--gc-max`: Maximum GC content (default: 0.6)
- `--homopolymer-max`: Maximum homopolymer length (default: 2)
- `--min-distance`: Minimum Hamming distance between sequences (default: 3)
- `--cpus`: Number of CPU cores to use (default: all available)
- `--paired`: Generate paired barcodes (doubles count, creates two files)
- `--seeds`: Seed sequence files (one sequence per line)
- `--output-dir`: Output directory (default: test)
- `--output-prefix`: Output filename prefix (default: barcodes)

**Output Files**:
- `{prefix}.txt` or `{prefix}_paired1.txt` & `{prefix}_paired2.txt`: Generated barcodes
- `generate_barcode_{timestamp}.log`: Detailed generation log

**Functions to add**:
- More elegant handling of paired output generation when paired seeds are detected
- Validation of seed sequences 
- Progress saving for really large target (>100K)

---

#### 2. `validate_barcodes.py` - Barcode Validation

**Purpose**: Validate existing barcode lists against quality filters with support for variable-length sequences.

**Validation Process**:
1. Load and parse input files (report lengths distribution)
2. Apply biological filters (GC content, homopolymer checks)
3. Validate distance constraints with early stopping on first violation
4. Generate detailed validation report

**Key Features**:
- Supports multiple input files (automatically concatenated) with variable lengths
- Uses Hamming distance for equal-length sequences, Levenshtein for mixed lengths
- Parallel processing for large datasets (>10,000 sequences)
- Early stopping of distance validation for efficiency
- Comprehensive reporting

**Basic Usage**:
```bash
# Validate a single file
python src/validate_barcodes.py --input barcodes.txt

# Validate multiple files
python src/validate_barcodes.py --input file1.txt file2.txt file3.txt

# Skip distance checking if biological filters fail
python src/validate_barcodes.py --input barcodes.txt --skip-distance
```

**Required Arguments**:
- `--input`: Input file(s) containing DNA barcodes (one per line)

**Optional Arguments**:
- `--gc-min`: Minimum GC content (default: 0.4)
- `--gc-max`: Maximum GC content (default: 0.6)
- `--homopolymer-max`: Maximum homopolymer length (default: 2)
- `--min-distance`: Minimum distance between sequences (default: 3)
- `--skip-distance`: Skip distance validation if biological filters fail (default: false)
- `--cpus`: Number of CPUs for parallel validation (default: all available)
- `--output-dir`: Output directory for logs and reports (default: test)

**Output Files**:
- `validation_report_{timestamp}.txt`: Detailed validation report
- `validate_barcode_{timestamp}.log`: Validation process log

**Functions to add**:
- Progress saving for really large input (i.e., when parallel processing is triggered)

---

### Utility Scripts

#### 3. `generate_random_sequences.py` - Test Data Generation

**Purpose**: Generate random DNA sequences for testing validation scripts.

**Usage**:
```bash
python src/utils/generate_random_sequences.py --count <num> --lengths <length1> [length2...] [--output <file>]
```

**Output**: Auto-generated filename in `test/` directory based on `count` and `length` if `--output` not specified.

---

#### 4. `memory_benchmark.py` - Performance Monitoring

**Purpose**: Monitor memory usage and performance of the main scripts.

**Usage**:
```bash
# General usage
python src/utils/memory_benchmark.py [--output-dir <dir>] <script_path> [script_args...]

# Examples
python src/utils/memory_benchmark.py src/generate_barcodes.py --count 1000 --length 12
python src/utils/memory_benchmark.py src/validate_barcodes.py --input barcodes.txt
```

**Output**: Memory usage report with peak memory consumption and execution time. Log saved to specified directory (default: `test/`).

## Additional Utilities

### Core Library Modules

#### `dna_utils.py`
DNA encoding/decoding utilities and validation functions:
- DNA bases encoded as 0-3 integers (A=0, T=1, G=2, C=3) for enhanced efficiency
- Validate command-line arguments and check feasibility
- Calculate and check for theoretical Hamming bounds for given length and min-distance constraints

#### `filter_utils.py`
High-performance filtering algorithms with Numba JIT compilation:
- Numba JIT compilation for speed
- GC content and homopolymer repeat filtering
- Hamming distance for equal-length sequences with early stopping 
- Levenshtein distance for variable-length sequences 
