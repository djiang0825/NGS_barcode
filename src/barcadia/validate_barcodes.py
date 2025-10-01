#!/usr/bin/env python3
"""
validate_barcodes.py

Validate if provided lists of NGS barcodes satisfy all quality filters (variable length sequences supported).

Program Overview:

1. Load and parse input file(s) and report lengths distribution
2. Check if sequences fail biological filters (GC content, homopolymer checks) - reports ALL violations per sequence
3. For sequences passing biological filters, apply intelligent algorithm selection for distance validation (unless --skip-distance flag is enabled):
   3a. Method selection logic (based on sequences passing biological filters):
       - Small barcode sets (<10K sequences): Pairwise sequential
       - Large barcode sets (≥10K sequences) with mixed lengths and/or min_distance > 4: Pairwise parallel (when multiple CPUs available, otherwise sequential)
       - Large equal-length barcode sets (≥10K sequences) with min_distance ≤ 4: Neighbor enumeration
   3b. Distance calculation: Hamming distance for equal-length sequences, Levenshtein for mixed lengths
   3c. Progress logging during validation (every 10 chunks for pairwise parallel, every 10K sequences for neighbor enumeration)
   3d. Early stopping on first distance violation with detailed reporting
4. Generate comprehensive validation report with violation details

Input: list(s) of NGS barcodes (one per line as .txt). Multiple files supported, concatenated automatically.

Output: validation report (validation_report_{timestamp}.txt) and validate_barcodes_{timestamp}.log file

Optional arguments:
--gc-min: minimum GC content (default: 0.4)
--gc-max: maximum GC content (default: 0.6)
--homopolymer-max: maximum allowed homopolymer repeat length (default: 2)
--min-distance: minimum edit distance between barcodes (default: 3)
--skip-distance: skip distance validation entirely (default: off)
--output-dir: output directory for validation logs and reports (default: test)
--cpus: number of CPUs to use for pairwise parallel distance validation (default: all available)

Required arguments:
--input: input file(s) containing NGS barcodes (one per line)
"""

import argparse
import logging
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime

# Import utility functions
from .config_utils import ExistingSequenceSet, decode_sequence, setup_logging
from .filter_utils import Filter, calculate_distance, check_gc_content_int, check_homopolymer_int, generate_hamming_neighbors, select_distance_method


@dataclass
class ValidationResult:
    """Result of barcode validation containing all validation data"""
    overall_valid: bool
    total_sequences: int
    biological_passed: int
    biological_failed: int
    biological_violations: list[tuple[int, str, str]]  # (line_num, sequence, reason)
    distance_skipped: bool
    distance_violation: tuple[int, int, str, str, int] | None  # (line1, line2, seq1, seq2, distance)
    validation_method: str
    features_checked: int


def validate_biological_filters(seq_array, gc_min, gc_max, homopolymer_max):
    """Check if sequence passes all biological filters and return all violations"""
    violations = []

    # Check GC content
    if not check_gc_content_int(seq_array, gc_min, gc_max):
        violations.append("GC content outside range")

    # Check homopolymer runs
    if not check_homopolymer_int(seq_array, homopolymer_max):
        violations.append("Homopolymer run too long")

    if violations:
        return False, "; ".join(violations)
    else:
        return True, "Passes all filters"


def log_violation_details(sequences, i, j, violation):
    """Log distance violation and create violation details with DNA strings for reporting"""
    logging.info(f"Early stopping: Found distance violation between sequences {violation[0] + 1} and {violation[1] + 1} (distance={violation[2]})")
    seq1_str = decode_sequence(sequences[i])
    seq2_str = decode_sequence(sequences[j])
    return (violation[0] + 1, violation[1] + 1, seq1_str, seq2_str, violation[2])


def validate_distances_neighbor_enum(sequences, min_distance):
    """Validate using neighbor enumeration - much faster for appropriate cases"""
    # Build hash set of all sequences for O(1) lookup
    sequence_set = set(tuple(seq) for seq in sequences)
    total_sequences = len(sequences)

    # Check each sequence for violations
    for i, seq in enumerate(sequences):
        # Progress logging every 10K sequences
        if i % 10_000 == 0 and i > 0:
            logging.info(f"Progress: {i:,}/{total_sequences:,} sequences processed ({i / total_sequences * 100:.1f}%)")

        seq_array = list(seq)  # Make mutable copy for neighbor generation

        # Generate all neighbors within min_distance
        for neighbor in generate_hamming_neighbors(seq_array, min_distance - 1):
            if neighbor in sequence_set and neighbor != tuple(seq):
                # Found a violation - get the index of the violating sequence
                for j, other_seq in enumerate(sequences):
                    if j != i and tuple(other_seq) == neighbor:
                        # Calculate actual distance for reporting
                        actual_distance = sum(a != b for a, b in zip(seq, other_seq, strict=False))
                        violation = (i, j, actual_distance)
                        violation_details = log_violation_details(sequences, i, j, violation)
                        sequences_processed = i + 1  # Number of sequences processed when violation found
                        return True, sequences_processed, violation_details

    # No violations found - processed all sequences
    return False, total_sequences, None


def generate_pair_chunk(start_idx, chunk_size, n):
    """Generate a chunk of pairs lazily starting from start_idx"""
    pairs_generated = 0
    current_idx = 0

    for i in range(n):
        for j in range(i + 1, n):
            if current_idx >= start_idx:
                if pairs_generated >= chunk_size:
                    return
                yield (i, j)
                pairs_generated += 1
            current_idx += 1


def validate_chunk(sequences_chunk, min_distance):
    """Worker function for distance validation (can be used sequentially or in parallel)"""
    pairs_checked = 0
    n = len(sequences_chunk)

    # Use the lazy pair generator for memory efficiency
    for i, j in generate_pair_chunk(0, n * (n - 1) // 2, n):
        distance = calculate_distance(sequences_chunk[i], sequences_chunk[j], min_distance)
        pairs_checked += 1

        if distance < min_distance:
            violation = (i, j, distance)
            violation_details = log_violation_details(sequences_chunk, i, j, violation)
            return True, pairs_checked, violation_details

    return False, pairs_checked, None


def validate_distances(sequences, min_distance, method, cpus, chunk_size):
    """Unified distance validation with method selection and parallel/sequential execution"""
    # Execute the chosen method
    if method == "neighbor_enumeration":
        # Use neighbor enumeration when set up is optimal (no parallelization involved)
        return validate_distances_neighbor_enum(sequences, min_distance)
    elif method == "pairwise_sequential" or cpus == 1:
        # Use sequential for small barcode sets or single CPU
        return validate_chunk(sequences, min_distance)
    else:  # method == "pairwise" and cpus > 1
        # Large dataset with multiple CPUs - always use parallel with pre-calculated chunk size
        chunks = [sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)]

        with ProcessPoolExecutor(max_workers=cpus) as executor:
            futures = [executor.submit(validate_chunk, chunk, min_distance) for chunk in chunks]

            # Process results with early stopping
            total_pairs_checked = 0
            for future in futures:
                early_stopped, pairs_checked, violation_info = future.result()
                total_pairs_checked += pairs_checked

                if early_stopped:
                    # Cancel remaining futures for early stopping
                    for f in futures:
                        f.cancel()
                    return True, total_pairs_checked, violation_info

            return False, total_pairs_checked, None


def validate_barcodes_core(sequences, filter_params, has_mixed_lengths=False, skip_distance=False, cpus=None) -> ValidationResult:
    """Main function to validate input barcode sets against biological filters and distance constraints"""
    if cpus is None:
        cpus = mp.cpu_count()
    start_time = time.time()

    logging.info("Starting barcode validation...")
    logging.info(f"Filter 1 (within-sequence), GC content: {filter_params.gc_min:.1%} - {filter_params.gc_max:.1%}")
    logging.info(f"Filter 2 (within-sequence), Max homopolymer repeat: {filter_params.homopolymer_max}")
    logging.info(f"Filter 3 (between-sequence), Minimum edit distance: {filter_params.min_distance}")

    # 1. Validate sequences for biological filters
    valid_sequences = []
    biological_violations = []

    for i, seq_array in enumerate(sequences):
        # Check biological filters
        is_valid, reason = validate_biological_filters(seq_array, filter_params.gc_min, filter_params.gc_max, filter_params.homopolymer_max)

        if is_valid:
            valid_sequences.append(seq_array)
        else:
            # Convert back to DNA string for reporting
            dna_string = decode_sequence(seq_array)
            biological_violations.append((i + 1, dna_string, reason))

    logging.info("Biological filter (GC content and homopolymer repeats) results:")
    logging.info(f"  Passed: {len(valid_sequences)} sequences")
    logging.info(f"  Failed: {len(biological_violations)} sequences")

    # 2. Validate sequences for distance constraints
    # Calculate total distance pairs for sequences that passed biological filters
    n = len(valid_sequences)
    total_pairs = n * (n - 1) // 2

    # Check if we should skip distance validation
    distance_skipped = False
    logging.info("Distance filter results:")
    if skip_distance:
        logging.info("Skipping distance validation (--skip-distance flag enabled)")
        early_stopped = False
        features_checked = 0
        distance_skipped = True
        validation_method = "skipped"
        violation_info = None
        logging.info("  Distance validation skipped")
    # If not, continue with distance validation
    else:
        logging.info("Validating distances for sequences that passed biological filters...")
        # Determine method and calculate chunk size for pairwise method
        method = select_distance_method(n, filter_params.min_distance, has_mixed_lengths)

        # Calculate chunk size for pairwise method with multiple CPUs
        chunk_size = None
        if method == "pairwise_sequential":
            logging.info("Using sequential pairwise distance checking (small barcode set for sequences that passed biological filters)")
        elif method == "pairwise" and cpus == 1:
            logging.info("Using sequential pairwise distance checking (1 CPU)")
        elif method == "pairwise" and cpus > 1:
            chunk_size = max(100000, total_pairs // (cpus * 10))
            logging.info(f"Using parallel pairwise distance checking (chunk size: {chunk_size})")

        # Execute validation
        early_stopped, features_checked, violation_info = validate_distances(valid_sequences, filter_params.min_distance, method, cpus, chunk_size)

        # Log results
        if method == "neighbor_enumeration":
            logging.info(f"  Total sequences (that passed biological filters): {n}")
            logging.info(f"  Sequences processed: {features_checked}")
        else:
            logging.info(f"  Total sequence pairs: {total_pairs:,} (sequences that passed biological filters)")
            logging.info(f"  Pairs checked: {features_checked:,}")

        validation_method = method

    overall_valid = len(valid_sequences) == len(sequences) and not early_stopped

    duration = time.time() - start_time

    logging.info("Validation complete!")
    logging.info(f"Overall validation: {'PASSED' if overall_valid else 'FAILED'}")
    logging.info(f"Total time: {duration:.2f} seconds")

    # Return structured result
    return ValidationResult(
        overall_valid=overall_valid,
        total_sequences=len(sequences),
        biological_passed=len(valid_sequences),
        biological_failed=len(biological_violations),
        biological_violations=biological_violations,
        distance_skipped=distance_skipped,
        distance_violation=violation_info,
        validation_method=validation_method,
        features_checked=features_checked
    )


def write_validation_report(result: ValidationResult, filter_params: Filter, args, log_filepath: str):
    """Write validation report to file"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = os.path.join(args.output_dir, f"validation_report_{timestamp}.txt")

    with open(report_file, "w") as f:
        f.write("Barcode Validation Report\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"Input file: {args.input}\n")
        f.write(f"Total sequences: {result.total_sequences}\n\n")
        f.write("Filter Settings:\n")
        f.write(f"  GC content: {filter_params.gc_min:.1%} - {filter_params.gc_max:.1%}\n")
        f.write(f"  Max homopolymer: {filter_params.homopolymer_max}\n")
        f.write(f"  Minimum distance: {filter_params.min_distance}\n\n")
        f.write(f"Biological filter passed: {result.biological_passed}\n")
        f.write(f"Biological filter failed: {result.biological_failed}\n")

        if result.distance_skipped:
            f.write("Distance validation: SKIPPED (--skip-distance flag enabled)\n")
        elif result.distance_violation:
            f.write("Distance validation: EARLY STOPPED (found first violation)\n")
            f.write(f"  Method used: {result.validation_method}\n")
            if result.validation_method == "neighbor_enumeration":
                f.write(f"  Sequences processed before stopping: {result.features_checked:,}\n")
            else:
                f.write(f"  Pairs checked before stopping: {result.features_checked:,}\n")
        else:
            f.write("Distance validation: PASSED (no violations found)\n")
            f.write(f"  Method used: {result.validation_method}\n")
            if result.validation_method == "neighbor_enumeration":
                f.write(f"  Total sequences processed: {result.features_checked:,}\n")
            else:
                f.write(f"  Total pairs checked: {result.features_checked:,}\n")
        f.write("\n")

        if result.biological_violations:
            f.write("Biological Filter (GC content and homopolymer) Violations:\n")
            f.write("-" * 30 + "\n")
            for line_num, seq, reason in result.biological_violations:
                f.write(f"Line {line_num}: {seq} - {reason}\n")
            f.write("\n")

        # Add distance violation details if available
        if result.distance_violation is not None:
            f.write("Distance Violations:\n")
            f.write("-" * 19 + "\n")
            seq1_line, seq2_line, seq1_str, seq2_str, distance = result.distance_violation
            f.write(f"Line {seq1_line}: {seq1_str} and Line {seq2_line}: {seq2_str} - distance {distance} (minimum required: {filter_params.min_distance})\n")
        f.write("\n")

    # Log file locations
    if log_filepath:
        logging.info(f"Log file: {log_filepath}")
    logging.info(f"Report file: {report_file}")
    
    # Print result
    if result.overall_valid:
        print("All barcodes are valid!")
    else:
        print("VALIDATION FAILED!")


def setup_argument_parser():
    """Setup and return the argument parser for barcode validation"""
    parser = argparse.ArgumentParser(
        description="Validate DNA barcodes against quality filters (GC content, homopolymer repeats, minimum distance)", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument("--input", type=str, required=True, nargs="+", help="Input file(s) containing DNA barcodes (one per line)")

    # Output arguments
    parser.add_argument("--output-dir", type=str, default="test", help="Output directory for validation logs and reports")

    # Filter arguments with defaults
    parser.add_argument("--gc-min", type=float, default=0.4, help="Minimum GC content (as fraction, e.g., 0.4 = 40%%)")
    parser.add_argument("--gc-max", type=float, default=0.6, help="Maximum GC content (as fraction, e.g., 0.6 = 60%%)")
    parser.add_argument("--homopolymer-max", type=int, default=2, help="Maximum allowed homopolymer repeat length")
    parser.add_argument("--min-distance", type=int, default=3, help="Minimum edit distance between sequences")

    # Performance arguments
    parser.add_argument("--cpus", type=int, default=mp.cpu_count(), help="Number of CPU cores to use for parallel distance validation")

    # Mode arguments
    parser.add_argument("--skip-distance", action="store_true", help="Skip distance validation entirely")

    return parser


def validate_validator_arguments(args, length_counts):
    """Validate validator-specific arguments (length, distance, homopolymer) and return has_mixed_lengths flag"""
    input_length = max(length_counts.keys())

    # Homopolymer repeat x max input length validation
    if args.homopolymer_max >= input_length:
        raise ValueError(f"Maximum homopolymer repeat length must be < max input length ({input_length}bp)")

    # Minimum distance x max input length validation
    if args.min_distance >= input_length:
        raise ValueError(f"Minimum distance must be < max input length ({input_length}bp)")

    # Check for mixed lengths
    has_mixed_lengths = len(length_counts) > 1

    return has_mixed_lengths


def main(argv=None):
    parser = setup_argument_parser()
    args = parser.parse_args(argv)
    log_filepath = setup_logging(args, "validate_barcodes")
    
    # Validate filter parameters immediately after parsing arguments
    filter_params = Filter(
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance
    )

    # Load input files using ExistingSequenceSet
    sequence_set = ExistingSequenceSet.from_input_files(args.input)

    # Validate validator-specific arguments and get mixed lengths flag
    has_mixed_lengths = validate_validator_arguments(args, sequence_set.length_counts)

    result = validate_barcodes_core(
        sequences=sequence_set.sequences,
        filter_params=filter_params,
        has_mixed_lengths=has_mixed_lengths,
        skip_distance=args.skip_distance,
        cpus=args.cpus
    )
    
    # Write report
    write_validation_report(result, filter_params, args, log_filepath)


if __name__ == "__main__":
    main(sys.argv[1:])
