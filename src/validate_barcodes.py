#!/usr/bin/env python3
"""
validate_barcodes.py

Validate if provided lists of NGS barcodes satisfy all quality filters (variable length sequences supported).

Algorithm overview:
1. Load and parse input file(s) and report lengths distribution
2. Check if sequences fail biological filters (GC content, homopolymer checks) - reports ALL violations per sequence
3. For sequences passing biological filters, apply intelligent algorithm selection for distance validation (unless --skip-distance flag is enabled):
   3a. Method selection logic (based on sequences passing biological filters):
       - Small datasets (<10K sequences): Sequential pairwise validation
       - Large mixed-length datasets: Parallel pairwise validation (neighbor enumeration complex for Levenshtein)
       - Large equal-length datasets: Efficiency-based choice between neighbor enumeration and pairwise
         * neighbor_enumeration if: (sequences × neighbors_per_seq) < (pairwise_operations × 0.8)
         * pairwise_parallel otherwise
   3b. Distance calculation: Hamming distance for equal-length sequences, Levenshtein for mixed lengths
   3c. Progress logging during validation (every 10 chunks for parallel, every 10K sequences for neighbor enumeration)
   3d. Early stopping on first distance violation with detailed reporting
4. Generate comprehensive validation report with violation details 

Input: list(s) of NGS barcodes (one per line as .txt). Multiple files supported, concatenated automatically.

Output: validation report (.txt) and .log file

Optional arguments:
--gc-min: minimum GC content (default: 0.4)
--gc-max: maximum GC content (default: 0.6)
--homopolymer-max: maximum allowed homopolymer length (default: 2)
--min-distance: minimum Hamming distance between sequences (default: 3)
--skip-distance: skip distance validation entirely (default: off)
--output-dir: output directory for validation logs and reports (default: test)
--cpus: number of CPUs to use for parallel distance validation (default: all available)

Required arguments:
--input: input file(s) containing NGS barcodes (one per line)
"""

import argparse
import logging
import os
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

# Import utility functions
from utils.dna_utils import DNA_BASES, encode_sequence, decode_sequence, validate_arguments
from utils.filter_utils import check_gc_content_int, check_homopolymer_int, calculate_distance, calculate_neighbor_count, generate_hamming_neighbors

def validate_sequence_biological(seq_array, gc_min, gc_max, homopolymer_max):
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

def choose_validation_method(n_sequences, sequences, min_distance):
    """Choose optimal validation method based on dataset characteristics"""
    # Calculate sequence lengths
    lengths = [len(seq) for seq in sequences]
    min_length = min(lengths) if lengths else 0
    max_length = max(lengths) if lengths else 0
    
    # Small datasets: always use pairwise
    if n_sequences < 10000:
        return "pairwise_sequential"
    
    # Large datasets: calculate crossover point
    pairwise_ops = n_sequences * (n_sequences - 1) // 2
    
    # For mixed lengths, neighbor enumeration is more complex, so bias toward pairwise
    if min_length != max_length:
        logging.info(f"Mixed sequence lengths detected ({min_length}-{max_length} bp), using pairwise validation")
        if pairwise_ops > 1_000_000_000:  # > 1B operations (100K sequences with 4.95B operations processed by pairwise_parallel took 10 hours!)
            logging.warning(f"⚠️  Large mixed-length dataset warning: {pairwise_ops:,} pairwise comparisons - this may take a long time if there's no violation in the dataset")
        return "pairwise_parallel"
    
    # Equal lengths: calculate neighbor enumeration efficiency
    # Use the actual sequence length since all sequences are equal length
    actual_length = len(sequences[0]) if sequences else 0
    neighbors_per_seq = calculate_neighbor_count(actual_length, min_distance - 1)
    neighbor_ops = n_sequences * neighbors_per_seq
    
    if neighbor_ops < pairwise_ops * 0.8:  # 20% buffer for neighbor enumeration overhead
        logging.info(f"Using neighbor enumeration: {neighbor_ops:,} operations vs {pairwise_ops:,} pairwise")
        return "neighbor_enumeration"
    else:
        logging.info(f"Using pairwise validation: neighbor enumeration would require {neighbor_ops:,} operations")
        if pairwise_ops > 1_000_000_000:  # > 1B operations (100K sequences with 4.95B operations processed by pairwise_parallel took 10 hours!)
            logging.warning(f"⚠️  Large equal-length dataset warning (pairwise selected over neighbor enumeration): {pairwise_ops:,} pairwise comparisons - this may take a long time if there's no violation in the dataset")
        return "pairwise_parallel"

def check_distance_violation(sequences, i, j, min_distance):
    """Check if a pair of sequences violates distance constraint"""
    distance = calculate_distance(sequences[i], sequences[j], min_distance)
    if distance < min_distance:
        return (i, j, distance)
    return None

def log_distance_violation(violation):
    """Log distance violation with consistent formatting"""
    logging.info(f"Early stopping: Found distance violation between sequences {violation[0]+1} and {violation[1]+1} (distance={violation[2]})")

def validate_distance_constraints(sequences, min_distance):
    """Check if all sequences meet minimum distance requirements - exits early on first violation"""
    pairs_checked = 0
    total_pairs = len(sequences) * (len(sequences) - 1) // 2
    
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            pairs_checked += 1
            
            # No progress logging for sequential validation (small datasets are fast)
            
            violation = check_distance_violation(sequences, i, j, min_distance)
            if violation is not None:
                # Early exit - we found a violation, no need to continue
                log_distance_violation(violation)
                # Return violation details: (seq1_line, seq2_line, seq1_str, seq2_str, distance)
                seq1_str = decode_sequence(sequences[i])
                seq2_str = decode_sequence(sequences[j])
                violation_details = (violation[0]+1, violation[1]+1, seq1_str, seq2_str, violation[2])
                return True, pairs_checked, violation_details
    
    return False, pairs_checked, None

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

def check_pair_chunk_worker(args):
    """Worker function for parallel distance validation"""
    start_idx, chunk_size, n, sequences, min_distance = args
    
    pairs_checked = 0
    
    # Generate pairs lazily for this chunk
    for i, j in generate_pair_chunk(start_idx, chunk_size, n):
        violation = check_distance_violation(sequences, i, j, min_distance)
        pairs_checked += 1
        
        if violation is not None:
            # Return violation with sequence details
            seq1_str = decode_sequence(sequences[i])
            seq2_str = decode_sequence(sequences[j])
            violation_details = (violation[0]+1, violation[1]+1, seq1_str, seq2_str, violation[2])
            return violation_details, pairs_checked
    
    return None, pairs_checked

def validate_distance_constraints_parallel(sequences, min_distance, cpus=None):
    """Parallel version of distance validation with early stopping (memory-efficient)"""
    n = len(sequences)
    cpus = cpus or mp.cpu_count()
    
    # Calculate chunk parameters
    total_pairs = n * (n - 1) // 2
    chunk_size = max(100000, total_pairs // (cpus * 10))
    
    # Generate chunk start indices (no pair data stored in memory)
    chunk_starts = list(range(0, total_pairs, chunk_size))
    
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        # Submit workers with lazy chunk parameters (not actual pair data)
        futures = [executor.submit(check_pair_chunk_worker, (start_idx, chunk_size, n, sequences, min_distance)) 
                  for start_idx in chunk_starts]
        
        # Process results with early stopping
        total_pairs_checked = 0
        chunk_count = 0
        
        for future in futures:
            violation_details, pairs_checked = future.result()
            total_pairs_checked += pairs_checked
            chunk_count += 1
            
            # Progress logging every 10 chunks for parallel validation
            if chunk_count % 10 == 0:
                logging.info(f"Progress: {total_pairs_checked:,}/{total_pairs:,} pairs checked "
                           f"({total_pairs_checked/total_pairs*100:.1f}%) - {chunk_count} chunks completed")
            
            if violation_details is not None:
                # Log the basic violation (keeping existing format)
                basic_violation = (violation_details[0]-1, violation_details[1]-1, violation_details[4])
                log_distance_violation(basic_violation)
                # Cancel remaining futures for early stopping
                for f in futures:
                    f.cancel()
                return True, total_pairs_checked, violation_details
        
        # If no violations found, we checked all pairs
        return False, total_pairs_checked, None

def validate_distance_constraints_neighbor_enumeration(sequences, min_distance):
    """Validate using neighbor enumeration - much faster for appropriate cases"""
    # Build hash set of all sequences for O(1) lookup
    sequence_set = set(tuple(seq) for seq in sequences)
    total_sequences = len(sequences)
    
    # Check each sequence for violations
    for i, seq in enumerate(sequences):
        # Progress logging every 10K sequences
        if i % 10_000 == 0 and i > 0:
            logging.info(f"Progress: {i:,}/{total_sequences:,} sequences processed "
                       f"({i/total_sequences*100:.1f}%)")
        
        seq_array = list(seq)  # Make mutable copy for neighbor generation
        
        # Generate all neighbors within min_distance
        for neighbor in generate_hamming_neighbors(seq_array, min_distance - 1):
            if neighbor in sequence_set and neighbor != tuple(seq):
                # Found a violation - get the index of the violating sequence
                for j, other_seq in enumerate(sequences):
                    if j != i and tuple(other_seq) == neighbor:
                        # Calculate actual distance for reporting
                        actual_distance = sum(a != b for a, b in zip(seq, other_seq))
                        violation = (i, j, actual_distance)
                        log_distance_violation(violation)
                        # Return violation details: (seq1_line, seq2_line, seq1_str, seq2_str, distance)
                        seq1_str = decode_sequence(seq)
                        seq2_str = decode_sequence(other_seq)
                        violation_details = (i+1, j+1, seq1_str, seq2_str, actual_distance)
                        sequences_processed = i + 1  # Number of sequences processed when violation found
                        return True, sequences_processed, violation_details
    
    # No violations found - processed all sequences
    return False, total_sequences, None

def validate_barcodes(input_files, gc_min, gc_max, homopolymer_max, min_distance, skip_distance=False, cpus=None):
    """Main validation function"""
    logging.info(f"Starting barcode validation...")
    logging.info(f"Input files: {input_files}")
    logging.info(f"Filter 1 (within-sequence), GC content: {gc_min:.1%} - {gc_max:.1%}")
    logging.info(f"Filter 2 (within-sequence), Max homopolymer: {homopolymer_max}")
    logging.info(f"Filter 3 (between-sequence), Minimum distance: {min_distance}")
    
    # Read sequences from all input files
    sequences = []
    line_num = 1
    
    # Convert single file to tuple for consistent handling
    if isinstance(input_files, str):
        input_files = (input_files,)
    
    for input_file in input_files:
        with open(input_file, 'r') as f:
            for line in f:
                seq = line.strip()
                if not seq:  # Skip empty lines
                    continue
                
                # Basic validation
                if not all(base in DNA_BASES for base in seq):
                    logging.error(f"Line {line_num}: Invalid DNA sequence '{seq}'")
                    continue
                
                sequences.append((line_num, seq))
                line_num += 1
    
    # Check sequence lengths
    lengths = [len(seq) for _, seq in sequences]
    length_counts = {}
    for length in lengths:
        length_counts[length] = length_counts.get(length, 0) + 1
    
    if len(length_counts) == 1:
        length_info = f"length {list(length_counts.keys())[0]}"
    else:
        length_breakdown = ", ".join([f"{count} at length {length}" for length, count in sorted(length_counts.items())])
        length_info = f"mixed lengths: {length_breakdown}"
    
    logging.info(f"Loaded {len(sequences)} sequences from file ({length_info})")
    
    # Convert to integer arrays for efficient processing
    seq_arrays = []
    valid_sequences = []
    biological_violations = []
    
    for line_num, seq in sequences:
        seq_array = encode_sequence(seq)
        
        # Check biological filters
        is_valid, reason = validate_sequence_biological(seq_array, gc_min, gc_max, homopolymer_max)
        
        if is_valid:
            seq_arrays.append(seq_array)
            valid_sequences.append((line_num, seq))
        else:
            biological_violations.append((line_num, seq, reason))
    
    logging.info(f"Biological filter results:")
    logging.info(f"  Passed: {len(valid_sequences)} sequences")
    logging.info(f"  Failed: {len(biological_violations)} sequences")
    
    # Calculate total pairs for sequences that passed biological filters
    n = len(seq_arrays)
    total_pairs = n * (n - 1) // 2
    
    # Check if we should skip distance validation
    distance_skipped = False
    if skip_distance:
        logging.info(f"Skipping distance validation (--skip-distance flag enabled)")
        early_stopped = False
        pairs_checked = 0
        distance_skipped = True
        validation_method = "skipped"
        violation_info = None
    else:
        # Choose optimal validation method based on dataset characteristics
        validation_method = choose_validation_method(n, seq_arrays, min_distance)
        logging.info(f"Distance validation method: {validation_method}")
        
        # Execute chosen validation method
        if validation_method == "neighbor_enumeration":
            early_stopped, pairs_checked, violation_info = validate_distance_constraints_neighbor_enumeration(seq_arrays, min_distance)
        elif validation_method == "pairwise_sequential":
            early_stopped, pairs_checked, violation_info = validate_distance_constraints(seq_arrays, min_distance)
        else:  # pairwise_parallel
            early_stopped, pairs_checked, violation_info = validate_distance_constraints_parallel(seq_arrays, min_distance, cpus)
    
    logging.info(f"Distance constraint results:")
    if validation_method == "neighbor_enumeration":
        logging.info(f"  Total sequences: {len(seq_arrays)}")
        logging.info(f"  Sequences processed: {pairs_checked}")
    elif validation_method == "skipped":
        logging.info(f"  Distance validation skipped")
    else:
        logging.info(f"  Total sequence pairs: {total_pairs:,} (sequences that passed biological filters)")
        logging.info(f"  Pairs checked: {pairs_checked:,}")
    
    # Summary
    total_violations = len(biological_violations) + (1 if early_stopped else 0)
    overall_valid = len(valid_sequences) == len(sequences) and not early_stopped
    
    # Check if early stopping occurred
    if validation_method == "neighbor_enumeration":
        # For neighbor enumeration, early stopping means we didn't process all sequences
        distance_early_stopped = early_stopped and pairs_checked < len(seq_arrays)
    else:
        # For pairwise methods, early stopping means we didn't check all pairs
        distance_early_stopped = early_stopped and pairs_checked < total_pairs
    
    logging.info(f"Validation complete!")
    logging.info(f"Overall validation: {'PASSED' if overall_valid else 'FAILED'}")
    
    return overall_valid, total_violations, biological_violations, valid_sequences, sequences, distance_early_stopped, distance_skipped, violation_info, validation_method, pairs_checked

def generate_validation_report(input_file, gc_min, gc_max, homopolymer_max, min_distance, 
                             biological_violations, valid_sequences, 
                             total_sequences, output_dir, early_stopped=False, distance_skipped=False, 
                             validation_method="unknown", pairs_checked=0, violation_info=None):
    """Generate detailed validation report"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = os.path.join(output_dir, f"validation_report_{timestamp}.txt")
    
    with open(report_file, 'w') as f:
        f.write("Barcode Validation Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Input file: {input_file}\n")
        f.write(f"Total sequences: {total_sequences}\n")
        f.write(f"Biological filter passed: {len(valid_sequences)}\n")
        f.write(f"Biological filter failed: {len(biological_violations)}\n")
        
        if distance_skipped:
            f.write(f"Distance validation: SKIPPED (--skip-distance flag enabled)\n")
        elif early_stopped:
            f.write(f"Distance validation: EARLY STOPPED (found first violation)\n")
            f.write(f"  Method used: {validation_method}\n")
            if validation_method == "neighbor_enumeration":
                f.write(f"  Sequences processed before stopping: {pairs_checked:,}\n")
            else:
                f.write(f"  Pairs checked before stopping: {pairs_checked:,}\n")
        else:
            f.write(f"Distance validation: PASSED (no violations found)\n")
            f.write(f"  Method used: {validation_method}\n")
            if validation_method == "neighbor_enumeration":
                f.write(f"  Total sequences processed: {pairs_checked:,}\n")
            else:
                f.write(f"  Total pairs checked: {pairs_checked:,}\n")
        f.write("\n")
        
        f.write("Filter Settings:\n")
        f.write(f"  GC content: {gc_min:.1%} - {gc_max:.1%}\n")
        f.write(f"  Max homopolymer: {homopolymer_max}\n")
        f.write(f"  Minimum distance: {min_distance}\n\n")
        
        if biological_violations:
            f.write("Biological Filter (GC content and homopolymer) Violations:\n")
            f.write("-" * 30 + "\n")
            for line_num, seq, reason in biological_violations:
                f.write(f"Line {line_num}: {seq} - {reason}\n")
            f.write("\n")
        
        # Add distance violation details if available
        if violation_info is not None:
            f.write("Distance Violations:\n")
            f.write("-" * 19 + "\n")
            seq1_line, seq2_line, seq1_str, seq2_str, distance = violation_info
            f.write(f"Line {seq1_line}: {seq1_str} and Line {seq2_line}: {seq2_str} - distance {distance} (minimum required: {min_distance})\n")
            f.write("\n")
    
    return report_file

def setup_argument_parser():
    """Setup and return the argument parser for barcode validation"""
    parser = argparse.ArgumentParser(
        description="Validate DNA barcodes against quality filters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--input', type=str, required=True, nargs='+',
                       help='Input file(s) containing DNA barcodes (one per line)')
    
    # Filter arguments
    parser.add_argument('--gc-min', type=float, default=0.4,
                       help='Minimum GC content (as fraction, e.g., 0.4 = 40%%)')
    parser.add_argument('--gc-max', type=float, default=0.6,
                       help='Maximum GC content (as fraction, e.g., 0.6 = 60%%)')
    parser.add_argument('--homopolymer-max', type=int, default=2,
                       help='Maximum allowed homopolymer length')
    parser.add_argument('--min-distance', type=int, default=3,
                       help='Minimum Hamming distance between sequences')
    
    # Performance arguments
    parser.add_argument('--skip-distance', action='store_true',
                       help='Skip distance validation entirely')
    parser.add_argument('--cpus', type=int, default=mp.cpu_count(),
                       help='Number of CPUs to use for parallel distance validation (default: all available)')
    
    # Output arguments
    parser.add_argument('--output-dir', type=str, default='test',
                       help='Output directory for validation logs and reports')
    
    return parser

def setup_logging(args):
    """Setup logging and create output directory. Returns log filepath."""
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup logging
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"validate_barcode_{timestamp}.log"
    log_filepath = os.path.join(args.output_dir, log_filename)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S',
        handlers=[
            logging.FileHandler(log_filepath),
            logging.StreamHandler()
        ]
    )
    
    return log_filepath

def validate_input_files(args):
    """Validate that input files exist"""
    for input_file in args.input:
        if not os.path.exists(input_file):
            raise ValueError(f"Input file does not exist: {input_file}")

def main():
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    log_filepath = setup_logging(args)
    validate_arguments(args)
    validate_input_files(args)
    
    # Run validation
    start_time = time.time()
    is_valid, total_violations, biological_violations, valid_sequences, sequences, early_stopped, distance_skipped, violation_info, validation_method, pairs_checked = validate_barcodes(
        input_files=args.input,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance,
        skip_distance=args.skip_distance,
        cpus=args.cpus
    )
    duration = time.time() - start_time
    
    # Generate report
    report_file = generate_validation_report(
        input_file=args.input,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance,
        biological_violations=biological_violations,
        valid_sequences=valid_sequences,
        total_sequences=len(sequences),
        output_dir=args.output_dir,
        early_stopped=early_stopped,
        distance_skipped=distance_skipped,
        validation_method=validation_method,
        pairs_checked=pairs_checked,
        violation_info=violation_info
    )
    
    # Log timing and file locations
    logging.info(f"Total time: {duration:.2f} seconds")
    logging.info(f"Log file: {log_filepath}")
    logging.info(f"Report file: {report_file}")
    
    # Final output
    if is_valid:
        print("✅ All barcodes are valid!")
    else:
        print(f"❌ VALIDATION FAILED: {total_violations} violations found")

if __name__ == "__main__":
    main()
