#!/usr/bin/env python3
"""
validate_barcodes.py

Validate if a provided list of DNA barcodes satisfies all quality filters.

This script checks if barcode sequences meet the same biological and distance
constraints used in the generation script.
"""

import argparse
import logging
import os
from datetime import datetime
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import threading

# Import utility functions
from utils.dna_utils import DNA_BASES, encode_sequence
from utils.filter_utils import  check_gc_content_int, check_homopolymer_int, calculate_distance



def validate_sequence(seq_array, gc_min, gc_max, homopolymer_max):
    """Check if sequence passes all biological filters and return reason if it fails"""
    # Check GC content
    if not check_gc_content_int(seq_array, gc_min, gc_max):
        return False, "GC content outside range"
    
    # Check homopolymer runs
    if not check_homopolymer_int(seq_array, homopolymer_max):
        return False, "Homopolymer run too long"
    
    return True, "Passes all filters"

def check_distance_chunk(chunk_data):
    """Check distance for a chunk of sequence pairs - returns immediately on first violation"""
    sequences, min_distance, start_idx = chunk_data
    pairs_checked = 0
    
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            pairs_checked += 1
            distance = calculate_distance(sequences[i], sequences[j], min_distance)
            if distance < min_distance:
                # Return immediately on first violation
                return (start_idx + i, start_idx + j, distance), pairs_checked, True
    
    return None, pairs_checked, False

def validate_distance_constraints_parallel(sequences, min_distance, n_cpus):
    """Check if all sequences meet minimum distance requirements with parallel processing and immediate early stopping"""
    if len(sequences) < n_cpus * 10:  # Skip parallelization for small datasets
        return validate_distance_constraints(sequences, min_distance)
    
    total_pairs = len(sequences) * (len(sequences) - 1) // 2
    pairs_checked = 0
    violations = []
    
    # Create chunks of sequence pairs for parallel processing
    chunk_size = max(100, len(sequences) // (n_cpus * 4))
    chunks = []
    
    for i in range(0, len(sequences), chunk_size):
        end_idx = min(i + chunk_size, len(sequences))
        chunk_sequences = sequences[i:end_idx]
        chunks.append((chunk_sequences, min_distance, i))
    
    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=n_cpus) as executor:
        futures = [executor.submit(check_distance_chunk, chunk) for chunk in chunks]
        
        for future in futures:
            result, chunk_pairs, early_stop = future.result()
            pairs_checked += chunk_pairs
            
            if early_stop:
                # Found violation - cancel remaining futures and return immediately
                violations = [result]  # Only store the first violation found
                for f in futures:
                    f.cancel()
                logging.info(f"Early stopping: Found distance violation between sequences {result[0]+1} and {result[1]+1} (distance={result[2]})")
                return violations, pairs_checked, total_pairs
    
    return violations, pairs_checked, total_pairs

def validate_distance_constraints(sequences, min_distance):
    """Check if all sequences meet minimum distance requirements - exits early on first violation"""
    violations = []
    total_pairs = len(sequences) * (len(sequences) - 1) // 2
    pairs_checked = 0
    
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            pairs_checked += 1
            distance = calculate_distance(sequences[i], sequences[j], min_distance)
            if distance < min_distance:
                violations.append((i, j, distance))
                # Early exit - we found a violation, no need to continue
                logging.info(f"Early stopping: Found distance violation between sequences {i+1} and {j+1} (distance={distance})")
                return violations, pairs_checked, total_pairs
    
    return violations, pairs_checked, total_pairs

def validate_barcodes(input_file, gc_min, gc_max, homopolymer_max, min_distance, output_dir, skip_distance=False, n_cpus=1):
    """Main validation function"""
    logging.info(f"Starting barcode validation...")
    logging.info(f"Input file: {input_file}")
    logging.info(f"Filter 1 (within-sequence), GC content: {gc_min:.1%} - {gc_max:.1%}")
    logging.info(f"Filter 2 (within-sequence), Max homopolymer: {homopolymer_max}")
    logging.info(f"Filter 3 (between-sequence), Minimum distance: {min_distance}")
    
    # Read sequences
    sequences = []
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            seq = line.strip()
            if not seq:  # Skip empty lines
                continue
            
            # Basic validation
            if not all(base in DNA_BASES for base in seq):
                logging.error(f"Line {line_num}: Invalid DNA sequence '{seq}'")
                continue
            
            sequences.append((line_num, seq))
    
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
        is_valid, reason = validate_sequence(seq_array, gc_min, gc_max, homopolymer_max)
        
        if is_valid:
            seq_arrays.append(seq_array)
            valid_sequences.append((line_num, seq))
        else:
            biological_violations.append((line_num, seq, reason))
    
    logging.info(f"Biological filter results:")
    logging.info(f"  Passed: {len(valid_sequences)} sequences")
    logging.info(f"  Failed: {len(biological_violations)} sequences")
    
    # Check if we should skip distance validation
    if skip_distance and len(biological_violations) > 0:
        logging.info(f"Skipping distance validation due to biological filter failures")
        distance_violations = []
        pairs_checked = 0
        total_pairs = 0
    else:
        # Check distance constraints (with early stopping and parallel processing)
        distance_violations, pairs_checked, total_pairs = validate_distance_constraints_parallel(seq_arrays, min_distance, n_cpus)
    
    logging.info(f"Distance constraint results:")
    logging.info(f"  Total sequence pairs: {total_pairs} (input passed biological filters)")
    logging.info(f"  Pairs checked: {pairs_checked} (early stopping on first violation)")
    
    # Summary
    total_violations = len(biological_violations) + len(distance_violations)
    overall_valid = len(valid_sequences) == len(sequences) and len(distance_violations) == 0
    
    # Check if early stopping occurred
    early_stopped = len(distance_violations) > 0 and pairs_checked < total_pairs
    
    logging.info(f"Validation complete!")
    logging.info(f"Overall validation: {'PASSED' if overall_valid else 'FAILED'}")
    
    return overall_valid, len(valid_sequences), total_violations, biological_violations, distance_violations, valid_sequences, sequences, early_stopped

def generate_validation_report(input_file, gc_min, gc_max, homopolymer_max, min_distance, 
                             biological_violations, distance_violations, valid_sequences, 
                             total_sequences, output_dir, early_stopped=False):
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
        
        if early_stopped:
            f.write(f"Distance validation: EARLY STOPPED (found first violation)\n")
        else:
            f.write(f"Distance violations: {len(distance_violations)} pairs\n\n")
        
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
        
        if distance_violations:
            f.write("Distance Constraint Violations:\n")
            f.write("-" * 35 + "\n")
            for i, j, distance in distance_violations:
                f.write(f"Sequences {i+1} and {j+1}: distance = {distance} (min required: {min_distance})\n")
            f.write("\n")
        
        if valid_sequences and not early_stopped:
            f.write("Valid Sequences (passed all filters):\n")
            f.write("-" * 35 + "\n")
            for line_num, seq in valid_sequences:
                f.write(f"Line {line_num}: {seq}\n")
        elif valid_sequences and early_stopped:
            f.write("Sequences that passed biological filters (distance validation incomplete):\n")
            f.write("-" * 55 + "\n")
            for line_num, seq in valid_sequences:
                f.write(f"Line {line_num}: {seq}\n")
    
    return report_file

def validate_arguments(args):
    """Validate command line arguments"""
    if args.gc_min < 0 or args.gc_max > 1 or args.gc_min >= args.gc_max:
        raise ValueError("GC content bounds must be: 0 ≤ gc_min < gc_max ≤ 1")
    
    if args.homopolymer_max < 1:
        raise ValueError("Maximum homopolymer length must be ≥ 1")
    
    if args.min_distance < 1:
        raise ValueError("Minimum distance must be ≥ 1")
    
    if not os.path.exists(args.input):
        raise ValueError(f"Input file does not exist: {args.input}")

def main():
    parser = argparse.ArgumentParser(
        description="Validate DNA barcodes against quality filters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--input', type=str, required=True,
                       help='Input file containing DNA barcodes (one per line)')
    
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
                       help='Skip distance validation if biological filters fail')
    parser.add_argument('--cpus', type=int, default=mp.cpu_count(),
                       help='Number of CPU cores to use for parallel distance checking')
    
    # Output arguments
    parser.add_argument('--output-dir', type=str, default='test',
                       help='Output directory for validation reports')
    
    args = parser.parse_args()
    
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
    
    # Validate arguments
    validate_arguments(args)
    
    # Run validation
    is_valid, valid_count, total_violations, biological_violations, distance_violations, valid_sequences, sequences, early_stopped = validate_barcodes(
        input_file=args.input,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance,
        output_dir=args.output_dir,
        skip_distance=args.skip_distance,
        n_cpus=args.cpus
    )
    
    # Generate detailed report
    report_file = generate_validation_report(
        input_file=args.input,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance,
        biological_violations=biological_violations,
        distance_violations=distance_violations,
        valid_sequences=valid_sequences,
        total_sequences=len(sequences),
        output_dir=args.output_dir,
        early_stopped=early_stopped
    )
    
    # Final output
    if is_valid:
        print(f"✅ VALIDATION PASSED: All sequences meet quality standards")
    else:
        print(f"❌ VALIDATION FAILED: {total_violations} violations found")
    
    print(f"Log file: {log_filepath}")
    print(f"Report file: {report_file}")

if __name__ == "__main__":
    main() 