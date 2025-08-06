#!/usr/bin/env python3
"""
generate_barcodes.py

Generate diverse DNA barcodes from scratch using iterative growth algorithm.

This script creates large sets of DNA sequences that satisfy minimum distance (between-sequence) and
biological quality constraints (within-sequence) using an efficient iterative growth approach.
"""

import numpy as np
import argparse
import logging
import time
import multiprocessing as mp
import os
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor

# Import utility functions
from utils.dna_utils import decode_sequence
from utils.filter_utils import hamming_distance_int, check_gc_content_int, check_homopolymer_int

def check_candidate_distance(candidate, existing_pool, min_distance):
    """Check if candidate sequence is sufficiently distant from all existing sequences"""
    for existing_seq in existing_pool:
        if hamming_distance_int(candidate, existing_seq, min_distance) < min_distance:
            return False
    return True

def passes_biological_filters_int(seq_array, gc_min, gc_max, homopolymer_max):
    """Check if sequence passes all biological quality filters (works with integer arrays)"""
    # Check GC content
    if not check_gc_content_int(seq_array, gc_min, gc_max):
        return False
    
    # Check homopolymer runs
    if not check_homopolymer_int(seq_array, homopolymer_max):
        return False
    
    return True

def generate_random_sequences(count, length, gc_min, gc_max, homopolymer_max):
    """Generate batch of random DNA sequences passing within-sequence biological filters"""
    sequences = []
    seen_sequences = set()
    
    # For reasonable batch sizes, the duplicate probability is low enough 
    # that we can optimize for speed over perfect duplicate detection
    max_attempts = count * 50  # Prevent infinite loops in saturated spaces
    attempts = 0
    
    while len(sequences) < count and attempts < max_attempts:
        attempts += 1
        # Generate integer array directly (A=0, T=1, G=2, C=3)
        seq_array = np.random.choice([0, 1, 2, 3], size=length).astype(np.int8)
        
        # Apply biological filters first (cheaper than duplicate check)
        if passes_biological_filters_int(seq_array, gc_min, gc_max, homopolymer_max):
            # Only check duplicates for sequences that pass bio filters
            seq_tuple = tuple(seq_array)
            if seq_tuple not in seen_sequences:
                seen_sequences.add(seq_tuple)
                sequences.append(seq_array)
    
    if len(sequences) < count:
        logging.warning(f"Could only generate {len(sequences)}/{count} unique sequences after {max_attempts} attempts")
    
    return sequences

def filter_chunk(candidates_chunk, existing_pool, min_distance):
    """Filter a chunk of candidates (helper function for parallel processing)"""
    valid_candidates = []
    for candidate in candidates_chunk:
        if check_candidate_distance(candidate, existing_pool, min_distance):
            valid_candidates.append(candidate)
    return valid_candidates

def filter_candidates_parallel(candidates, existing_pool, min_distance, n_cpus):
    """Parallel filtering of candidates using multiple processes"""
    # Skip parallelization for small batches (overhead not worth it)
    if len(candidates) < n_cpus * 10:
        return filter_chunk(candidates, existing_pool, min_distance)
    
    # Split candidates into chunks for parallel processing
    chunk_size = max(100, len(candidates) // (n_cpus * 4))
    chunks = [candidates[i:i+chunk_size] for i in range(0, len(candidates), chunk_size)]
    
    valid_candidates = []
    with ProcessPoolExecutor(max_workers=n_cpus) as executor:
        futures = [
            executor.submit(filter_chunk, chunk, existing_pool, min_distance)
            for chunk in chunks
        ]
        
        for future in futures:
            valid_candidates.extend(future.result())
    
    return valid_candidates

def generate_barcodes(target_count, length, gc_min, gc_max, homopolymer_max, min_distance, 
                     n_cpus, n_batches):
    """Main function to generate diverse barcode set using iterative growth"""
    # Calculate adaptive batch size based on number of batches
    # Ensure minimum effective batch size for reasonable success rates
    ideal_batch_size = target_count // n_batches
    min_effective_batch_size = max(50, target_count // 20)  # At least 50, or 5% of target
    batch_size = max(min_effective_batch_size, ideal_batch_size)
    
    logging.info(f"Starting barcode generation...")
    logging.info(f"Seed: None (building from scratch)")
    logging.info(f"Target count: {target_count} barcode sequences of length {length}")
    logging.info(f"Filter 1 (within-sequence), GC content: {gc_min:.1%} - {gc_max:.1%}")
    logging.info(f"Filter 2 (within-sequence), Max homopolymer: {homopolymer_max}")
    logging.info(f"Filter 3 (between-sequence), Minimum distance: {min_distance}")
    logging.info(f"CPUs: {n_cpus}, Initial batch size: {batch_size}")
    
    # Initialize empty pool
    selected_pool = []
    
    start_time = time.time()
    batch_num = 0
    total_generated = 0
    total_processed = 0
    
    while len(selected_pool) < target_count:
        batch_num += 1
        batch_start = time.time()
        
        # Generate random candidate sequences with biological filtering
        logging.info(f"Batch {batch_num}: (Step 1/3) Generating {batch_size} random candidates that pass within-sequence filters...")
        candidates = generate_random_sequences(batch_size, length, gc_min, gc_max, homopolymer_max)
        total_generated += len(candidates)
        
        # Filter candidates for distance constraints (parallel)
        logging.info(f"Batch {batch_num}: (Step 2/3) Filtering candidates for distance ≥{min_distance} with existing pool...")
        valid_candidates = filter_candidates_parallel(candidates, selected_pool, min_distance, n_cpus)
        
        # Within-batch distance checking
        newly_selected = []
        sequences_needed = target_count - len(selected_pool)
        
        logging.info(f"Batch {batch_num}: (Step 3/3) Filtering candidates for distance ≥{min_distance} within a batch...")
        for candidate in valid_candidates:
            valid = True
            # Check against sequences already accepted in THIS batch
            for new_seq in newly_selected:
                if hamming_distance_int(candidate, new_seq, min_distance) < min_distance:
                    valid = False
                    break
            if valid:
                newly_selected.append(candidate)
                # Early stopping if we have enough sequences
                if len(newly_selected) >= sequences_needed:
                    break
        
        # Add the verified batch to pool
        selected_pool.extend(newly_selected)
        total_processed += len(candidates)
        
        # Calculate metrics
        batch_time = time.time() - batch_start
        final_pass_rate = len(newly_selected) / len(candidates) * 100 if candidates else 0
        
        logging.info(f"Batch {batch_num}: Found {len(valid_candidates)} candidates, added {len(newly_selected)} sequences that satisfied all constraints")
        logging.info(f"  Final pass rate: {final_pass_rate:.1f}% ({batch_time:.1f}s)")
        logging.info(f"Progress: {len(selected_pool)}/{target_count} sequences "
                    f"({len(selected_pool)/target_count*100:.1f}%)")
        
        # Check if we've reached target
        if len(selected_pool) >= target_count:
            break
    
    total_time = time.time() - start_time
    overall_pass_rate = len(selected_pool) / total_processed * 100 if total_processed > 0 else 0
    
    logging.info(f"Generation complete!")
    logging.info(f"Final count: {len(selected_pool)} sequences")
    logging.info(f"Total candidates processed: {total_processed}")
    logging.info(f"Total time: {total_time:.1f} seconds")
    logging.info(f"Overall pass rate: {overall_pass_rate:.2f}%")
    logging.info(f"Sequences per second: {len(selected_pool)/total_time:.1f}")
    
    # Convert back to DNA strings
    return [decode_sequence(seq) for seq in selected_pool]

def validate_arguments(args):
    """Validate command line arguments and raise ValueError if invalid"""
    if args.gc_min < 0 or args.gc_max > 1 or args.gc_min >= args.gc_max:
        raise ValueError("GC content bounds must be: 0 ≤ gc_min < gc_max ≤ 1")
    
    if args.homopolymer_max < 1:
        raise ValueError("Maximum homopolymer length must be ≥ 1")
    
    if args.homopolymer_max >= args.length:
        raise ValueError("Maximum homopolymer length must be < sequence length")
    
    if args.min_distance < 1:
        raise ValueError("Minimum distance must be ≥ 1")
    
    if args.min_distance >= args.length:
        raise ValueError("Minimum distance must be < sequence length")
    
    if args.count <= 0:
        raise ValueError("Count must be > 0")
    
    if args.length <= 0:
        raise ValueError("Length must be > 0")

def main():
    parser = argparse.ArgumentParser(
        description="Generate diverse DNA barcodes from scratch",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--count', type=int, required=True,
                       help='Number of barcodes to generate')
    parser.add_argument('--length', type=int, required=True,
                       help='Length of each barcode sequence')
    
    # Output arguments
    parser.add_argument('--output-dir', type=str, default='test',
                       help='Output directory for barcodes and logs')
    parser.add_argument('--output-prefix', type=str, default='barcodes',
                       help='Output filename prefix (adds .txt automatically)')
    
    # Optional arguments with defaults
    parser.add_argument('--gc-min', type=float, default=0.4,
                       help='Minimum GC content (as fraction, e.g., 0.4 = 40%%)')
    parser.add_argument('--gc-max', type=float, default=0.6,
                       help='Maximum GC content (as fraction, e.g., 0.6 = 60%%)')
    parser.add_argument('--homopolymer-max', type=int, default=2,
                       help='Maximum allowed homopolymer length')
    parser.add_argument('--min-distance', type=int, default=3,
                       help='Minimum Hamming distance between sequences')
    
    # Performance arguments
    parser.add_argument('--cpus', type=int, default=mp.cpu_count(),
                       help='Number of CPU cores to use')
    parser.add_argument('--batches', type=int, default=10,
                       help='Number of batches to divide generation into (controls batch size)')
    
    # Output arguments
    parser.add_argument('--paired', action='store_true',
                       help='Generate paired barcodes (doubles target count, splits into two files)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup logging with file output
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"generate_barcode_{timestamp}.log"
    log_filepath = os.path.join(args.output_dir, log_filename)
    
    # Configure logging to both file and console
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
    
    # Adjust target count for paired mode
    actual_target_count = args.count * 2 if args.paired else args.count
    
    # Generate barcodes
    barcodes = generate_barcodes(
        target_count=actual_target_count,
        length=args.length,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance,
        n_cpus=args.cpus,
        n_batches=args.batches
    )
    
    # Write output
    if args.paired:
        # Split barcodes randomly into two files
        import random
        random.shuffle(barcodes)
        split_point = len(barcodes) // 2
        
        # Write first paired file
        paired1_filepath = os.path.join(args.output_dir, f"{args.output_prefix}_paired1.txt")
        with open(paired1_filepath, 'w') as f:
            for barcode in barcodes[:split_point]:
                f.write(barcode + '\n')
        
        # Write second paired file
        paired2_filepath = os.path.join(args.output_dir, f"{args.output_prefix}_paired2.txt")
        with open(paired2_filepath, 'w') as f:
            for barcode in barcodes[split_point:]:
                f.write(barcode + '\n')
        
        # Log file locations
        logging.info(f"Paired files written to:")
        logging.info(f"  {paired1_filepath} ({split_point} barcodes)")
        logging.info(f"  {paired2_filepath} ({len(barcodes) - split_point} barcodes)")
        logging.info(f"Log file: {log_filepath}")
        
        print(f"Successfully generated {len(barcodes)} barcodes (paired)")
    else:
        # Write single output file
        output_filepath = os.path.join(args.output_dir, f"{args.output_prefix}.txt")
        with open(output_filepath, 'w') as f:
            for barcode in barcodes:
                f.write(barcode + '\n')
        
        # Log file location
        logging.info(f"Output written to: {output_filepath}")
        logging.info(f"Log file: {log_filepath}")
        
        print(f"Successfully generated {len(barcodes)} barcodes")

if __name__ == "__main__":
    main()