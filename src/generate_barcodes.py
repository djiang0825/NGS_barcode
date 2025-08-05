#!/usr/bin/env python3
"""
generate_barcodes.py

Generate diverse DNA barcodes from scratch using iterative growth algorithm.

This script creates large sets of DNA sequences that satisfy minimum distance and
biological quality constraints using an efficient iterative growth approach.
"""

import numpy as np
import argparse
import logging
import time
import random
import multiprocessing as mp
from datetime import datetime
from numba import jit
from concurrent.futures import ProcessPoolExecutor

# DNA encoding constants
DNA_BASES = 'ATGC'
DNA_TO_INT = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
INT_TO_DNA = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}

def decode_sequence(seq_array):
    """Convert integer array back to DNA string"""
    return ''.join(INT_TO_DNA[base] for base in seq_array)

@jit(nopython=True)
def hamming_distance_int(seq1, seq2, min_distance):
    """Calculate Hamming distance with early stopping (assumes equal-length sequences)"""
    distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1
            if distance >= min_distance:
                return distance  # Early stopping
    return distance

def check_candidate_distance(candidate, existing_pool, min_distance):
    """Check if candidate sequence is sufficiently distant from all existing sequences"""
    for existing_seq in existing_pool:
        if hamming_distance_int(candidate, existing_seq, min_distance) < min_distance:
            return False
    return True

@jit(nopython=True)
def check_gc_content_int(seq_array, gc_min, gc_max):
    """Check if sequence passes GC content filter (works with integer arrays)"""
    # G=2, C=3 in our encoding - count them directly
    gc_count = 0
    for base in seq_array:
        if base == 2 or base == 3:  # G or C
            gc_count += 1
    gc_content = gc_count / len(seq_array)
    return gc_min <= gc_content <= gc_max

@jit(nopython=True)
def check_homopolymer_int(seq_array, max_length):
    """Check for homopolymer runs longer than max_length (works with integer arrays)"""
    if len(seq_array) == 0:
        return True
        
    current_base = seq_array[0]
    current_count = 1
    
    for base in seq_array[1:]:
        if base == current_base:
            current_count += 1
            if current_count > max_length:
                return False  # Fails check
        else:
            current_base = base
            current_count = 1
    
    return True  # Passes check

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
    """Generate batch of random DNA sequences passing biological filters"""
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
    
    logging.info(f"Starting barcode generation")
    logging.info(f"Target count: {target_count} sequences of length {length}")
    logging.info(f"GC content: {gc_min:.1%} - {gc_max:.1%}")
    logging.info(f"Max homopolymer: {homopolymer_max}")
    logging.info(f"Minimum distance: {min_distance}")
    logging.info(f"CPUs: {n_cpus}, Batches: {n_batches}, Batch size: {batch_size}")
    
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
        logging.info(f"Batch {batch_num}: Generating {batch_size} random candidates...")
        candidates = generate_random_sequences(batch_size, length, gc_min, gc_max, homopolymer_max)
        total_generated += len(candidates)
        
        # Filter candidates for distance constraints (parallel)
        logging.info(f"Batch {batch_num}: Filtering candidates for distance ≥{min_distance}...")
        valid_candidates = filter_candidates_parallel(candidates, selected_pool, min_distance, n_cpus)
        
        # Within-batch distance checking
        newly_selected = []
        sequences_needed = target_count - len(selected_pool)
        
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
        initial_success_rate = len(valid_candidates) / len(candidates) * 100 if candidates else 0
        final_success_rate = len(newly_selected) / len(candidates) * 100 if candidates else 0
        
        logging.info(f"Batch {batch_num}: Found {len(valid_candidates)} candidates, added {len(newly_selected)} sequences")
        logging.info(f"  Success rates: Initial {initial_success_rate:.1f}%, Final {final_success_rate:.1f}% ({batch_time:.1f}s)")
        logging.info(f"Progress: {len(selected_pool)}/{target_count} sequences "
                    f"({len(selected_pool)/target_count*100:.1f}%)")
        
        # Check if we've reached target
        if len(selected_pool) >= target_count:
            break
    
    total_time = time.time() - start_time
    overall_success_rate = len(selected_pool) / total_processed * 100 if total_processed > 0 else 0
    
    logging.info(f"Generation complete!")
    logging.info(f"Final count: {len(selected_pool)} sequences")
    logging.info(f"Total candidates processed: {total_processed}")
    logging.info(f"Total time: {total_time:.1f} seconds")
    logging.info(f"Overall success rate: {overall_success_rate:.2f}%")
    logging.info(f"Sequences per second: {len(selected_pool)/total_time:.1f}")
    
    # Convert back to DNA strings
    return [decode_sequence(seq) for seq in selected_pool]

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
    parser.add_argument('--output', type=str, required=True,
                       help='Output file path for generated barcodes')
    
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
    
    # Logging
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )
    
    # Validate arguments
    if args.gc_min < 0 or args.gc_max > 1 or args.gc_min >= args.gc_max:
        raise ValueError("GC content bounds must be: 0 ≤ gc_min < gc_max ≤ 1")
    
    if args.homopolymer_max < 1:
        raise ValueError("Maximum homopolymer length must be ≥ 1")
    
    if args.min_distance < 1:
        raise ValueError("Minimum distance must be ≥ 1")
    
    if args.count <= 0:
        raise ValueError("Count must be > 0")
    
    if args.length <= 0:
        raise ValueError("Length must be > 0")
    
    # Generate barcodes
    barcodes = generate_barcodes(
        target_count=args.count,
        length=args.length,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance,
        n_cpus=args.cpus,
        n_batches=args.batches
    )
    
    # Write output
    with open(args.output, 'w') as f:
        for barcode in barcodes:
            f.write(barcode + '\n')
    
    print(f"Successfully generated {len(barcodes)} barcodes")
    print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main()