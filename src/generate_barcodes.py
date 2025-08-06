#!/usr/bin/env python3
"""
generate_barcodes.py

Efficiently generate large set of NGS barcodes from scratch or from seed sequences using iterative growth algorithm (paired mode supported).

Algorithm overview (guarantees no duplicate sequences and satisfies >= min distance constraints in final pool):

1. Generate random sequences that pass within-sequence biological filters (i.e., GC content and homopolymer repeats)
2. Filter candidates to ensure >= min distance constraints compared to existing pool (parallel)
3. Filter candidates to ensure >= min distance constraints within a batch (sequential)
4. Add the verified batch to pool
5. Repeat until target count is reached

Input: none (or optionally, seed sequence files)

Output: barcode list (one per line as .txt) and .log file

Optional arguments:
--gc-min: minimum GC content (default: 0.4)
--gc-max: maximum GC content (default: 0.6)
--homopolymer-max: maximum allowed homopolymer length (default: 2)
--min-distance: minimum Hamming distance between sequences (default: 3)
--cpus: number of CPU cores to use (default: all available)
--paired: generate paired barcodes (doubles target count, splits into two files; default: False)
--seeds: seed sequence files (any number of files, one sequence per line as .txt; if not provided, will generate from scratch; default: None) 
NOTE: seed sequences are not validated. If necessary, please run validate_barcodes.py first to ensure they pass all the filters.
--output-dir: output directory for barcodes and logs (default: test)
--output-prefix: output filename prefix (adds .txt automatically, default: barcodes)

Required arguments:
--count: number of barcodes to generate
--length: length of each barcode sequence
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
from utils.dna_utils import decode_sequence, encode_sequence, validate_arguments
from utils.filter_utils import hamming_distance_int, check_gc_content_int, check_homopolymer_int, calculate_distance

def check_candidate_distance(candidate, existing_pool, seed_pool, min_distance):
    """Check if candidate sequence is sufficiently distant from all existing and seed sequences"""
    # Check against existing pool (all same length as candidate)
    for existing_seq in existing_pool:
        if hamming_distance_int(candidate, existing_seq, min_distance) < min_distance:
            return False
    
    # Check against seed pool (may have different lengths)
    for seed_seq in seed_pool:
        if calculate_distance(candidate, seed_seq, min_distance) < min_distance:
            return False
    
    return True

def passes_biological_filters_int(seq_array, gc_min, gc_max, homopolymer_max):
    """Check if sequence passes all biological quality filters (works with integer arrays)"""
    # Check GC content
    if not check_gc_content_int(seq_array, gc_min, gc_max):
        return False
    
    # Check homopolymer repeats
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

def load_seed_sequences(seed_files):
    """Load seed sequences from files and convert to integer arrays"""
    if not seed_files:
        return []
    
    seed_sequences = []
    
    for seed_file in seed_files:
        file_count = 0
        with open(seed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                seq = line.strip()
                if not seq:  # Skip empty lines
                    continue
                
                # Basic validation
                if not all(base in 'ATGC' for base in seq):
                    logging.warning(f"Seed file {seed_file}, line {line_num}: Invalid DNA sequence '{seq}', skipping")
                    continue
                
                seq_array = encode_sequence(seq)
                seed_sequences.append(seq_array)
                file_count += 1
        
        logging.info(f"Loaded {file_count} sequences from {seed_file}")
    
    # Calculate and log length distribution
    if seed_sequences:
        lengths = [len(seq) for seq in seed_sequences]
        length_counts = {}
        for length in lengths:
            length_counts[length] = length_counts.get(length, 0) + 1
        
        if len(length_counts) == 1:
            length_info = f"length {list(length_counts.keys())[0]}"
        else:
            length_breakdown = ", ".join([f"{count} at length {length}" for length, count in sorted(length_counts.items())])
            length_info = f"mixed lengths: {length_breakdown}"
        
        logging.info(f"Total loaded: {len(seed_sequences)} seed sequences from {len(seed_files)} files ({length_info})")
    else:
        logging.info(f"Total loaded: 0 seed sequences from {len(seed_files)} files")
    
    return seed_sequences

def filter_chunk(candidates_chunk, existing_pool, seed_pool, min_distance):
    """Filter a chunk of candidates (helper function for parallel processing)"""
    valid_candidates = []
    for candidate in candidates_chunk:
        if check_candidate_distance(candidate, existing_pool, seed_pool, min_distance):
            valid_candidates.append(candidate)
    return valid_candidates

def filter_candidates_parallel(candidates, existing_pool, seed_pool, min_distance, n_cpus):
    """Parallel filtering of candidates using multiple processes"""
    # Skip parallelization for small batches (overhead not worth it)
    if len(candidates) < n_cpus * 10:
        return filter_chunk(candidates, existing_pool, seed_pool, min_distance)
    
    # Split candidates into chunks for parallel processing
    chunk_size = max(100, len(candidates) // (n_cpus * 4))
    chunks = [candidates[i:i+chunk_size] for i in range(0, len(candidates), chunk_size)]
    
    valid_candidates = []
    with ProcessPoolExecutor(max_workers=n_cpus) as executor:
        futures = [
            executor.submit(filter_chunk, chunk, existing_pool, seed_pool, min_distance)
            for chunk in chunks
        ]
        
        for future in futures:
            valid_candidates.extend(future.result())
    
    return valid_candidates

def generate_barcodes(target_count, length, gc_min, gc_max, homopolymer_max, min_distance, 
                     n_cpus, seed_pool=None):
    """Main function to generate diverse barcode set using iterative growth"""
    if seed_pool:
        logging.info(f"Starting barcode generation...")
        logging.info(f"Seed: {len(seed_pool)} sequences loaded from seed files")
        logging.warning("Building from seed lists without validation. Assuming that seed lists pass all the filters. Please run validate_barcodes.py to ensure this if necessary.")
    else:
        logging.info(f"Starting barcode generation...")
        logging.info(f"Seed: None (building from scratch)")
    
    # Initialize pool with seeds if provided
    selected_pool = []
    if seed_pool:
        selected_pool = seed_pool.copy()
        logging.info(f"Initialized pool with {len(selected_pool)} seed sequences")
        # Adjust target count to account for seeds
        target_count += len(seed_pool)
        logging.info(f"Adjusted target count to {target_count} (including {len(seed_pool)} seed sequences)")
    
    # Calculate batch size after seed adjustment: cap at 10,000 for large datasets, use 10 batches for small datasets
    if target_count <= 10000:
        # Small datasets: use 10 batches
        batch_size = max(50, target_count // 10)
    else:
        # Large datasets: cap at 10,000 per batch
        batch_size = 10000
    
    logging.info(f"Target count: {target_count} barcode sequences of length {length}")
    logging.info(f"Filter 1 (within-sequence), GC content: {gc_min:.1%} - {gc_max:.1%}")
    logging.info(f"Filter 2 (within-sequence), Max homopolymer: {homopolymer_max}")
    logging.info(f"Filter 3 (between-sequence), Minimum distance: {min_distance}")
    logging.info(f"CPUs: {n_cpus}, Initial batch size: {batch_size}")
    
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
        valid_candidates = filter_candidates_parallel(candidates, selected_pool, seed_pool or [], min_distance, n_cpus)
        
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
    
    # Seed arguments
    parser.add_argument('--seeds', nargs='+', type=str, default=[],
                       help='Seed sequence files (any number of files, one sequence per line)')
    
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
    
    # Load seed sequences if provided
    seed_pool = None
    if args.seeds:
        seed_pool = load_seed_sequences(args.seeds)
    
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
        seed_pool=seed_pool
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