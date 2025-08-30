#!/usr/bin/env python3
"""
generate_barcodes.py

Generate high-performance DNA barcodes for NGS applications using optimized iterative growth algorithm (supports extension from seed sequences and paired mode for dual-indexing).

Program Overview:

1. Load seed sequence files as existing pool and report length distribution (will generate from scratch if no seeds are provided)
2. Validate intput arguments, including whether the requested count is possible for the barcode length and minimum distance given Hamming and Gilbert-Varshamov bounds (when no seeds provided or everything equal-length)
   * Target count is updated given whether seeds are provided and/or paired mode is on

3. Generate a batch of candidate unique random sequences that pass the within-sequence biological filters (i.e., GC content and homopolymer repeats) → step 1/3 of the overall filtering process
   * Batch size is capped at 10,000 for barcode sets larger than 10,000, otherwise enforces 10 batches of minimum size 50
4. Conduct two-step between-sequence distance filtering with optimized method selection (neighbor enumeration or pairwise):
   4a. Filter candidates against existing pool → step 2/3 of the overall filtering process
     * If pairwise, use parallel processing with a minimum chunk size of 100 when candidate size is large enough, otherwise sequential (determined adaptively)
   4b. Filter remaining candidates against each other within the current batch → step 3/3 of the overall filtering process
     * If pairwise, always proceed sequentially in this step to ensure that the final pool fully satisfies the minimum distance requirement

   Method selection (applies to both steps, decided once per generation):
   - Small barcode sets (<10K sequences counting seeds if seeds are present): Always use pairwise distance checking
   - Mixed-length sequences (within seeds and/or between seeds and new barcodes): Always use pairwise distance checking
   - Large equal-length barcode sets (no seeds or everything equal-length): Choose between neighbor enumeration vs pairwise based on min_distance
     * Neighbor enumeration: when min_distance <= 4 (limited number of neighbors to check)
     * Pairwise distance checking: when min_distance > 4 (large number of neighbors to check)
5. Add the verified batch (candidates that pass all 3 steps of filtering) to pool, repeat until target count is reached

6. Write outputs to files and log results

Input: none (or optionally, seed sequence files or paired seed files)

Output: barcode list (one per line as .txt; if paired mode: two files with suffixes _paired1.txt and _paired2.txt) and generate_barcodes_{timestamp}.log file

Optional arguments:
--gc-min: minimum GC content (default: 0.4)
--gc-max: maximum GC content (default: 0.6)
--homopolymer-max: maximum allowed homopolymer repeat length (default: 2)
--min-distance: minimum edit distance between barcodes (default: 3)
--cpus: number of CPU cores to use during the parallel filtering step (default: all available)
--seeds: seed sequence files (any number of .txt files with one sequence per line, multiple files will be concatenated automatically; if not provided, will generate from scratch; incompatible with --paired mode; default: None) 
--paired: generate paired barcodes (doubles target count, randomly splits output into two equal parts; incompatible with --seeds; default: off)
--paired-seed1: paired seed sequence file 1 (used only with --paired and --paired-seed2, only one file is accepted, all sequences must be same length and match count/length of --paired-seed2; default: None)
--paired-seed2: paired seed sequence file 2 (used only with --paired and --paired-seed1, only one file is accepted, all sequences must be same length and match count/length of --paired-seed1; default: None) 
--output-dir: output directory for barcodes and logs (default: test)
--output-prefix: output filename prefix (adds .txt automatically; if paired mode is on, adds _paired1.txt and _paired2.txt; default: barcodes)

Required arguments:
--count: number of barcodes or barcode pairs to generate
--length: length of barcodes or barcode pairs to generate

NOTE: seed sequences (paired or unpaired) are not validated against intended filters. If necessary, please run validate_barcodes.py first to ensure they pass all the filters before providing them as seeds here.
"""

import numpy as np
import argparse
import logging
import math
import time
import multiprocessing as mp
import os
import random
from concurrent.futures import ProcessPoolExecutor

# Import utility functions
from utils.config_utils import decode_sequence
from utils.filter_utils import hamming_distance_int, check_gc_content_int, check_homopolymer_int, calculate_distance, generate_hamming_neighbors, validate_filter_arguments, select_distance_method
from utils.config_utils import setup_logging, read_files

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
    
    # Use NumPy's newer, thread-safe random number generator
    rng = np.random.default_rng()
    
    # For reasonable batch sizes, the duplicate probability is low enough 
    # that we can optimize for speed over perfect duplicate detection
    max_attempts = count * 50  # Prevent infinite loops in saturated spaces
    attempts = 0
    
    while len(sequences) < count and attempts < max_attempts:
        attempts += 1
        # Generate integer array directly (A=0, T=1, G=2, C=3) with explicit dtype
        seq_array = rng.integers(0, 4, size=length, dtype=np.int8)
        
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
        # Check if candidate is sufficiently distant from all existing sequences
        is_valid = True
        for existing_seq in existing_pool:
            if calculate_distance(candidate, existing_seq, min_distance) < min_distance:
                is_valid = False
                break
        
        if is_valid:
            valid_candidates.append(candidate)
    return valid_candidates

def filter_candidates_neighbor_enum(candidates, existing_pool, min_distance):
    """Filter candidates using neighbor enumeration for equal-length sequences"""
    if not candidates:
        return []
    
    # Convert existing pool to hash set for O(1) lookup
    # Note: This function is only called when all sequences are guaranteed to be the same length
    existing_set = set(tuple(seq) for seq in existing_pool)
    
    valid_candidates = []
    for candidate in candidates:
        seq_array = list(candidate)  # Make mutable copy for neighbor generation
        
        # Generate all neighbors within min_distance-1 and check for collisions
        is_valid = True
        for neighbor in generate_hamming_neighbors(seq_array, min_distance - 1):
            if neighbor in existing_set:
                is_valid = False
                break
        
        if is_valid:
            valid_candidates.append(candidate)
    
    return valid_candidates

def filter_candidates(candidates, existing_pool, min_distance, n_cpus, method="pairwise"):
    """Enhanced filtering with method selection (neighbor enumeration or parallel pairwise)"""
    if not candidates:
        return []
    
    if method == "neighbor_enumeration":
        # Use neighbor enumeration when set up is optimal (no parallelization needed - it's fast enough)
        return filter_candidates_neighbor_enum(candidates, existing_pool, min_distance)
    else:
        # Calculate optimal chunk size and number of chunks
        chunk_size = max(100, len(candidates) // (n_cpus * 4)) # minimum chunk size is 100
        num_chunks = (len(candidates) + chunk_size - 1) // chunk_size  # Ceiling division
        
        # Use sequential if we'd only create 1 chunk (parallelization not worthwhile)
        if num_chunks < 2:
            return filter_chunk(candidates, existing_pool, min_distance)
        
        # Use parallel processing with multiple chunks
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

def filter_within_batch(valid_candidates, min_distance, method, sequences_needed):
    """Filter candidates within a batch using the specified method, stopping when we have enough sequences"""
    newly_selected = []
    
    if method == "neighbor_enumeration":
        # Use neighbor enumeration approach
        newly_selected_set = set()  # Hash set for O(1) neighbor lookups
        
        for candidate in valid_candidates:
            if len(newly_selected) >= sequences_needed:
                break
                
            seq_array = list(candidate)  # Make mutable copy for neighbor generation
            
            # Check if any neighbor of this candidate is already in newly_selected
            collision_found = False
            for neighbor in generate_hamming_neighbors(seq_array, min_distance - 1):
                if neighbor in newly_selected_set:
                    collision_found = True
                    break
            
            if not collision_found:
                newly_selected.append(candidate)
                newly_selected_set.add(tuple(candidate))
    else:
        # Use original pairwise approach
        for candidate in valid_candidates:
            if len(newly_selected) >= sequences_needed:
                break
                
            valid = True
            # Check against sequences already accepted in THIS batch
            for new_seq in newly_selected:
                if hamming_distance_int(candidate, new_seq, min_distance) < min_distance:
                    valid = False
                    break
            if valid:
                newly_selected.append(candidate)
    
    return newly_selected

def generate_barcodes(target_count, length, gc_min, gc_max, homopolymer_max, min_distance, 
                     n_cpus, seed_pool=None, is_paired=False, has_mixed_lengths=False):
    """Main function to generate diverse barcode set using iterative growth"""
    logging.info(f"Starting barcode generation...")

    # 1. First log mode information (paired vs non-paired)
    original_target = target_count
    if is_paired:
        target_count *= 2
        logging.info(f"Mode: Paired (target count {original_target} → {target_count}, doubled)")
    else:
        logging.info(f"Mode: Standard (target count {target_count})")
    
    # Store the target count before adding seeds
    target_before_seeds = target_count
    
    # 2. Initialize and log seed information
    selected_pool = []
    
    if seed_pool:
        seed_count = len(seed_pool)
        selected_pool = seed_pool.copy()
        
        # Adjust target count to account for seeds
        target_count += seed_count
        
        # Log seed initialization and target count adjustment
        if is_paired:
            logging.info(f"Initialized pool with {seed_count} paired seed sequences ({seed_count // 2} pairs)")
            logging.info(f"Adjusted target count: {target_before_seeds} → {target_count} (including {seed_count} paired seed sequences)")
        else:
            logging.info(f"Initialized pool with {seed_count} seed sequences")
            logging.info(f"Adjusted target count: {target_before_seeds} → {target_count} (including {seed_count} seed sequences)")
        logging.warning("Building from seed lists without validation. Assuming that seed lists pass all the filters. Please run validate_barcodes.py to ensure this if necessary.")
    else:
        # No seeds for either mode
        logging.info("Seeds: None (building from scratch)")
    
    # 3. Make one global decision about which distance method to use
    # Use the shared utility function to determine the method with pre-calculated has_mixed_lengths
    method = select_distance_method(target_count, min_distance, has_mixed_lengths)
    
    # Log the decision
    if method == "neighbor_enumeration":
        logging.info(f"Using neighbor enumeration for distance checking (target count: {target_count}, min distance: {min_distance})")
    else:
        if target_count < 10000:
            logging.info(f"Using pairwise distance checking (small barcode set: {target_count} < 10K)")
        elif has_mixed_lengths:
            logging.info(f"Using pairwise distance checking (large barcode set: {target_count} ≥ 10K), mixed-length sequences detected)")
        else:
            logging.info(f"Using pairwise distance checking (large barcode set: {target_count} ≥ 10K), min distance {min_distance} > 4)")
    
    # 4. Calculate batch size after seed adjustment: cap at 10,000 for large barcode sets, use 10 batches for small barcode sets (min batch size is 50)
    if target_count <= 10000:
        # Small barcode sets: use 10 batches
        batch_size = max(50, target_count // 10)
    else:
        # Large barcode sets: cap at 10,000 per batch
        batch_size = 10000
    
    logging.info(f"Target count: {target_count} barcode sequences of length {length}")
    logging.info(f"Filter 1 (within-sequence), GC content: {gc_min:.1%} - {gc_max:.1%}")
    logging.info(f"Filter 2 (within-sequence), Max homopolymer repeat length: {homopolymer_max}")
    logging.info(f"Filter 3 (between-sequence), Minimum edit distance: {min_distance}")
    logging.info(f"CPUs: {n_cpus}; batch size: {batch_size}")
    
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
        
        # Filter candidates for distance constraints
        logging.info(f"Batch {batch_num}: (Step 2/3) Filtering candidates for distance ≥{min_distance} with existing pool...")
        valid_candidates = filter_candidates(candidates, selected_pool, min_distance, n_cpus, method)
        
        # Within-batch distance checking
        logging.info(f"Batch {batch_num}: (Step 3/3) Filtering candidates for distance ≥{min_distance} within a batch...")
        sequences_needed = target_count - len(selected_pool)
        newly_selected = filter_within_batch(valid_candidates, min_distance, method, sequences_needed)
        
        # Add the verified batch to pool
        selected_pool.extend(newly_selected)
        total_processed += len(candidates)
        
        # Calculate metrics
        batch_time = time.time() - batch_start
        final_pass_rate = len(newly_selected) / len(candidates) * 100 if candidates else 0
        
        logging.info(f"Batch {batch_num}: Found {len(valid_candidates)} candidates, added {len(newly_selected)} sequences that satisfied all constraints")
        logging.info(f"  Final pass rate: {final_pass_rate:.1f}% ({batch_time:.1f}s)")
        logging.info(f"Progress: {len(selected_pool)}/{target_count} sequences ({len(selected_pool)/target_count*100:.1f}%)")
        
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

def write_barcode_outputs(barcodes, args, paired_seed1_pool, paired_seed2_pool, log_filepath):
    """Write barcode outputs to files and log results"""
    if args.paired:
        # Treat no seeds as empty seed pools for unified logic
        seed1_strings = [decode_sequence(seq) for seq in paired_seed1_pool] if paired_seed1_pool else []
        seed2_strings = [decode_sequence(seq) for seq in paired_seed2_pool] if paired_seed2_pool else []
        
        # Calculate how many seeds to skip from barcodes list
        num_seeds = len(seed1_strings) + len(seed2_strings)
        new_barcodes = barcodes[num_seeds:] if num_seeds > 0 else barcodes
        
        # Split new barcodes into two equal groups
        random.shuffle(new_barcodes)
        split_point = len(new_barcodes) // 2
        
        # Write first paired file (seed1 + first half of new barcodes)
        paired1_filepath = os.path.join(args.output_dir, f"{args.output_prefix}_paired1.txt")
        with open(paired1_filepath, 'w') as f:
            for barcode in seed1_strings + new_barcodes[:split_point]:
                f.write(barcode + '\n')
        
        # Write second paired file (seed2 + second half of new barcodes)
        paired2_filepath = os.path.join(args.output_dir, f"{args.output_prefix}_paired2.txt")
        with open(paired2_filepath, 'w') as f:
            for barcode in seed2_strings + new_barcodes[split_point:]:
                f.write(barcode + '\n')
        
        # Log file locations with unified logic
        logging.info(f"Paired files written to:")
        if num_seeds > 0:
            logging.info(f"  {paired1_filepath} ({len(seed1_strings)} seeds + {split_point} new = {len(seed1_strings) + split_point} total)")
            logging.info(f"  {paired2_filepath} ({len(seed2_strings)} seeds + {len(new_barcodes) - split_point} new = {len(seed2_strings) + len(new_barcodes) - split_point} total)")
            print(f"Successfully generated {len(barcodes)} barcodes (paired with seeds)")
        else:
            logging.info(f"  {paired1_filepath} ({split_point} barcodes)")
            logging.info(f"  {paired2_filepath} ({len(new_barcodes) - split_point} barcodes)")
            print(f"Successfully generated {len(barcodes)} barcodes (paired)")
        logging.info(f"Log file: {log_filepath}")
    else:
        # Write single output file
        output_filepath = os.path.join(args.output_dir, f"{args.output_prefix}.txt")
        with open(output_filepath, 'w') as f:
            for barcode in barcodes:
                f.write(barcode + '\n')
        
        # Log file location
        logging.info(f"Output written to: {output_filepath}")
        logging.info(f"Log file: {log_filepath}")
        
        # Print result with seed info if applicable
        if args.seeds:
            print(f"Successfully generated {len(barcodes)} barcodes (with seeds)")
        else:
            print(f"Successfully generated {len(barcodes)} barcodes")

def setup_argument_parser():
    """Setup and return the argument parser for barcode generation"""
    parser = argparse.ArgumentParser(
        description="Generate high-performance DNA barcodes for NGS applications (from scratch or by extending from provided seed sequences)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--count', type=int, required=True,
                       help='Number of barcodes or barcode pairs to generate')
    parser.add_argument('--length', type=int, required=True,
                       help='Length of barcodes or barcode pairs to generate')
    
    # Output arguments
    parser.add_argument('--output-dir', type=str, default='test',
                       help='Output directory for barcodes and logs')
    parser.add_argument('--output-prefix', type=str, default='barcodes',
                       help='Output filename prefix (adds .txt automatically; when paired mode is on, adds _paired1.txt and _paired2.txt)')
    
    # Filter arguments with defaults
    parser.add_argument('--gc-min', type=float, default=0.4,
                       help='Minimum GC content (as fraction, e.g., 0.4 = 40%%;)')
    parser.add_argument('--gc-max', type=float, default=0.6,
                       help='Maximum GC content (as fraction, e.g., 0.6 = 60%%;)')
    parser.add_argument('--homopolymer-max', type=int, default=2,
                       help='Maximum allowed homopolymer repeat length')
    parser.add_argument('--min-distance', type=int, default=3,
                       help='Minimum edit distance between barcodes')
    
    # Performance arguments
    parser.add_argument('--cpus', type=int, default=mp.cpu_count(),
                       help='Number of CPU cores to use during the parallel filtering step')
    
    # Seed arguments
    parser.add_argument('--seeds', nargs='+', type=str, default=[],
                       help='Seed sequence files with existing barcode lists to extend from (any number of .txt files, one sequence per line; multiple files will be concatenated automatically) - incompatible with --paired mode')
    
    # Mode arguments
    parser.add_argument('--paired', action='store_true',
                       help='Generate paired barcodes by doubling the target count and randomly splitting the output into two equal parts, saved as _paired1.txt and _paired2.txt.  - incompatible with --seeds')
    
    # Paired seed arguments (for paired mode only)
    parser.add_argument('--paired-seed1', type=str, default=None,
                       help='Seed sequence file for one side of the paired barcode set to extend from (single .txt, one sequence per line; requires matching file provided via --paired-seed2) - used only with --paired')
    parser.add_argument('--paired-seed2', type=str, default=None,
                       help='Seed sequence file for the other side of the paired barcode set to extend from (single .txt, one sequence per line; requires matching file provided via --paired-seed1) - used only with --paired')
    
    return parser

def validate_seed_arguments(args):
    """Validate seed-related argument combinations"""
    # 1. Validate paired seed arguments
    if args.paired:
        # In paired mode, either both paired seeds or no seeds at all
        if (args.paired_seed1 is not None) != (args.paired_seed2 is not None):
            raise ValueError("When using --paired with seeds, both --paired-seed1 and --paired-seed2 must be provided")
        elif args.seeds:
            raise ValueError("--seeds argument is incompatible with --paired mode. Use --paired-seed1 and --paired-seed2 instead")
    else:
        # In non-paired mode, paired seed arguments should not be used
        if args.paired_seed1 is not None or args.paired_seed2 is not None:
            raise ValueError("--paired-seed1 and --paired-seed2 can only be used with --paired mode")
    
    # 2. Validate that seed files exist
    if args.seeds:
        for seed_file in args.seeds:
            if not os.path.exists(seed_file):
                raise ValueError(f"Seed file does not exist: {seed_file}")
    
    elif args.paired_seed1 and not os.path.exists(args.paired_seed1):
        raise ValueError(f"Paired seed file 1 does not exist: {args.paired_seed1}")
    
    elif args.paired_seed2 and not os.path.exists(args.paired_seed2):
        raise ValueError(f"Paired seed file 2 does not exist: {args.paired_seed2}")

def load_and_validate_seeds(args):
    """Load seed sequences and validate paired seeds. Returns (seed_pool, paired_seed1_pool, paired_seed2_pool, length_counts)"""
    seed_pool = []
    paired_seed1_pool = []
    paired_seed2_pool = []
    length_counts = {} # will return this for later use
    
    if args.seeds:
        seed_pool, length_counts = read_files(args.seeds)
    elif args.paired and args.paired_seed1 and args.paired_seed2:
        # Load paired seeds separately and combine for generation
        paired_seed1_pool, seed1_length_counts = read_files([args.paired_seed1])
        paired_seed2_pool, seed2_length_counts = read_files([args.paired_seed2])
        
        # Check that paired seeds are truly paired (same count and all same length)
        # 1. Check that both files have the same number of sequences
        if len(paired_seed1_pool) != len(paired_seed2_pool):
            raise ValueError(f"Paired seed files must have the same number of sequences. "
                           f"Seed1: {len(paired_seed1_pool)} sequences, Seed2: {len(paired_seed2_pool)} sequences")
        
        # 2. Check that both files have sequences of the same length within the file
        elif len(seed1_length_counts) != 1:
            raise ValueError(f"All sequences in paired seed file 1 must be the same length. "
                           f"Found lengths: {sorted(seed1_length_counts.keys())}")
        elif len(seed2_length_counts) != 1:
            raise ValueError(f"All sequences in paired seed file 2 must be the same length. "
                           f"Found lengths: {sorted(seed2_length_counts.keys())}")
        
        # 3. Check that both files have sequences of the same length between the files
        elif list(seed1_length_counts.keys())[0] != list(seed2_length_counts.keys())[0]:
            raise ValueError(f"Paired seed files must have sequences of the same length. "
                           f"Seed1 length: {list(seed1_length_counts.keys())[0]}, Seed2 length: {list(seed2_length_counts.keys())[0]}")
        else:
            # All validations passed - combine both for generation pool
            seed_pool = paired_seed1_pool + paired_seed2_pool
            
            # Since paired seeds are validated to have the same length, just use seed1's length counts
            # and double the count since we have two files
            for length, count in seed1_length_counts.items():
                length_counts[length] = count * 2
    
    return seed_pool, paired_seed1_pool, paired_seed2_pool, length_counts

def calculate_hamming_bound(length, min_distance):
    """Calculate theoretical maximum number of sequences for given length and minimum distance"""
    
    # Hamming bound: M ≤ 4^n / V(n, t) where t = floor((d-1)/2)
    # V(n, t) is the volume of a Hamming sphere of radius t
    total_sequences = 4 ** length
    t = (min_distance - 1) // 2
    
    # Calculate volume of Hamming sphere: V(n, t) = sum(C(n,i) * 3^i) for i=0 to t
    sphere_volume = 0
    for i in range(t + 1):
        # C(n, i) = n! / (i! * (n-i)!)
        combinations = math.comb(length, i)
        sphere_volume += combinations * (3 ** i)  # 3 possible mutations per position
    
    # Hamming bound
    max_sequences = total_sequences // sphere_volume
    return max_sequences

def calculate_gv_bound(length, min_distance):
    """Calculate Gilbert-Varshamov bound (lower bound) for given length and minimum distance"""
    
    # GV bound: M ≥ 4^n / V(n, d-1) where d is the minimum distance
    # V(n, d-1) is the volume of a Hamming sphere of radius d-1
    total_sequences = 4 ** length
    
    # Calculate volume of Hamming sphere: V(n, d-1) = sum(C(n,i) * 3^i) for i=0 to d-1
    sphere_volume = 0
    for i in range(min_distance):
        # C(n, i) = n! / (i! * (n-i)!)
        combinations = math.comb(length, i)
        sphere_volume += combinations * (3 ** i)  # 3 possible mutations per position
    
    # GV bound (minimum possible)
    min_sequences = total_sequences // sphere_volume
    return min_sequences

def validate_generator_arguments(args, seed_pool=None, length_counts=None):
    """Validate generator-specific arguments (length, count, min_distance, homopolymer_max) and return has_mixed_lengths flag"""
    # 1. Length validation
    if args.length <= 0:
        raise ValueError("Length must be > 0")
    
    # Length warning for very long barcodes
    if args.length > 20:
        logging.warning(f"Requested barcode length ({args.length}) exceeds 20bp. Very long barcodes can be synthetically unstable and more error-prone, and this program is optimized for short sequences. It is recommended to generate barcodes of length 20bp or less for expected behavior.")
    
    # 2. Homopolymer repeat x barcode length validation
    if args.homopolymer_max >= args.length:
        raise ValueError("Maximum homopolymer repeat length must be < new barcode length")
    
    # 3. Minimum distance x barcode length validation
    # Determine effective length for distance validation and check for mixed lengths
    effective_length = args.length
    seed_length = None
    has_mixed_lengths = False # will return this flag for later use
    
    if seed_pool:
        # Use length_counts to determine seed lengths
        seed_length = max(length_counts.keys())
        effective_length = max(args.length, seed_length)
        
        # Distance validation for case with seeds
        if args.min_distance >= effective_length:
            raise ValueError(f"Minimum distance must be < {effective_length} (the longer of new barcode length {args.length} and max seed length {seed_length})")
        
        # Check for mixed lengths (within seeds or between seeds and target), will return this flag for later use
        if len(length_counts) > 1 or seed_length != args.length:
            has_mixed_lengths = True
    
    else:
        # Distance validation for case without seeds
        if args.min_distance >= effective_length:
            raise ValueError(f"Minimum distance must be < {effective_length} (new barcode length)")
    
    # 4. Count validation
    if args.count <= 0:
        raise ValueError("Count must be > 0")
    
    # 5. Count x barcode length x minimum distance validation given Hamming and Gilbert-Varshamov bounds
    logging.info(f"Performing Hamming and Gilbert-Varshamov bounds validation for barcodes of length {args.length} with minimum distance {args.min_distance} to see if the requested count is possible...")
    
    # Check bounds only for simple cases (no mixed lengths)
    if not has_mixed_lengths:
        # Calculate Hamming and Gilbert-Varshamov bounds
        max_possible = calculate_hamming_bound(args.length, args.min_distance)
        min_possible = calculate_gv_bound(args.length, args.min_distance)
        
        # Format large numbers for readability
        max_formatted = f"{max_possible:,}" if max_possible > 1000 else str(max_possible)
        min_formatted = f"{min_possible:,}" if min_possible > 1000 else str(min_possible)
        
        logging.info(f"Bounds for barcode length {args.length}, min distance {args.min_distance}: GV (lower) bound = {min_formatted}, Hamming (upper) bound = {max_formatted}")
        
        # Validate against bounds
        if not seed_pool:
            # Case 1: No seeds - simple check
            if args.count > max_possible:
                raise ValueError(f"Requested count ({args.count:,}) exceeds Hamming (upper) bound ({max_formatted}) for length {args.length}, min distance {args.min_distance}. Please reduce the requested count.")
            elif args.count > min_possible:
                logging.warning(f"Requested count ({args.count:,}) exceeds GV (lower) bound ({min_formatted}) - may take longer to generate and could potentially fail")
        else:
            # Case 2: With seeds - check total sequences
            total_sequences = args.count + len(seed_pool)
            if total_sequences > max_possible:
                raise ValueError(f"Total sequences ({total_sequences:,}) exceeds Hamming (upper) bound ({max_formatted}) for length {args.length}, min distance {args.min_distance}. Please reduce the requested count or seed size.")
            elif total_sequences > min_possible:
                logging.warning(f"Total sequences ({total_sequences:,}) exceeds GV (lower) bound ({min_formatted}) - may take longer to generate and could potentially fail")
    else:
        # Complex case: mixed lengths - skip bounds validation
        logging.warning(f"Seeds have mixed lengths or different length(s) from target. Hamming and Gilbert-Varshamov bounds validation skipped.")
    
    return has_mixed_lengths

def main():
    parser = setup_argument_parser()
    args = parser.parse_args()
    log_filepath = setup_logging(args, "generate_barcodes")
    validate_filter_arguments(args) # simple validation of filter arguments
        
    # Only load and validate seeds if any seed arguments are provided
    if args.seeds or args.paired_seed1 or args.paired_seed2:
        validate_seed_arguments(args)
        seed_pool, paired_seed1_pool, paired_seed2_pool, length_counts = load_and_validate_seeds(args)
    else:
        # Initialize seed-related variables as empty
        seed_pool = []
        paired_seed1_pool = []
        paired_seed2_pool = []
        length_counts = {}

    # Perform complex validation on generator-specific arguments and get the has_mixed_lengths flag
    has_mixed_lengths = validate_generator_arguments(args, seed_pool, length_counts)

    # Generate barcodes
    barcodes = generate_barcodes(
        target_count=args.count,
        length=args.length,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        homopolymer_max=args.homopolymer_max,
        min_distance=args.min_distance,
        n_cpus=args.cpus,
        seed_pool=seed_pool,
        is_paired=args.paired,
        has_mixed_lengths=has_mixed_lengths
    )
    
    # Write outputs
    write_barcode_outputs(barcodes, args, paired_seed1_pool, paired_seed2_pool, log_filepath)

if __name__ == "__main__":
    main()