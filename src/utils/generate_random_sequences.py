#!/usr/bin/env python3
"""
generate_random_sequences.py

Generate random DNA sequences for testing validation scripts (variable length sequences supported).

Input: none

Output: random DNA sequences (one per line as .txt)

Optional arguments:
--output: output file path (default: test/{count}_random_{min}to{max}_bp_sequences.txt)

Required arguments:
--count: number of sequences to generate
--lengths: possible lengths for sequences (e.g., 15 20 25)
"""

import argparse
import random
import os
import numpy as np

# Import DNA encoding utilities
from dna_utils import decode_sequence

def generate_random_sequence(length):
    """Generate a single random DNA sequence of given length as integer array"""
    return np.random.choice([0, 1, 2, 3], size=length).astype(np.int8)

def main():
    parser = argparse.ArgumentParser(
        description="Generate random DNA sequences for testing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--count', type=int, required=True,
                       help='Number of sequences to generate')
    parser.add_argument('--lengths', nargs='+', type=int, required=True,
                       help='Possible lengths for sequences (e.g., 15 20 25)')
    parser.add_argument('--output', type=str, 
                       help='Output file path (default: test/{count}_random_{min}to{max}_bp_sequences.txt)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.count <= 0:
        raise ValueError("Count must be > 0")
    
    for length in args.lengths:
        if length <= 0:
            raise ValueError(f"All lengths must be > 0, got {length}")
    
    # Check if count exceeds maximum possible sequences across all lengths
    total_max_possible = sum(4 ** length for length in args.lengths)
    if args.count > total_max_possible:
        raise ValueError(f"Count ({args.count}) exceeds maximum possible sequences for lengths {args.lengths} ({total_max_possible})")
    
    # Generate default output path if not specified
    if args.output is None:
        os.makedirs('test', exist_ok=True)
        min_length = min(args.lengths)
        max_length = max(args.lengths)
        
        # Handle single length vs range
        if min_length == max_length:
            length_range = str(min_length)
        else:
            length_range = f"{min_length}to{max_length}"
        
        args.output = f"test/{args.count}_random_{length_range}_bp_sequences.txt"
    
    # Generate sequences with random lengths
    sequences = []
    for _ in range(args.count):
        # Randomly choose a length from the provided options
        chosen_length = random.choice(args.lengths)
        seq_array = generate_random_sequence(chosen_length)
        sequences.append(seq_array)
    
    # Write to file (convert to DNA strings for output)
    with open(args.output, 'w') as f:
        for seq_array in sequences:
            dna_string = decode_sequence(seq_array)
            f.write(dna_string + '\n')
    
    # Count sequences by length
    length_counts = {}
    for seq_array in sequences:
        length = len(seq_array)
        length_counts[length] = length_counts.get(length, 0) + 1
    
    length_breakdown = ", ".join([f"{count} at length {length}" for length, count in sorted(length_counts.items())])
    
    print(f"Generated {args.count} random DNA sequences with lengths: {length_breakdown}")
    print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main() 