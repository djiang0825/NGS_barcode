#!/usr/bin/env python3
"""
generate_random_sequences.py

Generate random DNA sequences for testing validation scripts.
Simple utility that creates random sequences without any filtering.
"""

import argparse
import random
import os

# DNA bases
DNA_BASES = 'ATGC'

def generate_random_sequence(length):
    """Generate a single random DNA sequence of given length"""
    return ''.join(random.choice(DNA_BASES) for _ in range(length))

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
        seq = generate_random_sequence(chosen_length)
        sequences.append(seq)
    
    # Write to file
    with open(args.output, 'w') as f:
        for seq in sequences:
            f.write(seq + '\n')
    
    # Count sequences by length
    length_counts = {}
    for seq in sequences:
        length = len(seq)
        length_counts[length] = length_counts.get(length, 0) + 1
    
    length_breakdown = ", ".join([f"{count} at length {length}" for length, count in sorted(length_counts.items())])
    
    print(f"Generated {args.count} random DNA sequences with lengths: {length_breakdown}")
    print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main() 