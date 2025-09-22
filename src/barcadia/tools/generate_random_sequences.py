#!/usr/bin/env python3
"""
generate_random_sequences.py

Generate random DNA sequences for testing validation scripts (variable length sequences supported).

Example usage: python src/barcadia/tools/generate_random_sequences.py --count 10000 --lengths 12 13

Output: random DNA sequences (one per line as .txt)

Optional arguments:
--output: output file path (default: test/{count}_random_{min}to{max}bp_sequences.txt)

Required arguments:
--count: number of sequences to generate
--lengths: possible lengths for sequences
"""

import argparse
import os
import random
import sys

import numpy as np

# Add src directory to Python path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from config_utils import decode_sequence  # type: ignore


def generate_random_sequence(length):
    """Generate a single random DNA sequence of given length as integer array"""
    return np.random.randint(0, 4, size=length, dtype=np.int8)


def main():
    parser = argparse.ArgumentParser(
        description="Generate random DNA sequences for testing",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Example: python src/barcadia/tools/generate_random_sequences.py --count 10000 --lengths 12 13",
    )

    parser.add_argument("--count", type=int, required=True, help="Number of sequences to generate")
    parser.add_argument("--lengths", nargs="+", type=int, required=True, help="Possible lengths for sequences")
    parser.add_argument("--output", type=str, help="Output file path (default: test/{count}_random_{min}to{max}bp_sequences.txt)")

    args = parser.parse_args()

    # Validate arguments
    if args.count <= 0:
        raise ValueError("Count must be > 0")

    for length in args.lengths:
        if length <= 0:
            raise ValueError(f"All lengths must be > 0, got {length}")

    # Check if count exceeds maximum possible sequences across all lengths
    total_max_possible = sum(4**length for length in args.lengths)
    if args.count > total_max_possible:
        raise ValueError(f"Count ({args.count}) exceeds maximum possible sequences for lengths {args.lengths} ({total_max_possible})")

    # Generate default output path if not specified
    if args.output is None:
        os.makedirs("test", exist_ok=True)
        min_length = min(args.lengths)
        max_length = max(args.lengths)

        # Handle single length vs range
        length_range = str(min_length) if min_length == max_length else f"{min_length}to{max_length}"

        args.output = f"test/{args.count}_random_{length_range}bp_sequences.txt"

    # Generate sequences, write to file, and count lengths in one loop
    length_counts = {}
    with open(args.output, "w") as f:
        for _ in range(args.count):
            # Randomly choose a length from the provided options
            chosen_length = random.choice(args.lengths)
            seq_array = generate_random_sequence(chosen_length)

            # Write to file (convert to DNA string for output)
            dna_string = decode_sequence(seq_array)
            f.write(dna_string + "\n")

            # Count sequences by length
            length_counts[chosen_length] = length_counts.get(chosen_length, 0) + 1

    length_breakdown = ", ".join([f"{count} at length {length}" for length, count in sorted(length_counts.items())])

    print(f"Generated {args.count} random DNA sequences with lengths: {length_breakdown}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
