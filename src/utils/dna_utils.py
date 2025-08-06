#!/usr/bin/env python3
"""
dna_utils.py

DNA encoding and decoding utility functions for barcode generation and validation.
"""

import numpy as np

# DNA encoding constants
DNA_BASES = 'ATGC'
DNA_TO_INT = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
INT_TO_DNA = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}

def encode_sequence(dna_string):
    """Convert DNA string to integer array"""
    return np.array([DNA_TO_INT[base] for base in dna_string], dtype=np.int8)

def decode_sequence(seq_array):
    """Convert integer array back to DNA string"""
    return ''.join(INT_TO_DNA[base] for base in seq_array)

def validate_arguments(args):
    """Validate command line arguments and raise ValueError if invalid"""
    if args.gc_min < 0 or args.gc_max > 1 or args.gc_min >= args.gc_max:
        raise ValueError("GC content bounds must be: 0 ≤ gc_min < gc_max ≤ 1")
    
    if args.homopolymer_max < 1:
        raise ValueError("Maximum homopolymer length must be ≥ 1")
    
    if args.min_distance < 1:
        raise ValueError("Minimum distance must be ≥ 1")
    
    # Optional length validation (only if length attribute exists)
    if hasattr(args, 'length'):
        if args.homopolymer_max >= args.length:
            raise ValueError("Maximum homopolymer length must be < sequence length")
        
        if args.min_distance >= args.length:
            raise ValueError("Minimum distance must be < sequence length")
        
        if args.length <= 0:
            raise ValueError("Length must be > 0")
    
    # Optional count validation (only if count attribute exists)
    if hasattr(args, 'count') and args.count <= 0:
        raise ValueError("Count must be > 0") 