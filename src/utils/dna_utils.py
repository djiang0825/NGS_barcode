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