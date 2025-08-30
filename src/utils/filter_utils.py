#!/usr/bin/env python3
"""
filter_utils.py

Filter-related utility functions with Numba JIT compilation for efficient barcode generation and validation.
"""

from numba import jit
import numpy as np

# Simple validation of filter arguments used in both generation and validation
def validate_filter_arguments(args):
    """Validate filter-related command line arguments and raise ValueError if invalid"""
    if args.gc_min < 0 or args.gc_max > 1 or args.gc_min >= args.gc_max:
        raise ValueError("GC content bounds must be: 0 ≤ gc_min < gc_max ≤ 1")
    
    elif args.homopolymer_max < 1:
        raise ValueError("Maximum homopolymer repeat length must be ≥ 1")
    
    elif args.min_distance < 1:
        raise ValueError("Minimum edit distance must be ≥ 1")

# Biological filter functions
@jit(nopython=True, cache=True)
def check_gc_content_int(seq_array, gc_min, gc_max):
    """Check if sequence passes GC content filter (works with integer arrays)"""
    # G=2, C=3 in our encoding - count them directly
    gc_count = 0
    for base in seq_array:
        if base == 2 or base == 3:  # G or C
            gc_count += 1
    gc_content = gc_count / len(seq_array)
    return gc_min <= gc_content <= gc_max

@jit(nopython=True, cache=True)
def check_homopolymer_int(seq_array, homopolymer_max):
    """Check for homopolymer repeats longer than homopolymer_max (works with integer arrays)"""
    current_base = seq_array[0]
    current_count = 1
    
    for base in seq_array[1:]:
        if base == current_base:
            current_count += 1
            if current_count > homopolymer_max:
                return False  # Fails check
        else:
            current_base = base
            current_count = 1
    
    return True  # Passes check

# Distance calculation functions
@jit(nopython=True, cache=True)
def hamming_distance_int(seq1, seq2, min_distance):
    """Calculate Hamming distance with early stopping (assumes equal-length sequences, works with integer arrays)"""
    distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1
            if distance >= min_distance:
                return distance  # Early stopping
    return distance

@jit(nopython=True, cache=True)
def levenshtein_distance_int(seq1, seq2, min_distance):
    """Calculate Levenshtein distance with early stopping (assumes mixed-length sequences, works with integer arrays)"""
    if len(seq1) < len(seq2):
        return levenshtein_distance_int(seq2, seq1, min_distance)
    
    elif len(seq2) == 0:
        return len(seq1)
    
    # Use numpy arrays for better performance with numba
    previous_row = np.arange(len(seq2) + 1, dtype=np.int32)
    
    # Early stopping: if initial row already exceeds min_distance, return early
    if previous_row.min() >= min_distance:
        return min_distance
    
    for i in range(len(seq1)):
        current_row = np.zeros(len(seq2) + 1, dtype=np.int32)
        current_row[0] = i + 1
        for j in range(len(seq2)):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (seq1[i] != seq2[j])
            current_row[j + 1] = min(insertions, deletions, substitutions)
        
        # Early stopping: if minimum value in current row >= min_distance,
        # the final distance will be >= min_distance
        if current_row.min() >= min_distance:
            return min_distance
            
        previous_row = current_row
    
    return previous_row[-1]

def calculate_distance(seq1, seq2, min_distance):
    """Calculate distance between two sequences, using Hamming for equal length, Levenshtein otherwise"""
    if len(seq1) == len(seq2):
        return hamming_distance_int(seq1, seq2, min_distance)
    else:
        return levenshtein_distance_int(seq1, seq2, min_distance)

def select_distance_method(target_count, min_distance, has_mixed_lengths):
    """
    Determine which distance checking method to use based on barcode set characteristics.
    Returns: "pairwise" or "neighbor_enumeration"
    
    Rules:
    1. Small barcode sets (<10K sequences counting seeds if seeds are present): Always use pairwise
    2. Mixed-length sequences (within seeds and/or between seeds and new barcodes): Always use pairwise
    3. Minimum distance threshold: For min_distance <= 4, use neighbor enumeration, otherwise pairwise
    """
    # Rule 1: Small barcode sets always use pairwise
    if target_count < 10000:
        return "pairwise"
    
    # Rule 2: Mixed-length sequences always use pairwise
    elif has_mixed_lengths:
        return "pairwise"
    
    # Rule 3: For large minimum distance always use pairwise
    elif min_distance <= 4:
        return "neighbor_enumeration"
    else:
        return "pairwise"

def generate_hamming_neighbors(seq_array, max_distance, current_distance=0):
    """Generate all Hamming neighbors within max_distance of a sequence"""
    if current_distance == max_distance:
        yield tuple(seq_array)
        return
    
    # Yield current sequence if distance > 0
    if current_distance > 0:
        yield tuple(seq_array)
    
    # Generate neighbors by substitution
    for i in range(len(seq_array)):
        original_base = seq_array[i]
        for new_base in [0, 1, 2, 3]:  # A, T, G, C
            if new_base != original_base:
                seq_array[i] = new_base
                yield from generate_hamming_neighbors(seq_array, max_distance, current_distance + 1)
        seq_array[i] = original_base  # backtrack
