#!/usr/bin/env python3
"""
filter_utils.py

Core filter utility functions with Numba JIT compilation for efficient DNA barcode processing.
"""

from numba import jit
import numpy as np

# Distance calculation functions
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

@jit(nopython=True)
def levenshtein_distance(seq1, seq2, min_distance):
    """Calculate Levenshtein distance with early stopping (works with integer arrays)"""
    if len(seq1) < len(seq2):
        return levenshtein_distance(seq2, seq1, min_distance)
    
    if len(seq2) == 0:
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
        return levenshtein_distance(seq1, seq2, min_distance)

# Biological filter functions
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
    """Check for homopolymer repeats longer than max_length (works with integer arrays)"""
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

 