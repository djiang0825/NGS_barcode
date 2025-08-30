#!/usr/bin/env python3
"""
config_utils.py

Configuration and DNA-encoding/decoding utility functions for efficient barcode generation and validation.
"""

import os
import logging
import numpy as np
from datetime import datetime

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

def setup_logging(args, script_name):
    """Setup logging and create output directory. Returns log filepath."""
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup logging with file output
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{script_name}_{timestamp}.log"
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
    
    return log_filepath

def read_files(file_paths):
    """
    Read DNA sequences from files and convert to integer arrays.
    
    Args:
        file_paths: List of file paths or single file path
    
    Returns:
        tuple: (sequences, length_counts) where sequences are integer arrays
    """
    # Convert single file to list for consistent handling
    if isinstance(file_paths, str):
        file_paths = [file_paths]
    
    sequences = []
    length_counts = {}
    
    for file_path in file_paths:
        file_count = 0
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                seq = line.strip()
                if not seq:  # Skip empty lines
                    continue
                
                # Basic validation
                if not all(base in DNA_BASES for base in seq):
                    logging.warning(f"File {file_path}, line {line_num}: Invalid DNA sequence '{seq}', skipping")
                    continue
                
                # Convert to integer array for efficient processing
                seq_array = encode_sequence(seq)
                sequences.append(seq_array)
                
                # Count length while reading
                length = len(seq_array)
                length_counts[length] = length_counts.get(length, 0) + 1
                
                file_count += 1
        
        logging.info(f"Loaded {file_count} sequences from {file_path}")
    
    if not sequences:
        raise ValueError(f"File(s) are empty: {', '.join(file_paths)}")
    
    # Generate length info for logging
    if len(length_counts) == 1:
        length_info = f"length {list(length_counts.keys())[0]}"
    else:
        length_breakdown = ", ".join([f"{count} at length {length}" for length, count in sorted(length_counts.items())])
        length_info = f"mixed lengths: {length_breakdown}"
    
    logging.info(f"Total loaded: {len(sequences)} sequences from {len(file_paths)} file(s) ({length_info})")
    
    return sequences, length_counts
