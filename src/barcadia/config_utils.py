#!/usr/bin/env python3
"""
config_utils.py

Configuration and DNA-encoding/decoding utility functions for efficient barcode generation and validation.
"""

import logging
import os
from datetime import datetime

import numpy as np

# DNA encoding constants
DNA_BASES = "ATGC"
DNA_TO_INT = {"A": 0, "T": 1, "G": 2, "C": 3}
INT_TO_DNA = {0: "A", 1: "T", 2: "G", 3: "C"}


def encode_sequence(dna_string):
    """Convert DNA string to integer array"""
    return np.array([DNA_TO_INT[base] for base in dna_string], dtype=np.int8)


def decode_sequence(seq_array):
    """Convert integer array back to DNA string"""
    return "".join(INT_TO_DNA[base] for base in seq_array)


def setup_logging(args, script_name):
    """Setup logging and create output directory. Returns log filepath."""
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Setup logging with file output
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{script_name}_{timestamp}.log"
    log_filepath = os.path.join(args.output_dir, log_filename)

    # Configure logging to both file and console
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%H:%M:%S", handlers=[logging.FileHandler(log_filepath), logging.StreamHandler()])

    return log_filepath


class ExistingSequenceSet:
    """
    A class to manage existing DNA sequence sets with file operations and validation.

    This class consolidates file reading, existence checking, and sequence management
    for both generation and validation scripts.
    """

    def __init__(self, sequences=None, length_counts=None):
        """
        Initialize the sequence set.

        Args:
            sequences: List of integer arrays (encoded DNA sequences)
            length_counts: Dictionary mapping length to count
        """
        self.sequences = sequences or []
        self.length_counts = length_counts or {}

    def _read_files(self, file_paths):
        """
        Internal method to read DNA sequences from files and convert to integer arrays.
        Handles file existence checking and path normalization internally.

        Args:
            file_paths: List of file paths or single file path

        Returns:
            tuple: (sequences, length_counts) where sequences are integer arrays

        Raises:
            ValueError: If any file does not exist or files are empty
        """
        # Normalize file paths (convert single file to list)
        if isinstance(file_paths, str):
            file_paths = [file_paths]

        # Check that all files exist
        for file_path in file_paths:
            if not os.path.exists(file_path):
                raise ValueError(f"File does not exist: {file_path}")

        sequences = []
        length_counts = {}

        for file_path in file_paths:
            file_count = 0
            with open(file_path) as f:
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

    @classmethod
    def from_files(cls, file_paths):
        """
        Create ExistingSequenceSet from files (used by both validation and generation scripts).

        Args:
            file_paths: List of file paths or single file path

        Returns:
            ExistingSequenceSet: Instance with loaded sequences and length counts
        """
        instance = cls()
        sequences, length_counts = instance._read_files(file_paths)
        instance.sequences = sequences
        instance.length_counts = length_counts
        return instance

    @classmethod
    def from_input_files(cls, file_paths):
        """
        Create ExistingSequenceSet from input files (used by validation script).

        Args:
            file_paths: List of file paths or single file path

        Returns:
            ExistingSequenceSet: Instance with loaded sequences and length counts
        """
        return cls.from_files(file_paths)

    @classmethod
    def from_unpaired_seeds(cls, file_paths):
        """
        Create ExistingSequenceSet from unpaired seed files (used by generation script).

        Args:
            file_paths: List of file paths or single file path

        Returns:
            ExistingSequenceSet: Instance with loaded sequences and length counts
        """
        return cls.from_files(file_paths)

    @classmethod
    def from_paired_seeds(cls, file1, file2):
        """
        Create ExistingSequenceSet from paired seed files (used by generation script).

        Args:
            file1: Path to first paired seed file
            file2: Path to second paired seed file

        Returns:
            ExistingSequenceSet: Instance with combined sequences and length counts
        """
        instance = cls()

        # Load paired seeds separately
        paired_seed1_pool, seed1_length_counts = instance._read_files([file1])
        paired_seed2_pool, seed2_length_counts = instance._read_files([file2])

        # Validate paired seeds
        # 1. Check that both files have the same number of sequences
        if len(paired_seed1_pool) != len(paired_seed2_pool):
            raise ValueError(f"Paired seed files must have the same number of sequences. Seed1: {len(paired_seed1_pool)} sequences, Seed2: {len(paired_seed2_pool)} sequences")

        # 2. Check that both files have sequences of the same length within the file
        elif len(seed1_length_counts) != 1:
            raise ValueError(f"All sequences in paired seed file 1 must be the same length. Found lengths: {sorted(seed1_length_counts.keys())}")
        elif len(seed2_length_counts) != 1:
            raise ValueError(f"All sequences in paired seed file 2 must be the same length. Found lengths: {sorted(seed2_length_counts.keys())}")

        # 3. Check that both files have sequences of the same length between the files
        elif list(seed1_length_counts.keys())[0] != list(seed2_length_counts.keys())[0]:
            raise ValueError(f"Paired seed files must have sequences of the same length. Seed1 length: {list(seed1_length_counts.keys())[0]}, Seed2 length: {list(seed2_length_counts.keys())[0]}")
        else:
            # All validations passed - combine both for generation pool
            combined_sequences = paired_seed1_pool + paired_seed2_pool

            # Since paired seeds are validated to have the same length, just use seed1's length counts
            # and double the count since we have two files
            combined_length_counts = {}
            for length, count in seed1_length_counts.items():
                combined_length_counts[length] = count * 2

            instance.sequences = combined_sequences
            instance.length_counts = combined_length_counts

        return instance
