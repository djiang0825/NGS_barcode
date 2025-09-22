#!/usr/bin/env python3
"""
memory_benchmark.py

Simple memory benchmarking utility for tracking maximum memory usage of barcadia commands.

Example usage:
  python src/barcadia/tools/memory_benchmark.py barcadia generate --args
  python src/barcadia/tools/memory_benchmark.py barcadia validate --args

Output: memory usage report shown in terminal and memory_benchmark_{timestamp}.log file

Optional arguments:
--mem-output-dir: output directory for benchmark logs (default: test)

Required arguments:
command: the barcadia command to benchmark (e.g., generate or validate)
"""

import argparse
import logging
import os
import subprocess
import sys
import time

import psutil

# Add src directory to Python path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from config_utils import setup_logging  # type: ignore


def benchmark_command(command):
    """
    Benchmark memory usage of a barcadia command.

    Args:
        command: The command to run (e.g., ['barcadia', 'generate', '--count', '1000'])

    Returns:
        dict: Benchmark results
    """
    # Run the command as a subprocess and track its memory
    start_time = time.time()

    # Start the subprocess; let it inherit stdout/stderr so it can print directly to terminal
    process = subprocess.Popen(command)

    # Track memory usage of the subprocess
    peak_memory = 0
    try:
        while process.poll() is None:  # While process is still running
            try:
                # Get memory usage of the subprocess
                child_process = psutil.Process(process.pid)
                memory_bytes = child_process.memory_info().rss
                peak_memory = max(peak_memory, memory_bytes)
                time.sleep(0.1)  # Sample every 100ms
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                # Process might have finished or we don't have access
                break

        # Wait for process to complete
        process.wait()
        return_code = process.returncode
        stdout, stderr = None, None

    except Exception as e:
        logging.error(f"Error tracking memory: {e}")
        process.terminate()
        return None

    duration = time.time() - start_time
    peak_memory_mb = peak_memory / 1024 / 1024

    if return_code != 0:
        logging.error(f"Command failed with return code {return_code}")
        logging.error("Check the command's own log file for details")
        return None

    return {"command": " ".join(command), "duration_seconds": duration, "peak_memory_mb": peak_memory_mb, "return_code": return_code, "stdout": stdout, "stderr": stderr}


def main():
    """Command-line interface for benchmarking barcadia commands."""
    parser = argparse.ArgumentParser(
        description="Benchmark memory usage of barcadia commands",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n  python src/barcadia/tools/memory_benchmark.py barcadia validate --args\n  python src/barcadia/tools/memory_benchmark.py barcadia generate --args",
    )

    parser.add_argument("--mem-output-dir", type=str, default="test", help="Output directory for benchmark logs (default: test)")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="The barcadia command to benchmark (e.g., generate or validate)")

    args = parser.parse_args()

    if not args.command:
        parser.error("Please provide a barcadia command to benchmark")

    # Setup logging
    args.output_dir = args.mem_output_dir
    log_filepath = setup_logging(args, "memory_benchmark")

    # Log benchmark start
    command_str = " ".join(args.command)
    logging.info(f"Starting memory benchmark for: {command_str}")
    logging.info("-" * 50)

    # Run benchmark
    result = benchmark_command(args.command)

    if result:
        # Log benchmark results only
        logging.info("-" * 50)
        logging.info("BENCHMARK RESULTS:")
        logging.info(f"Command: {result['command']}")
        logging.info(f"Duration: {result['duration_seconds']:.2f} seconds")
        logging.info(f"Peak Memory: {result['peak_memory_mb']:.2f} MB")
        logging.info(f"Return Code: {result['return_code']}")
        logging.info("-" * 50)
        # Note: the benchmarked command prints its own log file path to the terminal.
        logging.info(f"Benchmark log file: {log_filepath}")
    else:
        logging.error("Benchmark failed!")


if __name__ == "__main__":
    main()
