#!/usr/bin/env python3
"""
memory_benchmark.py

Simple memory benchmarking utility for tracking maximum memory usage of processes.

Example usage: python src/utils/memory_benchmark.py src/{script_name}.py --args

Output: memory usage report shown in terminal and .log file
"""

import psutil
import subprocess
import sys
import time
import os
import logging
from datetime import datetime


def benchmark_script(script_path, *args):
    """
    Benchmark memory usage of a Python script.
    
    Args:
        script_path: Path to the Python script to run
        *args: Arguments to pass to the script
    
    Returns:
        dict: Benchmark results
    """
    script_name = os.path.basename(script_path)
    
    # Run the script as a subprocess and track its memory
    cmd = [sys.executable, script_path] + list(args)
    start_time = time.time()
    
    # Start the subprocess - capture output to find log file location
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
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
        stdout, stderr = process.communicate()
        return_code = process.returncode
        
    except Exception as e:
        logging.error(f"Error tracking memory: {e}")
        process.terminate()
        return None
    
    duration = time.time() - start_time
    peak_memory_mb = peak_memory / 1024 / 1024
    
    if return_code != 0:
        logging.error(f"Script failed with return code {return_code}")
        logging.error("Check the script's own log file for details")
        return None
    
    return {
        'script': script_name,
        'duration_seconds': duration,
        'peak_memory_mb': peak_memory_mb,
        'return_code': return_code,
        'stdout': stdout,
        'stderr': stderr
    }


def main():
    """Command-line interface for benchmarking scripts."""
    if len(sys.argv) < 2:
        print("Usage: python src/utils/memory_benchmark.py <script_path> [script_args...]")
        print("Example: python src/utils/memory_benchmark.py src/validate_barcodes.py --input test/barcodes.txt")
        sys.exit(1)
    
    script_path = sys.argv[1]
    script_args = sys.argv[2:]
    
    # Setup logging
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"memory_benchmark_{timestamp}.log"
    log_filepath = os.path.join("test", log_filename)
    os.makedirs("test", exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S',
        handlers=[
            logging.FileHandler(log_filepath),
            logging.StreamHandler()
        ]
    )
    
    # Log benchmark start
    logging.info(f"Starting memory benchmark for: {script_path}")
    logging.info(f"Arguments: {' '.join(script_args)}")
    logging.info("-" * 50)
    
    # Run benchmark
    result = benchmark_script(script_path, *script_args)
    
    if result:
        # Log benchmark results only
        logging.info("-" * 50)
        logging.info("BENCHMARK RESULTS:")
        logging.info(f"Script: {result['script']}")
        logging.info(f"Duration: {result['duration_seconds']:.2f} seconds")
        logging.info(f"Peak Memory: {result['peak_memory_mb']:.2f} MB")
        logging.info(f"Return Code: {result['return_code']}")
        logging.info("-" * 50)
        
        # Extract log file location from script output
        script_log_file = None
        
        # Check both stdout and stderr for log file information
        for output_stream in [result['stdout'], result['stderr']]:
            if output_stream:
                for line in output_stream.split('\n'):
                    if 'Log file:' in line:
                        script_log_file = line.split('Log file:')[-1].strip()
                        break
                if script_log_file:
                    break
        
        if script_log_file:
            logging.info(f"Script log file: {script_log_file}")
        else:
            logging.info("Script log file: Not found in output")
        
        logging.info(f"Benchmark log file: {log_filepath}")
    else:
        logging.error("Benchmark failed!")


if __name__ == "__main__":
    main() 