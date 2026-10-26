#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def extract_third_line(directory, source_filename, output_filename):
    """
    Walk through directory and subdirectories, find files with source_filename,
    extract the 3rd line, and save it to output_filename in the same directory.
    """
    for root, dirs, files in os.walk(directory):
        if source_filename in files:
            source_path = os.path.join(root, source_filename)
            output_path = os.path.join(root, output_filename)
            
            try:
                with open(source_path, 'r') as f:
                    lines = f.readlines()
                
                if len(lines) >= 3:
                    third_line = lines[2]  # 0-indexed, so 3rd line is index 2
                    with open(output_path, 'w') as f:
                        f.write(third_line)
                    print(f"Extracted 3rd line from {source_path} -> {output_path}")
                else:
                    print(f"Warning: {source_path} has fewer than 3 lines, skipping.")
            except Exception as e:
                print(f"Error processing {source_path}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract the 3rd line from files and save to a new file in the same directory"
    )
    parser.add_argument(
        "directory",
        help="Root directory to search"
    )
    parser.add_argument(
        "-s", "--source",
        required=True,
        help="Name of the source file to read from"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Name of the output file to write to"
    )
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a valid directory")
        exit(1)
    
    extract_third_line(args.directory, args.source, args.output)
