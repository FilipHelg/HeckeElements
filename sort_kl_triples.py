#!/usr/bin/env python3
"""
Sort kl_triple_matches CSV files by the first three columns (w, x, y).
Priority order: w > x > y (all in increasing order).
"""

import csv
import sys
from pathlib import Path


def sort_kl_triples(input_file):
    """Sort a kl_triple_matches CSV file and save to a new file."""
    
    input_path = Path(input_file)
    if not input_path.exists():
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)
    
    # Read the CSV file
    rows = []
    with open(input_path, 'r') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        for row in reader:
            rows.append(row)
    
    # Sort by w, x, y (all as integers)
    rows.sort(key=lambda row: (int(row['w']), int(row['x']), int(row['y'])))
    
    # Generate output filename
    output_path = input_path.parent / f"{input_path.stem}_sorted.csv"
    
    # Write sorted data to new file
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"Sorted CSV written to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sort_kl_triples.py <input_csv>")
        sys.exit(1)
    
    sort_kl_triples(sys.argv[1])
