#!/usr/bin/env python3
"""
Debug script to find mismatches in n=6 nonzero vs old checker data.
"""
import csv

n = 6
nonzero_file = f'kl_triple_matches_nonzero_n{n}.csv'
old_file = f'kl_triple_matches_n{n}.csv'

# Load nonzero triples
nonzero = set()
nonzero_rows = {}
try:
    with open(nonzero_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (int(row['w']), int(row['x']), int(row['y']))
            nonzero.add(key)
            nonzero_rows[key] = row
except FileNotFoundError as e:
    print(f"Error: {nonzero_file} not found: {e}")
    exit(1)

print(f"Loaded {len(nonzero)} nonzero triples from {nonzero_file}\n")

# Load old triples
old = set()
old_rows = {}
try:
    with open(old_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (int(row['w']), int(row['x']), int(row['y']))
            old.add(key)
            old_rows[key] = row
except FileNotFoundError as e:
    print(f"Error: {old_file} not found: {e}")
    exit(1)

print(f"Loaded {len(old)} old checker triples from {old_file}\n")

# Check for mismatches
found = nonzero & old
missing_from_old = nonzero - old

if missing_from_old:
    print(f"ERROR: Found {len(missing_from_old)} nonzero triples NOT in old results!")
    print("\nNonzero triples not in old results:")
    for key in sorted(missing_from_old):
        w, x, y = key
        nonzero_row = nonzero_rows[key]
        print(f"\n  Triple ({w}, {x}, {y})")
        print(f"    nonzero: {nonzero_row}")
        
        # Check if it's close in the old file
        close_matches = [k for k in old if k[0] == w and k[1] == x]
        if close_matches:
            print(f"    Similar triples with same (w, x) in old:")
            for k in close_matches:
                print(f"      {k}: {old_rows[k]}")
else:
    print(f"✓ All {len(nonzero)} nonzero triples are in old results.")

print(f"\nSummary:")
print(f"  Nonzero triples: {len(nonzero)}")
print(f"  Found in old: {len(found)}")
print(f"  Missing from old: {len(missing_from_old)}")
