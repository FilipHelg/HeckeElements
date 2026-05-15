#!/usr/bin/env python3
"""
Verify that nonzero triples are a subset of the old checker results.
"""
import csv
import sys

def validate_n(n):
    """Check that nonzero results for n are subset of old results."""
    nonzero_file = f'kl_triple_matches_nonzero_n{n}.csv'
    old_file = f'kl_triple_matches_n{n}.csv'
    
    # Load nonzero triples
    nonzero = set()
    try:
        with open(nonzero_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = (int(row['w']), int(row['x']), int(row['y']))
                nonzero.add(key)
    except FileNotFoundError:
        print(f"Error: {nonzero_file} not found")
        return False
    
    # Check if they're in the old results
    found = 0
    try:
        with open(old_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = (int(row['w']), int(row['x']), int(row['y']))
                if key in nonzero:
                    found += 1
    except FileNotFoundError:
        print(f"Error: {old_file} not found")
        return False
    
    is_subset = (found == len(nonzero))
    print(f"\nn={n}:")
    print(f"  Nonzero triples: {len(nonzero)}")
    print(f"  Found in old results: {found}")
    print(f"  Is subset: {is_subset}")
    
    if not is_subset:
        print(f"  WARNING: {len(nonzero) - found} nonzero triples are NOT in old results!")
        return False
    
    return True

if __name__ == "__main__":
    all_ok = True
    for n in [4, 5, 6]:
        if not validate_n(n):
            all_ok = False
    
    if all_ok:
        print("\n✓ All nonzero results are subsets of old checker results.")
        sys.exit(0)
    else:
        print("\n✗ Some nonzero results are not subsets!")
        sys.exit(1)
