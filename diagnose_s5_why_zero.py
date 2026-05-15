#!/usr/bin/env python3
"""
Diagnose why S5 produces zero triples.
Check which permutations are involutions,
and manually analyze products for these involutions.
"""

import struct
import sys

def read_dual_cache(path):
    """Load involutions from dual KL cache."""
    with open(path, 'rb') as f:
        magic = f.read(8)
        if magic != b'DKBULK1\x00':
            print(f"Warning: unexpected magic {magic}")
        
        version, n, count, inv_count = struct.unpack('<IIII', f.read(16))
        print(f"Dual cache: n={n}, count={count}, involution_count={inv_count}")
        
        # Read present flags
        present = []
        for i in range(count):
            present.append(f.read(1)[0])
        
        involutions = [i for i in range(count) if present[i]]
        print(f"Found {len(involutions)} involutions: {involutions}")
        return involutions

def perm_to_cycles(perm_list):
    """Convert permutation list to cycle notation."""
    n = len(perm_list)
    visited = [False] * n
    cycles = []
    
    for i in range(n):
        if not visited[i]:
            cycle = []
            j = i
            while not visited[j]:
                visited[j] = True
                cycle.append(j + 1)  # 1-indexed
                j = perm_list[j]
            if len(cycle) > 1:
                cycles.append(cycle)
    
    return cycles

def perm_inverse(perm):
    """Compute inverse of a permutation."""
    return [perm.index(i) for i in range(len(perm))]

def is_involution(perm):
    """Check if permutation is its own inverse."""
    return perm == perm_inverse(perm)

# For n=5, load Sn data
print("=== S5 Involution Analysis ===\n")

# Try to read dual cache
try:
    involutions = read_dual_cache('dual_kl_bulk_n5.bin')
except Exception as e:
    print(f"Error reading dual cache: {e}")
    sys.exit(1)

# Try to load S5.txt to get permutations
try:
    perms = []
    with open('S5.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                perm_str = line.split(',')[0] if ',' in line else line
                perm = [int(x) for x in perm_str.split()]
                perms.append(perm)
    
    print(f"\nLoaded {len(perms)} permutations from S5.txt")
    
    # Check which loaded perms are involutions
    print("\nInvolutions in loaded S5 data:")
    actual_involutions = []
    for i, perm in enumerate(perms):
        if is_involution(perm):
            actual_involutions.append(i)
            cycles = perm_to_cycles(perm)
            cycle_str = ' '.join([str(c) for c in cycles]) if cycles else '(1)(2)(3)(4)(5)'
            print(f"  w={i}: {perm} = {cycle_str}")
    
    print(f"\nTotal involutions in S5.txt: {len(actual_involutions)}")
    
    # Compare
    loaded_set = set(actual_involutions)
    cache_set = set(involutions)
    
    print(f"\nDual cache has {len(cache_set)} involutions")
    print(f"S5.txt has {len(loaded_set)} involutions")
    
    if loaded_set == cache_set:
        print("✓ Sets match!")
    else:
        print(f"✗ In cache but not S5.txt: {cache_set - loaded_set}")
        print(f"✗ In S5.txt but not cache: {loaded_set - cache_set}")

except Exception as e:
    print(f"Error loading S5.txt: {e}")
    sys.exit(1)

print("\n=== Analysis ===")
print(f"Theory predicts {len(actual_involutions)} involutions should have nonzero Kåhrström triples.")
print(f"Current code finds 0 nonzero triples.")
print("\nPossible issues:")
print("1. Theory expectation is wrong")
print("2. Involutions don't actually satisfy condition for S5")
print("3. Product computation is wrong (but S4 and S6 work)")
print("4. Left-cell detection is wrong")
