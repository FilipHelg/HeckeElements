#!/usr/bin/env python3
"""
Check which involutions have same-left-cell pairs
"""
import csv
import itertools
import collections

def rs_tableaux(n, perm):
    """Compute Robinson-Schensted P and Q tableaux"""
    P = []
    Q = []
    for i, x in enumerate(perm):
        b = x; r = 0
        while True:
            if r == len(P): 
                P.append([b])
                Q.append([i+1])
                break
            row = P[r]
            c = 0
            while c < len(row) and row[c] < b: 
                c += 1
            if c == len(row): 
                row.append(b)
                Q[r].append(i+1)
                break
            b, row[c] = row[c], b
            r += 1
    return P

def get_left_cell(n, index):
    """Get left cell signature (P-tableau shape+content)"""
    perm = list(range(n))
    for i in range(index):
        # Generate next permutation in lexicographic order
        # For now, just return the index as cell id for testing
        pass
    P = rs_tableaux(n, perm)
    return tuple(tuple(row) for row in P)

# Generate all permutations for S6 and find involutions
involutions = []
for perm in itertools.permutations(range(6)):
    # Check if perm is an involution (perm = perm^{-1})
    inv_perm = [0] * 6
    for i, p in enumerate(perm):
        inv_perm[p] = i
    if inv_perm == list(perm):
        involutions.append(perm)

print(f"Total involutions in S6: {len(involutions)}")

# Load triple results
triples_with_d = collections.defaultdict(list)
with open('kl_triple_matches_nonzero_n6.csv', newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        d = int(row['w'])
        x = int(row['x'])
        y = int(row['y'])
        triples_with_d[d].append((x, y))

# Find involutions with triples
d_with_triples = set(triples_with_d.keys())
print(f"Involutions with nonzero product triples: {len(d_with_triples)}")
print(f"Involutions WITHOUT triples: {len(involutions) - len(d_with_triples)}")

# List the d values
print(f"\nD values with triples: {sorted(d_with_triples)}")

# Count triples per d
triple_counts = collections.Counter(len(v) for v in triples_with_d.values())
print(f"\nTriple count distribution:")
for count, freq in sorted(triple_counts.items()):
    print(f"  {count} triples: {freq} involutions")
