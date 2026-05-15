#!/usr/bin/env python3
"""plot_heatmap_combined.py

Generate a heatmap of pairs (x,y) that satisfy BOTH:
1. Are left-cell equivalent (same P-tableau)
2. Are right-cell equivalent to elements in the right KL pre-order closure of d

Usage: python plot_heatmap_combined.py <n> <d>

Produces `heatmap_combined_n<n>_d<d>.png` in CWD.
"""
import sys
import math
import subprocess
import csv
import tempfile
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import pandas as pd
import numpy as np
from collections import defaultdict

sns.set(style="whitegrid")


def index_to_perm(n, index):
    """Convert Lehmer index to permutation in lexicographic order."""
    remaining = list(range(n))
    perm = []
    for k in range(n - 1, -1, -1):
        fact = math.factorial(k)
        digit = index // fact if k > 0 else 0
        index = index % fact if k > 0 else 0
        perm.append(remaining.pop(digit))
    return perm


def rs_tableaux_both(n, perm):
    """Compute Robinson-Schensted P and Q tableaux from a permutation."""
    P = []
    Q = []
    for i, x in enumerate(perm):
        b = x
        r = 0
        while True:
            if r == len(P):
                P.append([b])
                Q.append([i + 1])
                break
            row = P[r]
            c = 0
            while c < len(row) and row[c] < b:
                c += 1
            if c == len(row):
                row.append(b)
                Q[r].append(i + 1)
                break
            b, row[c] = row[c], b
            r += 1
    return P, Q


def left_cell_key_from_index(n, index):
    """Get the P-tableau (left cell signature) from a permutation index."""
    perm = index_to_perm(n, index)
    P, Q = rs_tableaux_both(n, perm)
    # Return P-tableau as the key for left cells
    return tuple(tuple(row) for row in P)


def right_cell_key_from_index(n, index):
    """Get the Q-tableau (right cell signature) from a permutation index."""
    perm = index_to_perm(n, index)
    P, Q = rs_tableaux_both(n, perm)
    # Return Q-tableau as the key for right cells
    return tuple(tuple(row) for row in Q)


def get_right_preorder_closure(n, d):
    """
    Call right_preorder_involutions to get all involutions d' <= d.
    Returns a set of involution indices.
    """
    try:
        result = subprocess.run(
            ['./right_preorder_involutions', str(n), str(d)],
            capture_output=True,
            text=True,
            timeout=300
        )
        if result.returncode != 0:
            print(f"Error calling right_preorder_involutions: {result.stderr}", file=sys.stderr)
            sys.exit(1)
        
        involutions = set()
        for line in result.stdout.strip().split('\n'):
            if line.startswith('#') or line.startswith('index'):
                continue
            if not line.strip():
                continue
            parts = line.split(',')
            if len(parts) >= 1:
                try:
                    inv_index = int(parts[0])
                    involutions.add(inv_index)
                except ValueError:
                    pass
        
        print(f"Found {len(involutions)} involutions d' <= d={d}", file=sys.stderr)
        return involutions
    
    except FileNotFoundError:
        print("Error: right_preorder_involutions not found. Please compile it first.", file=sys.stderr)
        sys.exit(1)


def build_cell_maps(n):
    """Build maps from element index to left and right cell keys."""
    total = math.factorial(n)
    left_cell_map = {}
    right_cell_map = {}
    for idx in range(total):
        left_cell_map[idx] = left_cell_key_from_index(n, idx)
        right_cell_map[idx] = right_cell_key_from_index(n, idx)
    return left_cell_map, right_cell_map


def main(n, d):
    if n < 1:
        raise SystemExit("n must be >= 1")
    
    total = math.factorial(n)
    
    if d < 0 or d >= total:
        raise SystemExit(f"d must be in range [0, {total-1}]")
    
    print(f"Computing right preorder closure for n={n}, d={d}...", file=sys.stderr)
    involution_closure = get_right_preorder_closure(n, d)
    
    print(f"Building left and right cell maps for S_{n}...", file=sys.stderr)
    left_cell_map, right_cell_map = build_cell_maps(n)
    
    # Collect right cells from elements in the closure involutions
    print(f"Processing elements in the involution closure...", file=sys.stderr)
    right_cells_from_closure = set()
    for inv_idx in involution_closure:
        right_cell = right_cell_map[inv_idx]
        right_cells_from_closure.add(right_cell)
    
    print(f"Unique right cells in closure: {len(right_cells_from_closure)}", file=sys.stderr)
    
    # Identify elements that are right-cell-equivalent to closure elements
    elements_in_right_closure = set()
    for idx in range(total):
        if right_cell_map[idx] in right_cells_from_closure:
            elements_in_right_closure.add(idx)
    
    print(f"Elements with right cells from closure: {len(elements_in_right_closure)} / {total}", file=sys.stderr)
    
    # Build the pair matrix
    # A pair (i, j) is colored black if:
    # 1. i and j are left-cell equivalent, AND
    # 2. Both i and j are right-cell equivalent to closure elements
    
    print(f"Building combined pair matrix...", file=sys.stderr)
    combined = pd.DataFrame(
        [[1 if (left_cell_map[i] == left_cell_map[j] and 
                i in elements_in_right_closure and 
                j in elements_in_right_closure) else 0
          for j in range(total)] for i in range(total)],
        index=range(total),
        columns=range(total),
    )
    
    # Count satisfied pairs
    black_pixels = int(combined.values.sum())
    print(f"Black pixels (satisfying both criteria): {black_pixels} / {total*total}", file=sys.stderr)
    
    # Save as pixel-exact image
    arr = combined.values.astype(np.uint8)
    # map 1 -> black (0), 0 -> white (255)
    gray = (255 - arr * 255).astype(np.uint8)
    img = np.stack([gray, gray, gray], axis=2)
    
    output_file = f'heatmap_combined_n{n}_d{d}.png'
    mpimg.imsave(output_file, img)
    
    print(f'Saved: {output_file}', file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python plot_heatmap_combined.py <n> <d>')
        sys.exit(1)
    try:
        n = int(sys.argv[1])
        d = int(sys.argv[2])
    except ValueError:
        print('n and d must be integers')
        sys.exit(1)
    main(n, d)
