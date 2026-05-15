#!/usr/bin/env python3
"""plot_heatmap_right_preorder.py

Generate a heatmap of pairs (x,y) where both x and y are right cell equivalent
to elements in the right KL pre-order closure of a given involution d.

Usage: python plot_heatmap_right_preorder.py <n> <d>

Produces `heatmap_right_preorder_n<n>_d<d>.png` in CWD.
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


def build_right_cell_map(n):
    """Build a map from element index to right cell key (Q-tableau)."""
    total = math.factorial(n)
    right_cell_map = {}
    for idx in range(total):
        right_cell_map[idx] = right_cell_key_from_index(n, idx)
    return right_cell_map


def main(n, d):
    if n < 1:
        raise SystemExit("n must be >= 1")
    
    total = math.factorial(n)
    
    if d < 0 or d >= total:
        raise SystemExit(f"d must be in range [0, {total-1}]")
    
    print(f"Computing right preorder closure for n={n}, d={d}...", file=sys.stderr)
    involution_closure = get_right_preorder_closure(n, d)
    
    print(f"Building right cell map for S_{n}...", file=sys.stderr)
    right_cell_map = build_right_cell_map(n)
    
    # Collect all right cells that correspond to elements in the closure involutions
    right_cells_in_closure = set()
    for idx in range(total):
        right_cell = right_cell_map[idx]
        # We want pairs where both x and y have right cells from elements
        # whose involution structure relates to our closure
        # For now, collect the right cells of all elements
        # and mark them as being in the closure if they appear in the involution closure
        right_cells_in_closure.add(right_cell)
    
    # Actually, we need to think about this differently.
    # The condition is: x is right-cell-equivalent to some element in involution closure
    # So we should collect right cells that contain elements from the involution closure
    
    print(f"Processing elements in the involution closure...", file=sys.stderr)
    right_cells_from_closure = set()
    for inv_idx in involution_closure:
        right_cell = right_cell_map[inv_idx]
        right_cells_from_closure.add(right_cell)
    
    # For each element, collect its right cell
    element_right_cell = {}
    for idx in range(total):
        element_right_cell[idx] = right_cell_map[idx]
    
    print(f"Elements with right cells from closure: "
          f"{sum(1 for idx in range(total) if element_right_cell[idx] in right_cells_from_closure)} / {total}",
          file=sys.stderr)
    
    # Now build the pair matrix
    # A pair (i, j) is colored black if:
    # - element i's right cell is in right_cells_from_closure, AND
    # - element j's right cell is in right_cells_from_closure
    
    print(f"Building pair matrix...", file=sys.stderr)
    in_closure_set = set()
    for idx in range(total):
        if element_right_cell[idx] in right_cells_from_closure:
            in_closure_set.add(idx)
    
    print(f"Elements in closure: {len(in_closure_set)} / {total}", file=sys.stderr)
    
    # Create the matrix
    same_closure = pd.DataFrame(
        [[1 if (i in in_closure_set and j in in_closure_set) else 0 
          for j in range(total)] for i in range(total)],
        index=range(total),
        columns=range(total),
    )
    
    # Save as pixel-exact image
    arr = same_closure.values.astype(np.uint8)
    # map 1 -> black (0), 0 -> white (255)
    gray = (255 - arr * 255).astype(np.uint8)
    img = np.stack([gray, gray, gray], axis=2)
    
    output_file = f'heatmap_right_preorder_n{n}_d{d}.png'
    mpimg.imsave(output_file, img)
    
    print(f'Saved: {output_file}', file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python plot_heatmap_right_preorder.py <n> <d>')
        sys.exit(1)
    try:
        n = int(sys.argv[1])
        d = int(sys.argv[2])
    except ValueError:
        print('n and d must be integers')
        sys.exit(1)
    main(n, d)
