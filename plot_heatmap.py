#!/usr/bin/env python3
"""plot_heatmap.py

Generate the ordinary and the ordered heatmap of left cells in S_n.

Usage: python plot_heatmap.py <n>

Produces `heatmap_left_cells.png` and `heatmap_left_cells_ordered.png` in CWD.
"""
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import pandas as pd
import numpy as np
from collections import defaultdict

sns.set(style="whitegrid")


def index_to_perm(n, index):
    remaining = list(range(n))
    perm = []
    for k in range(n - 1, -1, -1):
        fact = math.factorial(k)
        digit = index // fact if k > 0 else 0
        index = index % fact if k > 0 else 0
        perm.append(remaining.pop(digit))
    return perm


def rs_tableau_from_perm(perm):
    P = []
    for x in perm:
        bumped = x
        r = 0
        while True:
            if r == len(P):
                P.append([bumped])
                break
            row = P[r]
            c = 0
            while c < len(row) and row[c] < bumped:
                c += 1
            if c == len(row):
                row.append(bumped)
                break
            bumped, row[c] = row[c], bumped
            r += 1
    return tuple(tuple(row) for row in P)


def left_cell_key_from_index(n, index):
    return rs_tableau_from_perm(index_to_perm(n, index))


def build_left_cell_order(n):
    total = math.factorial(n)
    keys = [left_cell_key_from_index(n, idx) for idx in range(total)]
    groups = defaultdict(list)
    for idx, key in enumerate(keys):
        groups[key].append(idx)

    ordered_groups = sorted(groups.items(), key=lambda item: min(item[1]))
    order = []
    block_starts = []
    block_mins = []
    pos = 0
    for _, group in ordered_groups:
        group_sorted = sorted(group)
        block_starts.append(pos)
        block_mins.append(group_sorted[0])
        order.extend(group_sorted)
        pos += len(group_sorted)
    return order, keys, block_starts, block_mins


def apply_cell_start_ticks(ax, block_starts, block_mins):
    labels = [str(m) for m in block_mins]
    ax.set_xticks(block_starts)
    ax.set_yticks(block_starts)
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_yticklabels(labels, fontsize=6)


def main(n):
    if n < 1:
        raise SystemExit("n must be >= 1")
    total = math.factorial(n)

    order, keys, block_starts, block_mins = build_left_cell_order(n)

    same_left_cell = pd.DataFrame(
        [[1 if keys[i] == keys[j] else 0 for j in range(total)] for i in range(total)],
        index=range(total),
        columns=range(total),
    )

    # Save ordinary heatmap as a pixel-exact image (n! x n!) so each pair maps to one pixel
    arr = same_left_cell.values.astype(np.uint8)
    # map 1 -> black (0), 0 -> white (255)
    gray = (255 - arr * 255).astype(np.uint8)
    img = np.stack([gray, gray, gray], axis=2)
    mpimg.imsave('heatmap_left_cells.png', img)

    # Ordered heatmap (grouped by left cells) as pixel-exact image with boundaries
    ordered_arr = arr[np.ix_(order, order)]
    gray_o = (255 - ordered_arr * 255).astype(np.uint8)
    img_o = np.stack([gray_o, gray_o, gray_o], axis=2)
    # draw boundaries (lime) at block starts
    boundaries = block_starts[1:]
    lime = np.array([0, 255, 0], dtype=np.uint8)
    for b in boundaries:
        if 0 <= b < img_o.shape[0]:
            img_o[b:b+1, :, :] = lime
            img_o[:, b:b+1, :] = lime
    # also draw outer axes (row 0 / col 0)
    img_o[0:1, :, :] = lime
    img_o[:, 0:1, :] = lime
    mpimg.imsave('heatmap_left_cells_ordered.png', img_o)

    print('Saved: heatmap_left_cells.png, heatmap_left_cells_ordered.png')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python plot_heatmap.py <n>')
        sys.exit(1)
    try:
        n = int(sys.argv[1])
    except ValueError:
        print('n must be an integer')
        sys.exit(1)
    main(n)
