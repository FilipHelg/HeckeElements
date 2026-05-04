#!/usr/bin/env python3
# plot_kl_triples.py
import sys
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict

sns.set(style="whitegrid")
counts_cmap = LinearSegmentedColormap.from_list(
    "counts_blue_red",
    ["#e6f7ff", "#7fc8ff", "#ff9b9b", "#d7191c"]
)

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


def left_cell_boundaries(order, keys):
    boundaries = []
    if not order:
        return boundaries

    current_key = keys[order[0]]
    for pos, idx in enumerate(order):
        key = keys[idx]
        if key != current_key:
            boundaries.append(pos)
            current_key = key
    return boundaries


def apply_cell_start_ticks(ax, block_starts, block_mins):
    labels = [str(m) for m in block_mins]
    ax.set_xticks(block_starts)
    ax.set_yticks(block_starts)
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_yticklabels(labels, fontsize=6)


def rs_shape_from_perm(perm_str):
    perm = [int(c) for c in perm_str.strip()]
    P = []
    for x in perm:
        bumped = x
        r = 0
        while True:
            if r == len(P):
                P.append([bumped])
                break
            # find first element greater than bumped
            row = P[r]
            c = 0
            while c < len(row) and row[c] < bumped:
                c += 1
            if c == len(row):
                row.append(bumped)
                break
            # bump
            bumped, row[c] = row[c], bumped
            r += 1
    # return shape as tuple of row lengths
    return tuple(len(row) for row in P)

def main(csv_path):
    df = pd.read_csv(csv_path, dtype=str)
    # ensure numeric indices present
    for col in ['w','x','y']:
        if col in df.columns:
            df[col] = df[col].astype(int)
    n = len(df['x_perm'].iloc[0]) if len(df) else 0
    # compute shapes from permutation columns (use x_perm)
    df['x_shape'] = df['x_perm'].apply(rs_shape_from_perm)
    df['shape_str'] = df['x_shape'].apply(lambda s: "-".join(map(str,s)))
    # basic counts by involution w
    counts_w = df['w'].value_counts().sort_index()
    plt.figure(figsize=(10,4))
    counts_w.plot(kind='bar', width=0.8)
    plt.title('Matches per involution w')
    plt.xlabel('w (index)')
    plt.ylabel('Number of matching pairs (rows)')
    plt.tight_layout()
    plt.savefig('matches_per_w.png', dpi=200)

    # counts by RS-shape
    counts_shape = df['shape_str'].value_counts()
    plt.figure(figsize=(10,5))
    counts_shape.plot(kind='bar')
    plt.title('Matches per RS-shape (based on x_perm)')
    plt.xlabel('RS-shape')
    plt.ylabel('Number of matching pairs')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('matches_per_shape.png', dpi=200)

    # matches per x (how many rows involve each x)
    matches_per_x = pd.concat([df['x'], df['y']]).value_counts()
    plt.figure(figsize=(8,4))
    matches_per_x.hist(bins=50)
    plt.title('Histogram: matches per element (x or y)')
    plt.xlabel('Number of matches')
    plt.ylabel('Count of elements')
    plt.tight_layout()
    plt.savefig('hist_matches_per_element.png', dpi=200)

    # heatmap: w vs shape counts (pivot)
    pivot = df.pivot_table(index='w', columns='shape_str', values='x', aggfunc='count', fill_value=0)
    # limit columns to top shapes for readability
    top_shapes = counts_shape.index[:20]
    pivot = pivot.loc[:, pivot.columns.isin(top_shapes)]
    plt.figure(figsize=(12,6))
    sns.heatmap(pivot, cmap='viridis', norm=None)
    plt.title('Heatmap: counts by involution w (rows) vs RS-shape (cols)')
    plt.xlabel('RS-shape (top 20)')
    plt.ylabel('w')
    plt.tight_layout()
    plt.savefig('heatmap_w_vs_shape.png', dpi=200)

    # normalized match-rate by shape: observed matches / possible pairs in shape-buckets
    # derive bucket membership from unique indices in CSV (approximation)
    shape_members = defaultdict(set)
    for _, row in df.iterrows():
        shape_members[row['shape_str']].add((int(row['x'])))
        shape_members[row['shape_str']].add((int(row['y'])))
    observed = counts_shape.to_dict()
    normalized = {}
    for shape, obs in observed.items():
        m = len(shape_members.get(shape, []))
        possible = m * (m - 1) / 2 if m >= 2 else 1
        normalized[shape] = obs / possible if possible > 0 else 0
    norm_series = pd.Series(normalized).sort_values(ascending=False)
    plt.figure(figsize=(10,5))
    norm_series.plot(kind='bar')
    plt.title('Normalized match-rate by RS-shape (observed / possible pairs)')
    plt.xlabel('RS-shape')
    plt.ylabel('Normalized match-rate')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('normalized_rate_by_shape.png', dpi=200)

    # heatmap: counts for exact (x,y) pairs across all w
    pair_counts = df.groupby(['x', 'y']).size().unstack(fill_value=0)
    total = math.factorial(n) if n > 0 else 0
    full_index = list(range(0, total))
    pair_counts = pair_counts.reindex(index=full_index, columns=full_index, fill_value=0)
    pair_mask = (pair_counts == 0)
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    ax.set_facecolor('black')
    sns.heatmap(pair_counts, cmap=counts_cmap, mask=pair_mask)
    plt.title('Heatmap: counts of (x,y) pairs across all w')
    plt.xlabel('y')
    plt.ylabel('x')
    plt.tight_layout()
    plt.savefig('heatmap_xy_counts.png', dpi=200)

    # heatmap: same-left-cell relation for all permutations in S_n
    if n > 0:
        order, keys, block_starts, block_mins = build_left_cell_order(n)

        same_left_cell = pd.DataFrame(
            [[1 if keys[i] == keys[j] else 0 for j in range(total)] for i in range(total)],
            index=range(total),
            columns=range(total),
        )
        plt.figure(figsize=(10, 8))
        sns.heatmap(same_left_cell, cmap='Greys', cbar=False, square=True)
        plt.title('Same-left-cell relation in Lehmer order')
        plt.xlabel('y')
        plt.ylabel('x')
        plt.tight_layout()
        plt.savefig('heatmap_left_cells.png', dpi=200)

        ordered_same_left_cell = same_left_cell.iloc[order, order]
        boundaries = block_starts[1:]
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(ordered_same_left_cell, cmap='Greys', cbar=False, square=True)
        apply_cell_start_ticks(ax, block_starts, block_mins)
        ax.axhline(0, color='lime', linewidth=0.2, alpha=0.35)
        ax.axvline(0, color='lime', linewidth=0.2, alpha=0.35)
        for boundary in boundaries:
            ax.axhline(boundary, color='lime', linewidth=0.2, alpha=0.35)
            ax.axvline(boundary, color='lime', linewidth=0.2, alpha=0.35)
        plt.title('Same-left-cell relation ordered by left cells')
        plt.xlabel('y (labels = cell minima)')
        plt.ylabel('x (labels = cell minima)')
        plt.tight_layout()
        plt.savefig('heatmap_left_cells_ordered.png', dpi=200)

        ordered_pair_counts = pair_counts.iloc[order, order]
        ordered_pair_mask = (ordered_pair_counts == 0)
        plt.figure(figsize=(10, 8))
        ax = plt.gca()
        ax.set_facecolor('black')
        ax = sns.heatmap(ordered_pair_counts, cmap=counts_cmap, mask=ordered_pair_mask)
        apply_cell_start_ticks(ax, block_starts, block_mins)
        ax.axhline(0, color='lime', linewidth=0.2, alpha=0.35)
        ax.axvline(0, color='lime', linewidth=0.2, alpha=0.35)
        for boundary in boundaries:
            ax.axhline(boundary, color='lime', linewidth=0.2, alpha=0.35)
            ax.axvline(boundary, color='lime', linewidth=0.2, alpha=0.35)
        plt.title('Heatmap: counts of (x,y) pairs ordered by left cells')
        plt.xlabel('y (labels = cell minima)')
        plt.ylabel('x (labels = cell minima)')
        plt.tight_layout()
        plt.savefig('heatmap_xy_counts_ordered.png', dpi=200)

    print("Saved: matches_per_w.png, matches_per_shape.png, hist_matches_per_element.png, heatmap_w_vs_shape.png, normalized_rate_by_shape.png, heatmap_xy_counts.png, heatmap_left_cells.png, heatmap_left_cells_ordered.png, heatmap_xy_counts_ordered.png")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python plot_kl_triples.py <kl_triple_cs v>")
        sys.exit(1)
    main(sys.argv[1])