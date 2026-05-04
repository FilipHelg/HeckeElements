#!/usr/bin/env python3
import argparse
import csv
import itertools
from collections import Counter


def rs_shape(perm):
    p = []
    for x in perm:
        bumped = x
        r = 0
        while True:
            if r == len(p):
                p.append([bumped])
                break
            row = p[r]
            c = 0
            while c < len(row) and row[c] < bumped:
                c += 1
            if c == len(row):
                row.append(bumped)
                break
            bumped, row[c] = row[c], bumped
            r += 1
    return tuple(len(row) for row in p)


def inverse_tuple(perm):
    inv = [0] * len(perm)
    for i, v in enumerate(perm):
        inv[v] = i
    return tuple(inv)


def is_involution_str(perm_str):
    p = tuple(int(ch) for ch in perm_str)
    return inverse_tuple(p) == p


def all_involutions_perm_strings(n):
    out = set()
    for p in itertools.permutations(range(n)):
        t = tuple(p)
        if inverse_tuple(t) == t:
            out.add("".join(str(x) for x in t))
    return out


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def summarize(path, n):
    rows = read_rows(path)
    w_set = {r["w_perm"] for r in rows}
    w_shape_match = {
        r["w_perm"]
        for r in rows
        if rs_shape(tuple(int(ch) for ch in r["x_perm"]))
        == rs_shape(tuple(int(ch) for ch in r["w_perm"]))
    }

    by_shape = Counter(rs_shape(tuple(int(ch) for ch in w)) for w in sorted(w_set))
    by_shape_match = Counter(rs_shape(tuple(int(ch) for ch in w)) for w in sorted(w_shape_match))

    all_invol = all_involutions_perm_strings(n)

    print(f"FILE: {path}")
    print(f"rows: {len(rows)}")
    print(f"distinct w: {len(w_set)}")
    print(f"distinct w with shape(x)=shape(w): {len(w_shape_match)}")
    print(f"distinct w that are involutions: {sum(1 for w in w_set if is_involution_str(w))}")
    print(f"total involutions in S{n}: {len(all_invol)}")
    print(f"missing involutions from distinct w: {len(all_invol - w_set)}")
    print(f"missing involutions from shape-matched set: {len(all_invol - w_shape_match)}")
    print("shape distribution for distinct w:")
    for shape, cnt in sorted(by_shape.items()):
        print(f"  {shape}: {cnt}")
    print("shape distribution for shape-matched distinct w:")
    for shape, cnt in sorted(by_shape_match.items()):
        print(f"  {shape}: {cnt}")

    if len(all_invol - w_shape_match) <= 40:
        print("involutions excluded by shape-matched set:")
        for w in sorted(all_invol - w_shape_match):
            print(f"  {w}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze triple CSVs against candidate Conjecture 19E-style involution counts"
    )
    parser.add_argument("csv", nargs="+", help="CSV file(s) with columns w_perm,x_perm,y_perm")
    parser.add_argument("-n", type=int, default=6, help="n for S_n (default: 6)")
    args = parser.parse_args()

    for path in args.csv:
        summarize(path, args.n)
        print()


if __name__ == "__main__":
    main()
