#!/usr/bin/env python3
"""
Compare P-tableau and Q-tableau results for Conjecture 19E
"""
import csv
import collections

def load_results(filename):
    """Load CSV results and return (triples, distinct_d)"""
    triples = []
    with open(filename, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            d, x, y = int(row['w']), int(row['x']), int(row['y'])
            triples.append((d, x, y))
    
    # Count distinct d
    distinct_d = set(t[0] for t in triples)
    return triples, distinct_d

# Load both variants
p_triples, p_d = load_results('kl_triple_matches_nonzero_n6.csv')
q_triples, q_d = load_results('kl_triple_matches_q_tableau_n6.csv')

print("\n=== TABLEAU COMPARISON (n=6) ===\n")
print(f"P-TABLEAU (left cells):")
print(f"  Total triples: {len(p_triples)}")
print(f"  Distinct d: {len(p_d)}")

print(f"\nQ-TABLEAU (right cells):")
print(f"  Total triples: {len(q_triples)}")
print(f"  Distinct d: {len(q_d)}")

print(f"\n=== INTERSECTION ANALYSIS ===")
p_set = set(p_triples)
q_set = set(q_triples)

print(f"Triples in both: {len(p_set & q_set)}")
print(f"Triples only in P: {len(p_set - q_set)}")
print(f"Triples only in Q: {len(q_set - p_set)}")

if q_set - p_set:
    print(f"\nQ-only triples:")
    for t in sorted(q_set - p_set):
        print(f"  {t}")

d_only_p = p_d - q_d
d_only_q = q_d - p_d
d_both = p_d & q_d

print(f"\nInvolutions in both: {len(d_both)}")
print(f"Involutions only in P: {len(d_only_p)}")
print(f"Involutions only in Q: {len(d_only_q)}")

if d_both:
    print(f"\nD values in both tableaux: {sorted(d_both)}")

# Check article prediction
print(f"\n=== ARTICLE PREDICTION ===")
print(f"Article: 51 involutions should satisfy (76 - 25 = 51)")
print(f"Current P-tableau result: {len(p_d)} involutions with violations")
print(f"Expected violations: 51 or less")
print(f"Discrepancy: {len(p_d)} vs expected ~25")
print(f"Note: If article says 51 satisfy, then 25 violate = 25 expected d with triples")
print(f"      We found only {len(p_d)} - missing {25 - len(p_d)}")
