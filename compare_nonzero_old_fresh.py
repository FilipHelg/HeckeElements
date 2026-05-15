import csv
nonzero='kl_triple_matches_nonzero_n6.csv'
old='kl_triple_matches_n6.csv'
fresh='kl_triple_matches_n6_current.csv'
ns=set()
with open(nonzero) as f:
    r=csv.DictReader(f)
    for row in r:
        ns.add((int(row['w']),int(row['x']),int(row['y'])))
old_set=set()
with open(old) as f:
    r=csv.DictReader(f)
    for row in r:
        old_set.add((int(row['w']),int(row['x']),int(row['y'])))
found_in_old=[t for t in ns if t in old_set]
print('found in old:',len(found_in_old))
# check fresh
fresh_set=set()
with open(fresh) as f:
    r=csv.DictReader(f)
    for row in r:
        fresh_set.add((int(row['w']),int(row['x']),int(row['y'])))
found_in_fresh=[t for t in ns if t in fresh_set]
print('found in fresh:',len(found_in_fresh))
print('\nMissing nonzero triples (not in old, but in fresh):')
for t in sorted(ns):
    if t not in old_set and t in fresh_set:
        print(t)
