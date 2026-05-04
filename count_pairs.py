import csv
s=set()
with open('kl_triple_matches_nonzero_n4.csv','r', newline='') as f:
    r=csv.reader(f)
    next(r, None)
    for row in r:
        if len(row)>=3:
            x=int(row[1]); y=int(row[2])
            if x>y: x,y=y,x
            s.add((x,y))
print(len(s))
