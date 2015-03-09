#!/usr/bin/env python
# computes the sum of multiple input columns and 
# adds an additional column with the result.
# column numbers are supplied as comma-separated 
# list and one based. Negative column numbers are 
# used to signal subtraction instead of addition.

import sys
cols = [int(col) for col in sys.argv[1].split(",")]
signs = [(col >= 0)*2 -1 for col in cols]
cols = [abs(c)-1 for c in cols]


name = "sum"
if len(sys.argv) > 2:
    name = sys.argv[2]
    
for line in sys.stdin:
    if line.startswith("#"):
        print line.rstrip()+'\t%s' % name
        continue
    data = line.rstrip().split("\t")
    values = [float(data[col]) for col in cols]

    res = 0.
    for value,sign in zip(values,signs):
        res += value * sign

    data.append(str(res))
    print "\t".join(data)

