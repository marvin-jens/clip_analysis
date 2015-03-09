#!/usr/bin/env python
# Simple script to filter lines of the input for meeting a numeric 
# threshold on one, or one of multiple (comma-separated) column numbers.
# If threshold is required on each column, you can chain multiple 
# scorethresh.py calls with pipes.
# by default, the column value has to be larger or equal to the threshold.
# If the column number is negative, this changes to less or equal.

import sys
cols = [int(col) for col in sys.argv[1].split(",")]
flags = [(col >= 0)*2 -1 for col in cols]
cols = [abs(c)-1 for c in cols]
threshs = [float(th) for th in sys.argv[2].split(",")]*len(cols)

for line in sys.stdin:
    if line.startswith("#"):
        print line,
        continue

    data = line.split("\t") + [0,]
    values = [data[col] for col in cols]

    valid = False
    for value,flag,thresh in zip(values,flags,threshs):
        try:
            value = float(value)
        except ValueError:
            continue
        if value * flag >= thresh * flag:
            valid = True

    if valid:
        print line,

