#! /usr/bin/env python

# invoke with column nr to extract as first parameter followed by
# filenames. The files should all have the same number of rows

import sys

col = int(sys.argv[1])
res = {}

for file_name in sys.argv[2:]:
    for line_nr, line in enumerate(open(file_name)):
        res.setdefault(line_nr, []).append(line.strip().split('\t')[col-1])

for line_nr in sorted(res):
    print ("\t".join(res[line_nr]))
