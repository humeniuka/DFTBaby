#!/usr/bin/env python
"""
add numbers from all columns (which can be in different files)
"""
import sys
import numpy as np

if __name__ == "__main__":
    columns = []
    if len(sys.argv) < 2:
        print "Usage: python %s <files whose columns should be added>" % sys.argv[0]
        exit(-1)
    for data_file in sys.argv[1:]:
        data = np.loadtxt(data_file)
        if len(data.shape) == 1:
            # single column
            columns.append(data)
        else:
            nrow,ncol = data.shape
            for c in range(0, ncol):
                columns.append(data[:,c])
    columns = np.array(columns)
    # add all columns
    colsum = np.sum(columns,axis=0)

    for r in colsum:
        print "%s" % r

