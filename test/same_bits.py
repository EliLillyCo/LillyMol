#!/usr/bin/env python

import sys 
import pandas as pd

def same_bits(file1, file2):
    # Read both files into pandas DataFrames
    df1 = pd.read_csv(file1, sep=" ", index_col=0)
    df2 = pd.read_csv(file2, sep=" ", index_col=0)
    
    if len(df1.columns) != len(df2.columns):
        return 1
    
    for (idx1, row1), (idx2, row2) in zip(df1.iterrows(), df2.iterrows()):
        if idx1 != idx2:
            return 1
        if sorted(list(row1)) != sorted(list(row2)):
            return 1

    return 0

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Must specify two files\n")
        sys.exit(1)

    file1_path, file2_path = sys.argv[1], sys.argv[2]
    sys.exit(same_bits(file1_path, file2_path))
