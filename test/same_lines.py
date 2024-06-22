#!/usr/bin/env python

import sys
from collections import Counter

def same_lines(file1_path, file2_path):
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        # Count occurrences of each line in both files
        count1 = Counter(file1)
        count2 = Counter(file2)

    # Check if the counts are equal for each line
    if count1 == count2:
        return 0
    else:
        return 1

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Must specify two files\n")
        sys.exit(1)

    file1_path, file2_path = sys.argv[1], sys.argv[2]
    sys.exit(same_lines(file1_path, file2_path))
