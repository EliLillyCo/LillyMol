# marvin_pka
This tool parses the output from cxcalc (Chemaxon/Marvin) and generates a molecule with
assigned formal charges.

The output from
```
cxcalc pka file.smi > file.pka
```

might look like
```
id      apKa1   apKa2   bpKa1   bpKa2   atoms
1       16.76           0.07            3,2
2       3.55    5.69                    2,7
3       3.36            -1.17   -9.49   7,1,2
4                       10.45   4.50    7,1
5                       -0.22           9
6       11.47   15.28   0.69    0.07    6,3,1,8
7                       4.56            8
8       13.50   14.14   5.89    2.07    8,7,9,7
9       2.77    13.78   9.78    -3.80   8,6,3,6
```
where each record contains varying numbers of entries for most basic
and most acidic pK values.

Note that the last column is a list of the correspinding atom numbers
associated with the sites in the earlier columns. A complete workflow
might look like

```
cxcalc pka file.smi > file.pka
marvin_pka -M file.pka -P 7.4 -d 1.0 file.smi
```
where the resulting output will have formal charges applied if the
computed pK is more than 1 log unit away from the assumed pH of 7.4
in this case.
