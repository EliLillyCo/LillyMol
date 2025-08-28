# gfp_to_descriptors

The fingerprints in LillyMol have demonstrated their utility in model building
across a wide variety of targets. But for an SVMFP model, there are no individual
features, just a kernel function formed from the bits.

The ASCII representations used in .gfp files are obscure and can only be parsed with
gfp tools. This task is made more complex by the existence of multiple fingerprints in a file.
A typical record might look like
```
$SMI<S1C(C)(C)CC(=S)NC2=CC=CC=C12>
FPIW<bzzrfdy.2..2.6+.+.66......2.+2.E.U.0.V2.31.....+2..8U.6.2U6Y.U.....2U+...c.G86OU+.Y+Y.0+U..0.+U.+..G0.W2W+F..6A...2.W.K.I06E+3V..I+..+22..26s.......3E.2+.EEGE..M+.2kE..E..0.E00G6+0W.F2IU+.V.8.0.+U02.6Q6Y2...6E+U.6......6+MUEE....E.....+20.U4.+6.A.U2..2.3.EE62...6.......2E.6+M....EFE.03U+.6..2+2..8.1..6+0..6232k.82..+.E.+6..+.2.8U..+.2.+7.UE..1;2048;271;2048;271;1>
FPMK<k..E.+U3U+0UO24.+A4Ek+F4FEebU.E.3;192;42;192;42;1>
FPMK2<E..E.+U.U+0U.24.+.0E..F2+E.0U.E.3;192;22;192;22;1>
PCN<CHEMBL1522809>
MPR<1UQ00kM0.k6.2;64;16;64;16;1>
|
```
The only two fields that do not need decoding are '$SMI<' and 'PCN', the identifier. All
other data is encoded with [du_ascii](/src/Foundational/iwbits/du_bin2ascii.cc) from Daylight
Chemical Information Systems.

While the fingerprint above only involved binary fixed width fingerprints, gfp fingerprints
can also contain sparse, counted fingerprints, where only the 'on' bits are set.
```
$SMI<C(S)[C@@H](N)COCC>
NCEC3AYC<...9kE..0wM...j7...LTk21.E2..0Aq..+G2E4T7h61DfQZ.E2+.EhGb56QYbY65fiLy0uR5RU+.UA0D4lgI2KEFB+Bcf9YIQ2mVU2+.U7r6LgwTRqJQ7O+bZX0Cd62.U2+.Qr+.xHi8b7k.U2.....1>
PCN<CHEMBL3307163>
|
$SMI<S1C2=CC=C(C)C(=C2N=C1NNC(=O)CC)C>
NCEC3AYC<...9kE..0wM...j7...9ukA+.E6...ji...LV...3uY..0Al+E6+.E..IXg+bmPG.OGdGU6iPU.+.E6+0+aiLF9kOg.KsFKE3nrAE.2+.U2Mei2.4XwcO+oLLF.kmLnU.E6+.HdcpmEx9QqUG+WQs2nGyk.0.E2+HgtsE3GaObVIsOu.JS40dU20.E7V8ZX.Ob51H59Li1WBauH2.E2+.Nk9xameQ4HUeuuJu9eHc4U+.k2+nGzniBduCA1SjgisruuP2.6+.E5aSDY.tfOI1iUQFKXgNnrU.E60.T8jby1plH4UxooGsDh4As.+.E2+z9AH..2.....2>
PCN<CHEMBL1455862>
|
```
where we see that the number of bits associated with each molecule is variable.
Fixed width binary fingerprints are prefixed with 'FP', and non-colliding, sparse
counted forms with 'NC'.

If all fingerprints in a file are fixed width, it is very easy to determine
the number of columns in the output. If however there are sparse forms, the
whole file must be read first in order to determine how many features are present.

We observe that with most sparse fingerprints, there can be a large number of
bits set within a small file. For example, starting with 1000 randon Chembl
molecules, and running with all defaults.
```
gfp_make.sh -EC3:AYC file.smi > file.gfp
gfp_to_descriptors_multiple.sh -v /tmp/junk.gfp | tcount.sh -
```
yields 20k bits generated. Each column in the output will have
`1 + 20052` columns - the first column is the identifier.

Add in atom pair fingerprints
```
gfp_make.sh -MAP8HRY1 -EC3:AYC ...
```
and now our 1000 starting molecules generate 138k columns of output.
Note too that gfp_to_descriptors_multiple took 1 minutes 16 seconds in
order to generate a 266 MB output file.

This seems inefficient. We should investigate to see whether there
might be some optimisations that could be made to this code.

The header record of that file might look like
```
Name F0.4155077414 F0.4137615677 F0.4137602690 F0.4137373789 ...
```
where F0 means this bit came from the first fingerprint in the input file. Given
two fingerprints in the input there will also be features that look like
```
F1.4067627904 F1.3032426880 F1.1699724960 ...
```
Many of the bits identified are only present in a small number of molecules. We
can impose support requirements on what gets output. The `-x` and `-n` options
suppress bits that are rare, while the `-X` and `-N` options are used to suppress
bits that are overly common.

The `-x` and `-X` options express the support requirement in fractional form.
So in order to suppress features occuring less than than in 0.002 of the population
or more than 0.998 of the time, that invocation might look like
```
gfp_to_descriptors_multiple -x 0.002 -X 0.998 -v  file.gfp > file.dat
```
This completes in 0.3 seconds, reporting 
```
Will discard bits that occur in fewer than 0.002 molecules
Will discard bits that occur in more than 998 molecules
Found 118270 bits set in NCMAP8HRY1 fingerprint
Found 20052 bits set in NCEC3AYC fingerprint
100846 bits below threshold 2
0 bits above threshold 998
Will produce 37476 descriptors
```
We see that removing singleton features discards 100k of the 138k.

> ![!NOTE] with 1000 molecules, 0.001 corresponds
> to 1 molecule, but every bit occurs in at least 1 molecule, so nothing would be
> filtered.

The `-n` option can be used to specify a specific minimum number of occurrences
avoiding any unexpected problems with fractions and truncations.

As is typically found, the upper support requirement did not suppress any
molecules. Perhaps in a combinatorial library this might happen.

Another way of dealing with overly wide outputs is to fold the output to fixed width.
The `-d` options specifies the number of columns you want to produce. Note that this
will necessarily induce collisions, where multiple columns will be mapped to the
same output feature. It is upredictable if this will be problematic or not.

> [!NOTE] With gfp_to_descriptors_multiple all information about the fingerprint is
> captured and features get consistent names.

This means that if you separately generate two different files, they can be combined.
```
gfp_make.sh ... file1.smi > file1.gfp
gfp_make.sh ... file2.smi > file2.gfp

gfp_to_descriptors_multiple -n 2 file1.gfp > file1.dat
gfp_to_descriptors_multiple -n 2 file2.gfp > file2.dat
```
At this stage the two .dat files will likely have different column counts, and therefore
different columns. Some columns will be shared. No current LillyMol
tool combines these, but `iwcut` can be used to fetch all the features that are in
one of the files from another file.

Alternatively a tool like Julia can easily perform this operation.
```
using DataFrames

julia> d1 = DataFrame(id=["id1", "id2"], F1 = [1.1, 2.1], F2=[1.2, 2.2], F3=[1.3, 2.3])
2×4 DataFrame
 Row │ id      F1       F2       F3      
     │ String  Float64  Float64  Float64 
─────┼───────────────────────────────────
   1 │ id1         1.1      1.2      1.3
   2 │ id2         2.1      2.2      2.3

julia> d2 = DataFrame(id=["id3", "id4"], F1 = [3.1, 4.1], F3=[3.3, 4.3], F4=[3.4, 4.4])
2×4 DataFrame
 Row │ id      F1       F3       F4      
     │ String  Float64  Float64  Float64 
─────┼───────────────────────────────────
   1 │ id3         3.1      3.3      3.4
   2 │ id4         4.1      4.3      4.4

julia> coalesce.(vcat(d1, d2;cols=:union), 0)
4×5 DataFrame
 Row │ id      F1       F2       F3       F4      
     │ String  Float64  Float64  Float64  Float64 
─────┼────────────────────────────────────────────
   1 │ id1         1.1      1.2      1.3      0.0
   2 │ id2         2.1      2.2      2.3      0.0
   3 │ id3         3.1      0.0      3.3      3.4
   4 │ id4         4.1      0.0      4.3      4.4
```

For Pandas, something like
```
result = pd.concat([df, df1], axis=1, ignore_index=False)
```
might work. For Polars something like
```
result = pl.concat([df1, df2], how="diagonal")
```
might work.

The -D option to `iwcut` specifies a descriptor file. The descriptor names from that
file are extracted from the first line. Then the main input is processed and only
the features that were in the -D file are written.

So if 'file1.dat' contained
```
id F1 F2 F3
```
and 'file2.dat' contained
```
id F1 F3 F4
```
then
```
iwcut -D file1.dat -k . file2.dat
```
would yield
```
id F1 F3
```
The column 'F2' is not found in file2.dat so it cannot be written.



