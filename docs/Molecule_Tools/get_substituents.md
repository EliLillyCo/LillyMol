# get_substituents
## Problem Statement.
`get_substituents` is a tool for identifying the substituents defined
by a substructure query. This can be useful in *de-novo* molecule
generation where you may have a molecule with an undesirable substituent
and wish to explore replacements.
But in order to ensure maximum probability
of synthetic feasibility, replacement fragments should have been
observed previously in the same structural context. You might also want
to restrict attention to well precedented fragments, those that have
been commonly found in such a situation.

The tool described here with both identify substituents that have precedent
in a given atomic context, but can also provide information on
frequency of occurrence.

## Alternatives
This task can be accomplished via other means, but with greater complexity.

`dicer` can break molecules at any bond you specify - by over-riding the
default queries built into it. But dicer is slow and complex, and you then face
the difficult task of post-processing the output - since it writes both
sides of the broken bonds.

`tsubstructure` could be used (with the `-j` option) to label the atoms
either side of a bond of interest, but you then need to write a reaction
to break that bond and get the fragments. `tsubstructure` now has the
concept of a substituent, but it is complex, and post processing
would still be needed.

`trxn` would probably be the best choice as an alternative to this tool,
but it is more complex. You would need to identify the bonds to be broken via
a smarts, break those bonds, place isotopes, and then remove the
part that is not of interest. Figure out how to handle multiple
matches, etc... Then if you wanted populations of fragments
you would need to do that as a separate post-processing step
(unique_molecules -e).

This seems like a common task, and having a tool to do this specific
task seems worthwhile.

## HOTO
Here is the usage message from `get_substituents`.
```
Identify substituents that are defined by the atoms adjacent to matched query atoms
 -s <smarts>       identify the matched atoms
 -q <query>        identify the matched atoms
 -O <atom>         from which query atom does the substituent sprout (default 0)
 -m <natoms>       min fragment size
 -M <natoms>       max fragment size
 -r <nrings>       min number of rings in the fragment
 -R <nrings>       max number of rings in the fragment
 -Y <nsys>         max number of ring systems in the fragment
 -I <isotope>      place <isotope> at the join atoms
 -f <smarts>       fragment must contain smarts
 -F <query>        fragment must contain query specification
 -a                each of the fragment must contain queries must match (by default any)
 -n                suppress normal output (write summary via -S)
 -S <fname>        write summary data Dicer.DicerFragment proto, to <fname>
 -z i              ignore molecules not matching any query
 -v                verbose output
```
The point of substitution, the anchor atom that will be retained, must be
specified by the `-s` and/or `-q` options. If you were looking for aromatic
substituents, this would be the atom in the ring. By default, the matched
atom from which the tool looks for substituents, is the first matched atom
in the query, but this can be changed with the -O option.

The rest of the atoms in the query are considered **not** part of the
substituent. So while `c` is unambiguous, because it does not consider
ring bonds, `C` would examine the substituents of every atom connected
via a chain bond, maybe four of them. Use more query atoms to properly
define the atomic context.

For example, if you were looking for substituents that have appeared opposite
the `n` atom of pyridine, that could be
```
-s n1ccccc1 -O 3
```
which would of course be the same as
```
-s c1ccncc1
```
without the use of the `-O` option. But the smarts starting with `n` is
more efficient.

If the query is for chain atoms, make sure to specify all the fixed
atoms in the query. So, if you wanted to find all substituents found
after a C(CH3)(CH3) group `-s '[CD4](-[CH3])(-[CH3])-*'`, but note
that this query will match both left and right. This can be tricky,
just like any substructrure based matching.

### Annotation
You almost certainly want to specify an isotope for the attachment atom.
This will be the atom in the fragment that would join to the matched atom 
defined above. So, if you were looking for aromatic substituents, the
query atom above would be `a` and if the fragment was Fluorine, the
isotope would be on the `F` atom.

### Constraints
You can impose requirements on the substituents found.

* -m minimum number of atoms in the substituent
* -M maximum number of atoms in the substituent
* -f <smarts>   substructure query that must be present in the substituent
* -Q <query>    substructure query that must be present in the substituent
* -r <nrings>   min number of rings in the fragment
* -R <nrings>   max number of rings in the fragment
* -Y <nsys>     max number of ring systems in the fragment

Other constraints could be added.

By default at least one of the queries specified via `-f` or `-F` must match.
But if you add the `-a` option, then **each** of the must have queries must
match. Default match is **or**, but that can be changed to **and** with
the `-a` option.

So if you had an inclusion, and an exclusion, there are a couple of ways
that could be done. You want fragments with exactly one Chlorine atoms
but zero Bromine atoms.

```
-f '1Cl&&0Br'
```
or
```
-f '1Cl' -f '0Br' -a
```
If you wanted no ring atoms `-f '0[R]'`, or if you wanted at most one
ring system `-f '0[/IWfsid1R].[/IWfsid2R]` which means you want zero
occurrences of two atoms that are in different ring systems.

Of course all this processing is no different from what can be done with
`tsubstructure` so whether filtering is done here during generation or
as a post-processing step is arbitrary.

## Summary File
Frequently the purpose of this tool is to accumulate data about
the fragments present in a set of known molecules. A Dicer.DicerFragment
proto can be generated for each fragment, and written with the `-S` option.

For example if we run this on the first 500k molecules in Chembl
```
get_substituents -i do=500000 -s c -I 1 -z i -v -m 4 -M 15 -n -S Sfile chembl.smi
```
We are looking for aromatic substituents - something attached to an aromatic
carbon. Fragments must contain between 4 at 15 atoms, and we ignore any
molecules that do not contain an aromatic carbon. Isotope 1 is placed
on the attachment points. The output file `Sfile` looks like
```
O=C1N(CCCC12C[1NH]CC2)CCOC iso: ATT smi: "O=C1N(CCCC12C[1NH]CC2)CCOC" par: "CHEMBL3484591" nat: 15 n: 1 
O=C(N[1CH2]C)NCCCN(C)C iso: ATT smi: "O=C(N[1CH2]C)NCCCN(C)C" par: "CHEMBL3494056" nat: 12 n: 1 
ON=[1CH]N1CCN(CC=C)CC1 iso: ATT smi: "ON=[1CH]N1CCN(CC=C)CC1" par: "CHEMBL1861853" nat: 12 n: 1 
```
where the first column is the unique smiles of the fragment, followed by
a `text_format` representation of the proto.
This file is ready for loading into any key/value database if needed.

### Explanation
To understand the proto, take a look at the proto file in
`Molecule_Tools/dicer_fragments.proto`. For this discussion, the two
important fields are `par` which is the first molecule in
Chembl that exemplified the query, and `n` which is the number of
instances across all molecules examined.

The last column is the number of instances of this fragment in the input.
We can sort by that, and see the most common aromatic substituents.
```
COc1cc[1cH]cc1 iso: ATT smi: "COc1cc[1cH]cc1" par: "CHEMBL39210" nat: 8 n: 410 
O1CC[1NH]CC1 iso: ATT smi: "O1CC[1NH]CC1" par: "CHEMBL1506642" nat: 6 n: 507 
C[1CH](C)C iso: ATT smi: "C[1CH](C)C" par: "CHEMBL1904213" nat: 4 n: 594 
C[1SH](=O)=O iso: ATT smi: "C[1SH](=O)=O" par: "CHEMBL116841" nat: 4 n: 629 
O=C([1NH2])C iso: ATT smi: "O=C([1NH2])C" par: "CHEMBL3349336" nat: 4 n: 658 
O=[1CH]OC iso: ATT smi: "O=[1CH]OC" par: "CHEMBL1866040" nat: 4 n: 1365 
O=[1CH]OCC iso: ATT smi: "O=[1CH]OCC" par: "CHEMBL4175237" nat: 5 n: 1406 
N[1SH](=O)=O iso: ATT smi: "N[1SH](=O)=O" par: "CHEMBL3265256" nat: 4 n: 1705 
[1cH]1ccccc1 iso: ATT smi: "[1cH]1ccccc1" par: "CHEMBL38876" nat: 6 n: 2277 
F[1CH](F)F iso: ATT smi: "F[1CH](F)F" par: "CHEMBL503149" nat: 4 n: 2317 
```
Remember, by using the `-m 4` option combination, we do not see much
more common substituents such as `F` itself. Note too that we see overlapping
fragments, CHEMBL1866040 and CHEMBL4175237 as an example. This is of
course to be expected.

Or as a more complex query, what fragments have been found at the
c2 (between the nitrogens) of a benzimidazole. Drop the lower
atom count constraint, but drop things that have two different
ring systems (with the `-M` value set to 15, that ring system
constraint probably did very little.
```
get_substituents -s c1nc2ccccc2n1 -I 1 -z i -v -M 15 -n -S Sfile -f '0[/IWfsid1R].[/IWfsid2R]' chembl.smi
```
Now the most common fragments are
```
N1CC[1CH2]CC1 iso: ATT smi: "N1CC[1CH2]CC1" par: "CHEMBL159966" nat: 6 n: 19 
[1NH]1CCNCC1 iso: ATT smi: "[1NH]1CCNCC1" par: "CHEMBL292066" nat: 6 n: 22 
[1CH2]1CC1 iso: ATT smi: "[1CH2]1CC1" par: "CHEMBL172461" nat: 3 n: 23 
s1c[n][1cH]c1 iso: ATT smi: "s1c[n][1cH]c1" par: "CHEMBL625" nat: 5 n: 25 
s1[1cH]ccc1 iso: ATT smi: "s1[1cH]ccc1" par: "CHEMBL201221" nat: 5 n: 26 
o1[1cH]ccc1 iso: ATT smi: "o1[1cH]ccc1" par: "CHEMBL201094" nat: 5 n: 39 
[n]1c[1cH]ccc1 iso: ATT smi: "[n]1c[1cH]ccc1" par: "CHEMBL83103" nat: 6 n: 45 
[n]1cc[1cH]cc1 iso: ATT smi: "[n]1cc[1cH]cc1" par: "CHEMBL379376" nat: 6 n: 52 
[n]1[1cH]cccc1 iso: ATT smi: "[n]1[1cH]cccc1" par: "CHEMBL72683" nat: 6 n: 75 
[1cH]1ccccc1 iso: ATT smi: "[1cH]1ccccc1" par: "CHEMBL153535" nat: 6 n: 103 
```
We see that a phenyl ring is the most common benzimidazole substituent,
and then pyridine.

### Combinations
If, via multiple invocations, you generate multiple summary files, there is
a python tool that can combine them, so that you end up with count data 
that is the aggregate of all those individual files. I have used it,
but I am not recommending it for others. But for the brave, take a look
at `Molecule_Tools/dicer_fragments_collate.py` which was originally
written for dicer, but since this is the same proto, should work
for these files.

# Conclusion
`get_substituents` is designed to address a very specific problem, identification
of different substituents that have a high probability of synthetic success.
Given a summary file, you can filter by frequency of occurrence, thereby
increasing the probability of success. Of course more complex patterns for
the substitution point will also increase that probability, one could
easily imagine queries for substituted rings with various electron donating
or withdrawing groups. Or what substituents have been found on the other
side of an amide `-s 'C(=O)N' or the amide the other way round
`-s 'N-[CD3]=O`.

And of course, the no restriction substituent will also work
`-s '*'` which will match any atom. Perhaps that might be called
for sometimes. This will generate a large number of substituents,
although when I just tested it, if you impose a sensible upper
atom count (-M) the numbers and timing are good.

Note that currently the tool does **not** check the bond type
of the bond being broken. Perhaps it should... Or perhaps it
should have the ability to set the isotope based on the kind
of bond broken. Or...

Note that if you specifically wanted the case of a double bond,
something like `-s '[$([CD3]=*)](-*)-*'` should work - we are
looking for 3 connected carbon atoms that have a double bond. But
by using a recursive smarts, the =* atom is NOT part of the reported
matched atoms.  Quite complex matching conditions could be set up this
way. If we do that for the first 200k molecules in Chembl
 the most common fragments are
```
CN(=O)=O iso: ATT smi: "CN(=O)=O" par: "CHEMBL2251224" nat: 4 n: 142 
N#CN iso: ATT smi: "N#CN" par: "CHEMBL1984194" nat: 3 n: 166 
NC iso: ATT smi: "NC" par: "CHEMBL3228370" nat: 2 n: 171 
Cc1ccccc1 iso: ATT smi: "Cc1ccccc1" par: "CHEMBL1594090" nat: 7 n: 184 
S=C(N)NN iso: ATT smi: "S=C(N)NN" par: "CHEMBL500557" nat: 5 n: 189 
ON iso: ATT smi: "ON" par: "CHEMBL2139230" nat: 2 n: 776 
C iso: ATT smi: "C" par: "CHEMBL1529759" nat: 1 n: 1469 
N iso: ATT smi: "N" par: "CHEMBL225304" nat: 1 n: 3459 
S iso: ATT smi: "S" par: "CHEMBL256250" nat: 1 n: 4592 
O iso: ATT smi: "O" par: "CHEMBL1213530" nat: 1 n: 80233 
```
which is probably not surprising.

Or if we were looking for substituents that attach to a phenyl ring
via an oxygen `-s '[OD2]-c1ccccc1'` should do it. But note that the
atoms identified in the substituent will NOT contain the Oxygen atom.
Using `-s '[$(c-[OD2])]1ccccc1'` would remove the OD2 from the 
matched atoms, and would keep it in the fragment.

The isotopically labelled fragments generated should enable rapid
generation of plausible new molecules.

This tool is related to `substituent_identification` another tool
designed to deal with contextual replacements.
