# Pairwise Transformations.

We often encounter cases where a transformation is to be applied to a molecule
a limited number of times. Scenarios like this have been outlined under the
title of Positional Analogue Scanning.
[JMedChem](https://pubs.acs.org/doi/10.1021/acs.jmedchem.9b02092)
[PracticalCheminformatics](https://www.bing.com/ck/a?!&&p=daab85304f982088d17513f4540dc9d30d297fa3204d1c823d77c1452dde95d7JmltdHM9MTc0NTQ1MjgwMA&ptn=3&ver=2&hsh=4&fclid=01a16e41-0a03-6ce5-25b5-7bea0bae6db3&psq=positional+analogue+scanning&u=a1aHR0cHM6Ly9wcmFjdGljYWxjaGVtaW5mb3JtYXRpY3MuYmxvZ3Nwb3QuY29tLzIwMjAvMDQvcG9zaXRpb25hbC1hbmFsb2d1ZS1zY2FubmluZy5odG1sIzp-OnRleHQ9QnklMjBnZW5lcmF0aW5nJTIwYSUyMHNldCUyMG9mJTIwcG9zaXRpb25hbCUyMGFuYWxvZ3MlMkMlMjBvbmUsd2hlcmUlMjBQQVMlMjBjYW4lMjBiZSUyMGNvbWJpbmVkJTIwd2l0aCUyMGNvbXB1dGF0aW9uYWwlMjBhbmFseXNpcy4&ntb=1)

For example we might be changing a Methyl group to 
a Fluoro. The default behaviour in trxn is to generate a single product with
every methyl converted to Fluorine.
```
scaffold {
  id: 0
  smarts: "[CD1H3]"
  change_element {
    atom: 0
    element: "F"
  }
}
``` 
does this. If you add `-m each` to the trxn invocation it will generate
a separate product molecule for each query match site - each methyl.

So with trxn we can generate the case of all Flourinated or singly Fluorinated.
What if we wish to limit the number of instances processed to two. trxn cannot do this.
The easy way to do
that is to pipe the ouput of one invocation of trxn to another. The first invocation
does one transformation, and the second performs a second change.

Let's use a more general reaction that allows the methyl to be replaced by
any isotopically labelled sidechain.
```
scaffold {
  id: 0
  smarts: "[CD1H3]-*"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "[1]"
  join {
    a1: 1
    a2: 0
  }
}
This way any functional group in the sidechain file will replace the
existing methyl group.

```
trxn -d -m each -v -z i -P change.rxn -S - input.smi replacement.smi |\
        trxn -d -J wrscaf -m each -z i -P change.rxn - replacement.smi
```
The -d option suppresses duplicates within the same scaffold. It might be a 
good idea to add a second -d option so duplicates across the entire
set of molecules are suppressed.

The `-m each` directive says for each query match, generate a separate
product molecule.

Using `-z i` says just ignore molecules not matching the query - depends 
on how the input file has been prepared. Then the reaction file. 

The two command line arguments, 'input.smi' and 'replacement.smi' provide the scaffold and sidechain
molecules respectively.

The output is written to stdout `-S -`.

That output is then consumed by another trxn invocation. The only difference in
options is the addition of `-J wrscaf`. This means that as a molecule is read, 
the scaffold molecule, the singly changed molecule, is also written to the output,
the default 'trxn.smi' in this case. Note that maybe you want this option on the
first trxn invocation? It depends...

When doing something like this there are likely to be many instances of
multiple scaffold matches - that's the whole point of this exercise. Getting all
the informational messages about multiple scaffold matches may be less than helpful.
Add `-J nomshmsg` to the invocation to supress these messages.

Obviously this process could be extended to any number of levels, limited only
by the number of available functional groups being replaced.

Given a random set of 50k Chembl molecules tsubstructure reports methyl matches
```
Tested 50000 molecules, 36733 matched, 13267 did not. Fraction 0.73466
 13736 molecules had 1 hits
 10786 molecules had 2 hits
 6058 molecules had 3 hits
 3359 molecules had 4 hits
 1421 molecules had 5 hits
 706 molecules had 6 hits
 315 molecules had 7 hits
 201 molecules had 8 hits
 78 molecules had 9 hits
 43 molecules had 10 hits
 18 molecules had 11 hits
 11 molecules had 12 hits
 1 molecules had 13 hits
83547 reactive sites
```
We see the diminishing number of possibilities for large scale elaboration.

Using pipelined execution offers advantages on multi-core machines. For example
running the above on 50k random molecules, 36k of which have at least 1 methyl takes about 17
seconds to generate 68k product molecues. 

## Caution
You need to pay attention to how the second trxn invocation interacts with the first.
For example if the task was to add a methyl group to an existing methyl (producing ethyl,
or more generally extending a chain)
passing that ethyl group to the same reaction would then result in propyl. This might,
or might not be what you want. One easy way of handling these cases is the use of
isotopes. In the case of adding a methyl, you could make sure the added methyl
group has an isotope, and then the reaction would specify that it only reacts
with a non-isotopic existing methyl. Again, it all depends...
The last trxn invocation in the chain can remove the isotopes.

## Automation
In order to make this more straightforward the script
[consecutive_reactions](/contrib/bin/consecutive_reactions.rb) automates this task.

For example if you with to add two fluorines to aroamtic carbon sites in a molecule start
with a reaction which adds a single F atom: (id's omitted)
```
scaffold {
  smarts: "[cH]"
}
sidechain {
  reagent: "F"
  smarts: "F"
  join { a1: 0 a2: 0}
}
```
The invocation might look like
```
consecutive_reactions -P file.rxn -N 2 file.smi > new.smi
```
If instead you wanted to add multiple different sidechains, put those isotopically
labelled sidechains in a file and use this reaction
```
scaffold {
  smarts: "[cH]"
}
sidechain {
  smarts: "[1]"
  join { a1: 0 a2: 0}
}
```
This time the invocation might look like
```
consecutive_reactions -P file.rxn -N 2 file.smi sidechains.smi > new.smi
```
Make sure to make your query specific enough. The above is very generic. When used on
1000 random molecules, and with just 5 sidechains, it took 66 seconds and generated 556k products.

Use the -v option to see what options are being given to trxn. You can add more options
via `-trxn ... -trxn` of grab the default command generated here and configure to needs.

For example if you wished the second reaction to only add the second substituent to an atom in
the same ring system as the one onto which the first has been added, use the same reaction as
above for the first step. The second reacton might be:
```
scaffold {
  id: 0
  query {
    query {
      ring_system_specifier {
        base {
          environment: "[R]!@[1]"
          set_global_id: 1
        }
      }
      smarts: "[/IWgid1cH]"
    }
  }
}
sidechain {
  id: 1
  smarts: "[1]"
  join { a1: 0 a2: 0}
}
```
Provide this as the -P2 argument, which specifies a second reaction used for all
reactions after the first.
This time it takes 23 seconds to generate 232k molecules from the
starting 100 and 5 different sidechains. Beware of combinatorics.

If you wished to restrict the second substituent to the same ring, use `ring_specifier` rather
than the `ring_system_specifier` used above.

A word of explation on the reaction file.

Normally for such a complex query, it would be placed in a separate file and in the
proto specification we would use a `query_file: "fname"` construct. That also facilitates debugging
with tsubstructure. But it is also possible to inline the query as done here.

The ring system specifier looks for a ring system that contains a ring atom, bonded to an
atom with an isotope - something from the first reaction. All atoms in that ring system
are assigned global id 1.

The smarts then looks for an atom with global id 1 and which is an aromatic carbon
with an attached Hydrogen atom. This ensures that the new addition goes onto an atom
in the same ring system as where the first one was attached.

That invocation might look like
```
consecutive_reactions.sh -v -P add_sidechain.rxn -P2 add_to_same_ring.rxn \
        file.smi sidechain.smi
```

Again, this tool is just a convenience wrapper. All it does is to assemble
trxn commands, but what it does do is automate much of that, sometimes complex,
task.

## Python
In the contrib directory you will find LillyMol and RDKit implementations of a
simple Positional Analogue Scanning implementation - the RDKit implementation is
from Pat Walters's blog. The LillyMol implementation is about 3x faster. 

But interestingly that python implementation is faster than using trxn in a pipelined
mode - trxn cannot do such arbitrary n-way combinations of the matched atoms.

