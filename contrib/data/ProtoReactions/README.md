# Example Reactions
This directory contains examples of textproto reaction files.

As this collection was being assembled we added the
`single_bond`, `double_bond` and `triple_bond` directives,
as well as some other abbreviations.

We make no judgement about which form is better, that is matter
of taste

```
  single_bond: [2, 3]
```
or
```
make_bond {
  a1: 2
  a2: 3
  btype: SS_SINGLE_BOND
}
```
But note that the single_bond attribute is just a repeated value,
so any (even) number of values can be given.
```
single_bond: [2, 3, 2, 7]
```
makes single bonds between 2-3 and 2-7. The alternative would be
two separate `make_bond` messages. Brevity is valuable, but so
is clarity and precision.

Note too that the `single_bond` etc. directives function exactly
the same as a `make_bond` directive, which means they will add a bond
if the two atoms are not already bonded, or it will change the
bond type if the atoms are already bonded.

In many of the directories here there will be two .rxn files, one
with a suffix of .long.rxn. It will contain the expanded, and
backwards compatible, versions of the changes that have just been
added.

In each directory there will be run.sh script that shows how to
enumerate the reaction - assuming `trxn` is in your PATH. The example
reagents are from a Sigma catalogue and are all quite small.
These reagents were *not* carefully curated, so it is quite possible
that there might be incompatible functional groups mixed in. Preparing
reagents for a reaction is generally hard.

We have tried to use the minimal set of options in order to get
the reaction to perform as intended. Generally unique embeddings
only is desirable (-u) and sometime suppressing symmetry related
matches is important (-k). There do remain some deep seated
limitations on symmetry handling in certain situations.

Each of these reactions shows something interesting about performing
reactions with LillyMol and `trxn`. The Buchwald directory shows
a possibly interesting query for identifying a mono substituted
aryl halide. In practice it would be more likely that one would
want to have just one such group in the whole molecule in which case
`1[I,Br,Cl]-c' would be the (much simpler) smarts for that.
