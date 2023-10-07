# Generative Design

Among practitioners of Generative Design there are many approaches:

* Systemtic or stochastic enumeration using known reactions and reagents.
* Use of A/I generative models.
* Cheminformatics based approaches.
* Hyrbrid approaches, combining A/I and Cheminformatics.
* other
 
Numerous tools within LillyMol have been used in support of generative design.

## Enumeration
Tools like `trxn`, `make_these_molecules`, `multiple_reactions` or  `molecular_transformations` can be used
for enumeration. Which tool to use will depend on the circumstances.

## A/I Generative Models
One of the problems with such models is their tendency to produce unrealistic molecules.
Tools like `fileconv` can be used to quickly eliminate problems such as

* Valence Errors.
* Too few or too many atoms.
* Too few, too many rings or ring systems, or too large rings.
* Disallowed elements.

That can be followed by invocation of the Lilly Medchem Rules [githug](https://...), and
other custom filtering rules, using `tsubstructure`.

One fundamental problem of precedent based rule sets, such as the Lilly Medchem Rules
is that they are precedent based. This means that many truly outrageous structures are
not covered. For example, a bridged benzene, `c12ccc(C2)cc1` passes through the Medchem
Rules with no demerits. When the rules were being formulated, nobody suggested that
we needed to guard against such things. Nor will we be adding such things to that
query set. There are an unlimited number of ridiculous molecules.

One of the most effective tools used for 'taming' generated molecules is use of
`smi2rings` [smi2rings](/docs/Molecule_Tools/smi2rings.md). That tool operates in two 
phases - one to extract the known rings from various collections, and then when new
molecules are encountered, they are looked up in the existing database(s) of rings.
Any new molecule, that involves an unprecedented ring is discarded. Any ring
for which there are a very small number of examples is also discarded. There
is no correct threshold, just a matter of how much risk you are willing to assume.

The other effective way of dealing with likely hard-to-make molecules is to use
the synthetic precedent tool, which consists of the pair `iwecfp_database_load` and
`iwecfp_database_lookup`. Molecules containing unknown, or minimally precedented fragments
at radius 2 are likely to be harder to make. Again, there are no hard and fast
rules for what a good cutoff should be.

As an example of prefiltering, this command
```
fileconv -c 10 -C 40 -r 1 -R 4 -R RS:2 -m 7 -V -S - -v generated.smi | tsubstructure -s 'c12ccc(C2)cc1' -n - - | smi2rings_bdb -d rings_database.bdb -d LOOKUP -f 10 -n -v - | iwecfp_database_lookup -d lly.spc.bdb -w PSD -
```
Takes about 5 seconds to process 20k molecules.  Slower when the
databases are swapped out, faster when the databases have been RAM
cached.  In a real generative experiment there would be more queries
for `tsubstructure`, presumably a file of queries, and the Lilly
Medchem Rules would also be inserted.

Combined with custom scoring functions, and filtering like the above,
we have found that molecules produced by generative A/I can be quite
realistic.

## Cheminformatics Approaches
Using Cheminformatics as part of a generative design solution offers many advanages.
Generally I believe that the best approach will be a combination of approaches, with
each technology, Generative I/A and Cheminformatics used for those tasks where
they provide an advantage.

Cheminformatics can very directly harness known chemistry, leading to very efficient
structure explorations via such things as sidechain replacment, ring replacement,
and simple molecular transformations.

The tool `ring_replacement` looks in a database of known rings, and can replace
the rings in an input molecule with equivalent rings from what is known. This
can provide a very efficient way of generating significant numbers of 
closely related candidates.

Not just for generative design, replacing a sidechain is a common task in
Medicinal Chemistry. The tool `get_substituents` scans collections of molecules
identifying substituents. A common task is to identify all known substituents
on an aromatic ring. Substructure query and other filters can then be applied to
the fragments retrieved, and `trxn` used to join the new candidate sidechains
to the starting molecule. When substituents are idenified, their prevalence is
also recorded, so decisions can be made when rare groups are encountered.

Similarly, if the task is to identify a new linker, the tool `get_linkers`
scans sets of molecules and identifies linker groups. You will then need
to write a reaction to excise an existing linker and replace it with one
identified by `get_linkers`.

### molecular_variants
This is a relatively new tool that was originally set up to help identify
smaller variants of molecules. When dealing with reagent acquisition problems,
similarity may, or may not work as intended - because the molecules are
small, and we encounter the old 'small change to a small molecule is actually
a large change' problem. `molecular_variants` works with a set of 
isostere-like transformation reactions. Each molecule is transformed to
what is considered to be a closely related form (usually smaller) and if
that product is present, the starting molecule is discarded - given the
availability of the methyl and ethyl variants, you probably want the methyl.
