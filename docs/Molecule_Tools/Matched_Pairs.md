# Matched Pairs
LillyMol has been used for Molecular Matched Pairs analyses for many years. This has
typically used the complementary fragment capability available within [dicer](/docs/Molecule_Tools/dicer.md).
With complementary fragment generation enabled in `dicer`, each fragment is written along with the other
fragments that when combined, reconstruct the starting molecule. With a single cut fragment both the
fragment and the complementary fragment will be single fragment molecules. With a two cut fragment,
the fragment will have one fragment, but the complement will have two.

These fragments and complementary fragments can be stored in database systems and used for Matched Pairs
analyses.

In the simplest form, a single isotopic label is applied to each fragment as it is generated. Single
cut fragments can be recreated without ambiguity, but multi-cut fragments are not. For example if
'Fc1ncc(Cl)cc1' is diced into '[1cH]1ncc[1cH]cc1' and '[1FH].[1Cl]' there are two sites at which
the F and Cl atoms can be re-attached. This might be just fine.

# Subsituent Identification
Any discussion of Molecular Matched Pairs needs to mention the concept of a "local matched pair" as
an alternate way of thinking about this problem. [substituent_identification](/docs/Molecule_Tools/substituent_identification.md)
identifies substituents in a set of molecules, and stores the local context(s) in which fragments are found.
When a new molecule is examined, it can remove existing fragments, and retrieve fragments that have
been observed in that context, and suggest plausible new molecules.
