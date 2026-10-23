# Retrosynthesis
This directory contains reactions for building the infrastructure needed by
the retrosynthesis tool [retrosynthesis](/docs/Molecule_Tools/retrosynthesis.md).

The reactions in this directory perform preparatory reactions on reagent files and
apply isotopic labels.

For example, by convention, primary amines are assigned isotope 2. So the
[reaction](primary_amine.rxn) merely identifies the nitrogen atom and places
an isotope on that atom.

And acid on the other hand must remove the OH Oxygen and place isotope 1
on the Carbon atom.o

The following numeric values are used. Note that numbers can be re-used.

1. Carboxyllic Acid
2. Primary Amine
2. Secondary Amine
3. Boronic Acid
4. Aryl Halide
5. Aldehyde and Ketone for reductive amination,  =O is removed.
