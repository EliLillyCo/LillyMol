# Markonikov

In this case a Halide is added across an unsaturated Carbon-Carbon bond. The rule is that the Halogen
goes to the Carbon atom with the more Hydrogen atoms. This example shows how that can be achieved.

We use the OR construct in LillyMol smarts.
```
  smarts: "[CR0H3]!@=[CR0H{1-2}]||[CR0H2]!@=[CR0H1]||[CR0H>0]!@=[CR0H>0]"
```
This says to first fine a CH3 bonded to a Carbon that is either CH1 or CH2.

If that is not found, look for a CH2 bonded to a Carbon that has 1 Hydrogen.

If that is not found, then find a C=C bond where both have at least one Hydrogen -
and in this case, they will have the have equal numbers. This query will match
twice per bond identified.

When the reaction is run it is important to enumerate the multiple hits
separately
```
trxn -m each -P markovnikov.rxn alkene.smi halide.smi
```
