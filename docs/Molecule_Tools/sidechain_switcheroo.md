# sidechain_switheroo
This is a simple *de-novo* tool that can quickly generate large numbers
of molecules.

The tool works by first removing sidechains from molecules, and then re-attaching
those sidechains onto sites across all molecules in the input. The re-attachment
sites are the same as the sites from which sidechains were removed.

For this reason, it works by reading a set of input molecules into memory. Unlike
most LillyMol tools, it cannot work on a semi-infinite stream of molecules
read from stdin.

One or more queries are used to identify the bonds separating the 'scaffolds' from the
'sidechains'. By default this is `[aD3]-!@{a<10}[R0]`. This means that only
aromatic substituents are identified - by default. Also, substituents
are defined as having fewer than 10 atoms. Adjust to taste.

Any smarts can be specified via the -s option, with the first matched atom
assumed to be an atom in the scaffold, the second matched atom being an
atom in the sidechain, with the bond between them being broken.

During testing the following smarts were found to work well.
```
sidechain_switcheroo -s '[aD3]-!@{a<10}[R0]' -s '[CG0D3R0]-!@{a<10}[R0]' -h -z i -V -x 3 -v file.smi > new.smi
```
This breaks bonds between an aromatic atom, as well as substituents attached to 
three connected, fully saturated, non ring carbon atoms. Again, adjust to taste.

Note however that given just 10 random Chembl molecules as input, the command above
generated 5900 molecules in about 0.4 seconds. This tool needs to be used carefully
or else numbers may explode.  

Running 100 molecules, each with just 12 atoms, generates 2.25M molecules in 92
seconds. A further 126k were rejected for undesirable atom arrangements - see below.

Numbers can be reduced by lowering the -x option (see below) or reducing the
number of atoms allowed in a substituent - use a number lower than 10 in the
smarts above.

The `-h` option adds Hydrogen as a substituent, which means that some of the
molecules generated will have lost a substituent at a site where previously there
was a substituent.

The `-x` option is important. This specifies the maximum number of scaffold sites
that are simultaneously altered. This is one of the most significant sources of the 
combinatorics problems with this tool. By default, without the -x option, single
sidechains will be moved from one site to another - inter or intra molecular
transfers. The `-x` option specifies how many simultaneous transfers can be done.
The random 10 Chembl molecules used for testing, generate 14 unique sidechains, but
the number of ways one can select 3 from 14 is large (2184). Combine this with
the number of ways in which as many as 3 scaffold sites can be selected and
the combinatorics of this becomes obvious.

Output is straightforward, sidechains will be removed from one location and
moved to other sites - either in the same molecule or (more likely) in other molecules.
It will _not_ re-create starting molecules.
The tool is exhaustive and deterministic, so all possibilities will be generated.
Product molecules with adjacent heteroatoms, or Halogens attached to
aliphatic Carbon atoms are discarded. It is very likely that other filters should
be applied - pipe the output to tsubstructure or molecule_filter to gain speed from parallel
processing.

## Alternatives
Reactions can also accomplish sidechain switching. Consider the reaction
```
name: "same_ring"
comment: "identifies pairs of substituents, on the same ring, and swaps them."
scaffold {
  smarts: "[/IWrid1D3R1]-!@{a<10}[R0D<3].[/IWrid1D3R1]-!@{a<10}[R0D<3]"
  break_bond {
    a1: 0
    a2: 1
  }
  break_bond {
    a1: 2
    a2: 3
  }
  make_bond {
    a1: 0
    a2: 3
    btype: SS_SINGLE_BOND
  }
  make_bond {
    a1: 2
    a2: 1
    btype: SS_SINGLE_BOND
  }
}
```
This identifies pairs of substituents on a ring and swaps them. When run as
```
trxn -d -m each -P /path/to/reaction.rxn file.smi
```
it will generate a new product for each pair found. Note that if the
two substituents are the same, the starting molecule will be reproduced. The
duplicate suppression option `-d` will not detect that, it only considers
product molecules. The suppress symmetry related matches option to trxn, -k
may help, but will not solve this issue.

This can be extended to ring systems, rather than rings, by changing
/IWrid1 to /IWfsid1 in the smarts above.

It could also be made more 'adventurous' by removing the constraint
on what is attached. For example `[/IWrid1D3R1]-!@*` would consider
every exocyclic bond as defining a "substituent".

Note too that if implicit Hydrogens are converted to explicit, then
substituents will move to sites that previously did not have
a substituent
```
trxn -d -J exph -m each -P /path/to/reaction.rxn file.smi
```
Using the example above, rather than 500 products generated, we now get
2500.

Theoretically one could combine all molecules into a single, multi
fragment molecule,
```
mkfrag -R : file.smi > one.smi
```
and use this reaction (note the 3 atom limit on substituents)
```
name: "any_ring"
comment: "identifies pairs of substituents and swaps them."
scaffold {
  smarts: "[D3R1]-!@{a<4}[R0D<3].[D3R1]-!@{a<4}[R0D<3]"
  break_bond {
    a1: 0
    a2: 1
  }
  break_bond {
    a1: 2
    a2: 3
  }
  make_bond {
    a1: 0
    a2: 3
    btype: SS_SINGLE_BOND
  }
  make_bond {
    a1: 2
    a2: 1
    btype: SS_SINGLE_BOND
  }
}
```
which in the case of 200 random Chembl molecules, generates a 'molecule' with
3600 atoms -small. The reaction finds 48k matches to the query! It then takes 43 seconds
to write 48k product molecules.

Once the product is formed, regenerate individual molecules with mkfrag
```
mkfrag trxn.smi > split.smi
```

This approach does not scale - but for very small sets of molecules, may be useful.

# Swapping Rings.
Currently there is no LillyMol tool to do ring switching within a set of molecules.
We can however use the ring replacement tools to approximate this in various ways.

The first task is to use ring_extraction to identify and classify the rings in
the set of molecules. During testing this command proved interesting.
```
ring_extraction -X list=_RINGS.txt -v -S /tmp/RINGS -c file.smi
```
which generates the usual extracted rings data files, RINGS_*, as well as a
file containing a list of all the files generated, _RINGS.txt.

Then use ring_replacement, or ring_replacement_inexact to swap rings - and ring
systems.
```
ring_replacement_inexact -V -v -R F:_RINGS.txt -s '[R]' file.smi
```
we have requested that any ring/ring system be replaced with any/all of the
rings where replacement is feasible. Again, this can quickly generate quite
large numbers of molecules - many of which may not necessarily be reasonable
variants of the starting molecule.

A specialised tool to do this may be desirable - in order to control the numbers
of molecules generated, possibly by limiting the types of ring swaps that
could be attempted. And maybe separately handling 'terminal' and 'internal' rings.

# Swapping Fragments.
Use [safe_generate](/docs/Molecule_Tools/SAFE.md) although that tool is
stochastic. For systematic interiour fragment replacement, use dicer to 
identify plausible fragments, and write a reaction to remove the old
fragment and insert the new.
