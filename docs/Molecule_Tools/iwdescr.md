# iwdescr

`iwdescr` computes a broad set of interpretable 2D molecular descriptors. The
descriptors characterize molecular size, composition, connectivity, ring
systems, flexibility, polarity, hydrogen bonding, symmetry, and related
structural properties.

The tool generates approximately 280 descriptors when all optional descriptor
groups are enabled. These descriptors are designed to be fast to compute and
have been used extensively in QSAR models and descriptor-importance analyses.

## Quick start

```shell
iwdescr.sh file.smi > file.w
```

This reads molecules from `file.smi` and writes a space-separated descriptor
table to `file.w`.

Use the `iwdescr.sh` wrapper provided in `contrib/bin` for normal descriptor
generation. The wrapper supplies the standard preprocessing and descriptor selection defaults.
Invoking `iwdescr` directly is seldom recommended.

For larger datasets use iwdescr_parallel.sh
```
iwdescr_parallel.sh -thr 16 -S file file.smi 
```
which runs iwdescr.sh across 16 threads and creates "file.w". 

### Wrapper defaults

The wrapper removes explicit hydrogen atoms with `-g all` and enables all optional
descriptor groups with `-O all`.

It also uses `-u 0`, which writes `0` for any descriptor value that is not
computed because the calculation is not applicable to the molecule. This is a
general missing-value policy, not specific to any particular descriptor.


## Input and output

### Input

`iwdescr` uses the standard LillyMol molecule input mechanism and accepts any
structure file format supported by LillyMol. The input type is normally inferred
from the filename suffix and can also be specified using the standard LillyMol
input options.

Each molecule should have a name. If a molecule has no name, `iwdescr` assigns
an arbitrary identifier. If the name contains multiple whitespace-separated
tokens, only the first token is written to the descriptor output.

### Output

By default, `iwdescr` writes a space-separated descriptor table to standard
output. The first row contains the column names. Each subsequent row contains
the molecule identifier followed by its descriptor values.

Abbreviated output has the following form:

```text
Name w_natoms w_nrings w_amw ...
molecule1 12 1 180.20 ...
molecule2 18 2 254.29 ...
```


### Output separator

The default output separator is a space. Use `-B sep=VALUE` to select a
different separator.

For comma-separated output:

```shell
iwdescr.sh -B sep=, file.smi > file.csv
```

For tab-separated output:

```shell
iwdescr.sh -B sep=tab file.smi > file.tsv
```

## Selecting descriptors

Some descriptor groups are optional because they require additional computation
or are not needed in every workflow. Optional descriptor groups are controlled
with the `-O` option.

Use `-O all` to enable all optional groups:

```shell
iwdescr.sh -O all file.smi > file.w
```

Use `-O none` to disable all optional groups:

```shell
iwdescr.sh -O none file.smi > file.w
```

Options are processed from left to right, so individual groups can be enabled
after `-O none`:

```shell
iwdescr.sh -O none -O symm file.smi > file.w
```

An enabled group can also be disabled by prefixing its name with `-`:

```shell
iwdescr.sh -O all -O -symm file.smi > file.w
```

Descriptor groups can be enabled or disabled as units. Individual descriptors
within a group cannot be independently selected. Unwanted columns can instead
be removed from the resulting descriptor file, for example with `iwcut`.

The choice of optional groups can materially affect performance. In one
representative benchmark, `-O none` required approximately 21% of the runtime of
`-O all`. Enabling symmetry descriptors after `-O none` increased this to
approximately 37%. These figures illustrate relative cost only; actual runtime
depends on molecular complexity, hardware, build configuration, and enabled
descriptor groups.

Some options alter descriptor computation without adding columns. Existing option
behavior is retained for compatibility with established scripts.

Use `-O help` to display the available groups.

| Group | Descriptors enabled |
| ----- | ------------------- |
| `adjring` | Bonds adjacent to ring fusions |
| `bbr` | Bonds between ring systems |
| `charge` | Formal-charge descriptors |
| `chiral` | Enables more expensive chirality perception for `nchiral`; does not add columns |
| `complex` | Fused-ring complexity descriptors |
| `crowd` | Atomic crowding descriptors |
| `dm` | Distance-matrix and molecular-shape descriptors |
| `donacc` | Donor and acceptor descriptors |
| `hbond` | Legacy compatibility option associated with hydrogen-bond processing; does not independently add columns |
| `shbond` | Simplified hydrogen-bond descriptors |
| `lcc` | Long carbon-chain descriptors |
| `ncon` | Atomic connectivity-count descriptors |
| `pbond` | Polar-bond descriptors |
| `psa` | Novartis topological polar surface area |
| `psymm` | Partial-symmetry descriptors |
| `ramey` | Element-count and related descriptors |
| `rcj` | Ring-chain junction descriptors |
| `rfuse` | Ring-fusion descriptors |
| `rss` | Ring-substitution distance descriptors |
| `rssr` | Ring-substitution ratio descriptors |
| `satchain` | Saturated-chain descriptors |
| `spch` | Spinach descriptors |
| `spcgrp` | Specific functional-group descriptors |
| `symm` | Molecular-symmetry descriptors |
| `alogp` | ALogP |
| `xlogp` | XLogP |

The `nchiral` descriptor is always computed using fast chirality perception.
Enabling `-O chiral` applies more expensive chirality perception when calculating
that descriptor, but does not add another output column.

The `-b` option enables donor/acceptor distance descriptors and specifies the
minimum feature separation. These descriptors also require donor/acceptor
assignment. This historical interface is retained for compatibility with
established scripts.

## Descriptor conventions

A connection is a bond to another atom. For example, the carbon atom in methane has no connections, while each carbon atom in ethane has one connection. A chain atom is an atom that is not part of a ring.

Descriptor names beginning with `f` usually represent fractions, but their
denominators depend on the property being measured. Each definition identifies
the relevant denominator where it is not evident from the corresponding count.
If there are no instances of the numerator, the fraction is reported as zero.
Fractions are not evaluated with a zero denominator.

## Descriptors

Descriptors are grouped below by their primary chemical interpretation.
The `Group` column identifies the `-O` option that enables an optional
descriptor; `Always` indicates that the descriptor is always computed.

### Molecular size and composition

| Name | Definition | Group |
| ---- | ---------- | ----- |
| natoms | The number of atoms in the molecule. | Always |
| nelem | Number of different elements in the molecule. | Always |
| amw | Average molecular weight. | Always |
| platt | Total molecular connectivity. | Always |
| htroatom | Number of heteroatoms. | Always |
| htroaf | Fraction of atoms that are heteroatoms. | Always |
| hcount | Total number of hydrogen atoms. | Always |
| hperatom | Average number of hydrogens per heavy atom. | Always |
| halogen | The number of halogen atoms in the molecule. | Always |
| halogena | The number of atoms which have one or more halogens attached. | Always |
| bigatom | Number of atoms from beyond the 2nd row of the periodic table. | Always |
| fbigatom | Fraction of atoms from beyond the 2nd row of the periodic table. | Always |
| rmync | Number of carbon atoms. | `ramey` |
| rmynn | Number of nitrogen atoms. | `ramey` |
| rmyno | Number of oxygen atoms. | `ramey` |
| rmynf | Number of fluorine atoms. | `ramey` |
| rmyns | Number of sulfur atoms. | `ramey` |
| rmyncl | Number of chlorine atoms. | `ramey` |
| rmynbr | Number of bromine atoms. | `ramey` |
| rmyni | Number of iodine atoms. | `ramey` |
| heavy_halogen | Chlorine + bromine + iodine. | `ramey` |
| nrgnhlht | Number of non-ring, non-halogen heteroatoms. | Always |

### Atomic connectivity and bonding

| Name | Definition | Group |
| ---- | ---------- | ----- |
| ncon1 | Number of atoms with one connection. | `ncon` |
| ncon2 | Number of atoms with two connections. | `ncon` |
| ncon3 | Number of atoms with three connections. | `ncon` |
| ncon4 | Number of atoms with four connections. | `ncon` |
| fncon1 | Fraction of atoms with one connection. | `ncon` |
| fncon2 | Fraction of atoms with two connections. | `ncon` |
| fncon3 | Fraction of atoms with three connections. | `ncon` |
| fncon4 | Fraction of atoms with four connections. | `ncon` |
| frhc | Fraction of highly connected atoms in the molecule. Highly connected means two or more connections. | Always |
| avcon | Average number of connections per heavy atom. | Always |
| avchcon | Average connectivity of chain atoms. | Always |
| avalcon | Average connectivity of aliphatic atoms. | Always |
| mltbd | Number of non-aromatic, non-single bonds. | Always |
| fmltbd | Fraction of bonds that are non-aromatic and not single bonds. | Always |
| chmltbd | Number of non-aromatic, non-single bonds between chain atoms. | Always |
| fchmltbd | Fraction of all bonds that are non-aromatic, non-single bonds between chain atoms. | Always |
| rgmltbd | Number of non-aromatic, non-single ring bonds. | Always |
| frgmltbd | `rgmltbd` divided by the number of bonds in the molecule. | Always |
| dcca | Number of chain atoms with two connections. | Always |
| fdcca | Fraction of atoms that are chain atoms with two connections. | Always |
| ch2 | Number of CH2 groups. | Always |
| ch | Number of carbon atoms that have one or more hydrogens attached. | Always |
| d2sp3 | Number of fully saturated atoms with two connections. | Always |
| csp3 | Number of sp3 carbon atoms. | Always |
| fcsp3 | Fraction of atoms that are sp3 carbon. | Always |
| fccsp3 | Fraction of carbon atoms that are sp3. | Always |
| csp3_chain | Number of sp3 carbon atoms not in a ring. | Always |
| cd4ring | Number of carbon atoms with four connections in a ring. | Always |
| cd4chain | Number of carbon atoms with four connections outside rings. | Always |
| crowding | Crowding score: adjacent atoms that each have three or more connections contribute 1; such atoms separated by an atom with two connections contribute 0.5. | `crowd` |
| fcrowdng | Fraction of atoms that have more than two connections and at least one neighbour with more than two connections. | `crowd` |

### Rings and ring systems

Ring-size descriptors use the SSSR ring set. An isolated ring is not fused to
another ring. A fused ring system contains two or more fused rings. Terminal
(or peripheral) and internal (or junction) rings describe their connectivity
after the molecule is reduced to its ring scaffold.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| nrings | Number of rings in the molecule. | Always |
| ringatom | Number of atoms in rings. | Always |
| rngatmf | Fraction of atoms that are in rings. | Always |
| ringsys | Number of isolated rings plus fused ring systems. | Always |
| nrings3 | Number of SSSR rings containing 3 atoms. | Always |
| nrings4 | Number of SSSR rings containing 4 atoms. | Always |
| nrings5 | Number of SSSR rings containing 5 atoms. | Always |
| nrings6 | Number of SSSR rings containing 6 atoms. | Always |
| nrings7 | Number of SSSR rings containing 7 atoms. | Always |
| nrings8 | Number of SSSR rings containing 8 atoms. | Always |
| rng7atoms | Number of SSSR rings containing more than 7 atoms. | Always |
| srsz | Size, in atoms, of the smallest SSSR ring. | Always |
| lrsz | Size, in atoms, of the largest SSSR ring. | Always |
| lrsysz | Number of rings in the largest fused ring system; benzene gives 1 and naphthalene gives 2. | Always |
| mars | Maximum number of atoms in a ring system. | Always |
| nrsyscmr | Number of fused ring systems containing two or more rings. | Always |
| nnsssrng | Number of non-sssr rings. | Always |
| isolrc | Number of isolated, non-fused rings. | Always |
| isolhtrc | Number of isolated, non-fused heterocyclic rings. | Always |
| trmnlrng | Number of terminal or peripheral rings, having one connection to the remainder of the ring scaffold. | `spch` |
| intrnlrng | Number of internal or junction rings, having two or more connections to the remainder of the ring scaffold. | `spch` |
| alring | Number of aliphatic rings. | Always |
| arring | Number of aromatic rings. | Always |
| mhr | Maximum number of heteroatoms in a ring. | Always |
| mxhrf | Highest heteroatom fraction in a ring. | Always |
| mnhrf | Minimum heteroatom fraction in a ring. | Always |
| rhacnt | Number of heteroatoms in rings. | Always |
| rhaf | Fraction of ring atoms that are heteroatoms. | Always |
| frafus | Fraction of ring atoms that are involved in ring fusions. | Always |
| nspiro | Number of spiro joins. | `complex` |
| scaffoldbranches | Number of branches from the ring scaffold. | `spch` |

### Ring substitution and fusion

The `rssys1` through `rssys9` descriptors measure distances around the perimeter of fused ring systems. Strongly fused ring systems are excluded because the path between substitution points may not be uniquely defined.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| amrcj | Number of times a non-ring atom joins an aromatic ring. Includes singly connected atoms. | `rcj` |
| alrcj | Number of times a non-ring atom joins an aliphatic ring. Includes singly connected atoms. | `rcj` |
| rcj | Ring-chain join. | `rcj` |
| rchj | Ring-chain join to a heteroatom. | `rcj` |
| ringisol | For each non-fused ring, count the number of branches off the ring. Do not count terminal groups (like fluoro or methyl), although more complex terminal groups such as co2h and no2 are counted as branches off the ring rather than terminal groups. For each ring, take 1/branches and sum over all rings in the molecule. | Always |
| rsarom1 | Number of pairs of substituents separated by 1 bond along an unfused aromatic ring. | `rss` |
| rsarom2 | Number of pairs of substituents separated by 2 bonds along an unfused aromatic ring. | `rss` |
| rsarom3 | Number of pairs of substituents separated by 3 bonds along an unfused aromatic ring. | `rss` |
| rsaliph1 | Ring substituents one bond apart on an aliphatic ring. | `rss` |
| rsaliph2 | Ring substituents two bonds apart on an aliphatic ring. | `rss` |
| rsaliph3 | Ring substituents three bonds apart on an aliphatic ring. | `rss` |
| rsaliph4 | Ring substituents four or more bonds apart on an aliphatic ring. | `rss` |
| rssys1 | Number of pairs of ring-system substitution points separated by 1 bond around the ring-system perimeter. | `rss` |
| rssys2 | Number of pairs of ring-system substitution points separated by 2 bonds around the ring-system perimeter. | `rss` |
| rssys3 | Number of pairs of ring-system substitution points separated by 3 bonds around the ring-system perimeter. | `rss` |
| rssys4 | Number of pairs of ring-system substitution points separated by 4 bonds around the ring-system perimeter. | `rss` |
| rssys5 | Number of pairs of ring-system substitution points separated by 5 bonds around the ring-system perimeter. | `rss` |
| rssys6 | Number of pairs of ring-system substitution points separated by 6 bonds around the ring-system perimeter. | `rss` |
| rssys7 | Number of pairs of ring-system substitution points separated by 7 bonds around the ring-system perimeter. | `rss` |
| rssys8 | Number of pairs of ring-system substitution points separated by 8 bonds around the ring-system perimeter. | `rss` |
| rssys9 | Number of pairs of ring-system substitution points separated by 9 bonds around the ring-system perimeter. | `rss` |
| frsub | Fraction of ring atoms that are substituted outside the ring. | `rssr` |
| frssub | Fraction of ring atoms that have a single-atom substituent. | `rssr` |
| alorthoring | Number of ortho substituents on an aliphatic ring. | `rssr` |
| arorthoring | Number of ortho substituents on an aromatic ring. | `rssr` |
| bbr1 | One non-ring bond between two rings. | `bbr` |
| bbr2 | Two non-ring bonds between two rings. | `bbr` |
| bbr3 | Three non-ring bonds between two rings. | `bbr` |
| bbr4 | Four non-ring bonds between two rings. | `bbr` |
| bbr5 | Five non-ring bonds between two rings. | `bbr` |
| bbr6 | Six non-ring bonds between two rings. | `bbr` |
| sboradjf | Exocyclic single bonds (to terminal group) adjacent to a ring fusion. | `adjring` |
| dboradjf | Exocyclic double bonds (to terminal group) adjacent to a ring fusion. | `adjring` |
| nsfsdsys | Number of strongly fused ring systems found. A strongly fused system contains rings that share more than one bond between them. It includes only those rings that are strongly fused to one or more neighbours. | `complex` |
| rnginsfs | Rings in strongly fused systems. | `complex` |
| lgstrfsy | Largest strongly fused system size. | `complex` |
| htrcsfsy | Heterocycles in strongly fused systems. | `complex` |
| mxhtsfsy | Max heteroatoms in a strongly fused system. | `complex` |
| npfsdsys | Number of planar fused ring systems found. A planar fused system contains rings that share at most one bond between adjacent rings. | `complex` |
| rnginpfs | Rings in planar fused systems. | `complex` |
| lgplnfsy | Largest planar fused system size. | `complex` |
| htrcpfsy | Heterocycles in planar fused systems. | `complex` |
| mxhtpfsy | Max heteroatoms in a planar fused system. | `complex` |
| al5 | Aliphatic rings of size 5. | `rfuse` |
| al6 | Aliphatic rings of size 6. | `rfuse` |
| ar5 | Aromatic rings of size 5. | `rfuse` |
| ar6 | Aromatic rings of size 6. | `rfuse` |
| fsdrng5l5l | Number of fused pairs consisting of a 5-membered aliphatic ring and a 5-membered aliphatic ring. | `rfuse` |
| fsdrng5l5r | Number of fused pairs consisting of a 5-membered aliphatic ring and a 5-membered aromatic ring. | `rfuse` |
| fsdrng5l6l | Number of fused pairs consisting of a 5-membered aliphatic ring and a 6-membered aliphatic ring. | `rfuse` |
| fsdrng5l6r | Number of fused pairs consisting of a 5-membered aliphatic ring and a 6-membered aromatic ring. | `rfuse` |
| fsdrng5r5r | Number of fused pairs consisting of a 5-membered aromatic ring and a 5-membered aromatic ring. | `rfuse` |
| fsdrng5r6l | Number of fused pairs consisting of a 5-membered aromatic ring and a 6-membered aliphatic ring. | `rfuse` |
| fsdrng5r6r | Number of fused pairs consisting of a 5-membered aromatic ring and a 6-membered aromatic ring. | `rfuse` |
| fsdrng6l6l | Number of fused pairs consisting of a 6-membered aliphatic ring and a 6-membered aliphatic ring. | `rfuse` |
| fsdrng6l6r | Number of fused pairs consisting of a 6-membered aliphatic ring and a 6-membered aromatic ring. | `rfuse` |
| fsdrng6r6r | Number of fused pairs consisting of a 6-membered aromatic ring and a 6-membered aromatic ring. | `rfuse` |
| fsdrngalal | Number of fused pairs consisting of two aliphatic rings. | `rfuse` |
| fsdrngalar | Number of fused pairs consisting of an aliphatic ring and an aromatic ring. | `rfuse` |
| fsdrngarar | Number of fused pairs consisting of two aromatic rings. | `rfuse` |

### Aromaticity and unsaturation

An electron-rich area is a connected region of adjacent atoms with pi electrons.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| aroma | Aromatic atoms. | Always |
| aromha | Aromatic heteroatoms. | Always |
| aromc | Number of aromatic carbon atoms. | Always |
| aliphc | Aliphatic carbon count. | Always |
| fraromha | Number of aromatic heteroatoms divided by the number of aromatic atoms. | Always |
| aromdens | Number of aromatic atoms divided by `natoms`. | Always |
| atmpiele | Number of atoms with pi electrons. | Always |
| fratmpie | Fraction of atoms with pi electrons. | Always |
| unsatura | Number of non-aromatic atoms that are unsaturated. | Always |
| funsatura | Fraction of atoms that are unsaturated. | Always |
| nconjgsc | Number of conjugated sections in the molecule, including aromatic atoms. | Always |
| atincnjs | The total number of atoms in the conjugated sections. | Always |
| mxcnjscz | Largest number of atoms in a conjugated section. | Always |
| cinconjs | Number of carbon atoms in conjugated sections. | Always |
| erichsct | Number of separate electron-rich areas. | Always |
| aiercsct | Number of atoms in electron-rich areas. | Always |
| faiercst | Fraction of atoms in electron-rich areas. | Always |
| lercsct | Number of atoms in the largest electron-rich area. | Always |
| numcdb | Number of non-aromatic, doubly bonded carbon atoms. | Always |
| totdbsub | Total number of substituents on doubly bonded, non-aromatic carbon atoms. | Always |
| avcdbsub | Average substitution on doubly bonded, non-aromatic carbon atoms. | Always |
| dvinylb | Number of non-aromatic single bonds whose endpoint atoms are each incident to a multiple bond (`*=*-*=*`). | `pbond` |

### Flexibility and chains

A flexible chain is a chain region whose bonds are likely to be rotatable. "Spinach" refers to atoms outside the molecular scaffold. The term is used historically within LillyMol but is not universal.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| rotbond | Number of rotatable bonds in the molecule. | Always |
| frotbond | Fraction of bonds that are rotatable. | Always |
| nflxchn | The number of separate flexible chains in the molecule. | Always |
| atflxchn | Number of atoms involved in flexible chains. | Always |
| faflxchn | Number of atoms in flexible chains divided by `natoms`. | Always |
| fnflxchn | Number of atoms in flexible chains divided by the number of non-ring atoms. | Always |
| lflxchn | The longest flexible chain in the molecule. | Always |
| avflxchn | Average number of atoms in a flexible chain. | Always |
| rkentrpy | Entropic measure of flexibility - defined by radha. | Always |
| mxlencchain2 | Maximum length of an all [CD2] chain with no branching. | `lcc` |
| mxlencchain3 | Maximum length of an all [CD2] chain with at most [CD3] as a branch point. | `lcc` |
| nsatchain | Number of fully saturated chain regions. | `satchain` |
| mxsatchain | Number of atoms in the largest saturated chain region. | `satchain` |
| fsatchain | Number of atoms in saturated chain regions divided by `natoms`. | `satchain` |
| frspch | Fraction of atoms outside the molecular scaffold. | `spch` |
| spchtro | Number of heteroatoms outside the molecular scaffold. | `spch` |
| rbfrspch | Fraction of bonds outside the molecular scaffold that are rotatable. | `spch` |
| scafrotb | Number of rotatable bonds whose two endpoint atoms are in the molecular scaffold. | `spch` |
| spchrotb | Number of rotatable bonds not wholly contained in the molecular scaffold. This includes rotatable bonds wholly in spinach and rotatable bonds joining spinach to scaffold atoms. | `spch` |
| maxrotbgrp | Largest connected component of rotatable bonds, counted as bonds. Adjacent rotatable bonds are connected when they share an atom. | `spch` |
| nrnspch | Number of non-ring, non-spinach atoms. Chain atoms joining rings. | `spch` |
| fnrnspc | Fraction of non-spinach atoms that are non-ring. | `spch` |
| satspcha | Number of saturated atoms outside the molecular scaffold. | `spch` |
| unsatspcha | Number of unsaturated atoms outside the molecular scaffold. | `spch` |
| fsatspcha | Fraction of atoms outside the molecular scaffold that are saturated. | `spch` |
| rng2spch | Number of ring connections to atoms outside the molecular scaffold. | `spch` |
| rng2bridge | Ring connection to chain scaffold atom. | `spch` |

### Hydrogen bonding, polarity, and charge

The Bruns formal-charge and donor/acceptor assignments are based on work by
Robert F. Bruns at Lilly. For donor/acceptor distance descriptors, values are
reported as zero when a molecule does not contain the required donor or acceptor
features. A bond is considered non-polar when its endpoint atoms have the same
atomic number; otherwise it is considered polar.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| ohsh | Number of oxygen and sulfur atoms with an attached hydrogen atom. | Always |
| ro5_ohnh | Rule-of-Five count of hydrogen atoms attached to oxygen or nitrogen atoms. | Always |
| ro5_on | Rule-of-Five count of oxygen and nitrogen atoms. | Always |
| hacts | Hydrogen bond score (simplistic). | `shbond` |
| hdons | Hydrogen bond score (simplistic). | `shbond` |
| hduals | Hydrogen bond score (simplistic). | `shbond` |
| brunsacc | Number of hydrogen-bond acceptors identified by the Bruns assignment. | `donacc` |
| brunsdon | Number of hydrogen-bond donors identified by the Bruns assignment. | `donacc` |
| brnsdual | Number of sites identified by the Bruns assignment as both donor and acceptor. | `donacc` |
| brunshbdsum | `brunsacc + brunsdon - brnsdual`. | `donacc` |
| brunspos | Number of positive formal charges assigned by the Bruns method. | `charge` |
| brunsneg | Number of negative formal charges assigned by the Bruns method. | `charge` |
| formal_charge | Number of formally charged atoms (brunspos + brunsneg); this is not the net formal charge. | `charge` |
| nplus | Number of positive formal charges. | `donacc` |
| nminus | Number of negative formal charges. | `donacc` |
| aamind | Minimum bond separation between acceptors. | `donacc` and `-b` |
| aa2mind | Second minimum bond distance between acceptors. | `donacc` and `-b` |
| aaave | Average bond separation between acceptors. | `donacc` and `-b` |
| admind | Minimum bond separation between acceptor and donor. | `donacc` and `-b` |
| ad2mind | Second minimum bond separation between acceptor and donor. | `donacc` and `-b` |
| adave | Average number of bonds between acceptor and donor. | `donacc` and `-b` |
| ddmind | Minimum bond separation between donors. | `donacc` and `-b` |
| dd2mind | Second minimum bond separation between donors. | `donacc` and `-b` |
| ddave | Average bond separation between donors. | `donacc` and `-b` |
| pbcount | Number of bonds whose endpoint atoms have different atomic numbers. | `pbond` |
| frpbond | Number of polar bonds divided by the total number of bonds. | `pbond` |
| nonpbond | Number of bonds whose endpoint atoms have the same atomic number. | `pbond` |
| pbarom | Number of polar bonds within aromatic rings. | `pbond` |
| npbarom | Number of non-polar bonds within aromatic rings. | `pbond` |
| pbunset | Number of polar bonds in unsaturated regions. The historical name is a misspelling of `pbunsat`. | `pbond` |
| internalhbd | Number of possible internal hydrogen bonds, based on donor-acceptor bond separation and the rotatable bonds between each pair. | Always |
| nvrtspsa | Novartis topological polar surface area (J. Med. Chem. 2000, 43, 3714-3717). | `psa` |

### Topological distance and shape

Distances in this section are shortest through-bond distances, not spatial
distances. The eccentricity of an atom is its greatest shortest-path distance
to any other atom. The molecular diameter is the maximum atomic eccentricity,
and the molecular radius is the minimum atomic eccentricity.

A molecule may contain multiple longest-path instances. All such paths
contribute to descriptors derived from the longest path. Central atoms are
equidistant from the endpoints of a longest path; when the path has two central
atoms, both are processed. Shell descriptors measure topological atom density
as the distance from a starting atom is increased.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| mxdst | Molecular diameter: the maximum shortest through-bond distance between two atoms. | `dm` |
| fmxdst | `mxdst / natoms`; larger for elongated molecules and smaller for compact molecules. | `dm` |
| muldiam | Number of longest-path instances having length `mxdst`. | `dm` |
| rad | Molecular radius: the minimum atomic eccentricity. | `dm` |
| mulrad | Number of atoms whose eccentricity equals `rad`. | `dm` |
| ishape | Compactness measure: (maximum eccentricity - minimum eccentricity) / maximum eccentricity. | `dm` |
| weiner | Wiener index derived from the molecular distance matrix. | Always |
| harary | Harary index derived from the molecular distance matrix. | `dm` |
| maxdarom | Maximum shortest through-bond distance between aromatic atoms. | `dm` |
| maxdrng | Maximum shortest through-bond distance between ring atoms, which need not be in the same ring. | `dm` |
| maxdhtro | Maximum shortest through-bond distance between heteroatoms. | `dm` |
| maxdons | Maximum shortest through-bond distance between oxygen, nitrogen, or sulfur atoms. | `dm` |
| avebbtwn | Average shortest through-bond distance between atom pairs. | `dm` |
| normbbtwn | `avebbtwn / natoms`. | `dm` |
| compact | `1 - (mxdst / natoms)`; larger values indicate more compact molecules. | `dm` |
| nolp | `natoms - mxdst - 1`; a compactness measure despite the historical name. | `dm` |
| avsdlp | Average shortest distance of an atom from a longest path, accumulated over longest-path instances. | `dm` |
| mxsdlp | Maximum shortest distance of an atom from a longest path. | `dm` |
| mxsdlprl | `mxsdlp / mxdst`. | `dm` |
| mdallp | Mean shortest distance of atoms from the longest-path instances. | `dm` |
| fmdallp | `mdallp / mxdst`. | `dm` |
| fdiffallp | Measure of anisotopy between the atoms at the ends of the longest path. | `dm` |
| centre3 | Number of atoms within three bonds of a central atom. | `dm` |
| centre3h | Number of heteroatoms within three bonds of a central atom. | `dm` |
| stddcentre | Standard deviation of atom distances from the central atom or atoms. | `dm` |
| avdcentre | Average atom distance from the central atom or atoms. | `dm` |
| cntrdgncy | Number of equivalent central atoms. | `dm` |
| cntrdshell1 | Number of atoms within one bond of a central atom. | `dm` |
| cntrdshell2 | Number of atoms within two bonds of a central atom. | `dm` |
| cntrdshell3 | Number of atoms within three bonds of a central atom. | `dm` |
| aveshell1 | Average number of atoms in radius-one shells over all starting atoms. | `dm` |
| aveshell2 | Average number of atoms in radius-two shells over all starting atoms. | `dm` |
| aveshell3 | Average number of atoms in radius-three shells over all starting atoms. | `dm` |
| maxshell3 | Maximum number of atoms in any radius-three shell. | `dm` |
| mh3b | Maximum number of heteroatoms within three bonds of any starting atom. | `dm` |
| tg3 | Terminal groups separated by 3 bonds. | `dm` |
| tm | Terminal methyl groups. | `dm` |

### Symmetry

Molecular symmetry is determined during canonical atom ordering. Two or more
atoms participate in a symmetry relationship when they belong to the same
canonical symmetry class. Partial symmetry describes atoms whose local
environments are equivalent through one or more topological shells but become
different at a larger radius.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| symmatom | Number of atoms belonging to a canonical symmetry class containing at least two atoms. | `symm` |
| fsymmatom | `symmatom / natoms`. | `symm` |
| lsepsymatom | Maximum shortest through-bond distance between atoms in the same canonical symmetry class. | `symm` |
| flsepsymatom | `lsepsymatom / natoms`. | `symm` |
| maxsymmclass | Largest number of atoms in a canonical symmetry class. | `symm` |
| maxpsymd | Largest topological radius through which any pair of atoms remains partially symmetric. | `psymm` |
| fmaxpsymd | `maxpsymd / natoms`. | `psymm` |
| maxpsymdmean | Mean partial-symmetry radius over all atoms, including zero for atoms with no partial-symmetry relationship. | `psymm` |
| psymdnzero | Number of atoms with no partial-symmetry relationship. | `psymm` |

### Specific functional groups

An exocyclic atom is a non-ring atom directly bonded to a ring atom.

| Name | Definition | Group |
| ---- | ---------- | ----- |
| co2h | Number of carboxylic acid groups. | Always |
| amine | Very poor counting of amines - need to change this sometime. | `spcgrp` |
| pyridine | Number of pyridine-like aromatic nitrogen atoms without hydrogen. | `spcgrp` |
| pyrrole | Number of pyrrole-like aromatic nitrogen atoms bearing hydrogen. | `spcgrp` |
| nchiral | Number of explicitly specified chiral centres. | Always |
| excybond | Number of bonds between ring atoms and exocyclic atoms. | Always |
| excydbond | Number of double bonds between ring atoms and exocyclic atoms. | Always |
| excydscon | Number of exocyclic single bonds to singly connected atoms. | Always |
| excydsconh | Number of exocyclic single bonds to singly connected heteroatoms. | Always |
| excydscondon | Number of exocyclic bonds to heteroatoms bearing hydrogen. | Always |

### LogP and related properties

| Name | Definition | Group |
| ---- | ---------- | ----- |
| alogp | Wildman-Crippen atom-contribution estimate of logP (J. Chem. Inf. Comput. Sci. 1999, 39, 868-873). | `alogp` |
| xlogp | Wang-Fu-Lai atom-contribution estimate of logP (J. Chem. Inf. Comput. Sci. 1997, 37, 615-621). | `xlogp` |
| cmr | Historical atom-contribution estimate of molar refractivity; the original source is unknown. | Always |
| acmbe | Andrews-Craik-Martin atom-contribution binding-energy score; the original source has not been confirmed. | `charge` |
| obalance | Oxygen balance percentage: `-1600 * (2C + H/2 - O) / molecular weight`. | `ramey` |

## Descriptor names

### Prefix

Output descriptor names include a prefix that identifies their provenance. By
default, `iwdescr` uses `w_`, producing names such as `w_natoms` and `w_amw`.

Use `-B prefix=VALUE` to change the prefix. The supplied value is concatenated
directly with the internal descriptor name; no underscore or other separator is
added automatically.

```shell
-B prefix=Q
```

produces:

```text
Qnatoms Qamw
```

whereas:

```shell
-B prefix=Q_
```

produces:

```text
Q_natoms Q_amw
```

Use `-B prefix=none` to omit the prefix:

```text
natoms amw
```

Prefixes affect output column names and provide a convenient way to distinguish
descriptors generated by different tools or configurations.

### Name translation

Internal descriptor names are fixed, but selected names can be translated when
the output header is written. Supply a textproto translation file using:

```shell
-B namexref=names.textproto
```

Each entry maps an internal descriptor name to an output name:

```textproto
feature {
  computed_name: "amw"
  name: "molecular_weight"
}
```

Here, `computed_name` is the internal descriptor name and must be specified
without the output prefix. `name` is the replacement used in the output header.
The example changes `w_amw` to `w_molecular_weight`.

Any number of translations can be specified. Descriptors not present in the
translation file retain their internal names. Entries whose `computed_name`
does not correspond to a descriptor are ignored. The optional `description`
field is currently ignored.

Name translation is applied only when writing the output header. Internal
operations such as filtering continue to use the internal descriptor name.

Prefix selection and name translation can be combined:

```shell
-B prefix=none -B namexref=names.textproto
```

This produces `molecular_weight` rather than `w_amw`.

## Filtering

`iwdescr` can filter molecules using computed descriptor values. The `-F`
option writes a SMILES file containing the molecules that pass all filters; it
does not write their descriptor rows.

Given:

```text
C methane
CC ethane
CCC propane
C1CC1 cyclopropane
```

the command:

```shell
iwdescr.sh -F natoms.gt.2 file.smi
```

writes:

```text
CCC propane
C1CC1 cyclopropane
```

A filter has the form:

```text
descriptor.operator.value
```

| Operator | Meaning |
| -------- | ------- |
| `eq` | Equal |
| `ne` | Not equal |
| `lt` | Less than |
| `le` | Less than or equal |
| `gt` | Greater than |
| `ge` | Greater than or equal |

Multiple `-F` options are combined using logical AND:

```shell
iwdescr.sh -F natoms.ge.10 -F nrings.gt.0 file.smi
```

Filters use internal descriptor names. The output prefix is optional, and names
introduced through `-B namexref=...` are not recognized in filter expressions.

Only descriptors enabled for the current invocation can be used in filters. For
example, this fails because `mxdst` belongs to the disabled `dm` group:

```shell
iwdescr.sh -O none -F mxdst.lt.30 file.smi
```

Enable the required group before applying the filter:

```shell
iwdescr.sh -O none -O dm -F mxdst.lt.30 file.smi
```

### Filtering molecules versus descriptor rows

To retain descriptor rows rather than write SMILES, filter the generated
descriptor table with `dfilefilter`:

```shell
iwdescr.sh file.smi | dfilefilter -e 'w_natoms<40' - > passing.dat
```

| Task | Tool | Example |
| ---- | ---- | ------- |
| Filter molecules and write SMILES | `iwdescr -F` | `natoms.lt.40` |
| Filter descriptor rows | `dfilefilter` | `w_natoms<40` |

The expression syntax differs for historical reasons. `iwdescr -F` uses internal
descriptor names, whereas `dfilefilter` operates on the actual column names in
the generated table, including prefixes or translated names.

## Fingerprinting

`iwdescr` can encode numeric descriptors as a counted molecular fingerprint.
Each selected descriptor is assigned one or more fingerprint bits. Its numeric
value is discretized using a descriptor-specific range, and the resulting
bucket is stored as the bit count.

The `iwdescr.sh` wrapper supplies the range profile:

```text
${LILLYMOL_HOME}/data/chembl.ranges
```

This textproto file contains descriptor ranges derived from the central 99% of
values observed in a sample of one million ChEMBL molecules. Values outside a
descriptor's range are assigned to its nearest endpoint.

> Fingerprints depend on the range profile. Changing the profile can change the
> fingerprint generated for the same molecule. Use the same profile when
> generating fingerprints that will be compared with one another.

### Selecting fingerprint descriptors

Use `-G ALL` to include every active descriptor:

```shell
iwdescr.sh -G ALL file.smi
```

Use `-G BEST` to include a subset selected from previous calibration studies:

```shell
iwdescr.sh -G BEST file.smi
```

Individual descriptors can be selected by name:

```shell
iwdescr.sh -G w_natoms,w_nrings file.smi
```

Only active descriptors can be fingerprinted. For example, this fails because
`mxdst` belongs to the disabled `dm` group:

```shell
iwdescr.sh -O none -G mxdst file.smi
```

Enable the group before requesting the descriptor:

```shell
iwdescr.sh -O none -O dm -G mxdst file.smi
```

### Resolution and replication

By default, each descriptor range is divided into 10 buckets. Use `-G R=n` to
change the resolution:

```shell
iwdescr.sh -G R=20 -G BEST file.smi
```

Replicating a descriptor assigns its encoded value to multiple fingerprint bits.
A positive integer sets the default replication count:

```shell
iwdescr.sh -G 4 -G BEST file.smi
```

A descriptor-specific replication count can be supplied after a colon:

```shell
iwdescr.sh -G w_natoms:10,w_nrings:20 file.smi
```

This assigns 10 replicated bits to `w_natoms` and 20 to `w_nrings`.

### Fingerprint output modes

When reading a structure file, fingerprint output is written as TDT records
containing the SMILES, molecule identifier, and fingerprint.

Use `-G FILTER` when processing an existing TDT stream. In this mode, `iwdescr`
adds the fingerprint to each input record rather than creating a new SMILES and
identifier record.

The same fingerprint can be generated through `gfp_make.sh` using `-W`:

```shell
gfp_make.sh -W natoms,nrings file.smi
```

## Other options
The -B option controls some fine grain details of the computation. The usage message is
```
 -B quiet               turn off unclassified atom warnings from TPSA computation.
 -B ranges=<fname>      read fingerprint ranges from <fname>. Input is textproto.
 -B sep=<char>          output separator, default space. -B sep=, or -B sep=tab, or...
 -B namexref=<fname>    w.Feature textproto of name translations. See docs for info.
 -B prefix=<s>          descriptor name prefix, 'w_' by default. Use `-B prefix=none` for no prefix.
 -B dstats              report individual statistics for each descriptor.
 -B flush               flush output after each molecule processed.
 -B float               integer values written as floats, "10" becomes "10."
 -B rpt=<n>             during testing, report progress every <n> molecules tested.
 -B wpipe               include smiles as first column of output.
 -B rpipe               input is a descriptor file with a leading smiles column.
```

### -B quiet
By default alog and xlogp calculations can display warning messages about unclassified
atoms. The `-B quiet` option suppresses those warnings.

### -B ranges=\<fname\>
This was previously mentioned under fingerprinting. The file is a `w:Ranges`, defined
in [iwdescr](src/Molecule_Tools/iwdescr.proto) and might look like
```
range { name:"w_natoms" min:5 max:50 }
range { name:"w_nrings" min:0 max:10 }
```
### -B sep=
Discusses previously under output separator.

### -B namexref=\<fname\>
Discussed previously under name translations.

### -B prefix=\<s\>
Discussed previously under output prefix.

### -B dstats
At the end of a run, report statistics for each descriptor.

### -B flush
Flush the output after each molecule written. Usually not needed.

### -B float
All values are written in decimal form, so "1" will be written as "1.".


### -B rpipe -B wpipe
Experimental, see below.


### Compatibility

The `-S` option enables the RDKit-compatible Novartis TPSA behavior. The default LillyMol implementation follows the published method as closely as possible. Differences are usually small because they are confined to relatively rare atom cases. Use `-S` when you need closer alignment with RDKit conventions or when reproducing RDKit-based workflows. Either mode may perform better depending on the application.

### Testing
For any cheminformatics tool it is important that results be invariant
to the ordering of the atoms. With the `-T` option, iwdescr performs random
smiles tests and determines whether or not the results are invariant across
different atom ordering. 

#### -T \<n\>

Generate <n> random smiles variants for testing.

#### -T kg

By default, processing terminates on a failure, use `-T kg` to keep going.

#### -T rpt=\<n\>
Report progress every <n> molecules tested.

#### -T WRITE=\<fname\>
Write failing molecules to <fname>. A suffix of ".smi" will be added if not
specified.

Across 2.7M Chembl there are about 3000 molecules that show failures. Many of
these involve complex fused ring systems where the Smallest Set of Smallest
Rings paradigm used by LillyMol is problematic - cubane has only 5 SSSR
rings. Others are due to deep seated bugs or design flaws in LillyMol, which are being
investigated. The number of failures is unlikely to ever reach zero.

### Experimental pipelines

`-B wpipe` writes records with a leading SMILES column so the output
can be piped into another LillyMol tool that accepts `rpipe` input.

Example:

```shell
iwdescr.sh -B wpipe in.smi | another_tool -B rpipe -
```
The header record will have a label "smiles" so the tabular nature of the
file is preserved
```
smiles id natoms ...
C methane 1 ...
CC ethane 2  ...
```

This mode is experimental and may change.
