# iwdescr
`iwdescr` computes molecular descriptors. There is a strong emphasis on interpretable
descriptors, that can be quickly computed.

Basic usage
```
iwdescr file.smi > file.w
```
which will generate a tabular (space separated) file containing a couple of hundred
molecular descriptors.

There is also a shell wrapper in the contrib/bin directory that sets up some useful
defaults. It is strongly recommended that the wrapper be used always. Otherwise you will
see various missing columns.

## Descriptors.
The following descriptors are computed.

| ---- | ---------- |
| name | definition |
| ---- | ---------- |
| natoms | the number of atoms in the molecule |
| nrings | number of rings in the molecule |
| nelem | number of different elements in the molecule |
| amw | average molecular weight |
| ncon1 | number of singly connected atoms |
| ncon2 | number of doubly connected atoms |
| ncon3 | number of three connected atoms |
| ncon4 | number of four connected atoms |
| fncon1 | fraction of atoms that are singly connected |
| fncon2 | fraction of atoms that are doubly connected |
| fncon3 | fraction of atoms that are three connected |
| fncon4 | fraction of atoms that are four connected |
| frhc | fraction of highly connected atoms in the molecule. highly connected means 2 or more connections. |
| avcon | the average number of (heavy atom) connections in the molecule. |
| mltbd | number of non-aromatic, non-single bonds. |
| fmltbd | fraction of bonds that are non-aromatic and not single bonds |
| chmltbd | non-aromatic, non-single chain bonds |
| fchmltbd | fraction of bonds that are non-aromatic and not single |
| rgmltbd | ring bonds that are non aromatic and not single |
| frgmltbd | rgmltbd / number of bonds in the molecule |
| dcca | doubly connected chain atoms |
| fdcca | fraction of the atoms that are doubly connected chain atoms |
| mxdst | longst through bond distance in the molecule |
| fmxdst |  mxdst / matoms |
| muldiam | number of instances for which mxdst occurs |
| rad | shortest longest bond count between atoms |
| mulrad | number of atoms for which rad occurs |
| rotbond | number of rotatable bonds in the molecule |
| rhacnt | ring heteroatoms |
| rhaf | fraction of ring atoms that are heteroatoms |
| frafus | fraction of ring atoms that are involved in ring fusions |
| rngatmf | fraction of atoms that are in rings |
| aroma | aromatic atoms |
| aromha | aromatic heteroatoms |
| fraromha | aromatic heteroatoms divided by number of ring atoms |
| aromdens | fraction of the atoms that are aromatic |
| ch2 | number of ch2 groups |
| ch | number of carbon atoms that have one or more hydrogens attached |
| htroaf | fraction of atoms that are heteroatoms |
| ohsh | oxygen or sulphur with a hydrogen attached |
| co2h | carboxyllic acid |
| amine | very poor counting of amines - need to change this sometime |
| hacts: | a composite hydrogen bond acceptor score - radha assigned relative strengths to various groups |
| hdons: | a composite hydrogen bond donor score - radha |
| hduals: | a composite score from groups which can be both hydrogen bond donors and acceptors. |
| mhr | maximum number of heteroatoms in a ring |
| mxhrf | highest heteroatom fraction in a ring |
| mnhrf | minimum heteroatom fraction in a ring |
| lrsysz | largest ring size |
| mars | most atoms in a ring system |
| frspch | fraction of atoms in the spinach - chain atoms outside rings |
| nchiral | number of chiral centres in the molecule. includes both explicitly marked centres and carbon and sulphur atoms that are actually chiral without being explicitly marked. |
| amrcj | number of times a non-ring atom joins an aromatic ring. includes singly connected atoms. measure of substitutions on aromatic rings |
| alrcj | number of times a non-ring atom joins an aliphatic ring. includes singly connected atoms |
| isolhtrc | number of isolated (not fused) heterocyclic ring |
| rhacnt | ring (aromatic and aliphatic) ring heteroatom count |
| qtsubc | number of carbons bonded to 4 different heteroatoms. |
| halogen | the number of halogen atoms in the molecule |
| halogena | the number of atoms which have one or more halogens attached |
| nrgnhtht | non-ring, non halogen heteroatoms |
| numcdb | number of non-aromatic, doubly bonded carbon atoms |
| totdbsub | total number of substituents on doubly bonded, non aromatic carbon atoms |
| avcdbsub | average substitution on doubly bonded, non aromatic carbon atoms |
| nringsn | number of rings of size n |
| nrings3 | number of rings of size 3 |
| nrings4 | number of rings of size 4 |
| nrings5 | number of rings of size 5 |
| nrings6 | number of rings of size 6 |
| nrings7 | number of rings of size 7 |
| nrings8 | number of rings of size 8 |
| nspiro | number of spiro joins |
| nconjgsc | number of conjugated sections in the molecule. does not include aromatic rings. includes adjacent doubly bonded atoms. |
| atincnjs | the total number of atoms in the conjugated sections. |
| mxcnjscz | largest number of atoms in a conjugated section |
| cinconjs | number of carbon atoms in conjugated sections |
| nsfsdsys | number strongly fused ring systems found. a strongly fused system contains rings that share more than 1 bond between them. includes only those rings that are strongly fused to one or more neighbours. |
| rnginsfs | rings in strongly fused systems. |
| lgstrfsy | largest strongly fused system size |
| htrcsfsy | heterocycles in strongly fused systems |
| mxhtsfsy | max heteroatoms in a strongly fused system |
| npfsdsys | number planar fused ring systems found. a planar fused system contains rings that share at most 1 bond between adjacent rings. |
| rnginpfs | rings in planar fused systems. |
| lgplnfsy | largest planar fused system size |
| htrcpfsy | heterocycles in planar fused systems |
| mxhtpfsy | max heteroatoms in a planar fused system |
| fbigatom | fraction of atoms from beyond the 2nd row of the periodic table |
| bigatom | number of atoms from beyond the 2nd row of the periodic table |
| ringsys | number of ring systems |
| pbarom | number of aromatic polar bonds - different atom types bonded in aromatic ring - pyridine etc... |
| npbarom | number of non-polar (like atoms bonded) aromatic bonds |
| npnpbond | number of non-polar bonds (like atoms bonded) in the molecule |
| pbunset | number of (i should have called it pbunsat) - number of not-fully-saturated polar bonds |
| frpbond | fraction of polar bonds in the molecule |
| dvinylb | a single bond with two adjacent vinyl groups *=*-*=*, not in an aromatic ring |
| crowding | two adjacent atoms each with 3 or more connections contribute 1 if the two highly connected atoms are separated by a 2 connected spacer atom, that contributes 0.5 |
| fcrowding | just crowding / natoms |
| ringisol | for each non-fused ring, count the number of branches off the ring. don't count terminal groups (like fluoro, of methyl) - unfortunately more complex terminal things like co2h and no2 are counted as branches off the ring rather than terminal groups. for each ring take 1/branches, and sum for all rings in the molecule. |
| fratmpie | fraction of atoms with pi electrons |
| atmpiele | number of atoms with pi electrons |
| avalcon | average connectivity of aliphatic atoms |
| fcrowdng | fraction of atoms that are crowded. crowded atoms are those that have > 2 connections and also have one or more neighbours with > 2 connections |
| avchcon | average connectivity of non ring (chain) atoms |
| faiercst | fraction of atoms involved in electron rich areas of the molecule. these include aromatic rings |
| aiercsct | number of atoms in electron rich areas of the molecule. |
| erichsct | number of separate electron rich areas of the molecule |
| unsatura | number of non-aromatic atoms that are unsaturated |
| funsatura | fraction of atoms that are unsaturated |
| nvrtspsa | polar surface area from novartis, j. med chem 2000, 43, 3714-3717 |
| atflxchn | number of atoms involved in flexible chains |
| faflxchn | fraction of atoms that are in flexible chains |
| nflxchn | the number of separate flexible chains in the molecule |
| lflxchn | the longest flexible chain in the molecule |
| fnflxchn | fraction of non-ring atoms that are in flexible chains |
| rkentrpy | entropic measure of flexibility - defined by radha |
| spchtro | number of heteroatoms in the spinach |
| nrnspch | number of non-ring non-spinach atoms. chain atoms joining rings |
| fnrnspc | fraction of non-spinach atoms that are non-ring |
| rbfrspch | fraction of bonds in the spinach that are rotatable |
| aamind | minimum bond separation between acceptors |
| aa2mind | second minimum bond distance between acceptors |
| aaave | average bond separation between acceptors |
| admind | minimum bond separation between acceptor and donor |
| ad2mind | second minimum bond separation between acceptor and donor |
| adave | average number of bonds between acceptor and donor |
| ddmind | minimum bond separation between donors |
| dd2mind | second minimum bond separation between donors |
| ddave | average bond separation between donors |
| pbcount | number of polar bonds |
| ishape | molecular connectivity descriptor: max eccentricity - min eccentricity |
| dvinylb | atoms with two adjacent vinyl bonds *=*-*-*=* |
| aromha | number of aromatic heteroatoms |
| aliphc | aliphatic carbon count |
| rmync | number of carbon atoms |
| rmynn | number of nitrogen atoms |
| rmyno | number of oxygen atoms |
| rmynf | number of fluorine atoms |
| rmyns | number of sulphur atoms |
| rmyncl | number of chlorine atoms |
| rmynbr | number of bromine atoms |
| rmyni | number of iodine atoms |
| heavy_halogen | chlorine + bromine + iodine |
| aromc | number of aromatic carbon atoms |
| maxdarom | max bond separation between aromatic atoms |
| maxdrng | max bond separation between ring atoms (not necessarily in same ring) |
| maxdhtro | max bond separation between heteroatoms |
| maxdons | max bond separation between O, N or S atoms |
| lercsct | largest electron rich section - all adjacent atoms have pi electrons |
| avebbtwn | average bond separation between atoms in the molecule |
| excybond | number exocyclic bonds |
| rssys3 | substituents 3 bonds apart in a ring system |
| rssys2 | substituents 2 bonds apart in a ring system |
| bbr4 | rings separated by 4 bonds |
| tm | terminal methyl groups |
| nrgnhlht | non ring non halogen heteroatoms |
| rssys1 | ring system substituents 1 bond apart |
| rssys2 | ring system substituents 2 bond apart |
| rssys3 | ring system substituents 3 bond apart |
| rssys4 | ring system substituents 4 bond apart |
| rssys5 | ring system substituents 5 bond apart |
| rssys6 | ring system substituents 6 bond apart |
| rssys7 | ring system substituents 7 bond apart |
| rssys8 | ring system substituents 8 bond apart |
| rssys9 | ring system substituents 9 bond apart |
| srsz | smallest ring size |
| lrsz | largest ring size |
| trmnlrng | terminal ring (1 connection) |
| intrnlrng | internal ring (multiple connections) |
| isolrc | isolated (not fused) rings |
| hperatom | average number of hydrogens per heavy atom |
| hcount | total number of hydrogen atoms |
| htroatom | number of heteroatoms |
| csp3 | number of sp3 carbon atoms |
| fcsp3 | fraction of atoms that are sp3 carbon |
| fccsp3 | fraction of carbon atoms that are sp3 |
| pyridine | aromatic nitrogen, no hydrogen |
| pyrrole | aromatic nitrogen, with hydrogen |
| ringatom | atom in a ring |
| ro5_ohnh | number of oxygen and nitrogen atoms with hydrogens |
| ro5_on | number of oxygen and nitrogen atoms |
| nnsssrng | number of non-sssr rings |
| nolp | number of atoms with lone pairs |
| nonpbond | non polar bonds |
| pbcount | polar bonds |
| frpbond | fraction of bonds that are polar |
| rsarom1 | ring substitutents one bond apart (ortho) on an aromatic ring |
| rsarom2 | ring substitutents one bond apart (meta) on an aromatic ring |
| rsarom3 | ring substitutents one bond apart (para) on an aromatic ring |
| rsaliph1 | ring substitutents one bond apart on an aliphatic ring |
| rsaliph2 | ring substitutents two bonds apart on an aliphatic ring |
| rsaliph3 | ring substitutents three bonds apart on an aliphatic ring |
| rsaliph4 | ring substitutents four or more bonds apart on an aliphatic ring |
| acmbe | binding energy (arbitrary formula from old paper) |
| al5 | aliphatic rings of size 5 |
| al6 | aliphatic rings of size 6 |
| ar5 | aromatic rings of size 5 |
| ar6 | aromatic rings of size 6 |
| alring | number of aliphatic rings |
| arring | number of aromatic rings |
| cmr | molar refractivity |
| brunsacc | fred bruns: acceptor |
| brunsdon | fred bruns: donor |
| brnsdual | fred bruns: donor and acceptor |
| brunspos | fred bruns: likely positive charge |
| brunsneg | fred bruns: likely negative charge |
| formal_charge| sum of brunspos + brunsneg. Net formal charge |
| brunshbdsum | brunsacc + brunsdon - brnsdual |
| cd4ring | carbon atoms with four connections in a ring |
| cd4chain | carbon atoms with four connections not in a ring |
| csp3_chain | sp3 carbon atoms not in a ring |
| frsub | fraction of ring atoms that are subsituted outside the ring |
| frssub | fraction of ring atoms that have a single atom subsituent |
| alorthoring | number of ortho substituents on an aliphatic ring |
| arorthoring | number of ortho substituents on an aromatic ring |
| fsatspcha | fraction of spinach atoms that are saturated |
| satspcha | number of spinach atoms that are saturated |
| unsatspcha | number of unsaturated spinach atoms |
| fsdrng5l5l | number of 5 al fused 5 al rings |
| fsdrng5l5r | number of 5 al fused 5 ar rings |
| fsdrng5l6l | number of 5 al fused 6 al rings |
| fsdrng5l6r | number of 5 al fused 6 ar rings |
| fsdrng5r5r | number of 5 ar fused 5 ar rings |
| fsdrng5r6l | number of 5 ar fused 6 al rings |
| fsdrng5r6r | number of 5 ar fused 6 ar rings |
| fsdrng6l6l | number of 6 al fused 6 al rings |
| fsdrng6l6r | number of 6 al fused 6 ar rings |
| fsdrng6r6r | number of 6 ar fused 6 ar rings |
| fsdrngalal | number of al fused al rings |
| fsdrngalar | number of al fused ar rings |
| fsdrngarar | number of ar fused ar rings |
| centre3 | number of atoms within 3 bonds of the most central atom |
| centre3h | number of heteroatoms within 3 bonds of the most central atom |
| stddcentre | standard deviation of distances from the most central atom |
| avdcentre | average distance from the most central atom |
| cntrdgncy | degeneracy of the most central atom |
| cntrdshell1 | atoms within 1 bond of the most central atom |
| cntrdshell2 | atoms within 2 bonds of the most central atom |
| cntrdshell3 | atoms within 3 bonds of the most central atom |
| bbr1 | one non ring bond between two rings |
| bbr2 | two non ring bonds between two rings |
| bbr3 | three non ring bonds between two rings |
| bbr4 | four non ring bonds between two rings |
| bbr5 | five non ring bonds between two rings |
| bbr6 | six non ring bonds between two rings |
| avsdlp | average shortest distance to longest path |
| mxsdlp | maximum shortest distance to longest path |
| mxsdlprl | avsdlp / longest_path |
| cinconjs | carbon atoms in conjugated sections |
| compact | (1 - max_eccentricity) / natoms |
| nplus | number of positively charged nitrogen atoms (bruns) |
| nminus | number of negatively charged nitrogen atoms (bruns) |
| rcj | ring chain join |
| rchj | ring chain join - to a heteroatom |
| tg3 | terminal groups separated by 3 bonds |
| rng2bridge | ring connection to chain scaffold atom |
| rng2spch | ring connection to spinach |
| hacts | hydrogen bond score (simplistic) |
| hdons | hydrogen bond score (simplistic) |
| hduals | hydrogen bond score (simplistic) |
| mh3b | most heteroatoms within 3 bonds |
| nrsyscmr | ring systems containing multiple rings |
| sboradjf | exocyclic single bonds (to terminal group) adjacent to a ring fusion |
| dboradjf | exocyclic double bonds (to terminal group) adjacent to a ring fusion |
| normbbtwn | average distance between atoms / number of atoms |
| excydscondon | exocyclic bond to a heteroatom with a hydrogen (donor) |
| excydbond | number exocyclic double bonds |
| excydsconh | exocyclic single bond to a singly connected heteroatom |
| excydscon | exocyclic single bond to a singly connected atom |
| obalance | oxygen balance |
| mdallp | mean distance of an atom from the longest path |
| fmdallp | mdallp / longest path |
| fdiffallp | measure of anisotopy between the atoms at the ends of the longest path |
| rng7atoms | the number of rings with > 7 atoms |
| rng7atoms | the number of branches in the scaffold |
| aveshell1 | for all radius 1 shells, how many atoms are included |
| aveshell2 | for all radius 2 shells, how many atoms are included |
| aveshell3 | for all radius 3 shells, how many atoms are included |
| maxshell3 | across radiud 3 shells, the max number of atoms included |
| symmatom | the number of atoms involved in a symmetry relationship |
| fsymmatom | symmatom / natoms |
| lsepsymatom | longest through bond separation between symmetric atoms |
| flsepsymatom | lsepsymatom / natoms |
| maxsymmclass | max number of atoms in a symmetry class - CF3 is 3 |
| maxpsymd | max number of atoms in a partial symmetry relationship |
| fmaxpsymd | maxpsymd / natoms |   
| maxpsymdmean | the mean number of atoms in a partial symmetry relationship |
| psymdnumzero | number of atoms not involved in a partial symmetry relationship |
| alogp | local implementation of alogp |
| xlogp | local implementation of xlogp |
| -------- | ----------- |

All are mostly interpretable molecular descriptors that are fairly quick to
compute.

## Subsets
Some descriptors are more expensive than others to compute. Not all descriptors
are needed for any invocation. Internally, groups of descriptors get computed
in functions, and some of those functions can be turned off and on. This is controlled
via the `-O` option.

By specifying `-O all` all descriptors will be computed. Using `-O none` will turn
off all optional descriptors. This can make a big difference in run times. Running
```
iwdescr.sh -O all file.smi > file.w
```
takes 5.7 seconds to process 20k molecules, generating 264 columns of ouput. Running
```
iwdescr.sh -O none file.smi > file.w
```
takes 1.2 seconds and generated 100 features. Following a `-O none` directive,
other features can be turned on by adding subsequent `-O` options. For example
`-O none -O symm` adds in symmetry related descriptors. Run time jumps to 2.1 seconds
and now 105 features are produced. It is not possible to turn off individual features.
But it is a common workflow to do something like
```
iwdescr.sh ... file.smi | iwcut -d w_natoms,w_nrings -
```

## Descriptor Names
Descriptor names are historical. Once upon a time we faced an 8 character limit
to descriptor names, and the 'w_' prefix was desirable in order to keep track
of which descriptors came from where. But changing descriptor names is very
difficult. If you do not like the current names, take a look in `Molecule_Tools/iwdescr.proto`
and you will see how to create a name translation table, that can be used via the
`-B namexref=fname` construct.


## Filtering
`iwdescr` can be used for filtering, although it is also possible to use `dfilefilter`
to filter the results in a pipeline.

For example, given an input file
```
C methane
CC ethane
CCC propane
C1CC1 cyclopropane
```
the command
```
iwdescr -F w_natoms.gt.2
```
writes
```
CCC propane
C1CC1 cyclopropane
```
Multiple filters can be added and work in an 'and' fashion.

## Fingerprinting
Most of the features in iwdescr can be turned into a fingerprint. This is enabled
by data that is embedded in the source about the likely range of each feature. This
can also be be handled by an external file in order to make keeping this maintained
easier. This distribution comes with a textproto file containing 99% ranges according
to a random sample of 1M molecules from a recent Chembl. This profile is set as
a default in the `contrib` directory. Note that if this file is changed, fingerprints
generated will also change.

Fingerprinting is controlled via the `-G` option. The output from `iwdescr.sh -G help` is
```
 -G FILTER         work as a TDT filter
 -G RE=<n>         the number of buckets used in discretising the values
 -G ALL            use all descriptors to generate the fingperint
 -G BEST           from calibration runs, certain descriptors have been designated
                   as the 'best'. Use these designated to generate the fingperint
 -G d1,d2,d3...    specify individual descriptors to be fingerprinted
 <n>               generate <n> replicates for each bit.
 If any descriptor name is followed by a ':n', that feature will include 'n'
 replicates of that bit in the output. By default, all features get the same number
```
A typical usage might be
```
iwdescr.sh -G w_natoms:10 -G w_nrings:20 ...
```
which will produce a fingerprint based on two of the features computed. This is
also available in `gfp_make` via the `-W` option.
