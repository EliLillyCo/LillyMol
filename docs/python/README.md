# LillyMol Python

Python bindings for LillyMol, built with [nanobind](https://github.com/wjakob/nanobind).

These started as a prototyping convenience and have become fast enough to use in
real tasks. LillyMol is a C++ cheminformatics toolkit that has always been built around
speed, and most of that survives the trip into python, because almost every call
does real work in C++ rather than shuffling python objects. Where we have measured
against RDKit on the same task, LillyMol python has run 3 to 7 times faster - see
[Performance](#performance).

The bindings cover the Molecule, substructure searching, reactions, fingerprints
and similarity, chemical standardisation, and a growing number of the LillyMol
command line tools. Names and signatures are still occasionally adjusted, so
check the reference if something does not behave as you expect.

## Contents

| Document | What is in it |
| -------- | ------------- |
| [Build.md](Build.md) | Python packages you need, building, and wheel packaging |
| [LillyMolPython.md](LillyMolPython.md) | The API reference - Molecule, Atom, Bond, Set_of_Atoms and Ring methods, substructure searching, reactions, GFP |
| [Tools.md](Tools.md) | LillyMol command line tools exposed to python - QED, Lipinski, MedchemWizard, unique molecules, Position3D |
| [descriptors.md](descriptors.md) | Molecular descriptors, singly and in batches, numpy and pandas forms |
| [fingerprints.md](fingerprints.md) | GFP quickstart, similarity background, and legacy counted fingerprint helpers |
| [depiction.md](depiction.md) | Optional 2D coordinate generation and stereochemical wedge assignment |
| [Ring replacement](/docs/Molecule_Tools/ring_replacement.md#python) | Generate ring-replacement analogues from Python |
| [tsubstructure.md](tsubstructure.md) | Substructure searching as the `tsubstructure` tool does it |
| [structure_database.md](structure_database.md) | Looking molecules up in a structure database |
| [synthetic_precedent.md](synthetic_precedent.md) | Synthetic precedent from a database of known molecules |
| [substituent_identification.md](substituent_identification.md) | Python lookup/generation API for databases built by `substituent_identification` |
| [4-pyridone.py](4-pyridone.py), [reaction_example.py](reaction_example.py) | Worked examples |

## Installing and verifying

See [Build.md](Build.md) for the full story. The short version, for working
inside a LillyMol checkout:

```
cd ${LILLYMOL_HOME}/src
bazel build -c opt //nanobind:all
./copy_shared_libraries.sh ../lib
```

then run python through the wrapper, which sets `PYTHONPATH` and
`LD_LIBRARY_PATH` for you:

```
${LILLYMOL_HOME}/run_python.sh my_script.py
```

Build `-c opt`. A default `fastbuild` is enormously slower and is not what any of
the timings here describe.

Check it works:

```python
from lillymol import *

m = MolFromSmiles('c1ccccc1CC(=O)N')
print(m.natoms(), m.aromatic_ring_count(), m.unique_smiles())
# 10 1 O=C(N)Cc1ccccc1
```

## First steps

`MolFromSmiles` will look familiar to RDKit users, and returns `None` on failure.
Diagnostics go to stderr; nothing but the `None` arrives in python.

```python
from lillymol import *

m = MolFromSmiles('c1ccccc1')
bad = MolFromSmiles('c1ccc')          # None, with a complaint on stderr
mols = MolFromSmiles(['C', 'CC', 'C1CC1'])   # a list form also exists
```

Reading files is usually what you want. `MolReaderContext` infers the file type from
the name:

```python
import lillymol_io

with lillymol_io.MolReaderContext('/path/to/file.smi') as reader:
    for mol in reader:
        print(mol.name(), mol.natoms())
```

`MolReaderContext` can also apply common preprocessing and skipping as molecules
are read:

```python
with lillymol_io.MolReaderContext('/path/to/file.smi', largest_fragment=True,
                               remove_chirality=True,
                               remove_isotopes=True,
                               skip_invalid_valence=True,
                               skip_non_organic=True) as reader:
    for mol in reader:
        ...
    print(reader.molecules_skipped())
```

For explicit reuse, `lillymol_io.MoleculePreprocessing` offers the same operations
with `process(mol)` and `process_copy(mol)`.

SDF names can be taken from SD data tags directly from `MolReaderContext`:

```python
with lillymol_io.MolReaderContext('/path/to/file.sdf',
                               sdf_identifier='CHEMBL_ID') as reader:
    for mol in reader:
        print(mol.name())
```

For workflows that need the SD data, use JSON names:

```python
with lillymol_io.MolReaderContext('/path/to/file.sdf',
                               sdf_tags_to_json=True,
                               all_sdf_tags=True) as reader:
    for mol in reader:
        print(mol.name())
```

There is no `None` molecule from a reader. A connection table error stops the read
unless you have raised the allowed error count. Note also that in LillyMol a
molecule with a bad valence is still a valid molecule; use
`skip_invalid_valence=True` or call `valence_ok()` if you care.

`ReaderContext` remains available as a compatibility alias for `MolReaderContext`.
`ContextWriter` remains available as a compatibility alias for `MolWriterContext`.

Molecules are always mutable. There is no distinction between a molecule you may
edit and one you may not, which removes a whole category of ceremony, at the cost
of making it your job to know when a function you called has changed its argument.

## The lazy Molecule

**This is the one thing to understand before writing anything substantial.**

A LillyMol `Molecule` computes nothing it was not asked for. Ring membership,
fragment membership, aromaticity and canonical order are all derived properties,
computed on demand and thrown away when you edit the molecule. That is why
building a molecule from smiles and asking only for `natoms()` is nearly free.

The consequence is that `Atom`, `Bond` and `Ring` objects hold no pointer back to
their `Molecule`, and so cannot know anything molecule-derived. A freshly parsed
benzene will tell you its bonds are not aromatic and not in a ring:

```python
benzene = MolFromSmiles('c1ccccc1')
bond = next(iter(benzene[0]))
bond.is_aromatic(), bond.nrings()      # (False, 0)   <- not yet computed

benzene.compute_aromaticity_if_needed()
bond.is_aromatic(), bond.nrings()      # (True, 1)
```

Note what that first line does **not** do: it does not raise. You get `False` and
`0`, which are perfectly plausible answers, and wrong. Nothing will tell you.

So if you are going to interrogate `Atom`, `Bond` or `Ring` objects, force the
perception you need first. Which call to use:

| You need | Call | Also gives you |
| -------- | ---- | -------------- |
| `Bond.nrings()`, `Bond.IsInRing()` | `mol.ring_membership()` | nothing else |
| `Bond.is_aromatic()`, `Ring.is_aromatic()` | `mol.compute_aromaticity_if_needed()` | ring membership as well |

Both are cheap once done - a second call is a flag test. Two traps worth naming:

- **`mol.nrings()` is not a force call.** It returns the ring count and does not
  push ring membership down onto the bonds. Neither does `natoms()`,
  `number_fragments()`, or anything else that merely counts.
- Any molecule-level aromaticity question, `mol.is_aromatic(0)` for instance,
  forces aromaticity as a side effect, so code that happens to ask one first
  appears to work. Do not rely on it; say what you mean.

None of this applies if you ask the `Molecule` rather than its parts.
`mol.is_aromatic(atom)`, `mol.nrings(atom)`, `mol.ring_bond_count(atom)` and
friends always give the right answer, because the molecule perceives whatever it
needs. **Asking the Molecule is both safer and faster.** Which brings us to:

## Working efficiently

`Atom` and `Bond` objects are convenient and, in python, expensive - each one you
touch has to be wrapped in a python object. The same question asked of the
`Molecule` costs one call and no wrapper.

```python
m.atomic_number(3)      # prefer this
m[3].atomic_number()    # same answer, more work
```

This is not a micro-optimisation. In a real script that walks every atom's
neighbours, replacing bond iteration with `connections()` - one call returning a
list of atom numbers - took the neighbour-walking cost from **166 us to 9 us per
molecule**, and the whole script from 12.7 s to 8.9 s.

```python
# One call per atom, no wrapper objects. bond_type_between_atoms returns the
# BondType directly, again avoiding a python Bond.
for a in range(m.natoms()):
    for nbr in m.connections(a):
        if m.bond_type_between_atoms(a, nbr) == BondType.DOUBLE_BOND:
            ...
```

Use `mol.bonds()` when you genuinely want every bond once - iterating atoms and
then their bonds visits each bond twice. Use `bond_between_atoms(a1, a2)` when you
need a real `Bond` to ask about ring membership, bond type or direction, and
`bond_type_between_atoms(a1, a2)` when the type is all you want. Both return
`None` when the atoms are not bonded.

### Copying at the boundary

The costs above are all one thing: **what gets copied or wrapped when a call
crosses between python and C++.** It is worth knowing the going rate, because it
is usually larger than the chemistry.

| crossing | cost |
| -------- | ---- |
| a trivial bound call, `natoms()`, `atomic_number()` | ~150 ns, nearly all of it the crossing |
| constructing a python wrapper for one `Bond` | enough that iterating an atom's bonds cost 166 us per molecule against 9 us for `connections()` |
| copy constructing one `Molecule` | 5.4 us |
| parsing one smiles | 14.3 us |
| a substructure search on a drug sized molecule | 6 to 18 us |

A `Molecule` copy costs about as much as the search you are about to do on it, so
a call that copies its input has doubled the price before starting. This is not
hypothetical - it is the largest performance defect we have found in these
bindings. `TSubstructure`'s batch entry points took a `std::vector<Molecule>`,
which makes the binding layer copy construct every molecule in the list, and because
that happens with the GIL held it also cannot overlap between threads. The batch call
was slower than the python loop it exists to replace, and capped at under 3x on
32 cores. Changing the signature to `std::vector<Molecule*>` - same python
signature, callers pass the same list - took it from 99k to 570k molecules per
second.

**If you are extending these bindings**, that is the rule to internalise: a
parameter of type `std::vector<SomeBoundType>` copies every element, silently, and
under the GIL. Take `std::vector<SomeBoundType*>` unless you specifically want the
caller's objects protected from mutation. Two things follow from taking pointers,
both of which the existing bindings handle explicitly and yours should too:
the binding layer converts a python `None` to `nullptr`, so an unchecked list containing
`None` will dereference it and take the interpreter down; and anything the C++
code writes to a molecule now lands on the caller's object rather than on a
discarded copy. The second is usually what you wanted - it is why a batch
`label_matched_atoms` was impossible before and works now.

**If you are using them**, the practical form of the rule is to prefer calls that
answer a whole question at once, and to be suspicious of any batch call that is
not faster than the loop it replaces.

## Translating RDKit code, and writing new code with an LLM

Translating an existing RDKit program to LillyMol python is now a realistic thing
to hand to an LLM, and we have done it successfully on non-trivial programs. Two
things to keep in mind.

**LillyMol is deliberately narrower than RDKit.** RDKit is vast and tries to cover
everything; LillyMol concentrates on core cheminformatics and is unapologetic
about it. There is no depiction or 2D coordinate generation, no conformer
generation, no force fields, and nothing like RDKit's full descriptor catalogue.
If the program you are translating leans on any of that, the translation will stop
being a translation. Check what the program actually uses before starting.

**Do not expect identical numbers.** Both toolkits are correct; they answer
slightly different questions. There is no universally accepted definition of
aromaticity, and the same is true of several other everyday quantities. Real
divergences we have measured:

| Quantity | Why the two differ |
| -------- | ------------------ |
| Aromaticity | No agreed definition. Expect the great majority of rings to agree and a small tail not to |
| H-bond donors and acceptors | LillyMol's defaults are faithful to Lipinski's paper - donors are a count of N-H and O-H *hydrogens*, so an NH2 is two. RDKit counts matching atoms, and refines which heteroatoms qualify. `rdkit_num_h_donors()` and `rdkit_num_h_acceptors()` are provided for comparison |
| Rotatable bonds | Differing treatments of conjugated linkages, sulfonamides and terminal groups |
| Molecular weight | LillyMol has no table of isotopic masses, so `amw()` erases isotopic labels and warns. `[37C]OC` weighs the same as `COC` |
| Canonical smiles | Both are canonical, and they are different canonical forms. Never compare signature *strings* across toolkits - compare the groupings they induce |

If you are having an LLM write a tool from scratch, you have a choice. Writing it
against LillyMol alone is the simplest path. Writing it against both and comparing
is more work but tells you something worth knowing: where the two agree you can
believe the answer, where they disagree you have found a definition worth
understanding rather than an average worth taking, and you get a performance
comparison for free.

When comparing, compare the right thing. Canonical smiles differ between
toolkits by construction, so equality of output strings is the wrong test. Compare
the partition - do both tools put the same molecules in the same buckets, and find
the same number of distinct answers.

## Performance

Measured on the same input, same machine, `-c opt`:

| Task | RDKit | LillyMol python |
| ---- | ----- | --------------- |
| Multi-parameter optimisation score, 20k molecules | 9.13 s | 1.37 s |
| Position-resolved chemotype extraction, 25k molecules | 27.3 s | 8.4 s |

The second is graph traversal and canonicalisation heavy, the first is descriptor
heavy, which is roughly why the ratios differ. Neither number is a benchmark
suite, and your workload may sit anywhere.

The cost you can control is the number of times you cross the C++/python boundary.
A trivial bound call - `natoms()`, `atomic_number()` - costs on the order of 150 ns,
almost all of it the crossing rather than the work. Anything that hands back a
python object costs more. Favour calls that answer a whole question at once, and
prefer the `Molecule` over its parts.

## Threads, and substructure searching

Substructure searching releases the GIL for the duration of the C++ work, so
python threads can search concurrently. What you get depends entirely on which
entry point you use. Measured on 32 cores, one query, 4000 drug sized molecules
per thread. All three below are `TSubstructure`, so they are directly comparable.

| Approach | 1 thread | 4 threads | 16 threads |
| -------- | -------- | --------- | ---------- |
| `substructure_search(mol)`, one call per molecule | 54k mol/s | 181k (3.4x) | 199k (3.7x) |
| `substructure_search(list_of_molecules)` | **57k mol/s** | 218k (3.8x) | **570k (10.0x)** |
| `substructure_search(list_of_smiles)` | 45k mol/s | 163k (3.7x) | 517k (11.6x) |

**Pass a batch.** Either batch form is a good answer - a list of the molecules you
already have, or a list of smiles if that is what you are holding. Both do the
whole loop inside one C++ call with the GIL released, so threads genuinely
overlap.

**One call per molecule stops improving at about three or four cores.** Each call
releases the GIL, but the loop around it is python and reacquires the GIL every
iteration, so you end up bounded by how fast a single python loop can run. Adding
cores past that gains nothing. It is still the simplest thing to write, and fine
if its ceiling is fast enough for you.

`num_matches(list_of_molecules)` behaves the same way - 57k on one thread, 682k on
16 - and `label_matched_atoms(list_of_molecules)` now exists at all, which it could
not before.

The batch forms used to be the slowest option rather than the fastest, for the
reason set out in [Copying at the boundary](#copying-at-the-boundary). Worth
reading if you are ever tempted to reintroduce it.

One consequence of the batch forms no longer copying: if you have set
`set_reduce_to_largest_fragment` or `set_make_implicit_hydrogens_explicit`, those
transformations now apply to **your** molecules rather than to copies the call
throws away. That is what the single molecule calls have always done, and with the
default options neither transformation runs at all.

Whichever you use, the safety contract is the same, and the GIL was previously
hiding it:

- **Give each thread its own query.** `Substructure_Query` and `TSubstructure` are
  thread compatible, not thread safe; a search accumulates match state in the
  query. Building a query from smarts costs about 7 us, so per-thread queries are
  free in practice.
- **Give each thread its own molecules.** A search *writes* to the molecule it
  searches, forcing ring and aromaticity perception on it. Two threads searching
  one shared `Molecule` is a data race, which rules out the natural-looking
  "one molecule, many queries in parallel" pattern.
- **Settings that are process wide are shared by every thread.** The way aromatic
  bonds match query bonds, `set_aromatic_bonds_lose_kekule_identity`, is one of them.
  Set it once, before the threads start searching. See
  [Aromatic Bonds and Kekule Forms](tsubstructure.md#aromatic-bonds-and-kekule-forms).

## Reference

The API reference is [LillyMolPython.md](LillyMolPython.md). Every method named in
its class-method tables is checked against the nanobind module by
`//nanobind:lillymol_doc_test`, so a renamed binding fails a test rather than
quietly misleading a reader.

```shell
cd ${LILLYMOL_HOME}/src
bazel test -c opt //nanobind:lillymol_doc_test
```

That test does not check prose or examples. If you change behaviour, the
surrounding text still needs a human.
