# Molecular Descriptors

The `lillymol_tools` module exposes molecular descriptor calculations for use
directly in Python rather than via descriptor files.

`IWDescr` is the primary high-throughput descriptor engine. It computes the same
complete descriptor set as
[`iwdescr.sh`](https://github.com/IanAWatson/LillyMol/blob/master/docs/Molecule_Tools/iwdescr.md)
and returns NumPy arrays.

`MolecularDescriptors` is a small convenience wrapper around `IWDescr` for more
discoverable Python usage. It provides descriptor names, NumPy array output,
dictionary output for one molecule, and batch output.

The Python interface does not provide options for enabling or disabling
individual `IWDescr` descriptor families.

## Requirements

The environment variable `LILLYMOL_HOME` must point to the LillyMol installation
or source tree. `IWDescr` uses it to locate the standard charge and
donor/acceptor queries:

```
export LILLYMOL_HOME=/path/to/LillyMol
```

NumPy must also be installed:

```
pip install numpy
```

For DataFrame output, pandas is also required:

```
pip install pandas
```

Construction of `IWDescr` or `MolecularDescriptors` fails if `LILLYMOL_HOME` is
not defined or the standard query files cannot be loaded.

## Basic Usage

For exploratory Python code, `MolecularDescriptors` is usually the easiest
starting point:

```python
from lillymol import MolFromSmiles
from lillymol_tools import MolecularDescriptors

mol = MolFromSmiles("CCO ethanol")
descriptors = MolecularDescriptors()

values = descriptors.compute(mol)
print(values["amw"])
print(values["nrings"])
```

`compute()` returns a normal Python `dict` keyed by descriptor name.
`compute_array()` returns the same values as a one-dimensional NumPy array, and
`names()` gives the corresponding column order.

For code that is already array-oriented, use `IWDescr` directly:

```python
from lillymol import MolFromSmiles
from lillymol_tools import IWDescr

mol = MolFromSmiles("CCO ethanol")
if mol is None:
    raise ValueError("Invalid SMILES")

iwdescr = IWDescr()
names = iwdescr.feature_names()
values = iwdescr.process(mol)
```

`feature_names()` returns a `list[str]`. `process()` returns a one-dimensional
NumPy array with `dtype=float32`.

The names and values use the same ordering:

```python
name_to_column = {name: column for column, name in enumerate(names)}

amw = values[name_to_column["amw"]]
mxdst = values[name_to_column["mxdst"]]

print(f"AMW {amw:.2f}, longest path {int(mxdst)}")
```

The feature ordering is stable for the lifetime of an `IWDescr` object. Fetch
the names once and reuse the mapping for each molecule:

```python
iwdescr = IWDescr()
name_to_column = {
    name: column for column, name in enumerate(iwdescr.feature_names())
}

for mol in molecules:
    values = iwdescr.process(mol)
    print(values[name_to_column["amw"]])
```

All descriptors are computed by every call to `process()` or
`MolecularDescriptors.compute()`, even when the application uses only a few
values.

## Scalar Helpers

Common one-off descriptor calculations are also available as free functions from
`lillymol_tools`:

```python
from lillymol import MolFromSmiles
from lillymol_tools import alogp, molecular_weight, qed_score, tpsa, xlogp

mol = MolFromSmiles("CCO ethanol")

print(molecular_weight(mol))
print(alogp(mol))
print(xlogp(mol))
print(tpsa(mol))
print(qed_score(mol))
```

`alogp`, `xlogp`, `tpsa`, and `qed_score` return `None` if the calculation
fails. `qed_score` constructs and initialises a temporary QED calculator from
`LILLYMOL_HOME`; use the `QED` class directly when computing QED for many
molecules.

## Rotatable Bonds

The `RotatableBonds` helper can return either the count of rotatable bonds or
the atom pairs that define those bonds. Atom numbers are zero-based LillyMol atom
indices. The order of atoms within each pair is the order stored in the molecule
`Bond` object and should be treated as an implementation detail.

```python
from lillymol import MolFromSmiles, RotatableBonds, EXPENSIVE

mol = MolFromSmiles("CCCCC pentane")
rotb = RotatableBonds()
rotb.set_calculation_type(EXPENSIVE)

print(rotb.rotatable_bonds(mol))       # 2
print(rotb.rotatable_bond_atoms(mol))  # [(1, 2), (2, 3)]
```

`QUICK` uses the fast LillyMol definition. `EXPENSIVE` excludes additional
non-freely-rotatable linkages such as conjugated amide-like bonds.

## Batch Processing

For collections of molecules, `process_list()` is more convenient than calling
`process()` in a loop. It accepts a list of `Molecule` objects and returns
results for all of them in one call.

### NumPy Array (default)

```python
iwdescr = IWDescr()
molecules = [...]  # list of Molecule objects

X = iwdescr.process_list(molecules)
```

`X` is a two-dimensional NumPy array of `dtype=float32` with shape
`(n_molecules, n_descriptors)`. Rows correspond to molecules in the order they
were supplied; columns correspond to descriptors in the same order as
`feature_names()`. This layout is directly compatible with scikit-learn
estimators and other tools that expect a feature matrix.

### Pandas DataFrame

Pass `as_dataframe=True` to receive a `pandas.DataFrame` instead:

```python
df = iwdescr.process_list(molecules, as_dataframe=True)
```

The DataFrame uses descriptor names as column labels and molecule names
(from `Molecule.name()`) as the row index, so individual rows and columns can
be accessed by name:

```python
print(df["amw"])             # one descriptor across all molecules
print(df.loc["ethanol"])     # all descriptors for one molecule
```

pandas is only imported when `as_dataframe=True` is requested. If pandas is not
installed and this option is used, a standard `ImportError` is raised with a
clear message.

### Choosing between the two forms

Use `process_list()` (NumPy) when feeding descriptors directly into a
machine-learning pipeline, where a plain array is what the next step expects.

Use `process_list(..., as_dataframe=True)` when exploring results interactively,
writing to CSV, or working with downstream code that benefits from named columns
and a named index:

```python
# Quick exploration
df = iwdescr.process_list(molecules, as_dataframe=True)
df.describe()
df[["amw", "nrings"]].hist()
df.to_csv("descriptors.csv")
```

## Missing Values

Some descriptors are not defined for every molecule. Undefined values are
returned as NumPy `NaN` values:

```python
import numpy

values = iwdescr.process(mol)
defined = numpy.isfinite(values)
```

The same applies to entries in the array or DataFrame returned by
`process_list()`. Rows for molecules whose calculation fails will contain NaN
for the affected descriptors. In a DataFrame these show up naturally in pandas
summary statistics and can be handled with the usual `df.dropna()` /
`df.fillna()` methods.

## Descriptor Definitions

The descriptor names and their definitions are documented in the
[`iwdescr` descriptor table](https://github.com/IanAWatson/LillyMol/blob/master/docs/Molecule_Tools/iwdescr.md#descriptors).

## Molecule and Threading Behaviour

Descriptor calculation may modify the supplied `Molecule`. In particular,
legacy chirality calculations may remove chiral-centre information. Make a copy
before calling `process()` or `process_list()` if the original molecules must
remain unchanged.

An `IWDescr` instance must not be used concurrently from multiple threads.
Separate instances can be used independently.

## Performance

The Python interface is typically about 10–15% slower than descriptor
calculation through the C++ executable. This is relatively small overhead for a
Python binding because the descriptor calculation remains in C++ and each call
returns the complete result as one NumPy array.

For large collections, prefer `process_list()` over a Python-level loop calling
`process()`. Both execute the same C++ code per molecule, but `process_list()`
allocates the output array once and avoids the per-call Python/C++ boundary
overhead of the loop.

Create one `IWDescr` object and reuse it rather than constructing a new object
for every molecule or batch.

## Complete Example

See
[`contrib/examples/descriptor_computation.py`](https://github.com/IanAWatson/LillyMol/blob/master/contrib/examples/descriptor_computation.py)
for a command-line example that reads molecules and retrieves descriptors by
name.

A batch example using `process_list`:

```python
from lillymol import MolFromSmiles
from lillymol_tools import IWDescr

smiles_list = [
    ("CCO",           "ethanol"),
    ("c1ccccc1",      "benzene"),
    ("CC(=O)O",       "acetic_acid"),
    ("c1ccccc1O",     "phenol"),
]

molecules = [
    MolFromSmiles(f"{smi} {name}")
    for smi, name in smiles_list
]

or

molecules = slurp("/path/to/file.smi")

iwdescr = IWDescr()

# As a NumPy array for ML pipelines
X = iwdescr.process_list(molecules)
print(X.shape)   # (4, n_descriptors)

# As a DataFrame for exploration
df = iwdescr.process_list(molecules, as_dataframe=True)
print(df[["amw", "nrings", "aromring"]])
```
